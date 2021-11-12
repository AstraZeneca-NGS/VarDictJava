package com.astrazeneca.vardict.modules;

import com.astrazeneca.vardict.Configuration;
import com.astrazeneca.vardict.collection.VariationMap;
import com.astrazeneca.vardict.data.*;
import com.astrazeneca.vardict.data.scopedata.Scope;
import com.astrazeneca.vardict.data.scopedata.VariationData;
import com.astrazeneca.vardict.variations.Mate;
import com.astrazeneca.vardict.variations.Sclip;
import com.astrazeneca.vardict.variations.Variation;
import htsjdk.samtools.*;

import java.time.LocalDateTime;
import java.util.*;
import java.util.regex.Matcher;

import static com.astrazeneca.vardict.data.scopedata.GlobalReadOnlyScope.instance;
import static com.astrazeneca.vardict.Utils.*;
import static com.astrazeneca.vardict.modules.RecordPreprocessor.getChrName;
import static com.astrazeneca.vardict.modules.VariationRealigner.adjInsPos;
import static com.astrazeneca.vardict.variations.VariationUtils.*;
import static com.astrazeneca.vardict.data.Patterns.*;
import static java.lang.Math.abs;

/**
 * Utility methods for parsing CIGAR from SAMRecord.
 */
public class CigarParser implements Module<RecordPreprocessor, VariationData> {
    // Utils maps
    private Map<Integer, VariationMap<String, Variation>> nonInsertionVariants;
    private Map<Integer, VariationMap<String, Variation>> insertionVariants;
    private Map<Integer, Integer> refCoverage;
    private Map<Integer, Sclip> softClips3End; // soft clipped at 3'
    private Map<Integer, Sclip> softClips5End; // soft clipped at 5'
    private final Map<Integer, Map<String, Integer>> positionToInsertionCount;
    private final Map<Integer, Map<String, Integer>> mnp; // Keep track of MNPs
    private final Map<Integer, Map<String, Integer>> positionToDeletionCount;
    private final Map<String, int[]> spliceCount;
    private final SVStructures svStructures;

    private Region region;
    private Set<String> splice;
    private Reference reference;
    private int maxReadLength;

    private Cigar cigar;
    private int start;
    private int totalReads;
    private int duplicateReads;
    private int discordantCount;
    private boolean svflag;
    private int offset;
    private int cigarElementLength;
    private int readPositionIncludingSoftClipped; // keep track the read position, including softclipped
    private int readPositionExcludingSoftClipped; // keep track the position in the alignment, excluding softclipped

    /**
     * Creates new CigarParser
     * @param svflag duprate calculates different ways due to it needed in partial pipeline or not.
     */
    public CigarParser(boolean svflag) {
        this.positionToInsertionCount = new HashMap<>();
        this.mnp = new HashMap<>();
        this.positionToDeletionCount = new HashMap<>();
        this.spliceCount = new HashMap<>();
        this.svStructures = new SVStructures();
        this.svflag = svflag;
    }

    /**
     * Parse cigar and create variations for each SAMRecord.
     * @param scope contains RecordPreprocessor (non-filtered record only)
     * @return filled variation data (with insertions and non-insertion maps of variations)
     */
    @Override
    public Scope<VariationData> process(Scope<RecordPreprocessor> scope) {
        RecordPreprocessor processor = scope.data;
        initFromScope(scope);

        if (instance().conf.y) {
            System.err.printf("TIME: Start parsing SAM: " + LocalDateTime.now() + " %s %s:%d-%d\n",
                    scope.bam, region.chr, region.start, region.end);
        }

        SAMRecord record;
        while ((record = processor.nextRecord()) != null) {
            try {
                parseCigar(getChrName(scope.region), record);
            } catch (Exception exception) {
                printExceptionAndContinue(exception, "record", record.getReadName(), scope.region);
            }
        }

        if (instance().conf.y) {
            System.err.println("TIME: Finish parsing SAM: " + LocalDateTime.now() + " " + scope.bam + " "
                    + getChrName(scope.region) + ":" + region.start + "-" + region.end);
        }

        processor.close();

        if (instance().conf.outputSplicing) {
            for (Map.Entry<String, int[]> entry : spliceCount.entrySet()) {
                String intron = entry.getKey();
                int intronCount = entry.getValue()[0];
                System.out.printf("%s\t%s\t%s\t%s\n", instance().sample, region.chr, intron, intronCount);
            }
            return new Scope<>(
                    scope.bam,
                    scope.region,
                    scope.regionRef,
                    scope.referenceResource,
                    maxReadLength,
                    splice,
                    scope.out,
                    new VariationData());
        }

        double duprate;
        if (svflag) {
            duprate = (discordantCount + 1) / (totalReads - duplicateReads + 1) > 0.5 ? 0.0 : 1.0;
        } else {
            duprate = (instance().conf.removeDuplicatedReads && totalReads != 0)
                    ? Double.parseDouble(String.format("%.3f", ((double) duplicateReads) / totalReads))
                    : 0.0;
        }

        VariationData variationData = new VariationData();
        variationData.nonInsertionVariants = getNonInsertionVariants();
        variationData.insertionVariants = getInsertionVariants();
        variationData.positionToInsertionCount = getPositionToInsertionCount();
        variationData.positionToDeletionsCount = getPositionToDeletionCount();
        variationData.refCoverage = getRefCoverage();
        variationData.softClips5End = getSoftClips5End();
        variationData.softClips3End = getSoftClips3End();
        variationData.maxReadLength = getMaxReadLength();
        variationData.mnp = getMNP();
        variationData.spliceCount = getSpliceCount();
        variationData.svStructures = getSVStructures();
        variationData.duprate = duprate;

        return new Scope<>(
                scope.bam,
                scope.region,
                scope.regionRef,
                scope.referenceResource,
                maxReadLength,
                splice,
                scope.out,
                variationData);
    }

    private void initFromScope(Scope<RecordPreprocessor> scope) {
        this.insertionVariants = scope.data.insertionVariants;
        this.nonInsertionVariants = scope.data.nonInsertionVariants;
        this.refCoverage = scope.data.refCoverage;
        this.softClips3End = scope.data.softClips3End;
        this.softClips5End = scope.data.softClips5End;
        this.region = scope.region;
        this.splice = scope.splice;
        this.reference = scope.regionRef;
        this.maxReadLength = scope.maxReadLength;
        this.duplicateReads = scope.data.duplicateReads;
        this.totalReads = scope.data.totalReads;
    }

    public Map<Integer, VariationMap<String, Variation>> getNonInsertionVariants() {
        return nonInsertionVariants;
    }

    public Map<Integer, VariationMap<String, Variation>> getInsertionVariants() {
        return insertionVariants;
    }

    public Map<Integer, Integer> getRefCoverage() {
        return refCoverage;
    }

    public Map<Integer, Sclip> getSoftClips3End() {
        return softClips3End;
    }

    public Map<Integer, Sclip> getSoftClips5End() {
        return softClips5End;
    }

    public Map<Integer, Map<String, Integer>> getPositionToInsertionCount() {
        return positionToInsertionCount;
    }

    public Map<Integer, Map<String, Integer>> getMNP() {
        return mnp;
    }

    public Map<Integer, Map<String, Integer>> getPositionToDeletionCount() {
        return positionToDeletionCount;
    }

    public Map<String, int[]> getSpliceCount() {
        return spliceCount;
    }

    public SVStructures getSVStructures() {
        return svStructures;
    }

    public int getMaxReadLength() {
        return maxReadLength;
    }

    /**
     * Modify the CIGAR for potential mis-alignment for indels at the end of reads to softclipping
     * and let VarDict's algorithm to figure out indels
     * @param chrName name of chromosome from region
     * @param record SAMRecord to process
     */
    private void parseCigar(String chrName,
                           SAMRecord record) {
        String querySequence = record.getReadString();
        int mappingQuality = record.getMappingQuality();
        cigar = record.getCigar();
        Map<Integer, Character> ref = reference.referenceSequences;

        final int insertionDeletionLength = getInsertionDeletionLength(cigar);

        // Number of mismatches from NM tag (or nM for STAR).
        int totalNumberOfMismatches = 0;
        Integer numberOfMismatches_NM = record.getIntegerAttribute(SAMTag.NM.name());
        Integer numberOfMismatches_STAR = record.getIntegerAttribute("nM");
        if (numberOfMismatches_STAR != null) {
            numberOfMismatches_NM = numberOfMismatches_STAR;
        }

        // Don't use NM since it includes gaps, which can be from indels
        if (numberOfMismatches_NM != null) {
            // Edit distance - indels is the # of mismatches
            totalNumberOfMismatches = numberOfMismatches_NM - insertionDeletionLength;
            if (totalNumberOfMismatches > instance().conf.mismatch) {
                return;
            }
        } else { //Log the read if number of mismatches is not available
            if (instance().conf.y && !record.getCigarString().equals("*")) {
                System.err.println("No NM tag for mismatches. " + record.getSAMString());
            }
            // Skip unmapped reads (issue #56) and reads without alignment
            if (record.getReadUnmappedFlag() || record.getCigarString().equals(SAMRecord.NO_ALIGNMENT_CIGAR)) {
                return;
            }
        }

        String queryQuality = getBaseQualityString(record);
        final boolean isMateReferenceNameEqual = record.getReferenceName().equals(record.getMateReferenceName());
        final int numberOfMismatches = totalNumberOfMismatches;
        boolean direction = record.getReadNegativeStrandFlag();

        if (instance().ampliconBasedCalling != null) {
            if (parseCigarWithAmpCase(record, isMateReferenceNameEqual))  {
                return;
            }
        }

        final int position;
        readPositionIncludingSoftClipped = 0;
        readPositionExcludingSoftClipped = 0;

        if (instance().conf.performLocalRealignment) {
            // Modify the CIGAR for potential mis-alignment for indels at the end of reads to softclipping and let VarDict's
            // algorithm to figure out indels
            CigarModifier cigarModifier = new CigarModifier(record.getAlignmentStart(),
                    record.getCigarString(),
                    querySequence,
                    queryQuality,
                    reference,
                    insertionDeletionLength,
                    region,
                    maxReadLength);
            ModifiedCigar modifiedCigar = cigarModifier.modifyCigar();
            position = modifiedCigar.position;
            cigar = TextCigarCodec.decode(modifiedCigar.cigar);
            querySequence = modifiedCigar.querySequence;
            queryQuality = modifiedCigar.queryQuality;
        } else {
            position = record.getAlignmentStart();
            cigar = record.getCigar();
        }

        cleanupCigar(record);
        //adjusted start position
        start = position;
        offset = 0;

        //determine discordant reads
        if (!getMateReferenceName(record).equals("=")) {
            discordantCount++;
        }

        //Ignore reads that are softclipped at both ends and both greater than 10 bp
        if (BEGIN_dig_dig_S_ANY_dig_dig_S_END.matcher(cigar.toString()).find()) {
            return;
        }
        // Only match and insertion counts toward read length
        // For total length, including soft-clipped bases
        final int readLengthIncludeMatchingAndInsertions = getMatchInsertionLength(cigar); // The read length for matched bases
        if (instance().conf.minmatch != 0 && readLengthIncludeMatchingAndInsertions < instance().conf.minmatch) {
            return;
        }
        // The total length, including soft-clipped bases
        int totalLengthIncludingSoftClipped = getSoftClippedLength(cigar);
        if (totalLengthIncludingSoftClipped > maxReadLength) { // Determine the read length
            maxReadLength = totalLengthIncludingSoftClipped;
        }

        // If supplementary alignment is present
        if (!instance().conf.samfilter.equals("0") && record.getSupplementaryAlignmentFlag()) {
            return; // Ignore the supplementary for now so that it won't skew the coverage
        }

        // Skip sites that are not in region of interest in CRISPR mode
        if (skipSitesOutRegionOfInterest()) return;

        // true if mate is in forward forection
        final boolean mateDirection = (record.getFlags() & 0x20) != 0 ? false : true;

        if (record.getReadPairedFlag() && record.getMateUnmappedFlag()) { //Mate unmapped, potential insertion
            // To be implemented
        } else if(record.getMappingQuality() > 10 && !instance().conf.disableSV) {
            // Consider high mapping quality mates only
             prepareSVStructuresForAnalysis(record, queryQuality, numberOfMismatches, direction, mateDirection,
                     position, totalLengthIncludingSoftClipped);
        }
        int mateAlignmentStart = record.getMateAlignmentStart();

        processCigar:
        //Loop over CIGAR records
        for (int ci = 0; ci < cigar.numCigarElements(); ci++) {
            if (skipOverlappingReads(record, position, direction, mateAlignmentStart)) {
                break;
            }
            // Length of segment in CIGAR
            cigarElementLength = cigar.getCigarElement(ci).getLength();

            //Letter from CIGAR
            final CigarOperator operator = getCigarOperator(cigar, ci);
            switch (operator) {
                case N: // Not matched
                    processNotMatched();
                    continue;
                case S: // Soft clipping
                    processSoftClip(chrName, record, querySequence, mappingQuality, ref, queryQuality, numberOfMismatches,
                            direction, position, totalLengthIncludingSoftClipped, ci);
                    continue;
                case H: // Hard clipping - skip
                    offset = 0;
                    continue;
                case I: { // Insertion
                    offset = 0;
                    ci = processInsertion(querySequence, mappingQuality, ref, queryQuality, numberOfMismatches,
                            direction, position, readLengthIncludeMatchingAndInsertions, ci);
                    continue;
                }
                case D: { // Deletion
                    offset = 0;
                    ci = processDeletion(querySequence, mappingQuality, ref, queryQuality, numberOfMismatches,
                            direction, readLengthIncludeMatchingAndInsertions, ci);
                    continue;
                }
                default:
                    break;
            }
            //End branch on CIGAR segment type
            // Now dealing with matching part
            int nmoff = 0;
            int moffset = 0;
            //Loop over bases of CIGAR segment
            for (int i = offset; i < cigarElementLength; i++) {
                //flag to trim reads at opt_T bases from start or end (depending on direction)
                boolean trim = isTrimAtOptTBases(direction, totalLengthIncludingSoftClipped);

                //variation string. Initialize to first base of the read sequence
                final char ch1 = querySequence.charAt(readPositionIncludingSoftClipped);
                String s = String.valueOf(ch1);
                boolean startWithDeletion = false;
                //skip if base is unknown
                if (ch1 == 'N') {
                    if (instance().conf.includeNInTotalDepth) {
                        incCnt(refCoverage, start, 1);
                    }
                    start++;
                    readPositionIncludingSoftClipped++;
                    readPositionExcludingSoftClipped++;
                    continue;
                }

                //sum of qualities for bases
                double q = queryQuality.charAt(readPositionIncludingSoftClipped) - 33;
                //number of bases for quality calculation
                int qbases = 1;
                //number of bases in insertion for quality calculation
                int qibases = 0;
                // for more than one nucleotide mismatch
                StringBuilder ss = new StringBuilder();
                // More than one mismatches will only perform when all nucleotides have queryQuality > $GOODQ
                // Update: Forgo the queryQuality check. Will recover later

                /*
                Condition:
                1). start + 1 is in region of interest
                2). base index is not more than segment length
                3). reference sequence has base at start
                4). base at n in read is not equal to base in reference sequence
                 */
                while ((start + 1) >= region.start
                        && (start + 1) <= region.end && (i + 1) < cigarElementLength
                        && q >= instance().conf.goodq
                        && isHasAndNotEquals(ref, start, querySequence, readPositionIncludingSoftClipped)
                        && isNotEquals('N', ref.get(start))) {
                    //Require higher quality for MNV
                    if (queryQuality.charAt(readPositionIncludingSoftClipped + 1) - 33 < instance().conf.goodq + 5) {
                        break;
                    }
                    //Break if base is unknown in the read
                    char nuc = querySequence.charAt(readPositionIncludingSoftClipped + 1);
                    if (nuc == 'N') {
                        break;
                    }
                    if (isHasAndEquals('N', ref, start + 1)) {
                        break;
                    }

                    //Condition: base at n + 1 does not match reference base at start + 1
                    if (isNotEquals(ref.get(start + 1), nuc)) {
                        //append the base from read
                        ss.append(nuc);
                        //add quality to total sum
                        q += queryQuality.charAt(readPositionIncludingSoftClipped + 1) - 33;
                        //increase number of bases
                        qbases++;
                        //shift read position by 1
                        readPositionIncludingSoftClipped++;
                        readPositionExcludingSoftClipped++;
                        i++;
                        //shift reference position by 1
                        start++;
                        nmoff++;
                    } else { // look ahead to see whether there're more mismatches, if yes, add reverence to grow MNV
                        int ssn = 0;
                        for (int ssi = 1; ssi <= instance().conf.vext; ssi++) {
                            if (i + 1 + ssi >= cigarElementLength){
                                break;
                            }
                            if (readPositionIncludingSoftClipped + 1 + ssi < querySequence.length()
                                    && isHasAndNotEquals(querySequence.charAt(readPositionIncludingSoftClipped + 1 + ssi), ref, start + 1 + ssi)) {
                                ssn = ssi + 1;
                                break;
                            }
                        }
                        if (ssn == 0) {
                            break;
                        }
                        //Require higher quality for MNV
                        if (queryQuality.charAt(readPositionIncludingSoftClipped + ssn) - 33 < instance().conf.goodq + 5) {
                            break;
                        }
                        for (int ssi = 1; ssi <= ssn; ssi++) {
                            ss.append(querySequence.charAt(readPositionIncludingSoftClipped + ssi));
                            q += queryQuality.charAt(readPositionIncludingSoftClipped + ssi) - 33;
                            qbases++;
                        }
                        readPositionIncludingSoftClipped += ssn;
                        readPositionExcludingSoftClipped += ssn;
                        i += ssn;
                        start += ssn;
                    }
                }

                //If multi-base mismatch is found, append it to s after '&'
                if (ss.length() > 0) {
                    s += "&" + ss;
                }
                int ddlen = 0;

                /*
                Condition:
                1). index $i is no farther than $VEXT from end of CIGAR segment
                2). CIGAR string contains next entry
                3). next CIGAR entry is a deletion
                4). reference sequence contains a base at $start
                5). either multi-nucleotide mismatch is found or read base at $n does not match reference base
                6). read base at $n has good quality
                 */
                if (isCloserThenVextAndGoodBase(querySequence, ref, queryQuality, ci, i, ss, CigarOperator.D)) {
                    //loop until end of CIGAR segments
                    while (i + 1 < cigarElementLength) {
                        //append next base to s and add its quality to q
                        s += querySequence.charAt(readPositionIncludingSoftClipped + 1);
                        q += queryQuality.charAt(readPositionIncludingSoftClipped + 1) - 33;
                        //increase number of bases
                        qbases++;

                        //shift read and reference indices by 1
                        i++;
                        readPositionIncludingSoftClipped++;
                        readPositionExcludingSoftClipped++;
                        start++;
                    }

                    //remove '&' delimiter from s
                    s = s.replaceFirst("&", "");
                    //prepend s with deletion of length of next CIGAR segment + '&' delimiter
                    s = "-" + cigar.getCigarElement(ci + 1).getLength() + "&" + s;
                    startWithDeletion = true;
                    ddlen = cigar.getCigarElement(ci + 1).getLength();
                    ci += 1;

                    //If CIGAR has insertion two segments ahead
                    if (isTwoInsertionsAhead(ci)) {
                        //append '^' + next-next segment sequence
                        s += "^" + substr(querySequence, readPositionIncludingSoftClipped + 1, cigar.getCigarElement(ci + 1).getLength());

                        //Loop over next-next segment
                        int nextLen = cigar.getCigarElement(ci + 1).getLength();
                        for (int qi = 1; qi <= nextLen; qi++) {
                            //add base quality to total quality
                            q += queryQuality.charAt(readPositionIncludingSoftClipped + 1 + qi) - 33;
                            //increase number of insertion bases
                            qibases++;
                        }

                        //adjust read position by length of next-next segment
                        readPositionIncludingSoftClipped += nextLen;
                        readPositionExcludingSoftClipped += nextLen;
                        ci += 1;
                    }
                    if (isNextAfterNumMatched(ci, 1)) {
                        Offset tpl = findOffset(start + ddlen + 1,
                                readPositionIncludingSoftClipped + 1, cigar.getCigarElement(ci + 1).getLength(), querySequence,
                                queryQuality, ref, refCoverage);
                        int toffset = tpl.offset;
                        if (toffset != 0) {
                            moffset = toffset;
                            nmoff += tpl.offsetNumberOfMismatches;
                            s += "&" + tpl.sequence;
                            String tq = tpl.qualitySequence;
                            for (int qi = 0; qi < tq.length(); qi++) {
                                q += tq.charAt(qi) - 33;
                                qibases++;
                            }
                        }
                    }
                } else if (isCloserThenVextAndGoodBase(querySequence, ref, queryQuality, ci, i, ss, CigarOperator.I)) {
                    while (i + 1 < cigarElementLength) {
                        s += querySequence.charAt(readPositionIncludingSoftClipped + 1);
                        q += queryQuality.charAt(readPositionIncludingSoftClipped + 1) - 33;
                        //increase number of bases
                        qbases++;

                        //shift read and reference indices by 1
                        i++;
                        readPositionIncludingSoftClipped++;
                        readPositionExcludingSoftClipped++;
                        start++;
                    }
                    //remove '&' delimiter from s
                    s = s.replaceFirst("&", "");
                    int nextLen = cigar.getCigarElement(ci + 1).getLength();
                    s += substr(querySequence, readPositionIncludingSoftClipped + 1, nextLen);
                    s = substr(s, 0, nextLen) + "&" + s.substring(nextLen);
                    s = "+" + s;

                    //Loop over next-next segment
                    for (int qi = 1; qi <= nextLen; qi++) {
                        //add base quality to total quality
                        q += queryQuality.charAt(readPositionIncludingSoftClipped + 1 + qi) - 33;
                        //increase number of insertion bases
                        qibases++;
                    }
                    readPositionIncludingSoftClipped += nextLen;
                    readPositionExcludingSoftClipped += nextLen;
                    ci += 1;
                    qibases--;
                    qbases++; // need to add to set the correction insertion position
                }
                if (!trim) {
                    //If start - qbases + 1 is in region of interest
                    final int pos = start - qbases + 1;
                    if (pos >= region.start && pos <= region.end && !s.contains("N")) {
                        addVariationForMatchingPart(mappingQuality, numberOfMismatches,
                                direction, readLengthIncludeMatchingAndInsertions, nmoff, s, startWithDeletion,
                                q, qbases, qibases, ddlen, pos);
                    }
                }

                //If variation starts with a deletion ('-' character)
                if (startWithDeletion) {
                    start += ddlen;
                }

                //Shift reference position by 1 if CIGAR segment is not insertion
                if (operator != CigarOperator.I) {
                    start++;
                }
                //Shift read position by 1 if CIGAR segment is not deletion
                if (operator != CigarOperator.D) {
                    readPositionIncludingSoftClipped++;
                    readPositionExcludingSoftClipped++;
                }
                // Skip read if it is overlap
                if (skipOverlappingReads(record, position, direction, mateAlignmentStart)) {
                    break processCigar;
                }
            }
            if (moffset != 0) {
                offset = moffset;
                readPositionIncludingSoftClipped += moffset;
                start += moffset;
                readPositionExcludingSoftClipped += moffset;
            }
            if (start > region.end) { //end if reference position is outside region of interest
                break;
            }
        }
    }

    private boolean isTrimAtOptTBases(boolean direction, int totalLengthIncludingSoftClipped) {
        if (instance().conf.trimBasesAfter != 0) {
            if (!direction) {
                return readPositionIncludingSoftClipped > instance().conf.trimBasesAfter;
            } else {
                return totalLengthIncludingSoftClipped - readPositionIncludingSoftClipped > instance().conf.trimBasesAfter;
            }
        }
        return false;
    }

    private boolean skipSitesOutRegionOfInterest() {
        int cutSite = instance().conf.crisprCuttingSite;
        int filterBp = instance().conf.crisprFilteringBp;
        if (cutSite != 0) {
            //The total aligned length, excluding soft-clipped bases and insertions
            int rlen3 = sum(globalFind(ALIGNED_LENGTH_MD, cigar.toString()));

            if (filterBp != 0) {
                return !(cutSite - start > filterBp && start + rlen3 - cutSite > filterBp);
            }
        }
        return false;
    }

    /**
     * Process CIGAR deletions part. Will ignore indels next to introns and create Variations for deletions
     */
    private int processDeletion(String querySequence, int mappingQuality, Map<Integer, Character> ref,
                                String queryQuality, int numberOfMismatches, boolean direction,
                                int readLengthIncludeMatchingAndInsertions, int ci) {
        // Ignore deletions right after introns at exon edge in RNA-seq
        if (skipIndelNextToIntron(cigar, ci)) {
            readPositionExcludingSoftClipped += cigarElementLength;
            return ci;
        }

        // $s description string of deleted segment
        StringBuilder descStringOfDeletedElement = new StringBuilder("-").append(cigarElementLength);
        // $ss sequence to be appended if next segment is matched
        StringBuilder sequenceToAppendIfNextSegmentMatched = new StringBuilder();
        // $q1 quality of last base before deletion
        char qualityOfLastSegmentBeforeDel = queryQuality.charAt(readPositionIncludingSoftClipped - 1);
        // $q quality of this segment
        StringBuilder qualityOfSegment = new StringBuilder();

        // For multiple indels within $VEXT bp
        //offset for read position if next segment is matched
        int multoffs = 0;
        //offset for reference position if next segment is matched
        int multoffp = 0;
        int nmoff = 0;

        /*
        Condition:
        1). CIGAR string has next entry
        2). length of next CIGAR segment is less than conf.vext
        3). next segment is matched
        4). CIGAR string has one more entry after next one
        5). this entry is insertion or deletion
         */
        if (isInsertionOrDeletionWithNextMatched(ci)) {
            int mLen = cigar.getCigarElement(ci + 1).getLength();
            int indelLen = cigar.getCigarElement(ci + 2).getLength();

            int begin = readPositionIncludingSoftClipped;
            appendSegments(querySequence, queryQuality, ci, descStringOfDeletedElement, qualityOfSegment,
                    mLen, indelLen, begin, false);

            //add length of next segment to both read and reference offsets
            //add length of next-next segment to reference position (for insertion) or to read position(for deletion)
            multoffs += mLen + (cigar.getCigarElement(ci + 2).getOperator() == CigarOperator.D ? indelLen : 0);
            multoffp += mLen + (cigar.getCigarElement(ci + 2).getOperator() == CigarOperator.I ? indelLen : 0);
            if (isNextAfterNumMatched(ci, 3)) {
                int vsn = 0;
                int tn = readPositionIncludingSoftClipped + multoffp;
                int ts = start + multoffs + cigarElementLength;
                for (int vi = 0; vsn <= instance().conf.vext && vi < cigar.getCigarElement(ci + 3).getLength(); vi++) {
                    if (querySequence.charAt(tn + vi) == 'N') {
                        break;
                    }
                    if (queryQuality.charAt(tn + vi) - 33 < instance().conf.goodq) {
                        break;
                    }
                    if (isHasAndEquals('N', ref, ts + vi)) {
                        break;
                    }
                    Character refCh = ref.get(ts + vi);
                    if (refCh != null) {
                        if (isNotEquals(querySequence.charAt(tn + vi), refCh)) {
                            offset = vi + 1;
                            nmoff++;
                            vsn = 0;
                        } else {
                            vsn++;
                        }
                    }
                }
                if (offset != 0) {
                    sequenceToAppendIfNextSegmentMatched.append(substr(querySequence, tn, offset));
                    qualityOfSegment.append(substr(queryQuality, tn, offset));
                }
            }
            // skip next 2 CIGAR segments
            ci += 2;
        } else if (isNextInsertion(ci)) {
            /*
            Condition:
            1). CIGAR string has next entry
            2). next CIGAR segment is an insertion
             */
            int insLen = cigar.getCigarElement(ci + 1).getLength();
            //Append '^' + next segment (inserted)
            descStringOfDeletedElement.append("^").append(substr(querySequence, readPositionIncludingSoftClipped, insLen));
            //Append next segement to quality string
            qualityOfSegment.append(substr(queryQuality, readPositionIncludingSoftClipped, insLen));

            //Shift reference position by length of next segment
            //skip next CIGAR segment
            multoffp += insLen;
            if (isNextAfterNumMatched(ci, 2)) {
                int mLen = cigar.getCigarElement(ci + 2).getLength();
                int vsn = 0;
                int tn = readPositionIncludingSoftClipped + multoffp;
                int ts = start + cigarElementLength;
                for (int vi = 0; vsn <= instance().conf.vext && vi < mLen; vi++) {
                    char seqCh = querySequence.charAt(tn + vi);
                    if (seqCh == 'N') {
                        break;
                    }
                    if (queryQuality.charAt(tn + vi) - 33 < instance().conf.goodq) {
                        break;
                    }
                    Character refCh = ref.get(ts + vi);
                    if (refCh != null) {
                        if (isEquals('N', refCh)) {
                            break;
                        }
                        if (isNotEquals(seqCh, refCh)) {
                            offset = vi + 1;
                            nmoff++;
                            vsn = 0;
                        } else {
                            vsn++;
                        }
                    }
                }
                if (offset != 0) {
                    sequenceToAppendIfNextSegmentMatched.append(substr(querySequence, tn, offset));
                    qualityOfSegment.append(substr(queryQuality, tn, offset));
                }
            }
            ci += 1;
        } else {
            /*
            Condition:
            1). CIGAR string has next entry
            2). next CIGAR segment is matched
             */
            if (isNextMatched(ci)) {
                int mLen = cigar.getCigarElement(ci + 1).getLength();
                int vsn = 0;
                //Loop over next CIGAR segment (no more than conf.vext bases ahead)
                for (int vi = 0; vsn <= instance().conf.vext && vi < mLen; vi++) {
                    char seqCh = querySequence.charAt(readPositionIncludingSoftClipped + vi);
                    //If base is unknown, exit loop
                    if (seqCh == 'N') {
                        break;
                    }
                    //If base quality is less than $GOODQ, exit loop
                    if (queryQuality.charAt(readPositionIncludingSoftClipped + vi) - 33 < instance().conf.goodq) {
                        break;
                    }
                    //If reference sequence has base at this position and it matches read base, update offset
                    Character refCh = ref.get(start + cigarElementLength + vi);
                    if (refCh != null) {
                        if (isEquals('N', refCh)) {
                            break;
                        }
                        if (isNotEquals(seqCh, refCh)) {
                            offset = vi + 1;
                            nmoff++;
                            vsn = 0;
                        } else {
                            vsn++;
                        }
                    }
                }

                //If next CIGAR segment has good matching base
                if (offset != 0) {
                    //Append first offset bases of next segment to ss and q
                    sequenceToAppendIfNextSegmentMatched.append(substr(querySequence, readPositionIncludingSoftClipped, offset));
                    qualityOfSegment.append(substr(queryQuality, readPositionIncludingSoftClipped, offset));
                }
            }
        }
        //offset should be reset to 0 on every loop, so it is non-zero only if previous next CIGAR segment is matched
        // Append '&' and ss to s if next segment has good matching base
        if (offset > 0) {
            descStringOfDeletedElement.append("&").append(sequenceToAppendIfNextSegmentMatched);
        }

        //quality of first matched base after deletion
        //append best of $q1 and $q2
        if (readPositionIncludingSoftClipped + offset >= queryQuality.length()) {
            qualityOfSegment.append(qualityOfLastSegmentBeforeDel);
        } else {
            char qualityOfSegmentWithOffset = queryQuality.charAt(readPositionIncludingSoftClipped + offset); // $q2
            qualityOfSegment.append(qualityOfLastSegmentBeforeDel > qualityOfSegmentWithOffset
                    ? qualityOfLastSegmentBeforeDel
                    : qualityOfSegmentWithOffset);
        }

        //If reference position is inside region of interest
        if (start >= region.start && start <= region.end) {
            addVariationForDeletion(mappingQuality,
                    numberOfMismatches, direction, readLengthIncludeMatchingAndInsertions,
                    descStringOfDeletedElement, qualityOfSegment, nmoff);
        }

        //adjust reference position by offset + multoffs
        start += cigarElementLength + offset + multoffs;

        //adjust read position by m (CIGAR segment length) + offset + multoffp
        readPositionIncludingSoftClipped += offset + multoffp;
        readPositionExcludingSoftClipped += offset + multoffp;
        return ci;
    }

    private boolean isNextAfterNumMatched(int ci, int number) {
        return cigar.numCigarElements() > ci + number && cigar.getCigarElement(ci + number).getOperator() == CigarOperator.M;
    }

    private boolean isTwoInsertionsAhead(int ci) {
        return cigar.numCigarElements() > ci + 1 && cigar.getCigarElement(ci + 1).getOperator() == CigarOperator.I;
    }

    private boolean isNextInsertion(int ci) {
        return instance().conf.performLocalRealignment && cigar.numCigarElements() > ci + 1
                && cigar.getCigarElement(ci + 1).getOperator() == CigarOperator.I;
    }

    private boolean isNextMatched(int ci) {
        return instance().conf.performLocalRealignment && cigar.numCigarElements() > ci + 1
                && cigar.getCigarElement(ci + 1).getOperator() == CigarOperator.M;
    }

    private boolean isInsertionOrDeletionWithNextMatched(int ci) {
        return instance().conf.performLocalRealignment && cigar.numCigarElements() > ci + 2
                && cigar.getCigarElement(ci + 1).getLength() <= instance().conf.vext
                && cigar.getCigarElement(ci + 1).getOperator() == CigarOperator.M
                && (cigar.getCigarElement(ci + 2).getOperator() == CigarOperator.I
                || cigar.getCigarElement(ci + 2).getOperator() == CigarOperator.D)
                && cigar.getCigarElement(ci + 3).getOperator() != CigarOperator.I
                && cigar.getCigarElement(ci + 3).getOperator() != CigarOperator.D;
    }

    private boolean isCloserThenVextAndGoodBase(String querySequence, Map<Integer, Character> ref, String queryQuality,
                                                int ci, int i, StringBuilder ss, CigarOperator cigarOperator) {
        // Do not adjust complex if we have hard-clips after insertion/deletion
        if (cigar.numCigarElements() > ci + 2 && cigar.getCigarElement(ci + 2).getOperator() == CigarOperator.H) {
            return false;
        }
        return instance().conf.performLocalRealignment && cigarElementLength - i <= instance().conf.vext
                && cigar.numCigarElements() > ci + 1
                && cigar.getCigarElement(ci + 1).getOperator() == cigarOperator
                && ref.containsKey(start)
                && (ss.length() > 0 || isNotEquals(querySequence.charAt(readPositionIncludingSoftClipped), ref.get(start)))
                && queryQuality.charAt(readPositionIncludingSoftClipped) - 33 >= instance().conf.goodq;
    }

    /**
     * Process CIGAR insertion part. Will ignore indels next to introns and create Variations for insertions
     */
    private int processInsertion(String querySequence, int mappingQuality, Map<Integer, Character> ref,
                                 String queryQuality, int numberOfMismatches, boolean direction, int position,
                                 int readLengthIncludeMatchingAndInsertions, int ci) {
        // Ignore insertions right after introns at exon edge in RNA-seq
        if (skipIndelNextToIntron(cigar, ci)) {
            readPositionIncludingSoftClipped += cigarElementLength;
            return ci;
        }

        //Inserted segment of read sequence
        StringBuilder descStringOfInsertionSegment = new StringBuilder(substr(querySequence, readPositionIncludingSoftClipped, cigarElementLength));

        //Quality of this segment
        StringBuilder qualityString = new StringBuilder(substr(queryQuality, readPositionIncludingSoftClipped, cigarElementLength));
        //Sequence to be appended if next segment is matched
        String ss = "";

        // For multiple indels within 10bp
        //Offset for read position if next segment is matched
        int multoffs = 0;
        //offset for reference position if next segment is matched
        int multoffp = 0;
        int nmoff = 0;

        /*
        Condition:
        1). CIGAR string has next entry
        2). length of next CIGAR segment is less than conf.vext
        3). next segment is matched
        4). CIGAR string has one more entry after next one
        5). this entry is insertion or deletion
         */
        if (isInsertionOrDeletionWithNextMatched(ci)) {
            int mLen = cigar.getCigarElement(ci + 1).getLength();
            int indelLen = cigar.getCigarElement(ci + 2).getLength();

            int begin = readPositionIncludingSoftClipped + cigarElementLength;
            appendSegments(querySequence, queryQuality, ci, descStringOfInsertionSegment, qualityString,
                    mLen, indelLen, begin, true);
            //add length of next segment to both multoffs and multoffp
            //add length of next-next segment to multoffp (for insertion) or to multoffs (for deletion)
            multoffs += mLen + (cigar.getCigarElement(ci + 2).getOperator() == CigarOperator.D ? indelLen : 0);
            multoffp += mLen + (cigar.getCigarElement(ci + 2).getOperator() == CigarOperator.I ? indelLen : 0);

            int ci6 = cigar.numCigarElements() > ci + 3 ? cigar.getCigarElement(ci + 3).getLength() : 0;
            if (ci6 != 0 && cigar.getCigarElement(ci + 3).getOperator()  == CigarOperator.M) {
                Offset tpl = findOffset(start + multoffs,
                        readPositionIncludingSoftClipped + cigarElementLength + multoffp,
                        ci6, querySequence, queryQuality, ref, refCoverage);
                offset = tpl.offset;
                ss = tpl.sequence;
                qualityString.append(tpl.qualitySequence);
            }
            //skip 2 CIGAR segments
            ci += 2;
        } else {
            /*
            Condition:
            1). CIGAR string has next entry
            2). next CIGAR segment is matched
             */
            if (isNextMatched(ci)) {
                int vsn = 0;
                //Loop over next CIGAR segment (no more than conf.vext bases ahead)
                for (int vi = 0; vsn <= instance().conf.vext && vi < cigar.getCigarElement(ci + 1).getLength(); vi++) {
                    //If base is unknown, exit loop
                    if (querySequence.charAt(readPositionIncludingSoftClipped + cigarElementLength + vi) == 'N') {
                        break;
                    }
                    //If base quality is less than conf.goodq, exit loop
                    if (queryQuality.charAt(readPositionIncludingSoftClipped + cigarElementLength + vi) - 33 < instance().conf.goodq) {
                        break;
                    }
                    //If reference sequence has base at this position and it matches read base, update offset
                    if (ref.containsKey(start + vi)) {
                        if (isNotEquals(querySequence.charAt(readPositionIncludingSoftClipped + cigarElementLength + vi), ref.get(start + vi))) {
                            offset = vi + 1;
                            nmoff++;
                            vsn = 0;
                        } else {
                            vsn++;
                        }
                    }
                }
                if (offset != 0) { //If next CIGAR segment has good matching base
                    //Append first offset bases of next segment to ss and q
                    ss += substr(querySequence, readPositionIncludingSoftClipped + cigarElementLength, offset);
                    qualityString.append(substr(queryQuality, readPositionIncludingSoftClipped + cigarElementLength, offset));
                    //Increase coverage for positions corresponding to first offset bases of next segment
                    for (int osi = 0; osi < offset; osi++) {
                        incCnt(refCoverage, start + osi, 1);
                    }
                }
            }
        }

        // Offset should be reset to 0 on every loop, so it is non-zero only if previous part of
        // code was executed
        // Append '&' and $ss to s if next segment has good matching base
        if (offset > 0) {
            descStringOfInsertionSegment.append("&").append(ss);
        }

        //If start of the segment is within region of interest and the segment does not have unknown bases
        if (start - 1 >= region.start && start - 1 <= region.end && !descStringOfInsertionSegment.toString().contains("N")) {
            int insertionPosition = start - 1;
            Matcher mm = BEGIN_ATGC_END.matcher(descStringOfInsertionSegment);
            if (mm.find()) {
                BaseInsertion tpl = adjInsPos(start - 1, descStringOfInsertionSegment.toString(), ref);
                int adjustedInsertionPosition = tpl.baseInsert;
                StringBuilder adjustedDescStringOfInsertionSegment = new StringBuilder(tpl.insertionSequence);

                // check for the index, if it is negative, there can be deletions, will not be adjusted.
                if (readPositionIncludingSoftClipped - 1 - (start - 1 - adjustedInsertionPosition) > 0) {
                    insertionPosition = adjustedInsertionPosition;
                    descStringOfInsertionSegment = adjustedDescStringOfInsertionSegment;
                }
            }
            //add '+' + s to insertions at inspos
            incCnt(getOrElse(positionToInsertionCount, insertionPosition, new HashMap<>()), "+" + descStringOfInsertionSegment, 1);

            //add insertion to table of variations
            Variation hv = getVariation(insertionVariants, insertionPosition, "+" + descStringOfInsertionSegment); //variant structure for this insertion
            hv.incDir(direction);
            //add count
            hv.varsCount++;
            //minimum of positions from start of read and end of read
            int tp = readPositionExcludingSoftClipped < readLengthIncludeMatchingAndInsertions - readPositionExcludingSoftClipped
                    ? readPositionExcludingSoftClipped + 1
                    : readLengthIncludeMatchingAndInsertions - readPositionExcludingSoftClipped;

            //mean read quality of the segment
            double tmpq = 0;
            for (int i = 0; i < qualityString.length(); i++) {
                tmpq += qualityString.charAt(i) - 33;
            }
            tmpq = tmpq / qualityString.length();

            //pstd is a flag that is 1 if the variant is covered by at least 2 read segments with different positions
            if (!hv.pstd && hv.pp != 0 && tp != hv.pp) {
                hv.pstd = true;
            }
            //qstd is a flag that is 1 if the variant is covered by at least 2 segment reads with different qualities
            if (!hv.qstd && hv.pq != 0 && tmpq != hv.pq) {
                hv.qstd = true;
            }
            hv.meanPosition += tp;
            hv.meanQuality += tmpq;
            hv.meanMappingQuality += mappingQuality;
            hv.pp = tp;
            hv.pq = tmpq;
            if (tmpq >= instance().conf.goodq) {
                hv.highQualityReadsCount++;
            } else {
                hv.lowQualityReadsCount++;
            }
            hv.numberOfMismatches += numberOfMismatches - nmoff;

            // Adjust the reference count for insertion reads
            /*
            Condition:
            1). reference sequence has base for the position
            2). hash contains variant structure for the position
            3). read base at position n-1 matches reference at start-1
             */
            int indexOfInsertionInQuerySequence = readPositionIncludingSoftClipped - 1 - (start - 1 - insertionPosition);
            if (insertionPosition > position && isHasAndEquals(querySequence.charAt(indexOfInsertionInQuerySequence), ref, insertionPosition)) {
                Variation tv = getVariationMaybe(nonInsertionVariants, insertionPosition, querySequence.charAt(indexOfInsertionInQuerySequence));
                //Substract count.
                if (tv != null) {
                    subCnt(tv, direction, tp, queryQuality.charAt(indexOfInsertionInQuerySequence) - 33,
                            mappingQuality, numberOfMismatches - nmoff);
                }
            }
            // Adjust count if the insertion is at the edge so that the AF won't > 1
            /*
            Condition:
            1). looking at second segment in CIGAR string
            2). first segment is a soft-clipping or a hard-clipping
            */
            if (ci == 1 && (cigar.getCigarElement(0).getOperator() == CigarOperator.S
                    || cigar.getCigarElement(0).getOperator() == CigarOperator.H)) {
                //Add one more variant corresponding to base at start - 1 to hash
                Variation ttref = getVariation(nonInsertionVariants, insertionPosition, ref.get(insertionPosition).toString());
                ttref.incDir(direction);
                ttref.varsCount++;
                ttref.pstd = hv.pstd;
                ttref.qstd = hv.qstd;
                ttref.meanPosition += tp;
                ttref.meanQuality += tmpq;
                ttref.meanMappingQuality += mappingQuality;
                ttref.pp = tp;
                ttref.pq = tmpq;
                ttref.numberOfMismatches += numberOfMismatches - nmoff;
                incCnt(refCoverage, insertionPosition, 1);
            }
        }

        //adjust read position by m (CIGAR segment length) + offset + multoffp
        readPositionIncludingSoftClipped += cigarElementLength + offset + multoffp;
        readPositionExcludingSoftClipped += cigarElementLength + offset + multoffp;
        //adjust reference position by offset + multoffs
        start += offset + multoffs;
        return ci;
    }

    /**
     * Process CIGAR soft-clipped part. Will ignore large soft-clips and create Variations for mis-softclipping reads
     * due to alignment
     */
    private void processSoftClip(String chrName, SAMRecord record, String querySequence, int mappingQuality,
                                 Map<Integer, Character> ref, String queryQuality, int numberOfMismatches,
                                 boolean direction, int position, int totalLengthIncludingSoftClipped, int ci) {
        //First record in CIGAR
        if (ci == 0) { // 5' soft clipped
            // Ignore large soft clip due to chimeric reads in library construction
            if (!instance().conf.chimeric) {
                String saTagString = record.getStringAttribute(SAMTag.SA.name());

                if (cigarElementLength >= 20 && saTagString != null) {
                    if (isReadChimericWithSA(record, position, saTagString, direction, true)) {
                        readPositionIncludingSoftClipped += cigarElementLength;
                        offset = 0;
                        // Had to reset the start due to softclipping adjustment
                        start = position;

                        return;
                    }
                    //trying to detect chimeric reads even when there's no supplementary
                    // alignment from aligner
                } else if (cigarElementLength >= Configuration.SEED_1) {
                    Map<String, List<Integer>> referenceSeedMap = reference.seed;
                    String sequence = getReverseComplementedSequence(record, 0, cigarElementLength);
                    String reverseComplementedSeed = sequence.substring(0, Configuration.SEED_1);

                    if (referenceSeedMap.containsKey(reverseComplementedSeed)) {
                        List<Integer> positions = referenceSeedMap.get(reverseComplementedSeed);
                        if (positions.size() == 1 &&
                                abs(start - positions.get(0)) <  2 * maxReadLength) {
                            readPositionIncludingSoftClipped += cigarElementLength;
                            offset = 0;
                            // Had to reset the start due to softclipping adjustment
                            start = position;
                            if (instance().conf.y) {
                                System.err.println(sequence + " at 5' is a chimeric at "
                                        + start + " by SEED " + Configuration.SEED_1);
                            }
                            return;
                        }
                    }
                }
            }
            // Align softclipped but matched sequences due to mis-softclipping
            /*
            Conditions:
            1). segment length > 1
            2). start between 0 and chromosome length,
            3). reference genome is known at this position
            4). reference and read bases match
            5). read quality is more than 10
             */
            while (cigarElementLength - 1 >= 0 && start - 1 > 0 && start - 1 <= instance().chrLengths.get(chrName)
                    && isHasAndEquals(querySequence.charAt(cigarElementLength - 1), ref, start - 1)
                    && queryQuality.charAt(cigarElementLength - 1) - 33 > 10) {
                //create variant if it is not present
                Variation variation = getVariation(nonInsertionVariants, start - 1, ref.get(start - 1).toString());
                //add count
                addCnt(variation, direction, cigarElementLength, queryQuality.charAt(cigarElementLength - 1) - 33, mappingQuality, numberOfMismatches);
                //increase coverage
                incCnt(refCoverage, start - 1, 1);
                start--;
                cigarElementLength--;
            }
            if (cigarElementLength > 0) { //If there remains a soft-clipped sequence at the beginning (not everything was matched)
                int sumOfReadQualities = 0; // $q Sum of read qualities (to get mean quality)
                int numberOfHighQualityBases = 0; // $qn Number of quality figures
                int numberOfLowQualityBases = 0; // $lowcnt Number of low-quality bases in the sequence

                // Loop over remaining soft-clipped sequence
                for (int si = cigarElementLength - 1; si >= 0; si--) {
                    // Stop if unknown base (N - any of ATGC) is found
                    if (querySequence.charAt(si) == 'N') {
                        break;
                    }
                    // $tq - base quality
                    int baseQuality = queryQuality.charAt(si) - 33;
                    if (baseQuality <= 12)
                        numberOfLowQualityBases++;
                    //Stop if a low-quality base is found
                    if (numberOfLowQualityBases > 1) {
                        break;
                    }
                    sumOfReadQualities += baseQuality;
                    numberOfHighQualityBases++;
                }
                sclip5HighQualityProcessing(querySequence, mappingQuality, queryQuality, numberOfMismatches,
                        direction, sumOfReadQualities, numberOfHighQualityBases, numberOfLowQualityBases);
            }
            cigarElementLength = cigar.getCigarElement(ci).getLength();
        } else if (ci == cigar.numCigarElements() - 1) { // 3' soft clipped
            // Ignore large soft clip due to chimeric reads in library construction
            if (!instance().conf.chimeric) {
                String saTagString = record.getStringAttribute(SAMTag.SA.name());
                if (cigarElementLength >= 20 && saTagString != null) {
                    if (isReadChimericWithSA(record, position, saTagString, direction, false)) {
                        readPositionIncludingSoftClipped += cigarElementLength;
                        offset = 0;
                        // Had to reset the start due to softclipping adjustment
                        start = position;

                        return;
                    }
                } else if (cigarElementLength >= Configuration.SEED_1) {
                    Map<String, List<Integer>> referenceSeedMap = reference.seed;
                    String sequence = getReverseComplementedSequence(record, -cigarElementLength, cigarElementLength);
                    String reverseComplementedSeed = substr(sequence, -Configuration.SEED_1, Configuration.SEED_1);

                    if (referenceSeedMap.containsKey(reverseComplementedSeed)) {
                        List<Integer> positions = referenceSeedMap.get(reverseComplementedSeed);
                        if (positions.size() == 1 && abs(start - positions.get(0)) <  2 * maxReadLength) {
                            readPositionIncludingSoftClipped += cigarElementLength;
                            offset = 0;
                            // Had to reset the start due to softclipping adjustment
                            start = position;
                            if (instance().conf.y) {
                                System.err.println(sequence  + " at 3' is a chimeric at "
                                        + start + " by SEED " + Configuration.SEED_1);
                            }
                            return;
                        }
                    }
                }
            }
            /*
            Conditions:
            1). read position is less than sequence length
            2). reference base is defined for start
            3). reference base at start matches read base at n
            4). read quality is more than 10
             */
            while (readPositionIncludingSoftClipped < querySequence.length()
                    && isHasAndEquals(querySequence.charAt(readPositionIncludingSoftClipped), ref, start)
                    && queryQuality.charAt(readPositionIncludingSoftClipped) - 33 > 10) {
                //Initialize entry in $hash if not present
                Variation variation = getVariation(nonInsertionVariants, start, ref.get(start).toString());
                //Add count
                addCnt(variation, direction, totalLengthIncludingSoftClipped - readPositionExcludingSoftClipped,
                        queryQuality.charAt(readPositionIncludingSoftClipped) - 33, mappingQuality, numberOfMismatches);
                //Add coverage
                incCnt(refCoverage, start, 1);
                readPositionIncludingSoftClipped++;
                start++;
                cigarElementLength--;
                readPositionExcludingSoftClipped++;
            }
            //If there remains a soft-clipped sequence at the end (not everything was matched)
            if (querySequence.length() - readPositionIncludingSoftClipped > 0) {
                int sumOfReadQualities = 0; // $q Sum of read qualities (to get mean quality)
                int numberOfHighQualityBases = 0; // $qn Number of quality figures
                int numberOfLowQualityBases = 0; // $lowcnt Number of low-quality bases in the sequence
                for (int si = 0; si < cigarElementLength; si++) { //Loop over remaining soft-clipped sequence
                    //Stop if unknown base (N - any of ATGC) is found
                    if (querySequence.charAt(readPositionIncludingSoftClipped + si) == 'N') {
                        break;
                    }
                    // $tq - Base quality
                    int baseQuality = queryQuality.charAt(readPositionIncludingSoftClipped + si) - 33;
                    if (baseQuality <= 12) {
                        numberOfLowQualityBases++;
                    }
                    //Stop if a low-quality base is found
                    if (numberOfLowQualityBases > 1) {
                        break;
                    }
                    sumOfReadQualities += baseQuality;
                    numberOfHighQualityBases++;
                }
                sclip3HighQualityProcessing(querySequence, mappingQuality, queryQuality, numberOfMismatches,
                        direction, sumOfReadQualities, numberOfHighQualityBases, numberOfLowQualityBases);
            }
        }
        //Move read position by m (length of segment in CIGAR)
        readPositionIncludingSoftClipped += cigarElementLength;
        offset = 0;
        start = position; // Had to reset the start due to softclipping adjustment
    }

    /**
     * N in CIGAR - skipped region from reference
     * Skip the region and add string start-end to %SPLICE
     */
    private void processNotMatched() {
        String key = (start - 1) + "-" + (start + cigarElementLength - 1);
        splice.add(key);
        int[] cnt = spliceCount.get(key);
        if (cnt == null) {
            cnt = new int[] { 0 };
            spliceCount.put(key, cnt);
        }
        cnt[0]++;

        start += cigarElementLength;
        offset = 0;
        return;
    }

    /**
     * In amplicon mode try to check distances and overlap fraction of amplicon.
     * @param record SAMRecord to check
     * @param isMateReferenceNameEqual true if reard and mate reference names identical
     * @return true if read must be proccesed further
     */
    private boolean parseCigarWithAmpCase(SAMRecord record, boolean isMateReferenceNameEqual) {
        String[] split = instance().ampliconBasedCalling.split(":");
        // Distance to amplicon (specified in -a option)
        int distanceToAmplicon;
        // Overlap fraction (in -a option)
        double overlapFraction;
        try {
            distanceToAmplicon = toInt(split[0]);
            overlapFraction = Double.parseDouble(split.length > 1 ? split[1] : "");
        } catch (NumberFormatException e) {
            distanceToAmplicon = 10;
            overlapFraction = 0.95;
        }
        // $rlen3 holds sum of lengths of matched and deleted segments (The total aligned length,
        // excluding soft-clipped bases and insertions)
        int readLengthIncludeMatchedAndDeleted = getAlignedLength(cigar);
        int segstart = record.getAlignmentStart();
        int segend = segstart + readLengthIncludeMatchedAndDeleted - 1;

        //If read starts with soft-clipped sequence
        if (cigar.getCigarElement(0).getOperator() == CigarOperator.S) {
            //Ignore reads that overlap with region of interest by fraction less than ovlp
            int ts1 = segstart > region.start ? segstart : region.start;
            int te1 = segend < region.end ? segend : region.end;
            if (Math.abs(ts1 - te1) / (double)(segend - segstart) > overlapFraction == false) {
                return true;
            }
        } else if (cigar.getCigarElement(cigar.numCigarElements() - 1).getOperator() == CigarOperator.S) {
            //If read ends with contains soft-clipped sequence
            //Ignore reads that overlap with region of interest by fraction less than ovlp
            int ts1 = segstart > region.start ? segstart : region.start;
            int te1 = segend < region.end ? segend : region.end;
            if (Math.abs(te1 - ts1) / (double)(segend - segstart) > overlapFraction == false) {
                return true;
            }

        } else { // No soft-clipping
            // If RNEXT is identical to RNAME and TLEN is defined
            if (isMateReferenceNameEqual && record.getInferredInsertSize() != 0) {
                if (record.getInferredInsertSize() > 0) {
                    segend = segstart + record.getInferredInsertSize() - 1;
                } else {
                    segstart = record.getMateAlignmentStart();
                    segend = record.getMateAlignmentStart() - record.getInferredInsertSize() - 1;
                }
            }
            // No segment overlapping test since samtools should take care of it
            int ts1 = segstart > region.start ? segstart : region.start;
            int te1 = segend < region.end ? segend : region.end;
            // Ignore reads that are more than distance from region of interest
            // and overlap is less than overlap fraction from configuration (-a)
            if ((abs(segstart - region.start) > distanceToAmplicon || abs(segend - region.end) > distanceToAmplicon)
                    || abs((ts1 - te1) / (double)(segend - segstart)) <= overlapFraction) {
                return true;
            }
        }
        return false;
    }

    /**
     * Decrement variant counters.
     * @param variation reference variant to decrement    $vref
     * @param direction true if read has negative strand flag    $dir
     * @param readPosition read position    $rp
     * @param baseQuality base quality   $q
     * @param mappingBaseQuality bases's mapping quality    $Q
     * @param numberOfMismatches number of mismatches   $nm
     */
    void subCnt(Variation variation,
                              boolean direction,
                              int readPosition,
                              double baseQuality,
                              int mappingBaseQuality,
                              int numberOfMismatches) {
        // rp = read position, q = quality
        variation.varsCount--;
        variation.decDir(direction);
        variation.meanPosition -= readPosition;
        variation.meanQuality -= baseQuality;
        variation.meanMappingQuality -= mappingBaseQuality;
        variation.numberOfMismatches -= numberOfMismatches;
        if (baseQuality >= instance().conf.goodq) {
            variation.highQualityReadsCount--;
        } else {
            variation.lowQualityReadsCount--;
        }
    }

    /**
     * Increment variant counters.
     * @param variation variant to increment    $vref
     * @param direction if read has negative strand flag    $dir
     * @param readPosition read position    $rp
     * @param baseQuality base quality   $q
     * @param mappingBaseQuality bases's mapping quality    $Q
     * @param numberOfMismatches number of mismatches   $nm
     */
     void addCnt(Variation variation,
                              boolean direction,
                              int readPosition,
                              double baseQuality,
                              int mappingBaseQuality,
                              int numberOfMismatches) {
        variation.varsCount++;
        variation.incDir(direction);
        variation.meanPosition += readPosition;
        variation.meanQuality += baseQuality;
        variation.meanMappingQuality += mappingBaseQuality;
        variation.numberOfMismatches += numberOfMismatches;
        if (baseQuality >= instance().conf.goodq) {
            variation.highQualityReadsCount++;
        } else {
            variation.lowQualityReadsCount++;
        }
    }

    /**
     * Increase count for variation
     * @param counters map to increase count for variation
     * @param index index in <code>counters</code> map
     * @param descriptionString string with variations
     */
    private void increment(Map<Integer, Map<String, Integer>> counters,
                           int index,
                           String descriptionString) {
        Map<String, Integer> map = counters.get(index);
        if (map == null) {
            map = new HashMap<>();
            counters.put(index, map);
        }
        incCnt(map, descriptionString, 1);
    }

    boolean isBEGIN_ATGC_AMP_ATGCs_END(String sequence) {
        if (sequence.length() > 2) {
            char firstChar = sequence.charAt(0);
            char secondChar = sequence.charAt(1);
            if (secondChar == '&' && isATGC(firstChar)) {
                for (int i = 2; i < sequence.length(); i++) {
                    if (!isATGC(sequence.charAt(i))) {
                        return false;
                    }
                }
                return true;
            }
        }
        return false;
    }

    boolean isATGC(char ch) {
        switch (ch) {
            case 'A':
            case 'T':
            case 'G':
            case 'C':
                return true;

            default:
                return false;
        }
    }

    /**
     * Find closest mismatches to combine with indels.
     * @param referencePosition position for reference  $refp
     * @param readPosition position for base sequence   $readp
     * @param cigarLength length of CIGAR   $mlen
     * @param querySequence base sequence   $rdseq
     * @param queryQuality base quality sequence    $qstr
     * @param reference reference bases by position     $ref
     * @param refCoverage reference coverage    $cov
     * @return Offset object of (offset, querySequence's substring, queryQuality substring, number of mismatches)
     */
    Offset findOffset(int referencePosition,
                      int readPosition,
                      int cigarLength,
                      String querySequence,
                      String queryQuality,
                      Map<Integer, Character> reference,
                      Map<Integer, Integer> refCoverage) {
        int offset = 0;
        String ss = "";
        String q = "";
        int tnm = 0;
        int vsn = 0;
        for (int vi = 0; vsn <= instance().conf.vext && vi < cigarLength; vi++) {
            if (querySequence.charAt(readPosition + vi) == 'N') {
                break;
            }
            if (queryQuality.charAt(readPosition + vi) - 33 < instance().conf.goodq) {
                break;
            }
            Character refCh = reference.get(referencePosition + vi);
            if (refCh != null) {
                char ch = querySequence.charAt(readPosition + vi);
                if (isNotEquals(ch, refCh)) {
                    offset = vi + 1;
                    tnm++;
                    vsn = 0;
                } else {
                    vsn++;
                }
            }
        }
        if (offset > 0) {
            ss = substr(querySequence, readPosition, offset);
            q = substr(queryQuality, readPosition, offset);
            for (int osi = 0; osi < offset; osi++) {
                incCnt(refCoverage, referencePosition + osi, 1);
            }
        }
        return new Offset(offset, ss, q, tnm);
    }

    /**
     * Attempts to clean up a CIGAR string so that edge cases are avoided in the rest of the code.
     * Specifically this method will remove leading or trailing hard-clips, and then convert
     * leading or trailing insertions into soft-clips.
     * @param rec SAMRecord
     */
    void cleanupCigar(final SAMRecord rec) {
        if (rec.getCigar() != null) {
            final List<CigarElement> elems = new ArrayList<>(rec.getCigar().getCigarElements());

            // First leading elements
            {
                final ListIterator<CigarElement> iterator = elems.listIterator();
                boolean noMatchesYet = true;
                while (iterator.hasNext() && noMatchesYet) {
                    final CigarElement elem = iterator.next();
                    if (elem.getOperator() == CigarOperator.INSERTION) {
                        final CigarElement replacement = new CigarElement(elem.getLength(), CigarOperator.SOFT_CLIP);
                        iterator.set(replacement);
                    }
                    else if (elem.getOperator() == CigarOperator.HARD_CLIP) {
                        iterator.remove();
                    }
                    else if (elem.getOperator().consumesReadBases() && elem.getOperator().consumesReferenceBases()) {
                        noMatchesYet = false;
                    }
                }
            }
            // Then trailing elements
            {
                final ListIterator<CigarElement> iterator = elems.listIterator(elems.size());
                boolean noMatchesYet = true;
                while (iterator.hasPrevious() && noMatchesYet) {
                    final CigarElement elem = iterator.previous();
                    if (elem.getOperator() == CigarOperator.INSERTION) {
                        final CigarElement replacement = new CigarElement(elem.getLength(), CigarOperator.SOFT_CLIP);
                        iterator.set(replacement);
                    }
                    else if (elem.getOperator() == CigarOperator.HARD_CLIP) {
                        iterator.remove();
                    }
                    else if (elem.getOperator().consumesReadBases() && elem.getOperator().consumesReferenceBases()) {
                        noMatchesYet = false;
                    }
                }
            }
            // And lastly replace the cigar
            rec.setCigar(new Cigar(elems));
        }
    }

    static CigarOperator getCigarOperator(Cigar cigar,
                                                 int ci) {
        CigarOperator operator =  cigar.getCigarElement(ci).getOperator();
        // Treat insertions at the edge as soft-clipping
        if ((ci == 0 || ci == cigar.numCigarElements() - 1) && operator == CigarOperator.I) {
            operator = CigarOperator.S;
        }
        return operator;
    }

    /**
     * Method calculates the read length including deleted and matched segments (excluding soft-clipped and insertions)
     * @param readCigar current read CIGAR
     * @return length of aligned segment
     */
    static int getAlignedLength(Cigar readCigar) {
        int length = 0;
        for (CigarElement element : readCigar.getCigarElements()) {
            if (element.getOperator() == CigarOperator.M || element.getOperator() == CigarOperator.D) {
                length += element.getLength();
            }
        }
        return length;
    }

    /**
     * Method calculates the read length including matched, soft-clipped and insertions segments
     * @param readCigar current read CIGAR
     * @return length of segment
     */
    static int getSoftClippedLength(Cigar readCigar) {
        int length = 0;
        for (CigarElement element : readCigar.getCigarElements()) {
            CigarOperator operator = element.getOperator();
            if (operator == CigarOperator.M || operator == CigarOperator.I || operator == CigarOperator.S) {
                length += element.getLength();
            }
        }
        return length;
    }

    /**
     * Method calculates the read length including matched and insertions segments
     * @param readCigar current read CIGAR
     * @return length of segment
     */
    static int getMatchInsertionLength(Cigar readCigar) {
        int length = 0;
        for (CigarElement element : readCigar.getCigarElements()) {
            if (element.getOperator() == CigarOperator.M || element.getOperator() == CigarOperator.I) {
                length += element.getLength();
            }
        }
        return length;
    }

    /**
     * Method calculates the read length including deletions and insertions segments
     * @param readCigar current read CIGAR
     * @return length of segment
     */
    static int getInsertionDeletionLength(Cigar readCigar) {
        int length = 0;
        for (CigarElement element : readCigar.getCigarElements()) {
            if (element.getOperator() == CigarOperator.I || element.getOperator() == CigarOperator.D) {
                length += element.getLength();
            }
        }
        return length;
    }

    public static String getMateReferenceName(SAMRecord record) {
        if (record.getMateReferenceName() == null) {
            return "*";
        }

        if (record.getReferenceName().equals(record.getMateReferenceName())) {
            return "=";
        }
        return record.getMateReferenceName();
    }
    /**
     * Skip the insertions and deletions that are right after or before introns
     (they indicate of aligner problem)
     * @param cigar cigar of the read
     * @param ci current index of element in cigar
     * @return true if cigar operator must be skipped, false if not
     */
    private boolean skipIndelNextToIntron(Cigar cigar, int ci) {
        if ((cigar.numCigarElements() > ci + 1 && cigar.getCigarElement(ci + 1).getOperator() == CigarOperator.N)
                || (ci > 0 && cigar.getCigarElement(ci - 1).getOperator() == CigarOperator.N)) {
            return true;
        }
        return false;
    }

    /**
     * Creates variation for matching part of Cigar (M). If variation already exists, increment it's counters
     */
    private void addVariationForMatchingPart(int mappingQuality, int nm, boolean dir,
                                                    int rlen1, int nmoff, String s,
                                                    boolean startWithDelition, double q,
                                                    int qbases, int qibases, int ddlen, int pos) {
        Variation hv;
        if (s.startsWith("+")) {
            increment(positionToInsertionCount, pos, s);
            hv = getVariation(insertionVariants, pos, s); //reference to variant structure
        } else {
            if(isBEGIN_ATGC_AMP_ATGCs_END(s)) {
                //if s is one base followed by '&' and one or more bases
                //add variant record for s to mnp
                increment(mnp, pos, s);
            }
            hv = getVariation(nonInsertionVariants, pos, s); //reference to variant structure
        }
        hv.incDir(dir);

        //increment count
        hv.varsCount++;

        //minimum of positions from start of read and end of read
        int tp = readPositionExcludingSoftClipped < rlen1 - readPositionExcludingSoftClipped
                ? readPositionExcludingSoftClipped + 1
                : rlen1 - readPositionExcludingSoftClipped;

        //average quality of bases in the variation
        q = q / (qbases + qibases);

        //pstd is a flag that is 1 if the variant is covered by at least 2 read segments with different positions
        if (!hv.pstd && hv.pp != 0 && tp != hv.pp) {
            hv.pstd = true;
        }

        //qstd is a flag that is 1 if the variant is covered by at least 2 segment reads with different qualities
        if (!hv.qstd && hv.pq != 0 && q != hv.pq) {
            hv.qstd = true;
        }
        hv.meanPosition += tp;
        hv.meanQuality += q;
        hv.meanMappingQuality += mappingQuality;
        hv.pp = tp;
        hv.pq = q;
        hv.numberOfMismatches += nm - nmoff;
        if (q >= instance().conf.goodq) {
            hv.highQualityReadsCount++;
        } else {
            hv.lowQualityReadsCount++;
        }
        int shift = (s.startsWith("+") && s.contains("&")) ? 1 : 0;

        //increase coverage for bases covered by the variation
        for (int qi = 1; qi <= qbases - shift; qi++) {
            incCnt(refCoverage, start - qi + 1, 1);
        }

        //If variation starts with a deletion ('-' character)
        if (startWithDelition) {
            //add variation to deletions map
            increment(positionToDeletionCount, pos, s);

            //increase coverage for next CIGAR segment
            for (int qi = 1; qi < ddlen; qi++) {
                incCnt(refCoverage, start + qi, 1);
            }
        }
    }

    /**
     * Creates variation for deletion in Cigar (D). If variation already exists, increment it's counters
     */
    private void addVariationForDeletion(int mappingQuality, int nm, boolean dir, int rlen1,
                                         StringBuilder descStringOfDeletedElement, StringBuilder qualityOfSegment, int nmoff) {
        //add variant structure for deletion at this position
        Variation hv = getVariation(nonInsertionVariants, start, descStringOfDeletedElement.toString()); //variation structure
        //add record for deletion in deletions map
        increment(positionToDeletionCount, start, descStringOfDeletedElement.toString());
        hv.incDir(dir);
        //increase count
        hv.varsCount++;

        //minimum of positions from start of read and end of read
        int tp = readPositionExcludingSoftClipped < rlen1 - readPositionExcludingSoftClipped
                ? readPositionExcludingSoftClipped + 1
                : rlen1 - readPositionExcludingSoftClipped;

        //average quality of bases
        double tmpq = 0;

        for (int i = 0; i < qualityOfSegment.length(); i++) {
            tmpq += qualityOfSegment.charAt(i) - 33;
        }
        tmpq = tmpq / qualityOfSegment.length();

        //pstd is a flag that is 1 if the variant is covered by at least 2 read segments with different positions
        if (!hv.pstd && hv.pp != 0 && tp != hv.pp) {
            hv.pstd = true;
        }

        //qstd is a flag that is 1 if the variant is covered by at least 2 segment reads with different qualities
        if (!hv.qstd && hv.pq != 0 && tmpq != hv.pq) {
            hv.qstd = true;
        }
        hv.meanPosition += tp;
        hv.meanQuality += tmpq;
        hv.meanMappingQuality += mappingQuality;
        hv.pp = tp;
        hv.pq = tmpq;
        hv.numberOfMismatches += nm - nmoff;
        if (tmpq >= instance().conf.goodq) {
            hv.highQualityReadsCount++;
        } else {
            hv.lowQualityReadsCount++;
        }

        //increase coverage count for reference bases missing from the read
        for (int i = 0; i < cigarElementLength; i++) {
            incCnt(refCoverage, start + i, 1);
        }
    }

    /**
     * Check if read is chimeric and contains SA tag
     * @param record SAMRecord to check
     * @param position start position of read
     * @param saTagString content of SA tag
     * @param dir direction of the strand (true is reverse direction)
     * @param is5Side true if it is 5' side
     * @return true if record is chimeric with SA tag
     */
    private boolean isReadChimericWithSA(SAMRecord record, int position, String saTagString,
                                         boolean dir, boolean is5Side){
        String[] saTagArray = saTagString.split(",");
        String saChromosome = saTagArray[0];
        int saPosition = Integer.valueOf(saTagArray[1]);
        String saDirectionString = saTagArray[2];
        String saCigar = saTagArray[3];
        boolean saDirectionIsForward = saDirectionString.equals("+");
        Matcher mm = SA_CIGAR_D_S_3clip.matcher(saCigar);

        if (is5Side) {
            mm = SA_CIGAR_D_S_5clip.matcher(saCigar);
        }

        boolean isChimericWithSA = ((dir && saDirectionIsForward) || (!dir && !saDirectionIsForward))
                && saChromosome.equals(record.getReferenceName())
                && (abs(saPosition - position) < 2 * maxReadLength)
                && mm.find();

        if (instance().conf.y && isChimericWithSA) {
            System.err.println(record.getReadName() + " " + record.getReferenceName()
                    + " " + position + " " + record.getMappingQuality()
                    + " " + cigar.toString()
                    + " is ignored as chimeric with SA: " +
                    saPosition + "," + saDirectionString + "," + saCigar);
        }

        return isChimericWithSA;
    }

    /**
     * Append sequence for deletion or insertion cases to create description string and quality string.
     */
    private void appendSegments(String querySequence, String queryQuality, int ci,
                                StringBuilder descStringOfElement, StringBuilder qualitySegment,
                                int mLen, int indelLen, int begin, boolean isInsertion) {

        //begin is n + m for insertion and n for deletion
        //append to s '#' and part of read sequence corresponding to next CIGAR segment (matched one)
        descStringOfElement.append("#").append(substr(querySequence, begin, mLen));
        //append quality string of next matched segment from read
        qualitySegment.append(substr(queryQuality, begin, mLen));

        //if an insertion is two segments ahead, append '^' + part of sequence corresponding
        // to next-next segment otherwise (deletion) append '^' + length of a next-next segment
        descStringOfElement.append('^').append(cigar.getCigarElement(ci + 2).getOperator() == CigarOperator.I
                ? substr(querySequence, begin + mLen, indelLen)
                : indelLen);
        //if an insertion is two segments ahead, append part of quality string sequence
        // corresponding to next-next segment otherwise (deletion)
        // append first quality score of next segment or return empty string
        if (isInsertion) {
            qualitySegment.append(cigar.getCigarElement(ci + 2).getOperator() == CigarOperator.I
                    ? substr(queryQuality, begin + mLen, indelLen)
                    : queryQuality.charAt(begin + mLen));
        } else {
            qualitySegment.append(cigar.getCigarElement(ci + 2).getOperator() == CigarOperator.I
                    ? substr(queryQuality, begin + mLen, indelLen)
                    : "");
        }
    }

    /**
     * Skip overlapping reads to avoid double counts for coverage (alignment or second in pair flag is used)
     * @param record SAMRecord to check
     * @param position start position of read
     * @param dir direction of the strand (true is reverse direction)
     * @param mateAlignmentStart start position of mate
     * @return true if reads overlap and we must skip them
     */
    private boolean skipOverlappingReads(SAMRecord record, int position, boolean dir, int mateAlignmentStart) {
        if (instance().conf.uniqueModeAlignmentEnabled && isPairedAndSameChromosome(record)
                && !dir && start >= mateAlignmentStart) {
            return true;
        }
        if (instance().conf.uniqueModeSecondInPairEnabled && record.getSecondOfPairFlag() && isPairedAndSameChromosome(record)
                && isReadsOverlap(record, position, mateAlignmentStart)) {
            return true;
        }
        return false;
    }

    /**
     * Check if reads are overlapping. Two cases are considered: if position after mate start and if before.
     * @param record SAMRecord to check
     * @param position start position of read
     * @param mateAlignmentStart start position of mate
     * @return true if reads overlap
     */
    private boolean isReadsOverlap(SAMRecord record, int position, int mateAlignmentStart){
        if (position >= mateAlignmentStart) {
            return start >= mateAlignmentStart
                    && start <= (mateAlignmentStart + record.getCigar().getReferenceLength() - 1);
        }
        else {
            return start >= mateAlignmentStart
                    && record.getMateAlignmentStart() <= record.getAlignmentEnd();
        }
    }

    /**
     * Check if read is paired and has the same chromosome name for read and mate
     * @param record SAMRecord to check
     * @return if read is paired and it's chromosome name is identical to mate reference name
     */
    private boolean isPairedAndSameChromosome(SAMRecord record) {
        return record.getReadPairedFlag() && getMateReferenceName(record).equals("=");
    }

    /**
     * Process soft clip on 3' if it is has high quality reads
     */
    private void sclip3HighQualityProcessing(String querySequence, int mappingQuality, String queryQuality,
                                            int numberOfMismatches, boolean dir, int sumOfReadQualities,
                                            int numberOfHighQualityBases, int numberOfLowQualityBases) {
        //If we have at least 1 high-quality soft-clipped base of region of interest
        if (numberOfHighQualityBases >= 1 && numberOfHighQualityBases > numberOfLowQualityBases
                && start >= region.start && start <= region.end) {
            //add record to $sclip3
            Sclip sclip = softClips3End.get(start);
            if (sclip == null) {
                sclip = new Sclip();
                softClips3End.put(start, sclip);
            }
            for (int si = 0; si < numberOfHighQualityBases; si++) {
                Character ch = querySequence.charAt(readPositionIncludingSoftClipped + si);
                int idx = si;
                TreeMap<Character, Integer> cnts = sclip.nt.get(idx);
                if (cnts == null) {
                    cnts = new TreeMap<>();
                    sclip.nt.put(idx, cnts);

                }
                incCnt(cnts, ch, 1);
                Variation variation = getVariationFromSeq(sclip, idx, ch);
                addCnt(variation, dir, numberOfHighQualityBases - si,
                        queryQuality.charAt(readPositionIncludingSoftClipped + si) - 33, mappingQuality, numberOfMismatches);
            }
            addCnt(sclip, dir, cigarElementLength, sumOfReadQualities / (double) numberOfHighQualityBases, mappingQuality, numberOfMismatches);
        }
    }

    /**
     * Process soft clip on 5' if it is has high quality reads
     */
    private void sclip5HighQualityProcessing(String querySequence, int mappingQuality, String queryQuality,
                                             int numberOfMismatches, boolean dir, int sumOfReadQualities,
                                             int numberOfHighQualityBases, int numberOfLowQualityBases) {
        //If we have at least 1 high-quality soft-clipped base of region of interest
        if (numberOfHighQualityBases >= 1 && numberOfHighQualityBases > numberOfLowQualityBases && start >= region.start && start <= region.end) {
            //add record to $sclip5
            Sclip sclip = softClips5End.get(start);
            if (sclip == null) {
                sclip = new Sclip();
                softClips5End.put(start, sclip);
            }
            for (int si = cigarElementLength - 1; cigarElementLength - si <= numberOfHighQualityBases; si--) {
                Character ch = querySequence.charAt(si);
                int idx = cigarElementLength - 1 - si;
                TreeMap<Character, Integer> cnts = sclip.nt.get(idx);
                if (cnts == null) {
                    cnts = new TreeMap<>();
                    sclip.nt.put(idx, cnts);
                }
                incCnt(cnts, ch, 1);
                Variation seqVariation = getVariationFromSeq(sclip, idx, ch);
                addCnt(seqVariation, dir, si - (cigarElementLength - numberOfHighQualityBases),
                        queryQuality.charAt(si) - 33, mappingQuality, numberOfMismatches);
            }
            addCnt(sclip, dir, cigarElementLength, sumOfReadQualities / (double) numberOfHighQualityBases,
                    mappingQuality, numberOfMismatches);
        }
    }

    /**
     * Prepare and fill SV structures for deletions, duplications, inversions.
     */
    private void prepareSVStructuresForAnalysis(final SAMRecord record,
                                                final String queryQuality,
                                                final int numberOfMismatches,
                                                final boolean readDirection,
                                                final boolean mateDirection,
                                                final int start,
                                                final int totalLengthIncludingSoftClipped) {
        final int mateStart = record.getMateAlignmentStart();
        final int mend = mateStart + totalLengthIncludingSoftClipped;
        int end = start;
        List<String> msegs = globalFind(ALIGNED_LENGTH_MND, cigar.toString());
        end += sum(msegs);

        int soft5 = 0;

        jregex.Matcher matcher = BEGIN_NUM_S_OR_BEGIN_NUM_H.matcher(cigar.toString());
        if (matcher.find() && matcher.group(1) != null) {
            int tt = toInt(matcher.group(1));
            if (tt != 0 && queryQuality.charAt(tt - 1) - 33 > instance().conf.goodq) {
                soft5 = start;
            }
        }
        int soft3 = 0;

        matcher = END_NUM_S_OR_NUM_H.matcher(cigar.toString());
        if (matcher.find() && matcher.group(1) != null) {
            int tt = toInt(matcher.group(1));
            if (tt != 0 && queryQuality.charAt(queryQuality.length() - tt) - 33 > instance().conf.goodq) {
                soft3 = end;
            }
        }
        //Variables for direction will be simplier for calculations. Mate direction is true if it is forward,
        //read direction is true if it is reverse.
        final int readDirNum = readDirection ? -1 : 1;
        final int mateDirNum = mateDirection ? 1 : -1;
        //Minimal distance from start to end of SV
        final int MIN_D = 75;

        if (getMateReferenceName(record).equals("=")) {
            int mlen = record.getInferredInsertSize();
            if (record.getStringAttribute(SAMTag.MC.name()) != null
                    && MC_Z_NUM_S_ANY_NUM_S.matcher(record.getStringAttribute(SAMTag.MC.name())).find()) {
                // Ignore those with mates mapped with softcliping at both ends
            } else if (record.getIntegerAttribute(SAMTag.MQ.name()) != null
                    && record.getIntegerAttribute(SAMTag.MQ.name()) < 15) {
                // Ignore those with mate mapping quality less than 15
            } else if (readDirNum * mateDirNum == -1 && (mlen * readDirNum) > 0
                    && queryQuality.length() > Configuration.MINMAPBASE) {
                // deletion candidate
                mlen = mateStart > start ? mend - start : end - mateStart;
                if(abs(mlen) > instance().conf.INSSIZE + instance().conf.INSSTDAMT * instance().conf.INSSTD) {
                    if (readDirNum == 1) {
                        if (svStructures.svfdel.size() == 0
                                || start - svStructures.svdelfend > Configuration.MINSVCDIST * maxReadLength) {
                            Sclip sclip = new Sclip();
                            sclip.varsCount = 0;
                            svStructures.svfdel.add(sclip);
                        }
                        addSV(getLastSVStructure(svStructures.svfdel), start, end, mateStart,
                                mend, readDirNum, totalLengthIncludingSoftClipped, mlen, soft3, maxReadLength/2.0,
                                queryQuality.charAt(Configuration.MINMAPBASE) - 33, record.getMappingQuality(), numberOfMismatches);
                        svStructures.svdelfend = end;
                    } else {
                        if (svStructures.svrdel.size() == 0
                                || start - svStructures.svdelrend > Configuration.MINSVCDIST * maxReadLength) {
                            Sclip sclip = new Sclip();
                            sclip.varsCount = 0;
                            svStructures.svrdel.add(sclip);
                        }
                        addSV(getLastSVStructure(svStructures.svrdel), start, end, mateStart,
                                mend, readDirNum, totalLengthIncludingSoftClipped, mlen, soft5, maxReadLength/2.0,
                                queryQuality.charAt(Configuration.MINMAPBASE) - 33, record.getMappingQuality(), numberOfMismatches);
                        svStructures.svdelrend = end;
                    }

                    if (!svStructures.svfdel.isEmpty()
                            && abs(start - svStructures.svdelfend) <= Configuration.MINSVCDIST * maxReadLength ){
                        adddisccnt(getLastSVStructure(svStructures.svfdel));
                    }
                    if (!svStructures.svrdel.isEmpty()
                            && abs(start - svStructures.svdelrend) <= Configuration.MINSVCDIST * maxReadLength ){
                        adddisccnt(getLastSVStructure(svStructures.svrdel));
                    }
                    if (!svStructures.svfdup.isEmpty() && abs(start - svStructures.svdupfend) <= MIN_D){
                        adddisccnt(getLastSVStructure(svStructures.svfdup));
                    }
                    if (!svStructures.svrdup.isEmpty() && abs(start - svStructures.svduprend) <= MIN_D){
                        adddisccnt(getLastSVStructure(svStructures.svrdup));
                    }
                    if (!svStructures.svfinv5.isEmpty() && abs(start - svStructures.svinvfend5) <= MIN_D){
                        adddisccnt(getLastSVStructure(svStructures.svfinv5));
                    }
                    if (!svStructures.svrinv5.isEmpty() && abs(start - svStructures.svinvrend5) <= MIN_D){
                        adddisccnt(getLastSVStructure(svStructures.svrinv5));
                    }
                    if (!svStructures.svfinv3.isEmpty() && abs(start - svStructures.svinvfend3) <= MIN_D){
                        adddisccnt(getLastSVStructure(svStructures.svfinv3));
                    }
                    if (!svStructures.svrinv3.isEmpty() && abs(start - svStructures.svinvrend3) <= MIN_D){
                        adddisccnt(getLastSVStructure(svStructures.svrinv3));
                    }
                }
            } else if (readDirNum * mateDirNum == -1 && readDirNum * mlen < 0
                    && queryQuality.length() > Configuration.MINMAPBASE) {
                //duplication
                if (readDirNum == 1) {
                    if (svStructures.svfdup.size() == 0
                            || start - svStructures.svdupfend > Configuration.MINSVCDIST * maxReadLength) {
                        Sclip sclip = new Sclip();
                        sclip.varsCount = 0;
                        svStructures.svfdup.add(sclip);
                    }
                    addSV(getLastSVStructure(svStructures.svfdup), start, end, mateStart,
                            mend, readDirNum, totalLengthIncludingSoftClipped, mlen, soft3, maxReadLength/2.0,
                            queryQuality.charAt(Configuration.MINMAPBASE) - 33, record.getMappingQuality(), numberOfMismatches);
                    svStructures.svdupfend = end;
                } else {
                    if (svStructures.svrdup.size() == 0
                            || start - svStructures.svduprend > Configuration.MINSVCDIST * maxReadLength) {
                        Sclip sclip = new Sclip();
                        sclip.varsCount = 0;
                        svStructures.svrdup.add(sclip);
                    }
                    addSV(getLastSVStructure(svStructures.svrdup), start, end, mateStart,
                            mend, readDirNum, totalLengthIncludingSoftClipped, mlen, soft5, maxReadLength/2.0,
                            queryQuality.charAt(Configuration.MINMAPBASE) - 33, record.getMappingQuality(), numberOfMismatches);
                    svStructures.svduprend = end;
                }
                if (!svStructures.svfdup.isEmpty()
                        && abs(start - svStructures.svdupfend) <= Configuration.MINSVCDIST * maxReadLength ) {
                    getLastSVStructure(svStructures.svfdup).disc++;
                }
                if (!svStructures.svrdup.isEmpty()
                        && abs(start - svStructures.svduprend) <= Configuration.MINSVCDIST * maxReadLength ) {
                    getLastSVStructure(svStructures.svrdup).disc++;
                }
                if (!svStructures.svfdel.isEmpty() && abs(start - svStructures.svdelfend) <= MIN_D) {
                    adddisccnt(getLastSVStructure(svStructures.svfdel));
                }
                if (!svStructures.svrdel.isEmpty() && abs(start - svStructures.svdelrend) <= MIN_D) {
                    adddisccnt(getLastSVStructure(svStructures.svrdel));
                }
                if (!svStructures.svfinv5.isEmpty() && abs(start - svStructures.svinvfend5) <= MIN_D) {
                    adddisccnt(getLastSVStructure(svStructures.svfinv5));
                }
                if (!svStructures.svrinv5.isEmpty() && abs(start - svStructures.svinvrend5) <= MIN_D) {
                    adddisccnt(getLastSVStructure(svStructures.svrinv5));
                }
                if (!svStructures.svfinv3.isEmpty() && abs(start - svStructures.svinvfend3) <= MIN_D) {
                    adddisccnt(getLastSVStructure(svStructures.svfinv3));
                }
                if (!svStructures.svrinv3.isEmpty() && abs(start - svStructures.svinvrend3) <= MIN_D) {
                    adddisccnt(getLastSVStructure(svStructures.svrinv3));
                }
            } else if (readDirNum * mateDirNum == 1 && queryQuality.length() > Configuration.MINMAPBASE) { // Inversion
                if (readDirNum == 1 && mlen != 0 ) {
                    if (mlen < -3 * maxReadLength) {
                        if (svStructures.svfinv3.size() == 0
                                || start - svStructures.svinvfend3 > Configuration.MINSVCDIST * maxReadLength) {
                            Sclip sclip = new Sclip();
                            sclip.varsCount = 0;
                            svStructures.svfinv3.add(sclip);
                        }
                        addSV(getLastSVStructure(svStructures.svfinv3), start, end, mateStart,
                                mend, readDirNum, totalLengthIncludingSoftClipped, mlen, soft3, maxReadLength/2.0,
                                queryQuality.charAt(Configuration.MINMAPBASE) - 33, record.getMappingQuality(), numberOfMismatches);
                        svStructures.svinvfend3 = end;
                        getLastSVStructure(svStructures.svfinv3).disc++;
                    } else if (mlen > 3 * maxReadLength) {
                        if (svStructures.svfinv5.size() == 0
                                || start - svStructures.svinvfend5 > Configuration.MINSVCDIST * maxReadLength) {
                            Sclip sclip = new Sclip();
                            sclip.varsCount = 0;
                            svStructures.svfinv5.add(sclip);
                        }
                        addSV(getLastSVStructure(svStructures.svfinv5), start, end, mateStart,
                                mend, readDirNum, totalLengthIncludingSoftClipped, mlen, soft3, maxReadLength/2.0,
                                queryQuality.charAt(Configuration.MINMAPBASE) - 33, record.getMappingQuality(), numberOfMismatches);
                        svStructures.svinvfend5 = end;
                        getLastSVStructure(svStructures.svfinv5).disc++;
                    }
                } else if (mlen != 0) {
                    if (mlen < -3 * maxReadLength) {
                        if (svStructures.svrinv3.size() == 0
                                || start - svStructures.svinvrend3 > Configuration.MINSVCDIST * maxReadLength) {
                            Sclip sclip = new Sclip();
                            sclip.varsCount = 0;
                            svStructures.svrinv3.add(sclip);
                        }
                        addSV(getLastSVStructure(svStructures.svrinv3), start, end, mateStart,
                                mend, readDirNum, totalLengthIncludingSoftClipped, mlen, soft5, maxReadLength/2.0,
                                queryQuality.charAt(Configuration.MINMAPBASE) - 33, record.getMappingQuality(), numberOfMismatches);
                        svStructures.svinvrend3 = end;
                        getLastSVStructure(svStructures.svrinv3).disc++;

                    } else if (mlen > 3 * maxReadLength) {
                        if (svStructures.svrinv5.size() == 0
                                || start - svStructures.svinvrend5 > Configuration.MINSVCDIST * maxReadLength) {
                            Sclip sclip = new Sclip();
                            sclip.varsCount = 0;
                            svStructures.svrinv5.add(sclip);
                        }
                        addSV(getLastSVStructure(svStructures.svrinv5), start, end, mateStart,
                                mend, readDirNum, totalLengthIncludingSoftClipped, mlen, soft5, maxReadLength/2.0,
                                queryQuality.charAt(Configuration.MINMAPBASE) - 33, record.getMappingQuality(), numberOfMismatches);
                        svStructures.svinvrend5 = end;
                        getLastSVStructure(svStructures.svrinv5).disc++;
                    }
                }
                if (mlen != 0) {
                    if (!svStructures.svfdel.isEmpty() && (start - svStructures.svdelfend) <= MIN_D) {
                        adddisccnt(getLastSVStructure(svStructures.svfdel));
                    }
                    if (!svStructures.svrdel.isEmpty() && (start - svStructures.svdelrend) <= MIN_D) {
                        adddisccnt(getLastSVStructure(svStructures.svrdel));
                    }
                    if (!svStructures.svfdup.isEmpty() && (start - svStructures.svdupfend) <= MIN_D) {
                        adddisccnt(getLastSVStructure(svStructures.svfdup));
                    }
                    if (!svStructures.svrdup.isEmpty() && (start - svStructures.svduprend) <= MIN_D) {
                        adddisccnt(getLastSVStructure(svStructures.svrdup));
                    }
                }
            }
        } else if (queryQuality.length() > Configuration.MINMAPBASE){ // Inter-chr translocation
            // to be implemented
            final String mchr = getMateReferenceName(record);
            if (record.getStringAttribute(SAMTag.MC.name()) != null
                    && MC_Z_NUM_S_ANY_NUM_S.matcher(record.getStringAttribute(SAMTag.MC.name())).find()) {
                // Ignore those with mates mapped with softcliping at both ends
            } else if (record.getIntegerAttribute(SAMTag.MQ.name()) != null
                    && record.getIntegerAttribute(SAMTag.MQ.name()) < 15) {
                // Ignore those with mate mapping quality less than 15
            } else if (readDirNum == 1) {
                if (svStructures.svffus.get(mchr) == null
                        || start - svStructures.svfusfend.get(mchr) > Configuration.MINSVCDIST * maxReadLength) {
                    Sclip sclip = new Sclip();
                    sclip.varsCount = 0;
                    List<Sclip> sclips = svStructures.svffus.getOrDefault(mchr, new ArrayList<>());
                    sclips.add(sclip);
                    svStructures.svffus.put(mchr, sclips);
                }
                int svn = svStructures.svffus.get(mchr).size() - 1;
                addSV(svStructures.svffus.get(mchr).get(svn), start, end, mateStart,
                        mend, readDirNum, totalLengthIncludingSoftClipped, 0, soft3, maxReadLength/2.0,
                        queryQuality.charAt(Configuration.MINMAPBASE) - 33, record.getMappingQuality(), numberOfMismatches);
                svStructures.svfusfend.put(mchr, end);
                svStructures.svffus.get(mchr).get(svn).disc++;
            } else {
                if (svStructures.svrfus.get(mchr) == null
                        || start - svStructures.svfusrend.get(mchr) > Configuration.MINSVCDIST * maxReadLength) {
                    Sclip sclip = new Sclip();
                    sclip.varsCount = 0;
                    List<Sclip> sclips = svStructures.svrfus.getOrDefault(mchr, new ArrayList<>());
                    sclips.add(sclip);
                    svStructures.svrfus.put(mchr, sclips);
                }
                int svn = svStructures.svrfus.get(mchr).size() - 1;
                addSV(svStructures.svrfus.get(mchr).get(svn), start, end, mateStart,
                        mend, readDirNum, totalLengthIncludingSoftClipped, 0, soft5, maxReadLength/2.0,
                        queryQuality.charAt(Configuration.MINMAPBASE) - 33, record.getMappingQuality(), numberOfMismatches);
                svStructures.svfusrend.put(mchr, end);
                svStructures.svrfus.get(mchr).get(svn).disc++;
            }
            if (!svStructures.svfdel.isEmpty() && start - svStructures.svdelfend <= Configuration.MINSVPOS) {
                adddisccnt(getLastSVStructure(svStructures.svfdel));
            }
            if (!svStructures.svrdel.isEmpty() && start - svStructures.svdelrend <= Configuration.MINSVPOS) {
                adddisccnt(getLastSVStructure(svStructures.svrdel));
            }
            if (!svStructures.svfdup.isEmpty() && start - svStructures.svdupfend <= Configuration.MINSVPOS) {
                adddisccnt(getLastSVStructure(svStructures.svfdup));
            }
            if (!svStructures.svrdup.isEmpty() && start - svStructures.svduprend <= Configuration.MINSVPOS) {
                adddisccnt(getLastSVStructure(svStructures.svrdup));
            }
            if (!svStructures.svfinv5.isEmpty() && start - svStructures.svinvfend5 <= Configuration.MINSVPOS) {
                adddisccnt(getLastSVStructure(svStructures.svfinv5));
            }
            if (!svStructures.svrinv5.isEmpty() && start - svStructures.svinvrend5 <= Configuration.MINSVPOS) {
                adddisccnt(getLastSVStructure(svStructures.svrinv5));
            }
            if (!svStructures.svfinv3.isEmpty() && start - svStructures.svinvfend3 <= Configuration.MINSVPOS) {
                adddisccnt(getLastSVStructure(svStructures.svfinv3));
            }
            if (!svStructures.svrinv3.isEmpty() && start - svStructures.svinvrend3 <= Configuration.MINSVPOS) {
                adddisccnt(getLastSVStructure(svStructures.svrinv3));
            }
        }
    }
    /**
     * Add discordant count
     * @param svref
     */
    private void adddisccnt(Sclip svref) {
        svref.disc++;
    }

    private Sclip getLastSVStructure(List<Sclip> svStructure) {
        return svStructure.get(svStructure.size() - 1);
    }

    /**
     * Add possible structural variant to current structural variant structure and update
     * @param sdref current structural variant structure
     */
    private  void addSV(Sclip sdref,
                              int start_s,
                              int end_e,
                              int mateStart_ms,
                              int mateEnd_me,
                              int dir,
                              int rlen,
                              int mlen,
                              int softp,
                              double pmean_rp,
                              double qmean,
                              double Qmean,
                              double nm) {
        sdref.varsCount++;
        sdref.incDir(dir == 1 ? false : true);

        if (qmean >= instance().conf.goodq) {
            sdref.highQualityReadsCount++;
        } else {
            sdref.lowQualityReadsCount++;
        }

        if (sdref.start == 0 || sdref.start >= start_s) {
            sdref.start = start_s;
        }
        if (sdref.end == 0 || sdref.end <= end_e) {
            sdref.end = end_e;
        }
        sdref.mates.add(new Mate(mateStart_ms, mateEnd_me, mlen,
                start_s, end_e, pmean_rp, qmean, Qmean, nm));

        if (sdref.mstart == 0 || sdref.mstart >= mateStart_ms) {
            sdref.mstart = mateStart_ms;
        }
        if (sdref.mend == 0 || sdref.mend <= mateEnd_me ) {
            sdref.mend = mateStart_ms + rlen;
        }
        if (softp != 0) {
            if (dir == 1) {
                if (abs(softp - sdref.end) < 10 ) {
                    int softCount = sdref.soft.getOrDefault(softp, 0);
                    softCount++;
                    sdref.soft.put(softp, softCount);
                }
            } else {
                if (abs(softp - sdref.start) < 10 ) {
                    int softCount = sdref.soft.getOrDefault(softp, 0);
                    softCount++;
                    sdref.soft.put(softp, softCount);
                }
            }
        }
    }
    /**
     * Get quality string in a safe way, avoid exceptions due to quality scores beyond MAX_PHRED_SCORE
     * @param record current SAMRecord from the parsed BAM
     * @return base quality string contains only values less then MAX_PHRED_SCORE
     */
    String getBaseQualityString(SAMRecord record) {
        try {
            return record.getBaseQualityString();
        } catch(IllegalArgumentException iae) {
            // For consolidated reads the quality score can be outside the allowed ranges for proper character representation
            // In this case we saturate at the highest value
            System.err.println("WARNING: Cannot get encode qualities for SAM entry at " + record.getReferenceName() + ":" + record.getAlignmentStart() + ", record name: '" + record.getReadName()+"'. Qualities capped at " + SAMUtils.MAX_PHRED_SCORE);
            byte[] qs=record.getBaseQualities();
            char[] qsChar = new char[qs.length];
            for(int i=0 ; i < qs.length; i++) {
                byte q = qs[i];
                if((q < 0) || (q > SAMUtils.MAX_PHRED_SCORE)) qs[i] = SAMUtils.MAX_PHRED_SCORE;
                qsChar[i] = (char)(q + 33);
            }
            return new String(qsChar);
        }
    }

    class Offset {
        int offset; //$offset
        String sequence; //$ss
        String qualitySequence; //$q
        Integer offsetNumberOfMismatches; //$tnm

        public Offset(int offset, String sequence, String qualitySequence, Integer offsetNumberOfMismatches) {
            this.offset = offset;
            this.sequence = sequence;
            this.qualitySequence = qualitySequence;
            this.offsetNumberOfMismatches = offsetNumberOfMismatches;
        }
    }
}
