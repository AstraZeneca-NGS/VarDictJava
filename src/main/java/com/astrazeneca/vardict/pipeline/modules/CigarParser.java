package com.astrazeneca.vardict.pipeline.modules;

import com.astrazeneca.utils.Tuple;
import com.astrazeneca.vardict.Configuration;
import com.astrazeneca.vardict.pipeline.data.Flags;
import com.astrazeneca.vardict.Region;
import com.astrazeneca.vardict.pipeline.data.Scope;
import com.astrazeneca.vardict.pipeline.data.VariationData;
import com.astrazeneca.vardict.variations.SoftClip;
import com.astrazeneca.vardict.variations.Variation;
import com.astrazeneca.vardict.variations.VariationUtils;
import htsjdk.samtools.*;

import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Set;

import static com.astrazeneca.GlobalReadOnlyScope.instance;
import static com.astrazeneca.utils.Utils.*;
import static com.astrazeneca.vardict.pipeline.modules.CigarUtils.getAlignedLenght;
import static com.astrazeneca.vardict.pipeline.modules.RecordPreprocessor.getChrName;
import static com.astrazeneca.vardict.variations.VariationUtils.isHasAndEquals;
import static com.astrazeneca.vardict.variations.VariationUtils.isHasAndNotEquals;
import static com.astrazeneca.vardict.variations.VariationUtils.isNotEquals;
import static java.lang.Math.abs;

/**
 * Class for parsing CIGAR.
 */
public class CigarParser implements Module<RecordPreprocessor, VariationData> {
    // Utils maps
    private final Map<Integer, Map<String, Variation>> nonInsertionVariants;
    private final Map<Integer, Map<String, Variation>> insertionVariants;
    private final Map<Integer, Integer> refCoverage;
    private final Map<Integer, SoftClip> softClips3End; // soft clipped at 3'
    private final Map<Integer, SoftClip> softClips5End; // soft clipped at 5'
    private final Map<Integer, Map<String, Integer>> positionToInsertionCount;
    private final Map<Integer, Map<String, Integer>> mnp; // Keep track of MNPs
    private final Map<Integer, Map<String, Integer>> positionToDeletionCount;
    private final Map<String, int[]> spliceCount;

    private Region region;
    private Set<String> splice;
    private Map<Integer, Character> reference;
    private int maxReadLength;

    private Cigar cigar;
    private int start;
    private int offset;
    private int cigarElementLength;
    private int readPositionIncludingSoftClipped; // keep track the read position, including softclipped
    private int readPositionExcludingSoftClipped; // keep track the position in the alignment, excluding softclipped


    public CigarParser() {
        this.nonInsertionVariants = new HashMap<>();
        this.insertionVariants = new HashMap<>();
        this.refCoverage = new HashMap<>();
        this.softClips3End = new HashMap<>();
        this.softClips5End = new HashMap<>();
        this.positionToInsertionCount = new HashMap<>();
        this.mnp = new HashMap<>();
        this.positionToDeletionCount = new HashMap<>();
        this.spliceCount = new HashMap<>();
    }

    @Override
    public Scope<VariationData> process(Scope<RecordPreprocessor> scope) {
        SAMRecord record;
        RecordPreprocessor processor = scope.data;

        initFromScope(scope);

        while ((record = processor.next()) != null) {
            parseCigar(
                    getChrName(scope.region),
                    record,
                    record.getReadString(),
                    record.getMappingQuality(),
                    new Flags(record.getFlags()));
        }

        processor.close();

        // TODO Unexecutable option in prev. version
        /*if (instance().conf.outputSplicing) {
            for (Map.Entry<String, int[]> entry : spliceCount.entrySet()) {
                System.out.printf("%s\t%s\t%s\t%s\n", instance().sample, region.chr, entry.getKey(), entry.getValue()[0]);
            }
            return null;
        }*/

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

        return new Scope<>(
                scope.bam,
                scope.region,
                scope.regionRef,
                maxReadLength,
                splice,
                scope.out,
                variationData);
    }

    private void initFromScope(Scope<RecordPreprocessor> scope) {
        this.region = scope.region;
        this.splice = scope.splice;
        this.reference = scope.regionRef;
        this.maxReadLength = scope.maxReadLength;
    }

    public Map<Integer, Map<String, Variation>> getNonInsertionVariants() {
        return nonInsertionVariants;
    }

    public Map<Integer, Map<String, Variation>> getInsertionVariants() {
        return insertionVariants;
    }

    public Map<Integer, Integer> getRefCoverage() {
        return refCoverage;
    }

    public Map<Integer, SoftClip> getSoftClips3End() {
        return softClips3End;
    }

    public Map<Integer, SoftClip> getSoftClips5End() {
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

    public int getMaxReadLength() {
        return maxReadLength;
    }

    /**
     * Parse {@link SAMRecord}'s CIGAR.
     * @param chrName chromosome's name
     * @param record {@link SAMRecord}
     * @param querySequence base sequence
     * @param mappingQuality record's mapping quality
     * @param flag {@link SAMRecord}'s flag
     */
    // perl version: 685
    void parseCigar(String chrName,
                    SAMRecord record,
                    String querySequence,
                    int mappingQuality,
                    Flags flag) {
        cigar = record.getCigar();

        final int insertionDeletionLength = CigarUtils.getInsertionDeletionLength(cigar);
        final Integer editDistance = record.getIntegerAttribute(SAMTag.NM.name());
        if (processNMTag(record, insertionDeletionLength, editDistance)) {
            return;
        }

        final String queryQuality = record.getBaseQualityString();
        final boolean isMateReferenceNameEqual = record.getReferenceName().equals(record.getMateReferenceName());

        //number of mismatches
        final int numberOfMismatches = editDistance != null ? editDistance - insertionDeletionLength : 0;
        boolean direction = flag.isReverseStrand();

        //perl version: 649
        if (instance().ampliconBasedCalling != null) {
            if (parseCigarWithAmpCase(record, isMateReferenceNameEqual))  {
                return;
            }
        }
        //perl version: 672
        if (flag.isUnmappedMate()) {
            // to be implemented
        } else {
            if (isMateReferenceNameEqual) {
                //if 'SA' tag (supplementary alignment) is present
                if (record.getStringAttribute(SAMTag.SA.name()) != null) {
                    if (flag.isSupplementaryAlignment()) { // the supplementary alignment
                        return; // Ignore the supplementary for now so that it won't skew the coverage
                    }
                }
            }
        }
        final int position;
        readPositionIncludingSoftClipped = 0;
        readPositionExcludingSoftClipped = 0;

        if (instance().conf.performLocalRealignment) {
            // Modify the CIGAR for potential mis-alignment for indels at the end of reads to softclipping and let VarDict's
            // algorithm to figure out indels
            Tuple.Tuple2<Integer, String> modifiedCigar = CigarUtils.modifyCigar(
                    insertionDeletionLength,
                    reference,
                    record.getAlignmentStart(),
                    record.getCigarString(),
                    querySequence,
                    queryQuality
            );

            position = modifiedCigar._1;
            cigar = TextCigarCodec.decode(modifiedCigar._2);
        } else {
            position = record.getAlignmentStart();
            cigar = record.getCigar();
        }

        //adjusted start position
        start = position;
        offset = 0;
        // Only match and insertion counts toward read length
        // For total length, including soft-clipped bases
        int readLengthForMatchedBases = CigarUtils.getMatchInsertionLength(cigar); // The read length for matched bases
        if (instance().conf.minmatch != 0 && readLengthForMatchedBases < instance().conf.minmatch) {
            return;
        }
        int totalLengthIncludingSoftClipped = CigarUtils.getSoftClippedLength(cigar); // The total length, including soft-clipped bases
        if (totalLengthIncludingSoftClipped > maxReadLength) { // Determine the read length
            maxReadLength = totalLengthIncludingSoftClipped;
        }

        // perl version: 890
        //Loop over CIGAR records
        for (int cigarElementIndex = 0; cigarElementIndex < cigar.numCigarElements(); cigarElementIndex++) {

            //length of segment in CIGAR
            cigarElementLength = cigar.getCigarElement(cigarElementIndex).getLength();

            //letter from CIGAR
            final CigarOperator operator = CigarUtils.getCigarOperator(cigar, cigarElementIndex);

            switch (operator) {
                case N: //N in CIGAR - skipped region from reference
                    //Skip the region and add string start-end to %SPLICE
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
                    continue;
                case S:
                    //First record in CIGAR
                    if (cigarElementIndex == 0) { // 5' soft clipped
                        processFirstCigarElementWithSoftClipping(
                                chrName,
                                querySequence,
                                mappingQuality,
                                queryQuality,
                                numberOfMismatches,
                                direction,
                                cigarElementIndex);
                    } else if (cigarElementIndex == cigar.numCigarElements() - 1) { // 3' soft clipped
                        processLastCigarElementWithSoftClipping(
                                querySequence,
                                mappingQuality,
                                queryQuality,
                                numberOfMismatches,
                                direction,
                                totalLengthIncludingSoftClipped);
                    }
                    //move read position by cigarElementLength (length of segment in CIGAR)
                    readPositionIncludingSoftClipped += cigarElementLength;
                    offset = 0;
                    start = position; // had to reset the start due to softclipping adjustment
                    continue;
                case H: //Hard clipping - skip
                    offset = 0;
                    continue;
                case I: { //Insertion
                    offset = 0;
                    cigarElementIndex = processCigarInsertions(
                            querySequence,
                            mappingQuality,
                            queryQuality,
                            numberOfMismatches,
                            direction,
                            readLengthForMatchedBases,
                            cigarElementIndex);
                    continue;
                }
                case D: { //deletion
                    offset = 0;
                    cigarElementIndex = processCigarDeletions(
                            querySequence,
                            mappingQuality,
                            queryQuality,
                            numberOfMismatches,
                            direction,
                            readLengthForMatchedBases,
                            cigarElementIndex);
                    continue;
                }
                default:
                    break;
            }
            //End branch on CIGAR segment type

            // Now dealing with matching part
            cigarElementIndex = processCigarMatches(
                    querySequence,
                    mappingQuality,
                    queryQuality,
                    numberOfMismatches,
                    direction,
                    readLengthForMatchedBases,
                    totalLengthIncludingSoftClipped,
                    cigarElementIndex,
                    operator);
            if (start > region.end) { //end if reference position is outside region of interest
                break;
            }
        }
    }

    /**
     * Dealing with matching part of CIGAR.
     * @param querySequence base sequence
     * @param mappingQuality record's mapping quality
     * @param queryQuality base quality sequence
     * @param numberOfMismatches number of mismatches
     * @param direction is read has negative strand flag
     * @param readLengthForMatchedBases the read length for matched bases
     * @param totalLengthIncludingSoftClipped the total length, including soft-clipped bases
     * @param cigarElementIndex index number of CIGAR element
     * @param operator letter from CIGAR
     * @return index number of CIGAR element; incremented if CIGAR has insertion two segments ahead
     */
    // perl version: 1188
    private int processCigarMatches(String querySequence,
                                    int mappingQuality,
                                    String queryQuality,
                                    int numberOfMismatches,
                                    boolean direction,
                                    int readLengthForMatchedBases,
                                    int totalLengthIncludingSoftClipped,
                                    int cigarElementIndex,
                                    CigarOperator operator) {
        int nmoff = 0;
        int moffset = 0;
        //Loop over bases of CIGAR segment
        for (int i = offset; i < cigarElementLength; i++) {
            //flag to trim reads at opt_T bases from start or end (depending on direction)
            boolean trim = false;
            if (instance().conf.trimBasesAfter != 0) {
                if (direction == false) {
                    if (readPositionIncludingSoftClipped > instance().conf.trimBasesAfter) {
                        trim = true;
                    }
                } else {
                    if (totalLengthIncludingSoftClipped - readPositionIncludingSoftClipped > instance().conf.trimBasesAfter) {
                        trim = true;
                    }
                }
            }

            //variation string. Initialize to first base of the read sequence
            final char ch1 = querySequence.charAt(readPositionIncludingSoftClipped);
            String s = String.valueOf(ch1);
            boolean startWithDeletion = false;
            //skip if base is unknown
            if (ch1 == 'N') {
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
            4). base at descriptionString in read is not equal to base in reference sequence
             */
            while ((start + 1) >= region.start
                    && (start + 1) <= region.end && (i + 1) < cigarElementLength
                    && q >= instance().conf.goodq
                    && isHasAndNotEquals(reference, start, querySequence, readPositionIncludingSoftClipped)
                    && isNotEquals('N', reference.get(start))) {

                //Break if base is unknown in the read
                char nuc = querySequence.charAt(readPositionIncludingSoftClipped + 1);
                if (nuc == 'N') {
                    break;
                }
                if (isHasAndEquals('N', reference, start + 1)) {
                    break;
                }

                //Condition: base at descriptionString + 1 does not match reference base at start + 1
                if (isNotEquals(reference.get(start + 1), nuc)) {

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
                } else { //if bases match, exit loop
                    break;
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
            5). either multi-nucleotide mismatch is found or read base at $descriptionString does not match reference base
            6). read base at $descriptionString has good quality
             */
            if (instance().conf.performLocalRealignment && cigarElementLength - i <= instance().conf.vext
                    && cigar.numCigarElements() > cigarElementIndex + 1
                    && cigar.getCigarElement(cigarElementIndex + 1).getOperator() == CigarOperator.D
                    && reference.containsKey(start)
                    && (ss.length() > 0 || isNotEquals(querySequence.charAt(readPositionIncludingSoftClipped), reference.get(start)))
                    && queryQuality.charAt(readPositionIncludingSoftClipped) - 33 > instance().conf.goodq) {

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
                s = "-" + cigar.getCigarElement(cigarElementIndex + 1).getLength() + "&" + s;
                startWithDeletion = true;
                ddlen = cigar.getCigarElement(cigarElementIndex + 1).getLength();
                cigarElementIndex += 1;

                //If CIGAR has insertion two segments ahead
                if (cigar.numCigarElements() > cigarElementIndex + 1
                        && cigar.getCigarElement(cigarElementIndex + 1).getOperator() == CigarOperator.I) {
                    //append '^' + next-next segment sequence
                    s += "^" + substr(querySequence, readPositionIncludingSoftClipped + 1, cigar.getCigarElement(cigarElementIndex + 1).getLength());

                    //Loop over next-next segment
                    int nextLen = cigar.getCigarElement(cigarElementIndex + 1).getLength();
                    for (int qi = 1; qi <= nextLen; qi++) {
                        //add base quality to total quality
                        q += queryQuality.charAt(readPositionIncludingSoftClipped + 1 + qi) - 33;
                        //increase number of insertion bases
                        qibases++;
                    }

                    //adjust read position by length of next-next segment
                    readPositionIncludingSoftClipped += nextLen;
                    readPositionExcludingSoftClipped += nextLen;
                    cigarElementIndex += 1;
                }
                if (cigar.numCigarElements() > cigarElementIndex + 1
                        && cigar.getCigarElement(cigarElementIndex + 1).getOperator() == CigarOperator.M) {
                    Tuple.Tuple4<Integer, String, String, Integer> tpl =
                            CigarUtils.findOffset(
                                    start + ddlen + 1,
                                    readPositionIncludingSoftClipped + 1,
                                    cigar.getCigarElement(cigarElementIndex + 1).getLength(),
                                    querySequence,
                                    queryQuality,
                                    reference,
                                    refCoverage
                            );
                    int toffset = tpl._1;
                    if (toffset != 0) {
                        moffset = toffset;
                        nmoff += tpl._4;
                        s += "&" + tpl._2;
                        String tq = tpl._3;
                        for (int qi = 0; qi < tq.length(); qi++) {
                            q += tq.charAt(qi) - 33;
                            qibases++;
                        }
                    }
                }
            }
            if (trim == false) {
                //If start - qbases + 1 is in region of interest
                final int pos = start - qbases + 1;
                if (pos >= region.start && pos <= region.end) {
                    //add variation record for $s
                    Variation hv = VariationUtils.getVariation(nonInsertionVariants, pos, s); //reference to variant structure
                    hv.incCountForDir(direction);

                    if(CigarUtils.isBEGIN_ATGC_AMP_ATGCs_END(s)) {
                        //if s is one base followed by '&' and one or more bases
                        //add variant record for s to mnp
                        CigarUtils.increment(mnp, pos, s);
                    }

                    //increment count
                    ++hv.varsCount;

                    //minimum of positions from start of read and end of read
                    int tp = readPositionExcludingSoftClipped < readLengthForMatchedBases - readPositionExcludingSoftClipped
                            ? readPositionExcludingSoftClipped + 1
                            : readLengthForMatchedBases - readPositionExcludingSoftClipped;

                    //average quality of bases in the variation
                    q = q / (qbases + qibases);

                    //isAtLeastAt2Positions is a flag that is 1 if the variant is covered by at least 2 read segments with different positions
                    if (hv.isAtLeastAt2Positions == false && hv.previousPosition != 0 && tp != hv.previousPosition) {
                        hv.isAtLeastAt2Positions = true;
                    }

                    //hasAtLeast2DiffQualities is a flag that is 1 if the variant is covered by at least 2 segment reads with different qualities
                    if (hv.hasAtLeast2DiffQualities == false && hv.previousQuality != 0 && q != hv.previousQuality) {
                        hv.hasAtLeast2DiffQualities = true;
                    }
                    hv.meanPosition += tp;
                    hv.meanQuality += q;
                    hv.meanMappingQuality += mappingQuality;
                    hv.previousPosition = tp;
                    hv.previousQuality = q;
                    hv.numberOfMismatches += numberOfMismatches - nmoff;
                    if (q >= instance().conf.goodq) {
                        hv.highQualityReadsCount++;
                    } else {
                        hv.lowQualityReadsCount++;
                    }

                    //increase coverage for bases covered by the variation
                    for (int qi = 1; qi <= qbases; qi++) {
                        VariationUtils.incCnt(refCoverage, start - qi + 1, 1);
                    }

                    //If variation starts with a deletion ('-' character)
                    if (startWithDeletion) {
                        //add variation to deletions map
                        CigarUtils.increment(positionToDeletionCount, pos, s);

                        //increase coverage for next CIGAR segment
                        for (int qi = 1; qi < ddlen; qi++) {
                            VariationUtils.incCnt(refCoverage, start + qi, 1);
                        }
                    }
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
        }
        if (moffset != 0) {
            offset = moffset;
            readPositionIncludingSoftClipped += moffset;
            start += moffset;
            readPositionExcludingSoftClipped += moffset;
        }
        return cigarElementIndex;
    }

    /**
     * Dealing with CIGAR deletions.
     * @param querySequence base sequence
     * @param mappingQuality record's mapping quality
     * @param queryQuality base quality sequence
     * @param numberOfMismatches number of mismatches
     * @param direction is read has negative strand flag
     * @param readLengthForMatchedBases the read length for matched bases
     * @param cigarElementIndex index number of CIGAR element
     * @return index number of CIGAR element
     */
    // perl version: 1064
    private int processCigarDeletions(String querySequence,
                                      int mappingQuality,
                                      String queryQuality,
                                      int numberOfMismatches,
                                      boolean direction,
                                      int readLengthForMatchedBases,
                                      int cigarElementIndex) {
        //description string of deleted segment
        StringBuilder deletedSegmentBases = new StringBuilder("-").append(cigarElementLength);
        //sequence to be appended if next segment is matched
        StringBuilder additionalSequence = new StringBuilder();
        //quality of last base before deletion
        char lastBAseQuality = queryQuality.charAt(readPositionIncludingSoftClipped - 1);
        //quality of this segment
        StringBuilder deletedSegmentQuality = new StringBuilder();

        // For multiple indels within $VEXT bp
        //offset for read position if next segment is matched
        int offsetForReadPosition = 0;
        //offset for reference position if next segment is matched
        int offsetForReferencePosition = 0;
        int nmoff = 0;

                    /*
                    Condition:
                    1). CIGAR string has next entry
                    2). length of next CIGAR segment is less than instance().conf.vext
                    3). next segment is matched
                    4). CIGAR string has one more entry after next one
                    5). this entry is insertion or deletion
                     */
        if (instance().conf.performLocalRealignment && cigar.numCigarElements() > cigarElementIndex + 2
                && cigar.getCigarElement(cigarElementIndex + 1).getLength() <= instance().conf.vext
                && cigar.getCigarElement(cigarElementIndex + 1).getOperator() == CigarOperator.M
                && (cigar.getCigarElement(cigarElementIndex + 2).getOperator() == CigarOperator.I
                || cigar.getCigarElement(cigarElementIndex + 2).getOperator() == CigarOperator.D)) {

            int mLen = cigar.getCigarElement(cigarElementIndex + 1).getLength();
            int indelLen = cigar.getCigarElement(cigarElementIndex + 2).getLength();

            //append '#' + next matched segment from read
            deletedSegmentBases.append("#").append(substr(querySequence, readPositionIncludingSoftClipped, mLen));
            //append quality string of next matched segment from read
            deletedSegmentQuality.append(substr(queryQuality, readPositionIncludingSoftClipped, mLen));

            //if an insertion is two segments ahead, append '^' + part of sequence corresponding to next-next segment
            //otherwise (deletion) append '^' + length of a next-next segment
            deletedSegmentBases
                    .append('^')
                    .append(
                            cigar.getCigarElement(cigarElementIndex + 2).getOperator() == CigarOperator.I
                                    ? substr(querySequence, readPositionIncludingSoftClipped + mLen, indelLen)
                                    : indelLen
                    );
            //same for quality string
            deletedSegmentQuality
                    .append(
                            cigar.getCigarElement(cigarElementIndex + 2).getOperator() == CigarOperator.I
                                    ? substr(queryQuality, readPositionIncludingSoftClipped + mLen, indelLen)
                                    : ""
                    );

            //add length of next segment to both read and reference offsets
            //add length of next-next segment to reference position (for insertion) or to read position(for deletion)
            offsetForReadPosition += mLen + (cigar.getCigarElement(cigarElementIndex + 2).getOperator() == CigarOperator.D ? indelLen : 0);
            offsetForReferencePosition += mLen + (cigar.getCigarElement(cigarElementIndex + 2).getOperator() == CigarOperator.I ? indelLen : 0);
            if (cigar.numCigarElements() > cigarElementIndex + 3
                    && cigar.getCigarElement(cigarElementIndex + 3).getOperator() == CigarOperator.M) {
                int vsn = 0;
                int tn = readPositionIncludingSoftClipped + offsetForReferencePosition;
                int ts = start + offsetForReadPosition + cigarElementLength;
                for (int vi = 0; vsn <= instance().conf.vext && vi < cigar.getCigarElement(cigarElementIndex + 3).getLength(); vi++) {
                    if (querySequence.charAt(tn + vi) == 'N') {
                        break;
                    }
                    if (queryQuality.charAt(tn + vi) - 33 < instance().conf.goodq) {
                        break;
                    }
                    if (isHasAndEquals('N', reference, ts + vi)) {
                        break;
                    }
                    Character refCh = reference.get(ts + vi);
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
                    additionalSequence.append(substr(querySequence, tn, offset));
                    deletedSegmentQuality.append(substr(queryQuality, tn, offset));
                }
            }
            // skip next 2 CIGAR segments
            cigarElementIndex += 2;
        } else if (instance().conf.performLocalRealignment
                && cigar.numCigarElements() > cigarElementIndex + 1
                && cigar.getCigarElement(cigarElementIndex + 1).getOperator() == CigarOperator.I) {
            /*
            Condition:
            1). CIGAR string has next entry
            2). next CIGAR segment is an insertion
             */

            int insLen = cigar.getCigarElement(cigarElementIndex + 1).getLength();
            //Append '^' + next segment (inserted)
            deletedSegmentBases.append("^").append(substr(querySequence, readPositionIncludingSoftClipped, insLen));
            //Append next segement to quality string
            deletedSegmentQuality.append(substr(queryQuality, readPositionIncludingSoftClipped, insLen));

            //Shift reference position by length of next segment
            //skip next CIGAR segment
            offsetForReferencePosition += insLen;
            if (cigar.numCigarElements() > cigarElementIndex + 2
                    && cigar.getCigarElement(cigarElementIndex + 2).getOperator() == CigarOperator.M) {
                int mLen = cigar.getCigarElement(cigarElementIndex + 2).getLength();
                int vsn = 0;
                int tn = readPositionIncludingSoftClipped + offsetForReferencePosition;
                int ts = start + cigarElementLength;
                for (int vi = 0; vsn <= instance().conf.vext && vi < mLen; vi++) {
                    char seqCh = querySequence.charAt(tn + vi);
                    if (seqCh == 'N') {
                        break;
                    }
                    if (queryQuality.charAt(tn + vi) - 33 < instance().conf.goodq) {
                        break;
                    }
                    Character refCh = reference.get(ts + vi);
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
                    additionalSequence.append(substr(querySequence, tn, offset));
                    deletedSegmentQuality.append(substr(queryQuality, tn, offset));
                }
            }
            cigarElementIndex += 1;
        } else {
            /*
            Condition:
            1). CIGAR string has next entry
            2). next CIGAR segment is matched
             */
            if (instance().conf.performLocalRealignment && cigar.numCigarElements() > cigarElementIndex + 1
                    && cigar.getCigarElement(cigarElementIndex + 1).getOperator() == CigarOperator.M) {
                int mLen = cigar.getCigarElement(cigarElementIndex + 1).getLength();
                int vsn = 0;
                //Loop over next CIGAR segment (no more than instance().conf.vext bases ahead)
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
                    Character refCh = reference.get(start + cigarElementLength + vi);
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
                    //Append first offset bases of next segment to additionalSequence and deletedSegmentQuality
                    additionalSequence.append(substr(querySequence, readPositionIncludingSoftClipped, offset));
                    deletedSegmentQuality.append(substr(queryQuality, readPositionIncludingSoftClipped, offset));
                }
            }
        }
        //offset should be reset to 0 on every loop, so it is non-zero only if previous next CIGAR segment is matched
        // Append '&' and additionalSequence to deletedSegmentBases if next segment has good matching base
        if (offset > 0) {
            deletedSegmentBases.append("&").append(additionalSequence);
        }

        //quality of first matched base after deletion
        //append best of $lastBAseQuality and $q2
        if (readPositionIncludingSoftClipped + offset >= queryQuality.length()) {
            deletedSegmentQuality.append(lastBAseQuality);
        } else {
            char q2 = queryQuality.charAt(readPositionIncludingSoftClipped + offset);
            deletedSegmentQuality.append(lastBAseQuality > q2 ? lastBAseQuality : q2);
        }

        //If reference position is inside region of interest
        if (start >= region.start && start <= region.end) {
            //add variant structure for deletion at this position
            Variation hv = VariationUtils.getVariation(nonInsertionVariants, start, deletedSegmentBases.toString()); //variation structure
            //add record for deletion in deletions map
            CigarUtils.increment(positionToDeletionCount, start, deletedSegmentBases.toString());
            hv.incCountForDir(direction);
            //increase count
            hv.varsCount++;

            //minimum of positions from start of read and end of read
            int tp = readPositionExcludingSoftClipped < readLengthForMatchedBases - readPositionExcludingSoftClipped
                    ? readPositionExcludingSoftClipped + 1
                    : readLengthForMatchedBases - readPositionExcludingSoftClipped;

            //average quality of bases
            double tmpq = 0;

            for (int i = 0; i < deletedSegmentQuality.length(); i++) {
                tmpq += deletedSegmentQuality.charAt(i) - 33;
            }

            tmpq = tmpq / deletedSegmentQuality.length();

            //isAtLeastAt2Positions is a flag that is 1 if the variant is covered by at least 2 read segments with different positions
            if (hv.isAtLeastAt2Positions == false && hv.previousPosition != 0 && tp != hv.previousPosition) {
                hv.isAtLeastAt2Positions = true;
            }

            //hasAtLeast2DiffQualities is a flag that is 1 if the variant is covered by at least 2 segment reads with different qualities
            if (hv.hasAtLeast2DiffQualities == false && hv.previousQuality != 0 && tmpq != hv.previousQuality) {
                hv.hasAtLeast2DiffQualities = true;
            }
            hv.meanPosition += tp;
            hv.meanQuality += tmpq;
            hv.meanMappingQuality += mappingQuality;
            hv.previousPosition = tp;
            hv.previousQuality = tmpq;
            hv.numberOfMismatches += numberOfMismatches - nmoff;
            if (tmpq >= instance().conf.goodq) {
                hv.highQualityReadsCount++;
            } else {
                hv.lowQualityReadsCount++;
            }

            //increase coverage count for reference bases missing from the read
            for (int i = 0; i < cigarElementLength; i++) {
                VariationUtils.incCnt(refCoverage, start + i, 1);
            }
        }

        //adjust reference position by offset + offsetForReadPosition
        start += cigarElementLength + offset + offsetForReadPosition;

        //adjust read position by cigarElementLength (CIGAR segment length) + offset + offsetForReferencePosition
        readPositionIncludingSoftClipped += offset + offsetForReferencePosition;
        readPositionExcludingSoftClipped += offset + offsetForReferencePosition;
        return cigarElementIndex;
    }

    /**
     * Dealing with CIGAR insertions.
     * @param querySequence base sequence
     * @param mappingQuality record's mapping quality
     * @param queryQuality base quality sequence
     * @param numberOfMismatches number of mismatches
     * @param direction is read has negative strand flag
     * @param readLengthForMatchedBases the read length for matched bases
     * @param cigarElementIndex index number of CIGAR element
     * @return index number of CIGAR element
     */
    // perl version: 968
    private int processCigarInsertions(String querySequence,
                                       int mappingQuality,
                                       String queryQuality,
                                       int numberOfMismatches,
                                       boolean direction,
                                       int readLengthForMatchedBases,
                                       int cigarElementIndex) {
        //inserted segment of read sequence
        StringBuilder insertedSegmentBases = new StringBuilder(substr(querySequence, readPositionIncludingSoftClipped, cigarElementLength));
        //quality of this segment
        StringBuilder insertedSegmentQuality = new StringBuilder(substr(queryQuality, readPositionIncludingSoftClipped, cigarElementLength));
        //sequence to be appended if next segment is matched
        String additionalSequence = "";

        // For multiple indels within 10bp

        //offset for read position if next segment is matched
        int offsetForReadPosition = 0;
        //offset for reference position if next segment is matched
        int offsetForReferencePosition = 0;
        int nmoff = 0;

                    /*
                    Condition:
                    1). CIGAR string has next entry
                    2). length of next CIGAR segment is less than instance().conf.vext
                    3). next segment is matched
                    4). CIGAR string has one more entry after next one
                    5). this entry is insertion or deletion
                     */
        if (instance().conf.performLocalRealignment && cigar.numCigarElements() > cigarElementIndex + 2
                && cigar.getCigarElement(cigarElementIndex + 1).getLength() <= instance().conf.vext
                && cigar.getCigarElement(cigarElementIndex + 1).getOperator() == CigarOperator.M
                && (cigar.getCigarElement(cigarElementIndex + 2).getOperator()  == CigarOperator.I
                || cigar.getCigarElement(cigarElementIndex + 2).getOperator()  == CigarOperator.D)) {

            int mLen = cigar.getCigarElement(cigarElementIndex + 1).getLength();
            int indelLen = cigar.getCigarElement(cigarElementIndex + 2).getLength();
            CigarOperator indelOperator = cigar.getCigarElement(cigarElementIndex + 2).getOperator();
            //append to insertedSegmentBases '#' and part of read sequence corresponding to next CIGAR segment (matched one)
            insertedSegmentBases.append("#").append(substr(querySequence, readPositionIncludingSoftClipped + cigarElementLength, mLen));
            //append next segment quality to insertedSegmentQuality
            insertedSegmentQuality.append(substr(queryQuality, readPositionIncludingSoftClipped + cigarElementLength, mLen));

            //if an insertion is two segments ahead, append '^' + part of sequence corresponding to next-next segment
            //otherwise (deletion) append '^' + length of a next-next segment
            insertedSegmentBases
                    .append('^')
                    .append(
                            indelOperator == CigarOperator.I
                                    ? substr(querySequence, readPositionIncludingSoftClipped + cigarElementLength + mLen, indelLen)
                                    : indelLen
                    );
            //if an insertion is two segments ahead, append part of quality string sequence corresponding to next-next segment
            //otherwise (deletion) append first quality score of next segment
            insertedSegmentQuality
                    .append(
                            indelOperator == CigarOperator.I
                                    ? substr(queryQuality, readPositionIncludingSoftClipped + cigarElementLength + mLen, indelLen)
                                    : queryQuality.charAt(readPositionIncludingSoftClipped + cigarElementLength + mLen)
                    );

            //add length of next segment to both offsetForReadPosition and offsetForReferencePosition
            //add length of next-next segment to offsetForReferencePosition (for insertion) or to offsetForReadPosition (for deletion)
            offsetForReadPosition += mLen + (indelOperator == CigarOperator.D ? indelLen : 0);
            offsetForReferencePosition += mLen + (indelOperator == CigarOperator.I ? indelLen : 0);

            int ci6 = cigar.numCigarElements() > cigarElementIndex + 3 ? cigar.getCigarElement(cigarElementIndex + 3).getLength() : 0;
            if (ci6 != 0 && cigar.getCigarElement(cigarElementIndex + 3).getOperator()  == CigarOperator.M) {
                Tuple.Tuple4<Integer, String, String, Integer> tpl = CigarUtils.findOffset(
                        start + offsetForReadPosition,
                        readPositionIncludingSoftClipped + cigarElementLength + offsetForReferencePosition,
                        ci6,
                        querySequence,
                        queryQuality,
                        reference,
                        refCoverage
                );
                offset = tpl._1;
                additionalSequence = tpl._2;
                insertedSegmentQuality.append(tpl._3);
            }
            //skip 2 CIGAR segments
            cigarElementIndex += 2;
        } else {
            /*
            Condition:
            1). CIGAR string has next entry
            2). next CIGAR segment is matched
             */
            if (instance().conf.performLocalRealignment
                    && cigar.numCigarElements() > cigarElementIndex + 1
                    && cigar.getCigarElement(cigarElementIndex + 1).getOperator() == CigarOperator.M) {
                int vsn = 0;
                //Loop over next CIGAR segment (no more than instance().conf.vext bases ahead)
                for (int vi = 0; vsn <= instance().conf.vext && vi < cigar.getCigarElement(cigarElementIndex + 1).getLength(); vi++) {
                    //If base is unknown, exit loop
                    if (querySequence.charAt(readPositionIncludingSoftClipped + cigarElementLength + vi) == 'N') {
                        break;
                    }
                    //If base quality is less than instance().conf.goodq, exit loop
                    if (queryQuality.charAt(readPositionIncludingSoftClipped + cigarElementLength + vi) - 33 < instance().conf.goodq) {
                        break;
                    }
                    //If reference sequence has base at this position and it matches read base, update offset
                    if (reference.containsKey(start + vi)) {
                        if (isNotEquals(querySequence.charAt(readPositionIncludingSoftClipped + cigarElementLength + vi), reference.get(start + vi))) {
                            offset = vi + 1;
                            vsn = 0;
                        } else {
                            vsn++;
                        }
                    }
                }
                if (offset != 0) { //If next CIGAR segment has good matching base
                    //Append first offset bases of next segment to additionalSequence and insertedSegmentQuality
                    additionalSequence += substr(querySequence, readPositionIncludingSoftClipped + cigarElementLength, offset);
                    insertedSegmentQuality.append(substr(queryQuality, readPositionIncludingSoftClipped + cigarElementLength, offset));
                    //Increase coverage for positions corresponding to first offset bases of next segment
                    for (int osi = 0; osi < offset; osi++) {
                        VariationUtils.incCnt(refCoverage, start + osi, 1);
                    }
                }
            }
        }

        //offset should be reset to 0 on every loop, so it is non-zero only if previous part of code was executed
        //Append '&' and $additionalSequence to insertedSegmentBases if next segment has good matching base
        if (offset > 0) {
            insertedSegmentBases.append("&").append(additionalSequence);
        }

        //If start of the segment is within region of interest and the segment does not have unknown bases
        if (start - 1 >= region.start && start - 1 <= region.end && !insertedSegmentBases.toString().contains("N")) {
            //add '+' + insertedSegmentBases to insertions at start - 1
            VariationUtils.incCnt(getOrElse(positionToInsertionCount, start - 1, new HashMap<String, Integer>()), "+" + insertedSegmentBases, 1);
            //add insertion to table of variations
            Variation hv = VariationUtils.getVariation(insertionVariants, start - 1, "+" + insertedSegmentBases); //variant structure for this insertion
            hv.incCountForDir(direction);
            //add count
            hv.varsCount++;
            //minimum of positions from start of read and end of read
            int tp = readPositionExcludingSoftClipped < readLengthForMatchedBases - readPositionExcludingSoftClipped
                    ? readPositionExcludingSoftClipped + 1
                    : readLengthForMatchedBases - readPositionExcludingSoftClipped;

            //mean read quality of the segment
            double tmpq = 0;
            for (int i = 0; i < insertedSegmentQuality.length(); i++) {
                tmpq += insertedSegmentQuality.charAt(i) - 33;
            }
            tmpq = tmpq / insertedSegmentQuality.length();

            //isAtLeastAt2Positions is a flag that is 1 if the variant is covered by at least 2 read segments with different positions
            if (hv.isAtLeastAt2Positions == false && hv.previousPosition != 0 && tp != hv.previousPosition) {
                hv.isAtLeastAt2Positions = true;
            }
            //hasAtLeast2DiffQualities is a flag that is 1 if the variant is covered by at least 2 segment reads with different qualities
            if (hv.hasAtLeast2DiffQualities == false && hv.previousQuality != 0 && tmpq != hv.previousQuality) {
                hv.hasAtLeast2DiffQualities = true;
            }
            hv.meanPosition += tp;
            hv.meanQuality += tmpq;
            hv.meanMappingQuality += mappingQuality;
            hv.previousPosition = tp;
            hv.previousQuality = tmpq;
            if (tmpq >= instance().conf.goodq) {
                hv.highQualityReadsCount++;
            } else {
                hv.lowQualityReadsCount++;
            }
            hv.numberOfMismatches += numberOfMismatches - nmoff;

            //perl version: 1039
            // Adjust the reference count for insertion reads

            /*
            Condition:
            1). reference sequence has base for the position
            2). nonInsertionVariants contains variant structure for the position
            3). read base at position descriptionString-1 matches reference at start-1
             */
            if (VariationUtils.getVariationMaybe(nonInsertionVariants, start - 1, reference.get(start - 1)) != null
                    && isHasAndEquals(querySequence.charAt(readPositionIncludingSoftClipped - 1), reference, start - 1)) {

                // subCnt(getVariation(nonInsertionVariants, start - 1, reference.get(start - 1 ).toString()), dir, tp, tmpq,
                // meanMappingQuality, numberOfMismatches, instance().conf);
                Variation tv = VariationUtils.getVariation(nonInsertionVariants, start - 1, String.valueOf(querySequence.charAt(readPositionIncludingSoftClipped - 1)));
                //Substract count.
                CigarUtils.subCnt(tv, direction, tp, queryQuality.charAt(readPositionIncludingSoftClipped - 1) - 33, mappingQuality, numberOfMismatches);
            }
            // Adjust count if the insertion is at the edge so that the AF won't > 1
            /*
            Condition:
            1). looking at second segment in CIGAR string
            2). first segment is a soft-clipping or a hard-clipping
            */
            if (cigarElementIndex == 1 && (cigar.getCigarElement(0).getOperator() == CigarOperator.S || cigar.getCigarElement(0).getOperator() == CigarOperator.H)) {
                //Add one more variant corresponding to base at start - 1 to nonInsertionVariants
                Variation ttref = VariationUtils.getVariation(nonInsertionVariants, start - 1, reference.get(start - 1).toString());
                ttref.incCountForDir(direction);
                ttref.varsCount++;
                ttref.isAtLeastAt2Positions = hv.isAtLeastAt2Positions;
                ttref.hasAtLeast2DiffQualities = hv.hasAtLeast2DiffQualities;
                ttref.meanPosition += tp;
                ttref.meanQuality += tmpq;
                ttref.meanMappingQuality += mappingQuality;
                ttref.previousPosition = tp;
                ttref.previousQuality = tmpq;
                ttref.numberOfMismatches += numberOfMismatches - nmoff;
                VariationUtils.incCnt(refCoverage, start - 1, 1);
            }
        }

        //adjust read position by cigarElementLength (CIGAR segment length) + offset + offsetForReferencePosition
        readPositionIncludingSoftClipped += cigarElementLength + offset + offsetForReferencePosition;
        readPositionExcludingSoftClipped += cigarElementLength + offset + offsetForReferencePosition;
        //adjust reference position by offset + offsetForReadPosition
        start += offset + offsetForReadPosition;
        return cigarElementIndex;
    }

    /**
     * Dealing with last CIGAR element in case of soft-clipped reads (if CIGAR operator is S)
     * @param querySequence base sequence
     * @param mappingQuality record's mapping quality
     * @param queryQuality base quality sequence
     * @param numberOfMismatches number of mismatches
     * @param direction is read has negative strand flag
     * @param totalLengthIncludingSoftClipped the total length, including soft-clipped bases
     */
    // perl version: 932
    private void processLastCigarElementWithSoftClipping(String querySequence,
                                                         int mappingQuality,
                                                         String queryQuality,
                                                         int numberOfMismatches,
                                                         boolean direction,
                                                         int totalLengthIncludingSoftClipped) {
    /*
    Conditions:
    1). read position is less than sequence length
    2). reference base is defined for start
    3). reference base at start matches read base at descriptionString
    4). read quality is more than 10
     */
        while (readPositionIncludingSoftClipped < querySequence.length()
                && isHasAndEquals(querySequence.charAt(readPositionIncludingSoftClipped), reference, start)
                && queryQuality.charAt(readPositionIncludingSoftClipped) - 33 > 10) {
            //initialize entry in $nonInsertionVariants if not present
            Variation variation = VariationUtils.getVariation(nonInsertionVariants, start, reference.get(start).toString());
            //add count
            CigarUtils.addCnt(
                    variation,
                    direction,
                    totalLengthIncludingSoftClipped - readPositionExcludingSoftClipped,
                    queryQuality.charAt(readPositionIncludingSoftClipped) - 33,
                    mappingQuality,
                    numberOfMismatches
            );
            //add coverage
            VariationUtils.incCnt(refCoverage, start, 1);
            readPositionIncludingSoftClipped++;
            start++;
            cigarElementLength--;
            readPositionExcludingSoftClipped++;
        }
        if (querySequence.length() - readPositionIncludingSoftClipped > 0) { //If there remains a soft-clipped sequence at the end (not everything was matched)
            int q = 0; //sum of read qualities (to get mean quality)
            int qn = 0; //number of quality figures
            int lowqcnt = 0; //number of low-quality bases in the sequence
            for (int si = 0; si < cigarElementLength; si++) { //loop over remaining soft-clipped sequence
                if (querySequence.charAt(readPositionIncludingSoftClipped + si) == 'N') { //stop if unknown base (N - any of ATGC) is found
                    break;
                }
                int tq = queryQuality.charAt(readPositionIncludingSoftClipped + si) - 33; //base quality
                if (tq <= 12) {
                    lowqcnt++;
                }
                //Stop if a low-quality base is found
                if (lowqcnt > 1) {
                    break;
                }
                q += tq;
                qn++;
            }
            //If we have at least 1 high-quality soft-clipped base within instance().conf.buffer of region of interest
            if (qn >= 1 && qn > lowqcnt && start >= region.start - instance().conf.buffer && start <= region.end + instance().conf.buffer) {
                //add record to $softClips3End
                SoftClip softClip = softClips3End.get(start);
                if (softClip == null) {
                    softClip = new SoftClip();
                    softClips3End.put(start, softClip);
                }
                for (int si = 0; si < qn; si++) {
                    Character ch = querySequence.charAt(readPositionIncludingSoftClipped + si);
                    int idx = si;
                    Map<Character, Integer> cnts = softClip.nt.get(idx);
                    if (cnts == null) {
                        cnts = new HashMap<>();
                        softClip.nt.put(idx, cnts);
                    }
                    VariationUtils.incCnt(cnts, ch, 1);
                    Variation variation = CigarUtils.getVariationFromSeq(softClip, idx, ch);
                    CigarUtils.addCnt(
                            variation,
                            direction,
                            qn - si,
                            queryQuality.charAt(readPositionIncludingSoftClipped + si) - 33,
                            mappingQuality,
                            numberOfMismatches
                    );
                }
                CigarUtils.addCnt(
                        softClip,
                        direction,
                        cigarElementLength,
                        q / (double)qn,
                        mappingQuality,
                        numberOfMismatches
                );
            }
        }
    }

    /**
     * Dealing with first CIGAR element in case of soft-clipped reads (if CIGAR operator is S)
     * @param chrName name of chromosome
     * @param querySequence base sequence
     * @param mappingQuality record's mapping quality
     * @param queryQuality base quality sequence
     * @param numberOfMismatches number of mismatches
     * @param direction is read has negative strand flag
     * @param cigarElementIndex index number of CIGAR element
     */
    // perl version: 901
    private void processFirstCigarElementWithSoftClipping(String chrName,
                                                          String querySequence,
                                                          int mappingQuality,
                                                          String queryQuality,
                                                          int numberOfMismatches,
                                                          boolean direction,
                                                          int cigarElementIndex) {
        // align soft-clipped but matched sequences due to mis-soft-clipping
                        /*
                        Conditions:
                        1). segment length > 1
                        2). start between 0 and chromosome length,
                        3). reference genome is known at this position
                        4). reference and read bases match
                        5). read quality is more than 10
                         */
        while (cigarElementLength - 1 >= 0 && start - 1 > 0 && start - 1 <= instance().chrLengths.get(chrName)
                && isHasAndEquals(querySequence.charAt(cigarElementLength - 1), reference, start - 1)
                && queryQuality.charAt(cigarElementLength - 1) - 33 > 10) {
            //create variant if it is not present
            Variation variation = VariationUtils.getVariation(nonInsertionVariants, start - 1, reference.get(start - 1).toString());
            //add count
            CigarUtils.addCnt(
                    variation,
                    direction,
                    cigarElementLength,
                    queryQuality.charAt(cigarElementLength - 1) - 33,
                    mappingQuality,
                    numberOfMismatches
            );
            //increase coverage
            VariationUtils.incCnt(refCoverage, start - 1, 1);
            start--;
            cigarElementLength--;
        }
        if (cigarElementLength > 0) { //If there remains a soft-clipped sequence at the beginning (not everything was matched)
            int q = 0; //sum of read qualities (to get mean quality)
            int qn = 0; //number of quality figures
            int lowqcnt = 0; //number of low-quality bases in the sequence

            //loop over remaining soft-clipped sequence
            for (int si = cigarElementLength - 1; si >= 0; si--) {
                //stop if unknown base (N - any of ATGC) is found
                if (querySequence.charAt(si) == 'N') {
                    break;
                }
                //tq - base quality
                int tq = queryQuality.charAt(si) - 33;
                if (tq <= 12)
                    lowqcnt++;
                //Stop if a low-quality base is found
                if (lowqcnt > 1)
                    break;

                q += tq;
                qn++;
            }
            //If we have at least 1 high-quality soft-clipped base within instance().conf.buffer of region of interest
            if (qn >= 1 && qn > lowqcnt && start >= region.start - instance().conf.buffer && start <= region.end + instance().conf.buffer) {
                //add record to $softClips5End
                SoftClip softClip = softClips5End.get(start);
                if (softClip == null) {
                    softClip = new SoftClip();
                    softClips5End.put(start, softClip);
                }
                for (int si = cigarElementLength - 1; cigarElementLength - si <= qn; si--) {
                    Character ch = querySequence.charAt(si);
                    int idx = cigarElementLength - 1 - si;
                    Map<Character, Integer> cnts = softClip.nt.get(idx);
                    if (cnts == null) {
                        cnts = new LinkedHashMap<>();
                        softClip.nt.put(idx, cnts);
                    }
                    VariationUtils.incCnt(cnts, ch, 1);
                    Variation seqVariation = CigarUtils.getVariationFromSeq(softClip, idx, ch);
                    CigarUtils.addCnt(
                            seqVariation,
                            direction,
                            si - (cigarElementLength - qn),
                            queryQuality.charAt(si) - 33,
                            mappingQuality,
                            numberOfMismatches
                    );
                }
                CigarUtils.addCnt(
                        softClip,
                        direction,
                        cigarElementLength,
                        q / (double)qn, mappingQuality,
                        numberOfMismatches
                );
            }

        }
        cigarElementLength = cigar.getCigarElement(cigarElementIndex).getLength();
    }

    /**
     * Process CIGAR in amplicon based calling case.
     * @param record {@link SAMRecord}
     * @param isMateReferenceNameEqual true if reference name is equal to mate reference name
     * @return true if skipping records required.
     */
    //perl version: 649
    private boolean parseCigarWithAmpCase(SAMRecord record, boolean isMateReferenceNameEqual) {
        String[] split = instance().ampliconBasedCalling.split(":");
        //distance to amplicon (specified in -a option)
        int distanceToAmplicon;
        //overlap fraction (in -a option)
        double overlapFraction;
        try {
            distanceToAmplicon = toInt(split[0]);
            overlapFraction = Double.parseDouble(split.length > 1 ? split[1] : "");
        } catch (NumberFormatException e) {
            distanceToAmplicon = 10;
            overlapFraction = 0.95;
        }
        //totalAlignedLength holds sum of lengths of matched and deleted segments
        int totalAlignedLength = getAlignedLenght(cigar); // The total aligned length, excluding soft-clipped
        // bases and insertions
        int segmentStart = record.getAlignmentStart();
        int segmentEnd = segmentStart + totalAlignedLength - 1;

        if (cigar.getCigarElement(0).getOperator() == CigarOperator.S) { //If read starts with soft-clipped sequence
            //Ignore reads that overlap with region of interest by fraction less than overlapFraction
            int tmpSegmentStart = segmentStart > region.start ? segmentStart : region.start;
            int tmpSegmentEnd = segmentEnd < region.end ? segmentEnd : region.end;
            if (Math.abs(tmpSegmentStart - tmpSegmentEnd) / (double)(segmentEnd - segmentStart) > overlapFraction == false) {
                return true;
            }
        } else if (cigar.getCigarElement(cigar.numCigarElements() - 1).getOperator() == CigarOperator.S) { //If read ends with contains soft-clipped sequence
            //Ignore reads that overlap with region of interest by fraction less than overlapFraction
            int tmpSegmentStart = segmentStart > region.start ? segmentStart : region.start;
            int tmpSegmentEnd = segmentEnd < region.end ? segmentEnd : region.end;
            if (Math.abs(tmpSegmentEnd - tmpSegmentStart) / (double)(segmentEnd - segmentStart) > overlapFraction == false) {
                return true;
            }
        } else { //no soft-clipping
            //if RNEXT is identical to RNAME and TLEN is defined
            if (isMateReferenceNameEqual && record.getInferredInsertSize() != 0) {
                if (record.getInferredInsertSize() > 0) {
                    segmentEnd = segmentStart + record.getInferredInsertSize() - 1;
                } else {
                    segmentStart = record.getMateAlignmentStart();
                    segmentEnd = record.getMateAlignmentStart() - record.getInferredInsertSize() - 1;
                }
            }
            // No segment overlapping test since samtools should take care of it
            int tmpSegmentStart = segmentStart > region.start ? segmentStart : region.start;
            int tmpSegmentEnd = segmentEnd < region.end ? segmentEnd : region.end;
            //ignore reads that are more than dis from region of interest and overlap is less than ovlp
            if ((abs(segmentStart - region.start) > distanceToAmplicon || abs(segmentEnd - region.end) > distanceToAmplicon)
                    || abs((tmpSegmentStart - tmpSegmentEnd) / (double)(segmentEnd - segmentStart)) <= overlapFraction) {
                return true;
            }
        }
        return false;
    }

    /**
     * Skip the read if number of mismatches is not available or more than acceptable value {@link Configuration#mismatch}
     * @param record {@link SAMRecord}
     * @param insertionDeletionLength insertions and deletions length from CIGAR
     * @param editDistance value of NM tag from record
     * @return true if skipping reads required.
     */
    private boolean processNMTag(SAMRecord record, int insertionDeletionLength, Integer editDistance) {
        if (editDistance != null) { // number of mismatches. Don't use NM since it includes gaps, which can be from indels
            if (editDistance - insertionDeletionLength > instance().conf.mismatch) { // edit distance - indels is the # of mismatches
                return true;
            }
        } else { //Skip the read if number of mismatches is not available
            if (instance().conf.y && !record.getCigarString().equals("*")) {
                System.err.println("No NM tag for mismatches. " + record.getSAMString());
            }
            if (record.getReadUnmappedFlag() || record.getCigarString().equals(SAMRecord.NO_ALIGNMENT_CIGAR)) {
                return true;
            }
        }
        return false;
    }
}
