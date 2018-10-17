package com.astrazeneca.vardict.modules;

import com.astrazeneca.vardict.Configuration;
import com.astrazeneca.vardict.collection.Tuple;
import com.astrazeneca.vardict.collection.VariationMap;
import com.astrazeneca.vardict.data.*;
import com.astrazeneca.vardict.variations.Sclip;
import com.astrazeneca.vardict.variations.Variation;
import htsjdk.samtools.*;

import java.io.IOException;
import java.time.LocalDateTime;
import java.util.*;
import java.util.regex.Matcher;

import static com.astrazeneca.vardict.Utils.*;
import static com.astrazeneca.vardict.collection.Tuple.tuple;
import static com.astrazeneca.vardict.data.Patterns.*;
import static com.astrazeneca.vardict.modules.CigarUtils.*;
import static com.astrazeneca.vardict.modules.StructuralVariantsProcessor.*;
import static com.astrazeneca.vardict.modules.VariationRealigner.adjInsPos;
import static com.astrazeneca.vardict.variations.VariationUtils.*;
import static java.lang.Math.abs;
import static java.lang.String.format;

public class SAMFileParser {
    private static final Random RND = new Random(System.currentTimeMillis());
    /**
     * Construct a variant structure given a region and BAM files.
     * @param region region
     * @param bam BAM file name
     * @param chrs map of chromosome lengths
     * @param sample sample name
     * @param splice set of strings representing spliced regions
     * @param ampliconBasedCalling string of maximum_distance:minimum_overlap for amplicon based calling
     * @param rlen max read length
     * @param reference reference in a given region
     * @param conf Configuration
     * @return Tuple of (noninsertion variant  structure, insertion variant structure, coverage, maxmimum read length)
     * @throws IOException
     */
    public static Tuple.Tuple5<Map<Integer, VariationMap<String, Variation>>,
            Map<Integer, VariationMap<String, Variation>>,
            Map<Integer, Integer>, Integer, Double> parseSAM(Region region, String bam,
                                                             Map<String, Integer> chrs,
                                                             String sample, Set<String> splice,
                                                             String ampliconBasedCalling, int rlen,
                                                             Reference reference, Configuration conf,
                                                             Map<Integer, VariationMap<String, Variation>> hash,
                                                             Map<Integer, VariationMap<String, Variation>> iHash,
                                                             Map<Integer, Integer> cov,
                                                             Map<Integer, Sclip> sclip3, Map<Integer, Sclip> sclip5,
                                                             boolean svflag) throws IOException {

        String[] bams = bam.split(":");

        Map<Integer, Map<String, Integer>> ins = new HashMap<>();
        Map<Integer, Map<String, Integer>> mnp = new HashMap<>(); // Keep track of MNPs
        Map<Integer, Map<String, Integer>> dels5 = new HashMap<>();
        Map<String, int[]> spliceCnt = new HashMap<>();
        Map<Integer, Character> ref = reference.referenceSequences;
        String chr = region.chr;

        // Contains all structural variants ends for deletions, duplications, inversions, insertions
        // and structures for structural variant analysis
        SVStructures svStructures = new SVStructures();

        int totalreads = 0;
        int dupreads = 0;
        int disc = 0;

        //TODO: Remove it? -C deprecated
        if (conf.chromosomeNameIsNumber && chr.startsWith("chr")) { //remove prefix 'chr' if option -C is set
            chr = region.chr.substring("chr".length());
        }

        for (String bami : bams) {
            String samfilter = conf.samfilter == null || conf.samfilter.isEmpty() ? "" : conf.samfilter;
            if (conf.y) {
                System.err.printf("TIME: Start parsing SAM: " + LocalDateTime.now() + " %s %s:%d-%d\n",
                        bami, chr, region.start, region.end);
            }
            try (SamView reader =  new SamView(bami, samfilter, region, conf.validationStringency)) {
                // dup contains already seen reads. For each seen read dup contains either POS-RNEXT-PNEXT or POS-CIGAR
                // (if next segment in template is unmapped).
                Set<String> dup = new HashSet<>();
                // Position of first matching base (POS in SAM)
                int dupp = -1;
                SAMRecord record;
                while ((record = reader.read()) != null) {
                    if (conf.isDownsampling() && RND.nextDouble() <= conf.downsampling) {
                        continue;
                    }

                    final String querySequence = record.getReadString();
                    final Flags flag = new Flags(record.getFlags());
                    final int mappingQuality = record.getMappingQuality();

                    // Ignore low mapping quality reads
                    if (conf.hasMappingQuality() && mappingQuality < conf.mappingQuality) {
                        continue;
                    }

                    if (flag.isNotPrimaryAlignment() && conf.samfilter != null) {
                        continue;
                    }

                    if (querySequence.length() == 1 && querySequence.charAt(0) == '*') {
                        continue;
                    }
                    totalreads++;

                    final String mrnm = getMrnm(record);

                    // Filter duplicated reads if option -t is set
                    if (conf.removeDuplicatedReads) {
                        if (record.getAlignmentStart() != dupp) {
                            dup.clear();
                        }
                        if (record.getMateAlignmentStart() < 10) {
                            String dupKey = record.getAlignmentStart() + "-" + mrnm + "-" + record.getMateAlignmentStart();
                            if (dup.contains(dupKey)) {
                                dupreads++;
                                continue;
                            }
                            dup.add(dupKey);
                            dupp = record.getAlignmentStart();
                        } else if (flag.isUnmappedMate()) {
                            String dupKey = record.getAlignmentStart() + "-" + record.getCigarString();
                            if (dup.contains(dupKey)) {
                                dupreads++;
                                continue;
                            }
                            dup.add(dupKey);
                            dupp = record.getAlignmentStart();
                        }
                    }

                    final Cigar readCigar = record.getCigar();
                    final int indel = getInsertionDeletionLength(readCigar);

                    int tnm = 0;
                    Integer nmi = record.getIntegerAttribute(SAMTag.NM.name());

                    // Number of mismatches. Don't use NM since it includes gaps, which can be from indels
                    if (nmi != null) {
                        tnm = nmi - indel;
                        if (tnm > conf.mismatch) { // Edit distance - indels is the # of mismatches
                            continue;
                        }
                    } else { //Skip the read if number of mismatches is not available
                        if (conf.y && !record.getCigarString().equals("*")) {
                            System.err.println("No NM tag for mismatches. " + record.getSAMString());
                        }
                        // Skip unmapped reads (issue #56)
                        if (record.getReadUnmappedFlag() || record.getCigarString().equals(SAMRecord.NO_ALIGNMENT_CIGAR)) {
                            continue;
                        }
                    }

                    final String queryQuality = record.getBaseQualityString();
                    final boolean isMrnmEqual = record.getReferenceName().equals(record.getMateReferenceName());

                    // Number of mismatches
                    final int nm = tnm;
                    int n = 0; // Keep track the read position, including softclipped
                    int p = 0; // Keep track the position in the alignment, excluding softclipped
                    boolean dir = flag.isReverseStrand();
                    if (ampliconBasedCalling != null) {
                        String[] split = ampliconBasedCalling.split(":");
                        // Distance to amplicon (specified in -a option)
                        int dis;
                        // Overlap fraction (in -a option)
                        double ovlp;
                        try {
                            dis = toInt(split[0]);
                            ovlp = Double.parseDouble(split.length > 1 ? split[1] : "");
                        } catch (NumberFormatException e) {
                            dis = 10;
                            ovlp = 0.95;
                        }
                        // rlen3 holds sum of lengths of matched and deleted segments (The total aligned length,
                        // excluding soft-clipped bases and insertions)
                        int rlen3 = getAlignedLength(readCigar);
                        int segstart = record.getAlignmentStart();
                        int segend = segstart + rlen3 - 1;

                        //If read starts with soft-clipped sequence
                        if (readCigar.getCigarElement(0).getOperator() == CigarOperator.S) {
                            //Ignore reads that overlap with region of interest by fraction less than ovlp
                            int ts1 = segstart > region.start ? segstart : region.start;
                            int te1 = segend < region.end ? segend : region.end;
                            if (Math.abs(ts1 - te1) / (double)(segend - segstart) > ovlp == false) {
                                continue;
                            }
                        } else if (readCigar.getCigarElement(readCigar.numCigarElements() - 1).getOperator() == CigarOperator.S) {
                            //If read ends with contains soft-clipped sequence
                            //Ignore reads that overlap with region of interest by fraction less than ovlp
                            int ts1 = segstart > region.start ? segstart : region.start;
                            int te1 = segend < region.end ? segend : region.end;
                            if (Math.abs(te1 - ts1) / (double)(segend - segstart) > ovlp == false) {
                                continue;
                            }

                        } else { // No soft-clipping
                            // If RNEXT is identical to RNAME and TLEN is defined
                            if (isMrnmEqual && record.getInferredInsertSize() != 0) {
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
                            // Ignore reads that are more than distacnce from region of interest \
                            // and overlap is less than overlap fraction from configuration (-a)
                            if ((abs(segstart - region.start) > dis || abs(segend - region.end) > dis)
                                    || abs((ts1 - te1) / (double)(segend - segstart)) <= ovlp) {
                                continue;
                            }
                        }
                    }

                    final int position;
                    final Cigar cigar;
                    if (conf.performLocalRealignment) {
                        // Modify the CIGAR for potential mis-alignmet for indels at the end of reads to softclipping
                        // and let VarDict's algorithm to figure out indels
                        Tuple.Tuple2<Integer, String> mc = CigarUtils.modifyCigar(indel, ref, record.getAlignmentStart(),
                                record.getCigarString(), querySequence, queryQuality, conf.LOWQUAL);
                        position = mc._1;
                        cigar = TextCigarCodec.decode(mc._2);
                    } else {
                        position = record.getAlignmentStart();
                        cigar = record.getCigar();
                    }

                    cleanupCigar(record);

                    //determine discordant reads
                    boolean mdir = (record.getFlags() & 0x20) != 0 ? false : true;
                    if (!getMrnm(record).equals("=")) {
                        disc++;
                    }

                    if (BEGIN_dig_dig_S_ANY_dig_dig_S_END.matcher(cigar.toString()).find()) {
                        continue; //ignore reads that are softclipped at both ends and both greater than 10 bp
                    }
                    //adjusted start position
                    int start = position;
                    int offset = 0;
                    // Only match and insertion counts toward read length
                    // For total length, including soft-clipped bases
                    int rlen1 = getMatchInsertionLength(cigar); // The read length for matched bases
                    if (conf.minmatch != 0 && rlen1 < conf.minmatch) {
                        continue;
                    }
                    // The total length, including soft-clipped bases
                    int rlen2 = getSoftClippedLength(cigar);
                    if (rlen2 > rlen) { // Determine the read length
                        rlen = rlen2;
                    }

                    // If 'SA' tag (supplementary alignment) is present
                    if (conf.samfilter != null && record.getStringAttribute(SAMTag.SA.name()) != null) {
                        if (flag.isSupplementaryAlignment()) { // the supplementary alignment
                            continue; // Ignore the supplementary for now so that it won't skew the coverage
                        }
                    }

                    //TODO: Determine whether to filter a read in CRISPR mode

                    if (flag.isUnmappedMate()) { //Mate unmapped, potential insertion
                        // To be implemented
                    } else if(record.getMappingQuality() > 10 && !conf.disableSV) {
                        // Consider high mapping quality mates only
                        prepareSVStructuresForAnalysis(rlen, conf, svStructures, record, queryQuality,
                                nm, dir, cigar, mdir, start, rlen2);
                    }
                    int mateAlignmentStart = record.getMateAlignmentStart();

                    processCigar:
                    //Loop over CIGAR records
                    for (int ci = 0; ci < cigar.numCigarElements(); ci++) {
                        if (skipOverlappingReads(conf, record, position, dir, start, mateAlignmentStart)) {
                            break;
                        }
                        // Length of segment in CIGAR
                        int m = cigar.getCigarElement(ci).getLength();

                        //Letter from CIGAR
                        final CigarOperator operator = getCigarOperator(cigar, ci);

                        switch (operator) {
                            case N: //N in CIGAR - skipped region from reference
                                //Skip the region and add string start-end to %SPLICE
                                String key = (start - 1) + "-" + (start + m - 1);
                                splice.add(key);
                                int[] cnt = spliceCnt.get(key);
                                if (cnt == null) {
                                    cnt = new int[] { 0 };
                                    spliceCnt.put(key, cnt);
                                }
                                cnt[0]++;

                                start += m;
                                offset = 0;
                                continue;

                            case S:
                                //First record in CIGAR
                                if (ci == 0) { // 5' soft clipped
                                    // Ignore large soft clip due to chimeric reads in library construction
                                    if (!conf.chimeric) {
                                        String saTagString = record.getStringAttribute(SAMTag.SA.name());

                                        if (m >= 20 && saTagString != null) {
                                            if (isReadChimericWithSA(record, position, saTagString, rlen, dir, conf, true)) {
                                                n += m;
                                                offset = 0;
                                                // Had to reset the start due to softclipping adjustment
                                                start = position;

                                                continue;
                                            }
                                            //trying to detect chimeric reads even when there's no supplementary
                                            // alignment from aligner
                                        } else if (m >= Configuration.SEED_1) {
                                            Map<String, List<Integer>> referenceSeedMap = reference.seed;
                                            String sequence = getReverseComplementedSequence(record, 0, m);
                                            String reverseComplementedSeed = sequence.substring(0, Configuration.SEED_1);

                                            if (referenceSeedMap.containsKey(reverseComplementedSeed)) {
                                                List<Integer> positions = referenceSeedMap.get(reverseComplementedSeed);
                                                if (positions.size() == 1 &&
                                                        abs(start - positions.get(0)) <  2 * rlen) {
                                                    n += m;
                                                    offset = 0;
                                                    // Had to reset the start due to softclipping adjustment
                                                    start = position;
                                                    if (conf.y) {
                                                        System.err.println(sequence + " at 5' is a chimeric at "
                                                                + start + " by SEED " + Configuration.SEED_1);
                                                    }
                                                    continue;
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
                                    while (m - 1 >= 0 && start - 1 > 0 && start - 1 <= chrs.get(chr)
                                            && isHasAndEquals(querySequence.charAt(m - 1), ref, start - 1)
                                            && queryQuality.charAt(m - 1) - 33 > 10) {
                                        //create variant if it is not present
                                        Variation variation = getVariation(hash, start - 1, ref.get(start - 1).toString());
                                        //add count
                                        addCnt(variation, dir, m, queryQuality.charAt(m - 1) - 33, mappingQuality, nm, conf.goodq);
                                        //increase coverage
                                        incCnt(cov, start - 1, 1);
                                        start--;
                                        m--;
                                    }
                                    if (m > 0) { //If there remains a soft-clipped sequence at the beginning (not everything was matched)
                                        int q = 0; // Sum of read qualities (to get mean quality)
                                        int qn = 0; // Number of quality figures
                                        int lowqcnt = 0; // Number of low-quality bases in the sequence

                                        // Loop over remaining soft-clipped sequence
                                        for (int si = m - 1; si >= 0; si--) {
                                            // Stop if unknown base (N - any of ATGC) is found
                                            if (querySequence.charAt(si) == 'N') {
                                                break;
                                            }
                                            // Tq - base quality
                                            int tq = queryQuality.charAt(si) - 33;
                                            if (tq <= 12)
                                                lowqcnt++;
                                            //Stop if a low-quality base is found
                                            if (lowqcnt > 1) {
                                                break;
                                            }
                                            q += tq;
                                            qn++;
                                        }
                                        sclip5HighQualityProcessing(region, conf, sclip5, querySequence, mappingQuality,
                                                queryQuality, nm, dir, start, m, q, qn, lowqcnt);
                                    }
                                    m = cigar.getCigarElement(ci).getLength();
                                } else if (ci == cigar.numCigarElements() - 1) { // 3' soft clipped
                                    // Ignore large soft clip due to chimeric reads in library construction
                                    if (!conf.chimeric) {
                                        String saTagString = record.getStringAttribute(SAMTag.SA.name());
                                        if (m >= 20 && saTagString != null) {
                                            if (isReadChimericWithSA(record, position, saTagString, rlen,
                                                    dir, conf, false)) {
                                                n += m;
                                                offset = 0;
                                                // Had to reset the start due to softclipping adjustment
                                                start = position;

                                                continue;
                                            }
                                        } else if (m >= Configuration.SEED_1) {
                                            Map<String, List<Integer>> referenceSeedMap = reference.seed;
                                            String sequence = getReverseComplementedSequence(record, -m, m);
                                            String reverseComplementedSeed = substr(sequence, -Configuration.SEED_1, Configuration.SEED_1);

                                            if (referenceSeedMap.containsKey(reverseComplementedSeed)) {
                                                List<Integer> positions = referenceSeedMap.get(reverseComplementedSeed);
                                                if (positions.size() == 1 && abs(start - positions.get(0)) <  2 * rlen) {
                                                    n += m;
                                                    offset = 0;
                                                    // Had to reset the start due to softclipping adjustment
                                                    start = position;
                                                    if (conf.y) {
                                                        System.err.println(sequence  + " at 3' is a chimeric at "
                                                                + start + " by SEED " + Configuration.SEED_1);
                                                    }
                                                    continue;
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
                                    while (n < querySequence.length()
                                            && isHasAndEquals(querySequence.charAt(n), ref, start)
                                            && queryQuality.charAt(n) - 33 > 10) {
                                        //Initialize entry in $hash if not present
                                        Variation variation = getVariation(hash, start, ref.get(start).toString());
                                        //Add count
                                        addCnt(variation, dir, rlen2 - p, queryQuality.charAt(n) - 33, mappingQuality, nm, conf.goodq);
                                        //Add coverage
                                        incCnt(cov, start, 1);
                                        n++;
                                        start++;
                                        m--;
                                        p++;
                                    }
                                    //If there remains a soft-clipped sequence at the end (not everything was matched)
                                    if (querySequence.length() - n > 0) {
                                        int q = 0; //Sum of read qualities (to get mean quality)
                                        int qn = 0; //Number of quality figures
                                        int lowqcnt = 0; //Number of low-quality bases in the sequence
                                        for (int si = 0; si < m; si++) { //Loop over remaining soft-clipped sequence
                                            //Stop if unknown base (N - any of ATGC) is found
                                            if (querySequence.charAt(n + si) == 'N') {
                                                break;
                                            }
                                            int tq = queryQuality.charAt(n + si) - 33; //Base quality
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
                                        sclip3HighQualityProcessing(region, conf, sclip3, querySequence, mappingQuality,
                                                queryQuality, nm, n, dir, start, m, q, qn, lowqcnt);
                                    }
                                }
                                //Move read position by m (length of segment in CIGAR)
                                n += m;
                                offset = 0;
                                start = position; // Had to reset the start due to softclipping adjustment
                                continue;
                            case H: //Hard clipping - skip
                                offset = 0;
                                continue;
                            case I: { //Insertion
                                offset = 0;
                                //Inserted segment of read sequence
                                StringBuilder s = new StringBuilder(substr(querySequence, n, m));

                                //Quality of this segment
                                StringBuilder q = new StringBuilder(substr(queryQuality, n, m));
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
                                if (conf.performLocalRealignment && cigar.numCigarElements() > ci + 2
                                        && cigar.getCigarElement(ci + 1).getLength() <= conf.vext
                                        && cigar.getCigarElement(ci + 1).getOperator() == CigarOperator.M
                                        && (cigar.getCigarElement(ci + 2).getOperator()  == CigarOperator.I
                                        || cigar.getCigarElement(ci + 2).getOperator()  == CigarOperator.D)
                                        && cigar.getCigarElement(ci + 3).getOperator() != CigarOperator.I
                                        && cigar.getCigarElement(ci + 3).getOperator() != CigarOperator.D) {

                                    int mLen = cigar.getCigarElement(ci + 1).getLength();
                                    int indelLen = cigar.getCigarElement(ci + 2).getLength();

                                    int begin = n + m;
                                    appendSegments(querySequence, queryQuality, cigar, ci, s, q, mLen, indelLen, begin, true);
                                    //add length of next segment to both multoffs and multoffp
                                    //add length of next-next segment to multoffp (for insertion) or to multoffs (for deletion)
                                    multoffs += mLen + (cigar.getCigarElement(ci + 2).getOperator() == CigarOperator.D ? indelLen : 0);
                                    multoffp += mLen + (cigar.getCigarElement(ci + 2).getOperator() == CigarOperator.I ? indelLen : 0);

                                    int ci6 = cigar.numCigarElements() > ci + 3 ? cigar.getCigarElement(ci + 3).getLength() : 0;
                                    if (ci6 != 0 && cigar.getCigarElement(ci + 3).getOperator()  == CigarOperator.M) {
                                        Tuple.Tuple4<Integer, String, String, Integer> tpl = findOffset(start + multoffs,
                                                n + m + multoffp, ci6, querySequence, queryQuality, ref, cov, conf.vext, conf.goodq);
                                        offset = tpl._1;
                                        ss = tpl._2;
                                        q.append(tpl._3);
                                    }
                                    //skip 2 CIGAR segments
                                    ci += 2;
                                } else {
                                    /*
                                    Condition:
                                    1). CIGAR string has next entry
                                    2). next CIGAR segment is matched
                                     */
                                    if (conf.performLocalRealignment && cigar.numCigarElements() > ci + 1
                                            && cigar.getCigarElement(ci + 1).getOperator() == CigarOperator.M) {
                                        int vsn = 0;
                                        //Loop over next CIGAR segment (no more than conf.vext bases ahead)
                                        for (int vi = 0; vsn <= conf.vext && vi < cigar.getCigarElement(ci + 1).getLength(); vi++) {
                                            //If base is unknown, exit loop
                                            if (querySequence.charAt(n + m + vi) == 'N') {
                                                break;
                                            }
                                            //If base quality is less than conf.goodq, exit loop
                                            if (queryQuality.charAt(n + m + vi) - 33 < conf.goodq) {
                                                break;
                                            }
                                            //If reference sequence has base at this position and it matches read base, update offset
                                            if (ref.containsKey(start + vi)) {
                                                if (isNotEquals(querySequence.charAt(n + m + vi), ref.get(start + vi))) {
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
                                            ss += substr(querySequence, n + m, offset);
                                            q.append(substr(queryQuality, n + m, offset));
                                            //Increase coverage for positions corresponding to first offset bases of next segment
                                            for (int osi = 0; osi < offset; osi++) {
                                                incCnt(cov, start + osi, 1);
                                            }
                                        }
                                    }
                                }

                                // Offset should be reset to 0 on every loop, so it is non-zero only if previous part of
                                // code was executed
                                // Append '&' and $ss to s if next segment has good matching base
                                if (offset > 0) {
                                    s.append("&").append(ss);
                                }

                                //If start of the segment is within region of interest and the segment does not have unknown bases
                                if (start - 1 >= region.start && start - 1 <= region.end && !s.toString().contains("N")) {
                                    int inspos = start - 1;
                                    Matcher mm = BEGIN_ATGC_END.matcher(s);
                                    if (mm.find()) {
                                        Tuple.Tuple3<Integer, String, Integer> tpl = adjInsPos(start - 1, s.toString(), ref);
                                        inspos = tpl._1;
                                        s = new StringBuilder(tpl._2);
                                    }
                                    //add '+' + s to insertions at inspos
                                    incCnt(getOrElse(ins, inspos, new HashMap<>()), "+" + s, 1);

                                    //add insertion to table of variations
                                    Variation hv = getVariation(iHash, inspos, "+" + s); //variant structure for this insertion
                                    hv.incDir(dir);
                                    //add count
                                    hv.cnt++;
                                    //minimum of positions from start of read and end of read
                                    int tp = p < rlen1 - p ? p + 1 : rlen1 - p;

                                    //mean read quality of the segment
                                    double tmpq = 0;
                                    for (int i = 0; i < q.length(); i++) {
                                        tmpq += q.charAt(i) - 33;
                                    }
                                    tmpq = tmpq / q.length();

                                    //pstd is a flag that is 1 if the variant is covered by at least 2 read segments with different positions
                                    if (!hv.pstd && hv.pp != 0 && tp != hv.pp) {
                                        hv.pstd = true;
                                    }
                                    //qstd is a flag that is 1 if the variant is covered by at least 2 segment reads with different qualities
                                    if (!hv.qstd && hv.pq != 0 && tmpq != hv.pq) {
                                        hv.qstd = true;
                                    }
                                    hv.pmean += tp;
                                    hv.qmean += tmpq;
                                    hv.Qmean += mappingQuality;
                                    hv.pp = tp;
                                    hv.pq = tmpq;
                                    if (tmpq >= conf.goodq) {
                                        hv.hicnt++;
                                    } else {
                                        hv.locnt++;
                                    }
                                    hv.nm += nm - nmoff;

                                    // Adjust the reference count for insertion reads
                                    /*
                                    Condition:
                                    1). reference sequence has base for the position
                                    2). hash contains variant structure for the position
                                    3). read base at position n-1 matches reference at start-1
                                     */
                                    //#if (getVariationMaybe(hash, inspos, ref.get(start - 1)) != null
                                    // #&& isHasAndEquals(querySequence.charAt(n - 1 - (start - inspos)), ref, inspos)) {

                                    // subCnt(getVariation(hash, inspos, ref.get(inspos).toString()), dir, tp, tmpq,
                                    // Qmean, nm, conf);
                                    if (inspos > position) {
                                        Variation tv = getVariationMaybe(hash, inspos, querySequence.charAt(n - 1 - (start - 1 - inspos)));
                                        //Substract count.
                                        if (tv != null) {
                                            subCnt(tv, dir, tp, queryQuality.charAt(n - 1 - (start - 1 - inspos)) - 33,
                                                    mappingQuality, nm - nmoff, conf.goodq);
                                        }
                                    }
                                    // #}
                                    // Adjust count if the insertion is at the edge so that the AF won't > 1
                                    /*
                                    Condition:
                                    1). looking at second segment in CIGAR string
                                    2). first segment is a soft-clipping or a hard-clipping
                                    */
                                    if (ci == 1 && (cigar.getCigarElement(0).getOperator() == CigarOperator.S
                                            || cigar.getCigarElement(0).getOperator() == CigarOperator.H)) {
                                        //Add one more variant corresponding to base at start - 1 to hash
                                        Variation ttref = getVariation(hash, inspos, ref.get(inspos).toString());
                                        ttref.incDir(dir);
                                        ttref.cnt++;
                                        ttref.pstd = hv.pstd;
                                        ttref.qstd = hv.qstd;
                                        ttref.pmean += tp;
                                        ttref.qmean += tmpq;
                                        ttref.Qmean += mappingQuality;
                                        ttref.pp = tp;
                                        ttref.pq = tmpq;
                                        ttref.nm += nm - nmoff;
                                        incCnt(cov, inspos, 1);
                                    }
                                }

                                //adjust read position by m (CIGAR segment length) + offset + multoffp
                                n += m + offset + multoffp;
                                p += m + offset + multoffp;
                                //adjust reference position by offset + multoffs
                                start += offset + multoffs;
                                continue;
                            }
                            case D: { //deletion
                                offset = 0;
                                //description string of deleted segment
                                StringBuilder s = new StringBuilder("-").append(m);
                                //sequence to be appended if next segment is matched
                                StringBuilder ss = new StringBuilder();
                                //quality of last base before deletion
                                char q1 = queryQuality.charAt(n - 1);
                                //quality of this segment
                                StringBuilder q = new StringBuilder();

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
                                if (conf.performLocalRealignment && cigar.numCigarElements() > ci + 2
                                        && cigar.getCigarElement(ci + 1).getLength() <= conf.vext
                                        && cigar.getCigarElement(ci + 1).getOperator() == CigarOperator.M
                                        && (cigar.getCigarElement(ci + 2).getOperator() == CigarOperator.I
                                        || cigar.getCigarElement(ci + 2).getOperator() == CigarOperator.D)
                                        && cigar.getCigarElement(ci + 3).getOperator() != CigarOperator.I
                                        && cigar.getCigarElement(ci + 3).getOperator() != CigarOperator.D) {

                                    int mLen = cigar.getCigarElement(ci + 1).getLength();
                                    int indelLen = cigar.getCigarElement(ci + 2).getLength();

                                    int begin = n;
                                    appendSegments(querySequence, queryQuality, cigar, ci, s, q, mLen, indelLen, begin,
                                            false);

                                    //add length of next segment to both read and reference offsets
                                    //add length of next-next segment to reference position (for insertion) or to read position(for deletion)
                                    multoffs += mLen + (cigar.getCigarElement(ci + 2).getOperator() == CigarOperator.D ? indelLen : 0);
                                    multoffp += mLen + (cigar.getCigarElement(ci + 2).getOperator() == CigarOperator.I ? indelLen : 0);
                                    if (cigar.numCigarElements() > ci + 3
                                            && cigar.getCigarElement(ci + 3).getOperator() == CigarOperator.M) {
                                        int vsn = 0;
                                        int tn = n + multoffp;
                                        int ts = start + multoffs + m;
                                        for (int vi = 0; vsn <= conf.vext && vi < cigar.getCigarElement(ci + 3).getLength(); vi++) {
                                            if (querySequence.charAt(tn + vi) == 'N') {
                                                break;
                                            }
                                            if (queryQuality.charAt(tn + vi) - 33 < conf.goodq) {
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
                                            ss.append(substr(querySequence, tn, offset));
                                            q.append(substr(queryQuality, tn, offset));
                                        }
                                    }
                                    // skip next 2 CIGAR segments
                                    ci += 2;
                                } else if (conf.performLocalRealignment && cigar.numCigarElements() > ci + 1
                                        && cigar.getCigarElement(ci + 1).getOperator() == CigarOperator.I) {
                                    /*
                                    Condition:
                                    1). CIGAR string has next entry
                                    2). next CIGAR segment is an insertion
                                     */

                                    int insLen = cigar.getCigarElement(ci + 1).getLength();
                                    //Append '^' + next segment (inserted)
                                    s.append("^").append(substr(querySequence, n, insLen));
                                    //Append next segement to quality string
                                    q.append(substr(queryQuality, n, insLen));

                                    //Shift reference position by length of next segment
                                    //skip next CIGAR segment
                                    multoffp += insLen;
                                    if (cigar.numCigarElements() > ci + 2
                                            && cigar.getCigarElement(ci + 2).getOperator() == CigarOperator.M) {
                                        int mLen = cigar.getCigarElement(ci + 2).getLength();
                                        int vsn = 0;
                                        int tn = n + multoffp;
                                        int ts = start + m;
                                        for (int vi = 0; vsn <= conf.vext && vi < mLen; vi++) {
                                            char seqCh = querySequence.charAt(tn + vi);
                                            if (seqCh == 'N') {
                                                break;
                                            }
                                            if (queryQuality.charAt(tn + vi) - 33 < conf.goodq) {
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
                                            ss.append(substr(querySequence, tn, offset));
                                            q.append(substr(queryQuality, tn, offset));
                                        }
                                    }
                                    ci += 1;
                                } else {
                                    /*
                                    Condition:
                                    1). CIGAR string has next entry
                                    2). next CIGAR segment is matched
                                     */
                                    if (conf.performLocalRealignment && cigar.numCigarElements() > ci + 1
                                            && cigar.getCigarElement(ci + 1).getOperator() == CigarOperator.M) {
                                        int mLen = cigar.getCigarElement(ci + 1).getLength();
                                        int vsn = 0;
                                        //Loop over next CIGAR segment (no more than conf.vext bases ahead)
                                        for (int vi = 0; vsn <= conf.vext && vi < mLen; vi++) {
                                            char seqCh = querySequence.charAt(n + vi);
                                            //If base is unknown, exit loop
                                            if (seqCh == 'N') {
                                                break;
                                            }
                                            //If base quality is less than $GOODQ, exit loop
                                            if (queryQuality.charAt(n + vi) - 33 < conf.goodq) {
                                                break;
                                            }
                                            //If reference sequence has base at this position and it matches read base, update offset
                                            Character refCh = ref.get(start + m + vi);
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
                                            ss.append(substr(querySequence, n, offset));
                                            q.append(substr(queryQuality, n, offset));
                                        }
                                    }
                                }
                                //offset should be reset to 0 on every loop, so it is non-zero only if previous next CIGAR segment is matched
                                // Append '&' and ss to s if next segment has good matching base
                                if (offset > 0) {
                                    s.append("&").append(ss);
                                }

                                //quality of first matched base after deletion
                                //append best of $q1 and $q2
                                if (n + offset >= queryQuality.length()) {
                                    q.append(q1);
                                } else {
                                    char q2 = queryQuality.charAt(n + offset);
                                    q.append(q1 > q2 ? q1 : q2);
                                }

                                //If reference position is inside region of interest
                                if (start >= region.start && start <= region.end) {
                                    addVariationForDeletion(conf, hash, cov, dels5, mappingQuality,
                                            nm, p, dir, start, rlen1, m, s, q, nmoff);
                                }

                                //adjust reference position by offset + multoffs
                                start += m + offset + multoffs;

                                //adjust read position by m (CIGAR segment length) + offset + multoffp
                                n += offset + multoffp;
                                p += offset + multoffp;
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
                        for (int i = offset; i < m; i++) {
                            //flag to trim reads at opt_T bases from start or end (depending on direction)
                            boolean trim = false;
                            if (conf.trimBasesAfter != 0) {
                                if (!dir) {
                                    if (n > conf.trimBasesAfter) {
                                        trim = true;
                                    }
                                } else {
                                    if (rlen2 - n > conf.trimBasesAfter) {
                                        trim = true;
                                    }
                                }
                            }

                            //variation string. Initialize to first base of the read sequence
                            final char ch1 = querySequence.charAt(n);
                            String s = String.valueOf(ch1);
                            boolean startWithDeletion = false;
                            //skip if base is unknown
                            if (ch1 == 'N') {
                                if (conf.includeNInTotalDepth) {
                                    incCnt(cov, start, 1);
                                }
                                start++;
                                n++;
                                p++;
                                continue;
                            }

                            //sum of qualities for bases
                            double q = queryQuality.charAt(n) - 33;
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
                                    && (start + 1) <= region.end && (i + 1) < m
                                    && q >= conf.goodq
                                    && isHasAndNotEquals(ref, start, querySequence, n)
                                    && isNotEquals('N', ref.get(start))) {

                                //Break if base is unknown in the read
                                char nuc = querySequence.charAt(n + 1);
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
                                    q += queryQuality.charAt(n + 1) - 33;
                                    //increase number of bases
                                    qbases++;
                                    //shift read position by 1
                                    n++;
                                    p++;
                                    i++;
                                    //shift reference position by 1
                                    start++;
                                    nmoff++;
                                } else { // look ahead to see whether there're more mismatches, if yes, add reverence to grow MNV
                                    int ssn = 0;
                                    for (int ssi = 1; ssi <= conf.vext; ssi++) {
                                        if (i + 1 + ssi >= m){
                                            break;
                                        }
                                        if (n + 1 + ssi < querySequence.length()
                                                && isHasAndNotEquals(querySequence.charAt(n + 1 + ssi), ref, start + 1 + ssi)) {
                                            ssn = ssi + 1;
                                            break;
                                        }
                                    }
                                    if (ssn == 0) {
                                        break;
                                    }
                                    for (int ssi = 1; ssi <= ssn; ssi++) {
                                        ss.append(querySequence.charAt(n + ssi));
                                        q += queryQuality.charAt(n + ssi) - 33;
                                        qbases++;
                                    }
                                    n += ssn;
                                    p += ssn;
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
                            if (conf.performLocalRealignment && m - i <= conf.vext
                                    && cigar.numCigarElements() > ci + 1
                                    && cigar.getCigarElement(ci + 1).getOperator() == CigarOperator.D
                                    && ref.containsKey(start)
                                    && (ss.length() > 0 || isNotEquals(querySequence.charAt(n), ref.get(start)))
                                    && queryQuality.charAt(n) - 33 > conf.goodq) {

                                //loop until end of CIGAR segments
                                while (i + 1 < m) {
                                    //append next base to s and add its quality to q
                                    s += querySequence.charAt(n + 1);
                                    q += queryQuality.charAt(n + 1) - 33;
                                    //increase number of bases
                                    qbases++;

                                    //shift read and reference indices by 1
                                    i++;
                                    n++;
                                    p++;
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
                                if (cigar.numCigarElements() > ci + 1
                                        && cigar.getCigarElement(ci + 1).getOperator() == CigarOperator.I) {
                                    //append '^' + next-next segment sequence
                                    s += "^" + substr(querySequence, n + 1, cigar.getCigarElement(ci + 1).getLength());

                                    //Loop over next-next segment
                                    int nextLen = cigar.getCigarElement(ci + 1).getLength();
                                    for (int qi = 1; qi <= nextLen; qi++) {
                                        //add base quality to total quality
                                        q += queryQuality.charAt(n + 1 + qi) - 33;
                                        //increase number of insertion bases
                                        qibases++;
                                    }

                                    //adjust read position by length of next-next segment
                                    n += nextLen;
                                    p += nextLen;
                                    ci += 1;
                                }
                                if (cigar.numCigarElements() > ci + 1
                                        && cigar.getCigarElement(ci + 1).getOperator() == CigarOperator.M) {
                                    Tuple.Tuple4<Integer, String, String, Integer> tpl = findOffset(start + ddlen + 1,
                                            n + 1, cigar.getCigarElement(ci + 1).getLength(), querySequence,
                                            queryQuality, ref, cov, conf.vext, conf.goodq);
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
                            if (!trim) {
                                //If start - qbases + 1 is in region of interest
                                final int pos = start - qbases + 1;
                                if (pos >= region.start && pos <= region.end) {
                                    addVariationForMatchingPart(conf, hash, cov, mnp, dels5, mappingQuality, nm, p,
                                            dir, start, rlen1, nmoff, s, startWithDeletion, q, qbases, qibases, ddlen, pos);
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
                                n++;
                                p++;
                            }
                            // Skip read if it is overlap
                            if (skipOverlappingReads(conf, record, position, dir, start, mateAlignmentStart)) {
                                break processCigar;
                            }
                        }
                        if (moffset != 0) {
                            offset = moffset;
                            n += moffset;
                            start += moffset;
                            p += moffset;
                        }
                        if (start > region.end) { //end if reference position is outside region of interest
                            break;
                        }
                    }
                }
            }
        }

        if (conf.y) {
            System.err.println("TIME: Finish parsing SAM: " + LocalDateTime.now() + " " + bam + " "
                    + chr + ":" + region.start + "-" + region.end);
        }

        if (conf.outputSplicing) {
            for (Map.Entry<String, int[]> entry : spliceCnt.entrySet()) {
                System.out.printf("%s\t%s\t%s\t%s\n", sample, region.chr, entry.getKey(), entry.getValue()[0]);
            }
            return null;
        }

        if (svflag) {
            return tuple(hash, iHash, cov, rlen, (disc + 1)/(totalreads - dupreads + 1) > 0.5 ? 0.0 : 1.0);
        }

        if (conf.y) {
            System.err.println("TIME: Start realign: " + LocalDateTime.now());
        }

        if(!conf.disableSV) {
            filterAllSVStructures(rlen, conf, svStructures);
        }
        adjMNP(hash, mnp, cov, ref, sclip3, sclip5, conf);

        if (conf.performLocalRealignment) {
            VariationRealigner.realignIndels(region, chrs, rlen, reference, conf, bam, hash, iHash, cov, sclip3,
                    sclip5, sample, splice, ampliconBasedCalling, ins, dels5,
                    svStructures);
        }

        if(!conf.disableSV) {
            findAllSVs(region, bam, chrs, sample, splice, ampliconBasedCalling, rlen,
                    reference, conf, hash, iHash, cov, sclip3, sclip5, svStructures);
        }
        if (conf.y) {
            outputClipping(sclip5, sclip3, conf);
            System.err.println("TIME: Finish realign:" + LocalDateTime.now());
        }

        return tuple(hash, iHash, cov, rlen,
                (conf.removeDuplicatedReads && totalreads != 0) ? Double.parseDouble(String.format("%.3f", ((double) dupreads)/totalreads)) : 0.0);
    }

    private static void addVariationForMatchingPart(Configuration conf,
                                                    Map<Integer, VariationMap<String, Variation>> hash,
                                                    Map<Integer, Integer> cov,
                                                    Map<Integer, Map<String, Integer>> mnp,
                                                    Map<Integer, Map<String, Integer>> dels5,
                                                    int mappingQuality,
                                                    int nm, int p, boolean dir, int start,
                                                    int rlen1, int nmoff, String s,
                                                    boolean startWithDelition, double q,
                                                    int qbases, int qibases, int ddlen, int pos) {
        //add variation record for $s
        Variation hv = getVariation(hash, pos, s); //reference to variant structure
        hv.incDir(dir);

        if(isBEGIN_ATGC_AMP_ATGCs_END(s)) {
            //if s is one base followed by '&' and one or more bases
            //add variant record for s to mnp
            increment(mnp, pos, s);
        }

        //increment count
        hv.cnt++;

        //minimum of positions from start of read and end of read
        int tp = p < rlen1 - p ? p + 1 : rlen1 - p;

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
        hv.pmean += tp;
        hv.qmean += q;
        hv.Qmean += mappingQuality;
        hv.pp = tp;
        hv.pq = q;
        hv.nm += nm - nmoff;
        if (q >= conf.goodq) {
            hv.hicnt++;
        } else {
            hv.locnt++;
        }

        //increase coverage for bases covered by the variation
        for (int qi = 1; qi <= qbases; qi++) {
            incCnt(cov, start - qi + 1, 1);
        }

        //If variation starts with a deletion ('-' character)
        if (startWithDelition) {
            //add variation to deletions map
            increment(dels5, pos, s);

            //increase coverage for next CIGAR segment
            for (int qi = 1; qi < ddlen; qi++) {
                incCnt(cov, start + qi, 1);
            }
        }
    }

    private static void addVariationForDeletion(Configuration conf,
                                                Map<Integer, VariationMap<String, Variation>> hash,
                                                Map<Integer, Integer> cov,
                                                Map<Integer, Map<String, Integer>> dels5,
                                                int mappingQuality, int nm, int p,
                                                boolean dir, int start, int rlen1,
                                                int m, StringBuilder s, StringBuilder q, int nmoff) {
        //add variant structure for deletion at this position
        Variation hv = getVariation(hash, start, s.toString()); //variation structure
        //add record for deletion in deletions map
        increment(dels5, start, s.toString());
        hv.incDir(dir);
        //increase count
        hv.cnt++;

        //minimum of positions from start of read and end of read
        int tp = p < rlen1 - p ? p + 1 : rlen1 - p;

        //average quality of bases
        double tmpq = 0;

        for (int i = 0; i < q.length(); i++) {
            tmpq += q.charAt(i) - 33;
        }
        tmpq = tmpq / q.length();

        //pstd is a flag that is 1 if the variant is covered by at least 2 read segments with different positions
        if (!hv.pstd && hv.pp != 0 && tp != hv.pp) {
            hv.pstd = true;
        }

        //qstd is a flag that is 1 if the variant is covered by at least 2 segment reads with different qualities
        if (!hv.qstd && hv.pq != 0 && tmpq != hv.pq) {
            hv.qstd = true;
        }
        hv.pmean += tp;
        hv.qmean += tmpq;
        hv.Qmean += mappingQuality;
        hv.pp = tp;
        hv.pq = tmpq;
        hv.nm += nm - nmoff;
        if (tmpq >= conf.goodq) {
            hv.hicnt++;
        } else {
            hv.locnt++;
        }

        //increase coverage count for reference bases missing from the read
        for (int i = 0; i < m; i++) {
            incCnt(cov, start + i, 1);
        }
    }

    private static boolean isReadChimericWithSA(SAMRecord record, int position, String saTagString, int rlen,
                                                boolean dir, Configuration conf, boolean is5Side){
        String[] saTagArray = saTagString.split(",");
        String saChromosome = saTagArray[0];
        int saPosition = Integer.valueOf(saTagArray[1]);
        String saDirectionString = saTagArray[2];
        String saCigar = saTagArray[3];
        boolean saDirectionIsForward = saDirectionString.equals("+");
        Matcher mm = SA_CIGAR_D_S_3clip.matcher(saCigar);;

        if (is5Side) {
            mm = SA_CIGAR_D_S_5clip.matcher(saCigar);
        }

        boolean isChimericWithSA = ((dir && saDirectionIsForward) || (!dir && !saDirectionIsForward))
                && saChromosome.equals(record.getReferenceName())
                && (abs(saPosition - position) < 2 * rlen)
                && mm.find();

        if (conf.y && isChimericWithSA) {
            System.err.println(record.getReadName() + " " + record.getReferenceName()
                    + " " + position + " " + record.getMappingQuality()
                    + " " + record.getCigarString()
                    + " is ignored as chimeric with SA: " +
                    saPosition + "," + saDirectionString + "," + saCigar);
        }

        return isChimericWithSA;
    }

    private static void appendSegments(String querySequence, String queryQuality, Cigar cigar, int ci,
                                       StringBuilder s, StringBuilder q, int mLen, int indelLen, int begin,
                                       boolean isInsertion) {

        //begin is n + m for insertion and n for deletion
        //append to s '#' and part of read sequence corresponding to next CIGAR segment (matched one)
        s.append("#").append(substr(querySequence, begin, mLen));
        //append quality string of next matched segment from read
        q.append(substr(queryQuality, begin, mLen));

        //if an insertion is two segments ahead, append '^' + part of sequence corresponding
        // to next-next segment otherwise (deletion) append '^' + length of a next-next segment
        s.append('^').append(cigar.getCigarElement(ci + 2).getOperator() == CigarOperator.I
                ? substr(querySequence, begin + mLen, indelLen)
                : indelLen);
        //if an insertion is two segments ahead, append part of quality string sequence
        // corresponding to next-next segment otherwise (deletion)
        // append first quality score of next segment or return empty string
        if (isInsertion) {
            q.append(cigar.getCigarElement(ci + 2).getOperator() == CigarOperator.I
                    ? substr(queryQuality, begin + mLen, indelLen)
                    : queryQuality.charAt(begin + mLen));
        } else {
            q.append(cigar.getCigarElement(ci + 2).getOperator() == CigarOperator.I
                    ? substr(queryQuality, begin + mLen, indelLen)
                    : "");
        }
    }

    private static boolean skipOverlappingReads(Configuration conf, SAMRecord record, int position, boolean dir, int start, int mateAlignmentStart) {
        if (conf.uniqueModeAlignmentEnabled && isPairedAndSameChromosome(record)
                && !dir && start >= mateAlignmentStart) {
            return true;
        }
        if (conf.uniqueModeSecondInPairEnabled && record.getSecondOfPairFlag() && isPairedAndSameChromosome(record)
                && isReadsOverlap(record, start, position, mateAlignmentStart)) {
            return true;
        }
        return false;
    }

    private static boolean isReadsOverlap(SAMRecord record, int start, int position, int mateAlignmentStart){
        if (position >= mateAlignmentStart) {
            return start >= mateAlignmentStart
                    && start <= (mateAlignmentStart + record.getCigar().getReferenceLength() - 1);
        }
        else {
            return start >= mateAlignmentStart
                    && record.getMateAlignmentStart() <= record.getAlignmentEnd();
        }
    }

    private static boolean isPairedAndSameChromosome(SAMRecord record) {
        return record.getReadPairedFlag() && getMrnm(record).equals("=");
    }

    private static void sclip3HighQualityProcessing(Region region, Configuration conf, Map<Integer, Sclip> sclip3,
                                                    String querySequence, int mappingQuality, String queryQuality,
                                                    int nm, int n, boolean dir, int start, int m, int q,
                                                    int qn, int lowqcnt) {
        //If we have at least 1 high-quality soft-clipped base of region of interest
        if (qn >= 1 && qn > lowqcnt && start >= region.start && start <= region.end) {
            //add record to $sclip3
            Sclip sclip = sclip3.get(start);
            if (sclip == null) {
                sclip = new Sclip();
                sclip3.put(start, sclip);
            }
            for (int si = 0; si < qn; si++) {
                Character ch = querySequence.charAt(n + si);
                int idx = si;
                Map<Character, Integer> cnts = sclip.nt.get(idx);
                if (cnts == null) {
                    cnts = new HashMap<>();
                    sclip.nt.put(idx, cnts);

                }
                incCnt(cnts, ch, 1);
                Variation variation = getVariationFromSeq(sclip, idx, ch);
                addCnt(variation, dir, qn - si, queryQuality.charAt(n + si) - 33, mappingQuality, nm, conf.goodq);
            }
            addCnt(sclip, dir, m, q / (double) qn, mappingQuality, nm, conf.goodq);
        }
    }

    private static void sclip5HighQualityProcessing(Region region, Configuration conf, Map<Integer, Sclip> sclip5,
                                                    String querySequence, int mappingQuality, String queryQuality,
                                                    int nm, boolean dir, int start, int m, int q,
                                                    int qn, int lowqcnt) {
        //If we have at least 1 high-quality soft-clipped base of region of interest
        if (qn >= 1 && qn > lowqcnt && start >= region.start && start <= region.end) {
            //add record to $sclip5
            Sclip sclip = sclip5.get(start);
            if (sclip == null) {
                sclip = new Sclip();
                sclip5.put(start, sclip);
            }
            for (int si = m - 1; m - si <= qn; si--) {
                Character ch = querySequence.charAt(si);
                int idx = m - 1 - si;
                Map<Character, Integer> cnts = sclip.nt.get(idx);
                if (cnts == null) {
                    cnts = new HashMap<>();
                    sclip.nt.put(idx, cnts);
                }
                incCnt(cnts, ch, 1);
                Variation seqVariation = getVariationFromSeq(sclip, idx, ch);
                addCnt(seqVariation, dir, si - (m - qn), queryQuality.charAt(si) - 33, mappingQuality, nm, conf.goodq);
            }
            addCnt(sclip, dir, m, q / (double) qn, mappingQuality, nm, conf.goodq);
        }
    }


    /**
     * Adjust MNP when there're breakpoints within MNP (multi-nucleotide polymorphism)
     * @param hash map of variants produced by parseSAM method
     * @param mnp map of MNP variants
     * @param cov coverage
     * @param ref map of reference bases
     * @param sclip5 map of 5' softclips
     * @param sclip3 map of 3' softclips
     * @param conf configuration
     */
    public static void adjMNP(Map<Integer, VariationMap<String, Variation>> hash,
                       Map<Integer, Map<String, Integer>> mnp,
                       Map<Integer, Integer> cov, Map<Integer, Character> ref,
                       Map<Integer, Sclip> sclip3, Map<Integer, Sclip> sclip5, Configuration conf) {

        for (Map.Entry<Integer, Map<String, Integer>> entry : mnp.entrySet()) {
            final Integer p = entry.getKey();
            Map<String, Integer> v = entry.getValue();

            for (Map.Entry<String, Integer> en : v.entrySet()) {
                final String vn = en.getKey();
                final Map<String, Variation> hashP = hash.get(p);
                if (hashP == null) {
                    continue;
                }
                final Variation vref = hashP.get(vn);
                if (vref == null ) { // The variant is likely already been used by indel realignment
                    continue;
                }
                if (conf.y) {
                    System.err.printf("  AdjMnt: %d %s %d\n", p, vn, vref.cnt);
                }

                final String mnt = vn.replaceFirst("&", "");
                for (int i = 0; i < mnt.length() - 1; i++) {
                    String left = substr(mnt, 0, i + 1);
                    if (left.length() > 1) {
                        StringBuilder sb = new StringBuilder(left);
                        sb.insert(1, "&");
                        left = sb.toString();
                    }

                    String right = substr(mnt, -(mnt.length() - i - 1));
                    if (right.length() > 1) {
                        StringBuilder sb = new StringBuilder(right);
                        sb.insert(1, "&");
                        right = sb.toString();
                    }
                    {
                        Variation tref = hashP.get(left);
                        if (tref != null) {
                            if (tref.cnt <= 0) {
                                continue;
                            }
                            if (tref.cnt < vref.cnt && tref.pmean / tref.cnt <= i + 1) {
                                if (conf.y) {
                                    System.err.printf("    AdjMnt Left: %s %s Left: %s Cnt: %s\n",p, vn, left, tref.cnt);
                                }
                                adjCnt(vref, tref, conf);
                                hashP.remove(left);
                            }
                        }
                    }
                    if (hash.containsKey(p + i + 1)) {
                        Variation tref = hash.get(p + i + 1).get(right);
                        if (tref != null) {
                            if (tref.cnt < 0) {
                                continue;
                            }
                            // #&& tref.pmean / tref.cnt <= mnt.length() - i - 1)
                            if (tref.cnt < vref.cnt) {
                                if (conf.y) {
                                    System.err.printf("    AdjMnt Right: %s %s Right: %s Cnt: %s\n", p, vn, right, tref.cnt);
                                }
                                adjCnt(vref, tref, conf);
                                incCnt(cov, p, tref.cnt);
                                hash.get(p + i + 1).remove(right);
                            }
                        }
                    }
                }
                if (sclip3.containsKey(p)) {
                    final Sclip sc3v = sclip3.get(p);
                    if (!sc3v.used) {
                        final String seq = findconseq(sc3v, conf, 0);
                        if (seq.startsWith(mnt)) {
                            if(seq.length() == mnt.length()
                                    || ismatchref(seq.substring(mnt.length()), ref, p + mnt.length(), 1, conf.y)) {
                                adjCnt(hash.get(p).get(vn), sc3v, conf);
                                incCnt(cov, p, sc3v.cnt);
                                sc3v.used = true;
                            }
                        }
                    }
                }
                if (sclip5.containsKey(p + mnt.length())) {
                    final Sclip sc5v = sclip5.get(p + mnt.length());
                    if (!sc5v.used) {
                        String seq =  findconseq(sc5v, conf, 0);
                        if (!seq.isEmpty() && seq.length() >= mnt.length()) {
                            seq =  new StringBuffer(seq).reverse().toString();
                            if (seq.endsWith(mnt)) {
                                if (seq.length() == mnt.length()
                                        || ismatchref(seq.substring(0, seq.length() - mnt.length()), ref, p - 1, -1, conf.y)) {
                                    adjCnt(hash.get(p).get(vn), sc5v, conf);
                                    incCnt(cov, p, sc5v.cnt);
                                    sc5v.used = true;
                                }
                            }
                        }

                    }
                }
            }
        }
    }

    public static boolean ismatchref(String seq, Map<Integer, Character> ref, int p, int dir, boolean debugLog) {
        int MM = 3;
        return ismatchref(seq, ref, p, dir, debugLog, MM);
    }

    /**
     * Utility method for adjustMNP method
     * @param sequence subsequence consensus sequence in soft-clipped reads  $seq
     * @param position key for MNP map     $p
     */
    public static boolean ismatchref(String sequence,
                              Map<Integer, Character> ref,
                              int position,
                              int dir,
                              boolean debugLog,
                              int MM) {
        if (debugLog) {
            System.err.println(format("      Matching REF %s %s %s %s", sequence, position, dir, MM));
        }

        int mm = 0;
        for (int n = 0; n < sequence.length(); n++) {
            final Character refCh = ref.get(position + dir * n);
            if (refCh == null) {
                return false;
            }
            if (charAt(sequence, dir == 1 ? n : dir * n - 1) != refCh) {
                mm++;
            }
        }
        return mm <= MM && mm / (double)sequence.length() < 0.15;
    }

    public static String getMrnm(SAMRecord record) {
        if (record.getMateReferenceName() == null) {
            return "*";
        }

        if (record.getReferenceName().equals(record.getMateReferenceName())) {
            return "=";
        }
        return record.getMateReferenceName();
    }

}
