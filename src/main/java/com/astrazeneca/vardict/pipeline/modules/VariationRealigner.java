package com.astrazeneca.vardict.pipeline.modules;

import com.astrazeneca.utils.Utils;
import com.astrazeneca.vardict.Region;
import com.astrazeneca.vardict.pipeline.data.SamView;
import com.astrazeneca.vardict.pipeline.data.RealignedVariationData;
import com.astrazeneca.vardict.pipeline.data.Scope;
import com.astrazeneca.vardict.pipeline.data.VariationData;
import com.astrazeneca.vardict.variations.SoftClip;
import com.astrazeneca.vardict.variations.Variation;
import htsjdk.samtools.SAMRecord;
import com.astrazeneca.utils.Tuple;

import java.io.IOException;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static com.astrazeneca.GlobalReadOnlyScope.*;
import static com.astrazeneca.vardict.pipeline.modules.CigarUtils.getAlignedLenght;
import static com.astrazeneca.vardict.pipeline.modules.PatternsUtils.AMP_ATGC;
import static com.astrazeneca.vardict.pipeline.modules.PatternsUtils.CARET_ATGNC;
import static com.astrazeneca.vardict.pipeline.modules.RecordPreprocessor.getChrName;
import static com.astrazeneca.vardict.variations.VariationUtils.*;
import static java.util.Collections.singletonMap;
import static com.astrazeneca.utils.Tuple.tuple;
import static com.astrazeneca.utils.Utils.charAt;
import static com.astrazeneca.utils.Utils.substr;
import static com.astrazeneca.utils.Utils.toInt;

/**
* Class for the realignment process for variations (insertions and deletions).
* Also realigns lager variations.
*/
public class VariationRealigner implements Module<VariationData, RealignedVariationData> {

    //Longest continued mismatches typical aligned at the end
    private static final int MAX_MISMATCH_ALLOWED = 3;

    private final static Pattern B_A7 = Pattern.compile("^.AAAAAAA");
    private final static Pattern B_T7 = Pattern.compile("^.TTTTTTT");
    private final static Pattern B_A8 = Pattern.compile("^.AAAAAAAA");
    private final static Pattern B_T8 = Pattern.compile("^.TTTTTTTT");

    private static final Pattern BEGIN_PLUS_ATGC = Pattern.compile("^\\+([ATGC]+)");
    private static final Pattern MINUS_NUMBER_AMP_ATGCs_END = Pattern.compile("(-\\d+)&[ATGC]+$");
    private static final Pattern HASH_ATGC = Pattern.compile("#([ATGC]+)");
    private static final Pattern BEGIN_MINUS_NUMBER = Pattern.compile("^-(\\d+)");
    private static final Pattern BEGIN_MINUS_NUMBER_ANY = Pattern.compile("^-\\d+(.*)");
    private static final Pattern UP_NUMBER_END = Pattern.compile("\\^(\\d+)$");
    private static final Pattern ATGSs_AMP_ATGSs_END = Pattern.compile("(\\+[ATGC]+)&[ATGC]+$");

    private static final Comparator<Object[]> REALIGNDEL_COMPARATOR = (o1, o2) -> {
        int x1 = (Integer)o1[2];
        int x2 = (Integer)o2[2];
        int f = Integer.compare(x2, x1);
        if (f != 0)
            return f;
        x1 = (Integer)o1[0];
        x2 = (Integer)o2[0];
        f =  Integer.compare(x1, x2);
        if (f != 0)
            return f;
        String s1 = (String)o1[1];
        String s2 = (String)o2[1];
        return s2.compareTo(s1);

    };

    private static final Comparator<Map.Entry<Integer, SoftClip>> COMP2 = (o1, o2) -> {
        int f = Integer.compare(o2.getValue().varsCount, o1.getValue().varsCount);
        if (f != 0)
            return f;
        return Integer.compare(o1.getKey(), o2.getKey());
    };

    private static final Comparator<Tuple.Tuple3<Integer, SoftClip, Integer>> COMP3 = (o1, o2) -> Integer.compare(o2._3, o1._3);


    private Map<Integer, Map<String, Variation>> nonInsertionVariants;
    private Map<Integer, Map<String, Variation>> insertionVariants;
    private Map<Integer, Map<String, Integer>> positionToInsertionCount;
    private Map<Integer, Map<String, Integer>> positionToDeletionsCount;
    private Map<Integer, Integer> refCoverage;
    private Map<Integer, SoftClip> softClips5End;
    private Map<Integer, SoftClip> softClips3End;
    private Map<Integer, Character> reference;
    private String chr;
    private int maxReadLength;
    private String[] bams;
    private Map<Integer, Map<String, Integer>> mnp;


    @Override
    public Scope<RealignedVariationData> process(Scope<VariationData> scope) {

        initFromScope(scope);

        if (instance().conf.performLocalRealignment) {
            try {
                realign();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }

        SAMFileParser.adjustMNP(nonInsertionVariants, mnp, refCoverage, reference, softClips3End, softClips5End);

        return new Scope<>(
                scope,
                new RealignedVariationData(nonInsertionVariants, insertionVariants, refCoverage, maxReadLength));
    }

    private void initFromScope(Scope<VariationData> scope) {
        this.nonInsertionVariants = scope.data.nonInsertionVariants;
        this.insertionVariants = scope.data.insertionVariants;
        this.positionToInsertionCount = scope.data.positionToInsertionCount;
        this.positionToDeletionsCount = scope.data.positionToDeletionsCount;
        this.refCoverage = scope.data.refCoverage;
        this.softClips5End = scope.data.softClips5End;
        this.softClips3End = scope.data.softClips3End;
        this.reference = scope.regionRef;
        this.chr = getChrName(scope.region);
        this.maxReadLength = scope.maxReadLength;
        this.bams = scope.bam.split(":");
        this.mnp = scope.data.mnp;
    }

    /**
    * The main public method for interaction with this class, this method starts realign logic for a given variation.
    * @throws IOException
    */
    public void realign() throws IOException {
        if (instance().conf.y)
            System.err.println("Start realignDeletions");
        realignDeletions();
        if (instance().conf.y)
            System.err.println("Start realignInsertions");
        realignInsertions();
        if (instance().conf.y)
            System.err.println("Start realignLagerDeletions");
        realignLagerDeletions();
        if (instance().conf.y)
            System.err.println("Start realignLagerInsertions");
        realignLagerInsertions();
        if (instance().conf.y)
            System.err.println("Start realignInsertionsLagerThan30");
        realignInsertionsLagerThan30();
    }

    /**
     * Realign deletions if already present in alignment, uses positionToDeletionsCount
     * @throws IOException
     */
    //perl version: 2383
    void realignDeletions() throws IOException {
        realignDeletions(positionToDeletionsCount);
    }

    /**
     * Realign deletions if already present in alignment
     * @param positionToIndelsCount - map with position as a key and map with information about deletions of this position
     *                        and count of this position as a value
     * @throws IOException
     */
    //perl version: 2383
    void realignDeletions(Map<Integer, Map<String, Integer>> positionToIndelsCount) throws IOException {
        List<Object[]> listPositionToIndelsCount = makeListWithPositionToVariationCountInfo(positionToIndelsCount);
        //perl version: 2396
        for (Object[] positionToIndelCount : listPositionToIndelsCount) {
            Integer position = (Integer)positionToIndelCount[0];
            String indel = (String)positionToIndelCount[1];
            Integer indelCount = (Integer)positionToIndelCount[2];

            int deletionLen = 0;
            Matcher mtch = BEGIN_MINUS_NUMBER.matcher(indel);
            if (mtch.find()) {
                deletionLen = toInt(mtch.group(1));
            }
            mtch = UP_NUMBER_END.matcher(indel);
            if (mtch.find()) {
                deletionLen += toInt(mtch.group(1));
            }
            String extrains = "";
            mtch = CARET_ATGNC.matcher(indel);
            if (mtch.find()) {
                extrains = mtch.group(1);
            }
            String extra = "";
            mtch = BEGIN_MINUS_NUMBER_ANY.matcher(indel);
            if (mtch.find()) {
                extra = mtch.group(1).replaceAll("\\^|&|#", "");
            }

            Variation vref = getVariation(nonInsertionVariants, position, indel);

            processRealignment(position, indelCount, deletionLen,  0, deletionLen + extra.length() - extrains.length() - 1, extra, vref, true);

            int endPosition = position + deletionLen + extra.length() - extrains.length();
            Variation h = getVariationMaybe(nonInsertionVariants, position, reference.get(position));
            if (bams != null && bams.length > 0
                    && endPosition - position >= 5
                    && endPosition - position < maxReadLength - 10
                    && h != null && h.varsCount != 0
                    && noPassingReads(chr, position, endPosition, bams)
                    && vref.varsCount > 2 * h.varsCount * (1 - (endPosition - position) / (double) maxReadLength)) {

                adjCnt(vref, h, h);
            }

        }
        updateStatistics(listPositionToIndelsCount, true);
    }

    private void processRealignment(Integer position, Integer indelCount, int indelLen, int indelLen2,
                                    int shift, String extra, Variation vref, boolean isDeletion) {
    /*
    5' flanking region is a region of DNA that is adjacent to the 5' end of the gene.
    The 5' flanking region contains the promoter, and may contain enhancers or other protein binding sites.
    It is the region of DNA that is not transcribed into RNA.

    3' flanking region is a region of DNA which is NOT copied into the mature mRNA, but which is present adjacent to 3' end of the gene.
    It was originally thought that the 3' flanking DNA was not transcribed at all, but it was discovered to be transcribed into RNA, but
    quickly removed during processing of the primary transcript to form the mature mRNA. The 3' flanking region often contains sequences
    which affect the formation of the 3' end of the message. It may also contain enhancers or other sites to which proteins may bind.

    e.g.:
     flankingSeq5End                flankingSeq3End
     ---------------                ---------------
     ACGATCCAGTACGTAACGACGAGCAGCGCGTGACTGATCGTAGCTA
     */
        int flankingSeqStart5End = get5EndFlankingSeqStart(position, isDeletion ? indelLen : indelLen2 - 1);
        int flankingSeqEnd3End = get3EndFlankingSeqEnd(position, isDeletion ? 2 * indelLen : indelLen2);

        String flankingSeq5End = joinRef(reference, flankingSeqStart5End, position + (isDeletion ? - 1 : 0)) + extra;
        String flankingSeq3End = extra + joinRef(reference, position + shift + 1, flankingSeqEnd3End);

        MismatchResult mismatchesOn3End = findMismatchesOn35End(reference, position + (isDeletion ? 0 : 1), flankingSeq3End, softClips3End, true);
        MismatchResult mismatchesOn5End = findMismatchesOn35End(reference, position + shift, flankingSeq5End, softClips5End, false);

        processIndels(position, indelCount, vref, indelLen, mismatchesOn3End, mismatchesOn5End, isDeletion);

        if (mismatchesOn3End.mismatchPosition != 0 && mismatchesOn3End.mismatches.size() == 1 && nonInsertionVariants.containsKey(mismatchesOn3End.mismatchPosition)) {
            nonInsertionVariants.get(mismatchesOn3End.mismatchPosition).remove(mismatchesOn3End.mismatchVariant);
        }
        if (mismatchesOn5End.mismatchPosition != 0 && mismatchesOn5End.mismatches.size() == 1 && nonInsertionVariants.containsKey(mismatchesOn5End.mismatchPosition)) {
            nonInsertionVariants.get(mismatchesOn5End.mismatchPosition).remove(mismatchesOn5End.mismatchVariant);
        }

        processIndelsSoftClip(softClips5End, position, vref, flankingSeq5End, mismatchesOn5End.softClipPositions, false, isDeletion);
        processIndelsSoftClip(softClips3End, position, vref, flankingSeq3End, mismatchesOn3End.softClipPositions, true, isDeletion);
    }

    private void updateStatistics(List<Object[]> listPositionToIndelsCount, boolean isDeletion) {
        for (int i = 0; i < listPositionToIndelsCount.size(); i++) {
            int j = isDeletion ? listPositionToIndelsCount.size() - 1 - i : i;
            Object[] positionToIndelsCount = listPositionToIndelsCount.get(j);
            Integer position = (Integer)positionToIndelsCount[0];
            String variation = (String)positionToIndelsCount[1];
            Map<Integer, Map<String, Variation>> variantsToProcess = isDeletion ? nonInsertionVariants : insertionVariants;
            if (!variantsToProcess.containsKey(position)) {
                continue;
            }
            Variation vref = variantsToProcess.get(position).get(variation);
            if (vref == null) {
                continue;
            }
            Matcher matcher = (isDeletion ? MINUS_NUMBER_AMP_ATGCs_END : ATGSs_AMP_ATGSs_END).matcher(variation);
            if (matcher.find()) {
                String tn = matcher.group(1);
                Variation tref = variantsToProcess.get(position).get(tn);
                if (tref != null) {
                    if (vref.varsCount < tref.varsCount) {
                        adjCnt(tref, vref, isDeletion ? null : getVariationMaybe(nonInsertionVariants, position, reference.get(position)));
                        variantsToProcess.get(position).remove(variation);
                    }
                }
            }
        }
    }

    /**
     * Processes insertion or deletion variants
     * @param position start position
     * @param count indel count
     * @param vref variant
     * @param indelLength indel length
     * @param mismatchesOn3End number of nucleotides to adjust
     * @param mismatchesOn5End number of nucleotides to adjust
     * @param isDeletion flag for deletion processing
     */
    //perl version: 2680
    private void processIndels(Integer position, Integer count, Variation vref, int indelLength,
                               MismatchResult mismatchesOn3End, MismatchResult mismatchesOn5End, boolean isDeletion) {
        // TODO logging should be added
        List<Tuple.Tuple3<String, Integer, Integer>> listOfMismatches = new ArrayList<>(mismatchesOn3End.mismatches);
        listOfMismatches.addAll(mismatchesOn5End.mismatches);

        for (Tuple.Tuple3<String, Integer, Integer> tuple : listOfMismatches) {
            //mismatches nucleotide
            String mismatch = tuple._1;
            //start position of clip that contains mismatches
            Integer mismatchPosition = tuple._2;
            //end (3 or 5)
            Integer mismatchEndType = tuple._3;
            if (mismatch.length() > 1) {
                mismatch = mismatch.charAt(0) + "&" + mismatch.substring(1);
            }
            if (!nonInsertionVariants.containsKey(mismatchPosition)) {
                continue;
            }
            Variation tv = nonInsertionVariants.get(mismatchPosition).get(mismatch);
            if (!isRelevantVariant(count, indelLength, mismatchesOn3End.numberOfMismatches, mismatchesOn5End.numberOfMismatches, mismatchEndType, tv)) {
                continue;
            }
            // Adjust reference varsCount so that AF won't > 1
            if (mismatchPosition > position && mismatchEndType == 5) {
                double f = (isDeletion && tv.meanPosition != 0) ? (mismatchPosition - position) / (tv.meanPosition / (double)tv.varsCount) : 1;
                if (f > 1) {
                    f = 1;
                }
                incCnt(refCoverage, position, (int)(tv.varsCount * f));
                if (isDeletion) {
                    adjRefCnt(tv, getVariationMaybe(nonInsertionVariants, position, reference.get(position)), indelLength);
                }
            }

            Variation lref = null;
            if (mismatchPosition > position && mismatchEndType == 3 && nonInsertionVariants.containsKey(position) && reference.containsKey(position)
                    && nonInsertionVariants.get(position).containsKey(reference.get(position).toString())) {
                lref = nonInsertionVariants.get(position).get(reference.get(position).toString());
            }

            adjCnt(vref, tv, lref);
            nonInsertionVariants.get(mismatchPosition).remove(mismatch);
            if (nonInsertionVariants.get(mismatchPosition).isEmpty()) {
                nonInsertionVariants.remove(mismatchPosition);
            }
        }
    }

    /**
     * Realign insertion (from state of the instance) if already present in alignment, uses positionToInsertionCount
     */
    void realignInsertions() {
        realignInsertions(positionToInsertionCount);
    }

    /**
     * Realign insertions
     *  @param positionToIndelsCount - map with position as a key and map with information about insertions of this position
     *                        and count of this position as a value
     */
    //perl version: 2650
    void realignInsertions(Map<Integer, Map<String, Integer>> positionToIndelsCount) {
        List<Object[]> listPositionToIndelsCount = makeListWithPositionToVariationCountInfo(positionToIndelsCount);
        //perl version: 2662
        for (Object[] positionToIndelCount : listPositionToIndelsCount) {
            Integer position = (Integer)positionToIndelCount[0];
            String indel = (String)positionToIndelCount[1];
            Integer indelCount = (Integer)positionToIndelCount[2];
            //perl version: 2665
            String insert;
            Matcher mtch = BEGIN_PLUS_ATGC.matcher(indel);
            if (mtch.find()) {
                insert = mtch.group(1);
            } else {
                continue;
            }
            //perl version: 2668
            String extrains = "";
            mtch = AMP_ATGC.matcher(indel);
            if (mtch.find()) {
                extrains = mtch.group(1);
            }
            String complexVariantMatch = ""; // the match part for a complex variant
            mtch = HASH_ATGC.matcher(indel);
            if (mtch.find()) {
                complexVariantMatch = mtch.group(1);
            }

            //perl version: 2671
            int newDeletion = 0; // the adjacent deletion
            try {
                newDeletion = toInt(indel);
            } catch (NumberFormatException ignore) {}

            String extra = indel.replaceFirst("^\\+", "")
                    .replaceFirst("&", "")
                    .replaceFirst("#", "")
                    .replaceFirst("\\^\\d+$", "")
                    .replaceFirst("\\^", "");

            Variation vref = getVariation(insertionVariants, position, indel);

            processRealignment(position, indelCount, insert.length(), indel.length(), extrains.length() + complexVariantMatch.length() + newDeletion, extra, vref, false);
        }
        //perl version: 2734
        updateStatistics(listPositionToIndelsCount, false);
    }

    /**
     * Processes insertion or deletion softClips
     * @param softClips softClips to process
     * @param position start position
     * @param vref variant
     * @param sequence sequence
     * @param softClipEndPositions list of softClip end positions
     * @param isForward +1 or -1, isForward use revert strand: 1 -> 3' 5'
     * @param isDeletion flag for deletion processing
     */
    private void processIndelsSoftClip(Map<Integer, SoftClip> softClips, Integer position, Variation vref,
                                       String sequence, List<Integer> softClipEndPositions, boolean isForward, boolean isDeletion) {
        //union of four methods: perl version 2684, perl version 2698, perl version 2704, perl version 2718
        //isForward use revert strand: 1 -> 3' 5'
        //TODO decide what we are going to do with sys.err
        for (Integer endPosition : softClipEndPositions) {
            SoftClip tv = softClips.get(endPosition);
            if (tv != null && !tv.used) {
                String consensusSequence = findconseq(tv);
                String cutSequence = isForward ? substr(sequence, endPosition - position - (isDeletion ? 0 : 1)) : sequence;
                if (!consensusSequence.isEmpty() && isMatch(consensusSequence, cutSequence, isForward ? 1 : -1, instance().conf.y)) {
                    if (isForward ? endPosition <= position : endPosition > position) {
                        incCnt(refCoverage, position, tv.varsCount);
                    }
                    Variation lref= null;
                    if (isForward) {
                        if (isDeletion) {
                            lref = endPosition <= position ? null : getVariationMaybe(nonInsertionVariants, position, reference.get(position));
                        } else lref = getVariant(position, endPosition);
                    }
                    adjCnt(vref, tv, lref);
                    tv.used = true;
                }
            }
        }
    }

    /**
     * Realign large deletions that are not present in alignment
     * @throws IOException
     */
    //perl version: 2013
    void realignLagerDeletions() throws IOException {
        //TODO can we figure out a solution with out using a direction
        processLagerDeletionFromOneEnd(softClips5End, -1);
        processLagerDeletionFromOneEnd(softClips3End, 1);
        if (instance().conf.y) {
            System.err.println("  Done: Realignlgdel\n");
        }
    }

    /**
     * Realign large deletions from 5' or 3'
     * @param softClips - softClips from 5' or 3'
     * @param direction - direction, for 3' it is equal to 1, for 5' it is equal to -1
     * @throws IOException
     */
    private void processLagerDeletionFromOneEnd(Map<Integer, SoftClip> softClips, int direction) throws IOException {
        ArrayList<Map.Entry<Integer, SoftClip>> sortedSoftClips = new ArrayList<>(softClips.entrySet());
        sortedSoftClips.sort(COMP2);

        for (Map.Entry<Integer, SoftClip> softClip : sortedSoftClips) {
            int position = softClip.getKey();
            final SoftClip softClipVariation = softClip.getValue();

            if (softClipVariation.used) {
                continue;
            }
            String consensusSequence = findconseq(softClipVariation);
            if (consensusSequence.isEmpty()) {
                continue;
            }
            //TODO We have to use also B_A8 and B_T8 in this case, have we?
            if (B_A7.matcher(consensusSequence).find() || B_T7.matcher(consensusSequence).find()) {
                continue;
            }
            if (consensusSequence.length() < 7) {
                continue;
            }
            if (isLowComplexSeq(consensusSequence)) {
                continue;
            }
            int breakpoint = findBreakpoint(consensusSequence, position + 5 * direction, direction);
            final int deletionLen = direction > 0 ? breakpoint - position : position - breakpoint;
            if (breakpoint == 0) {
                continue;
            }

            StringBuilder extra = new StringBuilder();
            int extraCurrentSize = 0;
            int refPosition = direction > 0 ? breakpoint + extraCurrentSize : breakpoint - extraCurrentSize - 1;
            while (extraCurrentSize < consensusSequence.length()
                    && isNotEquals(consensusSequence.charAt(extraCurrentSize), reference.get(refPosition))) {
                extra.append(consensusSequence.charAt(extraCurrentSize));
                extraCurrentSize++;
                // move to right or left, depends on direction
                refPosition = direction > 0 ? refPosition + 1 : refPosition - 1;
            }

            final String gt;
            if (direction > 0) {
                breakpoint = position; // Set it to 5'
                if (extra.length() > 0) {
                    gt = "-" + deletionLen + "&" + extra;
                } else {
                    gt = "-" + deletionLen;
                    while (Utils.isEquals(reference.get(breakpoint - 1), reference.get(breakpoint + deletionLen - 1))) {
                        breakpoint--;
                    }
                }
            } else {
                if (extra.length() > 0) {
                    extra = extra.reverse();
                    gt = "-" + deletionLen + "&" + extra;
                    breakpoint -= extra.length();
                } else {
                    gt = String.valueOf(-deletionLen);
                }
            }
            if (instance().conf.y) {
                System.err.printf("  Found Realignlgdel: %s %s %s %s %s %s\n", breakpoint, gt,
                        direction > 0 ? " 3' " : " 5' ", position, consensusSequence, softClipVariation.varsCount);
            }
            Variation tv = getVariation(nonInsertionVariants, breakpoint, gt);
            tv.hasAtLeast2DiffQualities = true; // more accurate implementation later
            tv.isAtLeastAt2Positions = true; // more accurate implementation later
            for (int tp = breakpoint; tp < breakpoint + deletionLen; tp++) {
                incCnt(refCoverage, tp, softClipVariation.varsCount);
            }
            if (direction > 0) {
                softClipVariation.meanPosition += deletionLen * softClipVariation.varsCount;
            }
            adjCnt(tv, softClipVariation);
            softClipVariation.used = true;

            if (direction < 0) {
                // Work on the softclipped read at 3'
                processSoftclippedReadsOn3End(breakpoint, deletionLen, tv);
            }

            Map<Integer, Map<String, Integer>> deletions = singletonMap(breakpoint, singletonMap(gt, tv.varsCount));
            realignDeletions(deletions);
            if (instance().conf.y) {
                System.err.printf("  Found lgdel: %s %s $p %s %s %s\n\n", breakpoint, gt,
                        direction > 0 ? "3'" : "5'",  position, tv.varsCount);
            }
        }
    }

    /**
     * Work on the softclipped read at 3' during the process of working on realigning lager deletion
     * @param breakpoint breakpoint
     * @param deletionLen deletion
     * @param tv variant
     */
    private void processSoftclippedReadsOn3End(int breakpoint, int deletionLen, Variation tv) {
        int n = 0;
        while (reference.containsKey(breakpoint + n)
                && reference.containsKey(breakpoint + deletionLen + n)
                && Utils.isEquals(reference.get(breakpoint + n), reference.get(breakpoint + deletionLen + n))) {
            n++;
        }
        int softClipPosition = breakpoint + n;
        StringBuilder str = new StringBuilder();
        int mismatchCount = 0;
        while (mismatchCount <= MAX_MISMATCH_ALLOWED
                && reference.containsKey(breakpoint + n)
                && reference.containsKey(breakpoint + deletionLen + n)
                && isNotEquals(reference.get(breakpoint + n), reference.get(breakpoint + deletionLen + n))) {
            str.append(reference.get(breakpoint + deletionLen + n));
            n++;
            mismatchCount++;
        }
        if (str.length() == 1) {
            while (reference.containsKey(breakpoint + n)
                    && reference.containsKey(breakpoint + deletionLen + n)
                    && Utils.isEquals(reference.get(breakpoint + n), reference.get(breakpoint + deletionLen + n))) {
                n++;
            }
            softClipPosition = breakpoint + n;
        }
        if (softClips3End.containsKey(softClipPosition) && !softClips3End.get(softClipPosition).used) {
            SoftClip softClip = softClips3End.get(softClipPosition);
            if (softClipPosition > breakpoint) {
                adjCnt(tv, softClip, getVariationMaybe(nonInsertionVariants, breakpoint, reference.get(breakpoint)));
            } else {
                adjCnt(tv, softClip);
            }

            if (softClipPosition == breakpoint) {
                for (int tp = breakpoint; tp < breakpoint + deletionLen; tp++) {
                    incCnt(refCoverage, tp, softClips3End.get(softClipPosition).varsCount);
                }
            }
            for (int ip = breakpoint + 1; ip < softClipPosition; ip++) {
                Variation vv = getVariation(nonInsertionVariants, ip, reference.get(deletionLen + ip).toString());
                vv.rmCnt(softClip);
                if (vv.varsCount == 0) {
                    nonInsertionVariants.get(ip).remove(reference.get(deletionLen + ip).toString());
                }
                if (nonInsertionVariants.get(ip).size() == 0) {
                    nonInsertionVariants.remove(ip);
                }
            }
            softClip.used = true;
        }
    }

    /**
     * Realign large insertions that are not present in alignment
     * @throws IOException
     */
    //perl version: 1912
    void realignLagerInsertions() throws IOException {
        processLagerInsertionFromOneEnd(softClips5End, -1);
        processLagerInsertionFromOneEnd(softClips3End, 1);
    }

    /**
     * Realign large insertions from one end that are not present in alignment
     * @param softClips softClips to process
     * @param direction direction to move (+1 or -1)
     * @throws IOException
     */
    private void processLagerInsertionFromOneEnd(Map<Integer, SoftClip> softClips, int direction) throws IOException {
        ArrayList<Map.Entry<Integer, SoftClip>> sortedSoftClips = new ArrayList<>(softClips.entrySet());
        sortedSoftClips.sort(COMP2);
        for (Map.Entry<Integer, SoftClip> softClip : sortedSoftClips) {
            int position = softClip.getKey();
            final SoftClip softClipVariation = softClip.getValue();
            //already been used in
            if (softClipVariation.used) {
                continue;
            }
            String consensusSequence = findconseq(softClipVariation);
            if (consensusSequence.isEmpty()) {
                continue;
            }
            if (instance().conf.y) {
                System.err.println("  Working lgins: 5: " + position + " " + consensusSequence);
            }
            if (direction > 0
                    ? B_A7.matcher(consensusSequence).find() || B_T7.matcher(consensusSequence).find()
                    : B_A8.matcher(consensusSequence).find() || B_T8.matcher(consensusSequence).find()) {
                continue;
            }
            if (consensusSequence.length() < 12) {
                continue;
            }
            if (isLowComplexSeq(consensusSequence)) {
                continue;
            }
            Tuple.Tuple3<Integer, String, Integer> insertionBounds = findInsertionBounds(consensusSequence, position, reference, direction, chr);
            final int firstBound = insertionBounds._1;
            final String insertion = insertionBounds._2;
            if (firstBound == 0) {
                continue;
            }
            if (instance().conf.y) {
                System.err.printf("  Found candidate lgins from %s : %s +%s %s %s\n", direction < 0 ? "5" : "3", firstBound, insertion, position, consensusSequence);
            }
            final Variation iref = getVariation(insertionVariants, firstBound, "+" + insertion);
            iref.isAtLeastAt2Positions = true;
            iref.hasAtLeast2DiffQualities = true;
            final Variation lref = direction > 0
                    ? getVariationMaybe(nonInsertionVariants, firstBound, reference.get(firstBound))
                    : null;
            adjCnt(iref, softClipVariation, lref);
            boolean rpflag = true; // A flag to indicate whether an insertion is a repeat
            for (int i = 0; i < insertion.length(); i++) {
                if (!Utils.isEquals(reference.get(firstBound + 1 + i), insertion.charAt(i))) {
                    rpflag = false;
                    break;
                }
            }
            if (direction < 0) {
                incCnt(refCoverage, firstBound, softClipVariation.varsCount);
            }
            int insertionLen = insertion.length();
            if (insertion.indexOf('&') != -1) {
                insertionLen--;
            }
            int seqLen = softClipVariation.seq.lastKey() + 1;
            for (int ii = direction > 0 ? insertionLen : insertionLen + 1; ii < seqLen; ii++) {
                int pii = direction > 0 ? position + ii - insertionLen :  firstBound - ii + insertionLen;
                if (!softClipVariation.seq.containsKey(ii)) {
                    continue;
                }
                for (Map.Entry<Character, Variation> ent : softClipVariation.seq.get(ii).entrySet()) {
                    Character tnt = ent.getKey();
                    Variation tv = ent.getValue();
                    Variation tvr = getVariation(nonInsertionVariants, pii, tnt.toString());
                    adjCnt(tvr, tv);
                    tvr.isAtLeastAt2Positions = true;
                    tvr.hasAtLeast2DiffQualities = true;
                    incCnt(refCoverage, pii, tv.varsCount);
                }
            }
            softClipVariation.used = true;
            Map<Integer, Map<String, Integer>> tins = singletonMap(firstBound, singletonMap("+" + insertion, iref.varsCount));
            realignInsertions(tins);
            Variation mref = getVariationMaybe(nonInsertionVariants, firstBound, reference.get(firstBound));
            if (rpflag && bams.length > 0 && insertion.length() >= 5
                    && insertion.length() < maxReadLength - 10
                    && mref != null && mref.varsCount != 0
                    && noPassingReads(chr, firstBound, firstBound + insertion.length(), bams)
                    && iref.varsCount > 2 * mref.varsCount) {
                adjCnt(iref, mref, mref);
            }
        }
    }

    /**
     * This will try to realign large insertions (typically larger than 30bp)
     * @throws IOException
     */
    //perl version: 1770
    void realignInsertionsLagerThan30() throws IOException {

        List<Tuple.Tuple3<Integer, SoftClip, Integer>> tmp5 = new ArrayList<>();
        for (Map.Entry<Integer, SoftClip> ent5 : softClips5End.entrySet()) {
            //ent5: (position, soft clipping 5' variant)
            tmp5.add(tuple(ent5.getKey(), ent5.getValue(), ent5.getValue().varsCount));
        }
        tmp5.sort(COMP3);

        List<Tuple.Tuple3<Integer, SoftClip, Integer>> tmp3 = new ArrayList<>();
        for (Map.Entry<Integer, SoftClip> ent3 : softClips3End.entrySet()) {
            tmp3.add(tuple(ent3.getKey(), ent3.getValue(), ent3.getValue().varsCount));
        }
        tmp3.sort(COMP3);

        for (Tuple.Tuple3<Integer, SoftClip, Integer> t5 : tmp5) {
            final int p5 = t5._1;
            final SoftClip sc5v = t5._2;
            final int cnt5 = t5._3;
            if (sc5v.used) {
                continue;
            }

            for (Tuple.Tuple3<Integer, SoftClip, Integer> t3 : tmp3) {
                final int p3 = t3._1;
                final SoftClip sc3v = t3._2;
                final int cnt3 = t3._3;
                if (sc3v.used) {
                    continue;
                }
                if (p5 - p3 > maxReadLength / 1.5) {
                    continue;
                }
                if (p3 - p5 > maxReadLength - 10) { // if they're too far away, don't even try
                    continue;
                }
                final String seq5 = findconseq(sc5v);
                final String seq3 = findconseq(sc3v);
                //next until at least one of consensus sequences has length > 10
                if (seq5.length() <= 10 || seq3.length() <= 10) {
                    continue;
                }
                if (instance().conf.y) {
                    System.err.printf("  Working lgins30: %s %s 3: %s %s 5: %s %s\n",
                            p3, p5, seq3, cnt3, new StringBuilder(seq5).reverse(), cnt5);
                }
                Tuple.Tuple3<Integer, Integer, Integer> tpl = find35match(seq5, seq3);
                int bp5 = tpl._1;
                int bp3 = tpl._2;
                //length of match
                int score = tpl._3;

                if (score == 0) {
                    continue;
                }
                //add matched part of seq3 (left clip)
                String ins = bp3 > 1 ? substr(seq3, 0, -bp3 + 1) : seq3;
                if (bp5 > 0) { //add not matched part of seq5 (right clip)
                    ins += new StringBuilder(substr(seq5, 0, bp5)).reverse();
                }
                if (isLowComplexSeq(ins)) {
                    if (instance().conf.y) {
                        System.err.println("  Discard low complex insertion found " + ins + ".");
                    }
                    continue;
                }
                int bi;
                Variation vref;
                if (instance().conf.y) {
                    System.err.printf("  Found candidate lgins30: %s %s %s\n", p3, p5, ins);
                }
                if (p5 > p3) {
                    String partOfRef = joinRef(reference, p5, p5 + seq3.length() - ins.length() + 2);
                    if (seq3.length() > ins.length() && !isMatch(substr(seq3, ins.length()), partOfRef, 1, instance().conf.y)) {
                        continue;
                    }
                    partOfRef = joinRef(reference, p3 - (seq5.length() - ins.length()) - 2, p3 - 1);
                    if (seq5.length() > ins.length() && !isMatch(substr(seq5, ins.length()), partOfRef, -1, instance().conf.y)) {
                        continue;
                    }
                    if (instance().conf.y) {
                        System.err.printf("  Found lgins30 complex: %s %s %s %s\n", p3, p5, ins.length(), ins);
                    }
                    String tmp = joinRef(reference, p3, p5 - 1);
                    if (tmp.length() > ins.length()) { // deletion is longer
                        ins = (p3 - p5) + "^" + ins;
                        bi = p3;
                        vref = getVariation(nonInsertionVariants, p3, ins);
                    } else if (tmp.length() < ins.length()) {
                        ins = substr(ins, 0, ins.length() - tmp.length()) + "&" + substr(ins, p3 - p5);
                        ins = "+" + ins;
                        bi = p3 - 1;
                        vref = getVariation(insertionVariants, bi, ins);
                    } else { // long MNP
                        ins = "-" + ins.length() + "^" + ins;
                        bi = p3;
                        vref = getVariation(nonInsertionVariants, p3, ins);
                    }
                } else {
                    String partOfRef = joinRef(reference, p5, p5 + seq3.length() - ins.length() + 2);
                    if (seq3.length() > ins.length() && !isMatch(substr(seq3, ins.length()), partOfRef, 1, instance().conf.y)) {
                        continue;
                    }
                    partOfRef = joinRef(reference, p3 - (seq5.length() - ins.length()) - 2, p3 - 1);
                    if (seq5.length() > ins.length() && !isMatch(substr(seq5, ins.length()), partOfRef, -1, instance().conf.y)) {
                        continue;
                    }
                    String tmp = ins.length() > p3 - p5 ? joinRef(reference, p5, p3)
                            : joinRef(reference, p5, p5 + (p3 - p5 - ins.length()) / 2); // Tandem duplication
                    if (instance().conf.y) {
                        System.err.printf("Found lgins30: %s %s %s %s + %s\n", p3, p5, ins.length(), tmp, ins);
                    }
                    ins = "+" + tmp + ins;
                    bi = p5 - 1;
                    vref = getVariation(insertionVariants, bi, ins);
                }
                sc3v.used = true;
                sc5v.used = true;
                vref.isAtLeastAt2Positions = true;
                vref.hasAtLeast2DiffQualities = true;
                incCnt(refCoverage, bi, sc5v.varsCount);
                if (instance().conf.y) {
                    System.err.printf(" lgins30 Found: '%s' %s %s %s\n", ins, bi, bp3, bp5);
                }

                if (ins.startsWith("+")) {
                    Variation mvref = getVariationMaybe(nonInsertionVariants, bi, reference.get(bi));
                    adjCnt(vref, sc3v, mvref);
                    adjCnt(vref, sc5v);
                    if (bams != null && bams.length > 0
                            && p3 - p5 >= 5 && p3 - p5 > maxReadLength - 10
                            && mvref != null && mvref.varsCount != 0
                            && vref.varsCount > 2 * mvref.varsCount
                            && noPassingReads(chr, p5, p3, bams)) {
                        adjCnt(vref, mvref, mvref);
                    }
                    Map<Integer, Map<String, Integer>> tins = new HashMap<>();
                    Map<String, Integer> map = new HashMap<>();
                    map.put(ins, vref.varsCount);
                    tins.put(bi, map);
                    realignInsertions(tins);
                } else if (ins.startsWith("-")) {
                    adjCnt(vref, sc3v, getVariationMaybe(nonInsertionVariants, bi, reference.get(bi)));
                    adjCnt(vref, sc5v);
                    Map<Integer, Map<String, Integer>> tdel = new HashMap<>();
                    Map<String, Integer> map = new HashMap<>();
                    map.put(ins, vref.varsCount);
                    tdel.put(bi, map);
                    realignDeletions(tdel);
                } else {
                    adjCnt(vref, sc3v);
                    adjCnt(vref, sc5v);
                }
                break;
            }

        }
        if (instance().conf.y) {
            System.err.println("Done: lgins30\n");
        }
    }

    /**
     * Checks that the variant has good parameters and meets to the definition "relevant"
     */
    static boolean isRelevantVariant(Integer variationCount, int varLength, int nm3, int nm5, Integer me, Variation tv) {
        if (tv == null) {
            return false;
        }
        if (tv.varsCount == 0) {
            return false;
        }
        if (tv.meanQuality / tv.varsCount < instance().conf.goodq) {
            return false;
        }
        if (tv.meanPosition / tv.varsCount > (me == 3 ? nm3 + 4 : nm5 + 4)) { // opt_k;
            return false;
        }
        if (tv.varsCount >= variationCount + varLength || tv.varsCount / variationCount >= 8) {
            return false;
        }
        return true;
    }

    /**
     * check whether there're reads supporting wild type in deletions
     * Only for indels that have micro-homology
     * @param chr chromosome name
     * @param s start position
     * @param e end position
     * @param bams BAM file list
     * @return true if any read was found in chr:s-e
     * @throws IOException
     */
    //perl version: 2495
    static boolean noPassingReads(String chr, int s, int e, String[] bams) throws IOException {
        int cnt = 0;
        int midcnt = 0; // Reads end in the middle
        int dlen = e - s;
        String dlenqr = dlen + "D";
        Region region = new Region(chr, s, e, "");
        for (String bam : bams) {
            try (SamView reader = new SamView(bam, "", region, instance().conf.validationStringency)) {
                SAMRecord record;
                while ((record = reader.read()) != null) {
                    if (record.getCigarString().contains(dlenqr)) {
                        continue;
                    }
                    int rs = record.getAlignmentStart();
                    int rlen = getAlignedLenght(record.getCigar()); // The total aligned length, excluding soft-clipped bases and
                    // insertions
                    int re = rs + rlen;
                    if (re > e + 2 && rs < s - 2) {
                        cnt++;
                    }
                    if (rs < s - 2 && re > s && re < e) {
                        midcnt++;
                    }
                }

            }
        }
        if (instance().conf.y) {
            System.err.printf("    Passing Read CNT: %s %s %s %s %s\n", cnt, chr, s, e, midcnt);
        }
        return cnt <= 0;
    }

    /**
     * Find if sequences match
     * @param seq1 first sequence
     * @param seq2 second sequence
     * @param dir direction of seq2 (1 or -1)
     * @param debugLog print debug message if true
     * @return true if seq1 matches seq2 with no more than 2 mismathes
     */
    //perl version: 2751
    static boolean isMatch(String seq1, String seq2, int dir, boolean debugLog) {
        if (debugLog) {
            System.err.printf("    Matching %s %s %s\n", seq1, seq2, dir);
        }
        seq2 = seq2.replaceAll("#|\\^", "");
        int mm = 0;
        for (int n = 0; n < seq1.length() && n < seq2.length(); n++) {
            if (seq1.charAt(n) != substr(seq2, dir * n - (dir == -1 ? 1 : 0), 1).charAt(0)) {
                mm++;
            }
        }

        return (mm <= 2 && mm / (double)seq1.length() < 0.15);
    }

    /**
     * Given a variant sequence, find the mismatches and potential softclipping positions
     * @param ref map of reference bases
     * @param p position
     * @param sequence sequence
     * @param sclip map of 3' or 5' softclips
     * @param is3End flag for 3' softclips processing
     * @return mismatchResult - the class that incorporates mismatches - the list of [mismatches, softClipPositions, 3],
     * softClipPositions - list of soft clip positions, nm - number of nucleotides to adjust
     */
    //Union of two methods. Perl version: 2572, 2612
    static MismatchResult findMismatchesOn35End(Map<Integer, Character> ref, int p, String sequence, Map<Integer, SoftClip> sclip, boolean is3End) {
        int direction = (is3End) ? 1 : -1;
        String seq = sequence.replaceAll("#|\\^", "");
        List<Tuple.Tuple3<String, Integer, Integer>> mismatches = new ArrayList<>(); // mismatches, mismatch positions, 5 or 3 ends
        int n = 0;
        int numberOfMismatches = 0;
        int mcnt = 0;
        List<Integer> softClipPositions = new ArrayList<>();
        StringBuilder str = new StringBuilder();
        if (is3End) {
            while (n < seq.length() && Utils.isEquals(ref.get(p + direction * n), charAt(seq, direction * n))) {
                n++;
            }
            softClipPositions.add(p + n);
        } else softClipPositions.add(p + 1);

        while (is3End ?
                n < seq.length() && isNotEquals(ref.get(p + n), charAt(seq,n)) && mcnt <= MAX_MISMATCH_ALLOWED
                : isHasAndNotEquals(charAt(seq, -1 - n), ref, p - n) && mcnt < MAX_MISMATCH_ALLOWED) {
            str.append(charAt(seq, getNIndex(is3End, n)));
            mismatches.add(tuple(str.toString(), p + direction * n, (is3End) ? 3 : 5));
            n++;
            mcnt++;
        }
        // Adjust clipping position if only one mismatch
        int mismatchPosition = 0;
        Character misnt = null;
        if (str.length() == 1) {
            while ((!is3End || n < seq.length()) && isHasAndEquals(charAt(seq, getNIndex(is3End, n)), ref, p + direction * n)) {
                n++;
                numberOfMismatches++;
            }
            int currentPos = p + direction * n;
            if (numberOfMismatches > 1) {
                int n2 = 0;
                while ((is3End ? n + n2 + 1 < seq.length() : -1 - n - 1 - n2 >= 0) &&
                        isHasAndEquals(charAt(seq, getNIndex(is3End, n) + direction * (n2 + 1)), ref, p + direction * (n + 1 + n2))) {
                    n2++;
                }
                if (n2 > 2 && (!is3End || n + n2 + 1 < seq.length())) {
                    softClipPositions.add(currentPos + direction * n2);
                    mismatchPosition = currentPos;
                    misnt = charAt(seq, getNIndex(is3End, n));
                    if (sclip.containsKey(currentPos + direction * n2)) {
                        sclip.get(currentPos + direction * n2).used = true;
                    }
                    numberOfMismatches += n2;
                } else {
                    softClipPositions.add(currentPos);
                    if (sclip.containsKey(currentPos)) {
                        sclip.get(currentPos).used = true;
                    }
                }
            }
        }
        return new MismatchResult(mismatches, softClipPositions, numberOfMismatches, mismatchPosition, misnt == null ? "" : misnt.toString());
    }

    private static int getNIndex(boolean is3End, int n) {
        return is3End ? n : -1 - n;
    }

    /**
     * Test whether the two soft-clipped reads match
     * Returns the breakpoints in 5 and 3 prime soft-clipped reads
     * @param seq5 consensus sequence 5' strand
     * @param seq3 consensus sequence 3' strand
     * @return Tuple of (start position of 5' strand, start position of 3' strand, max match length)
     */
    //perl version: 1883
    static Tuple.Tuple3<Integer, Integer, Integer> find35match(String seq5, String seq3) {
        int max = 0;
        int b3 = 0;
        int b5 = 0;
        for (int i = 0; i < seq5.length() - 8; i++) {
            for (int j = 1; j < seq3.length() - 8; j++) {
                int nm = 0;
                int n = 0;
                while (n + j <= seq3.length() && i + n <= seq5.length()) {
                    if (!substr(seq3, -j - n, 1).equals(substr(seq5, i + n, 1))) {
                        nm++;
                    }
                    if (nm > MAX_MISMATCH_ALLOWED) {
                        break;
                    }
                    n++;
                }
                if (n - nm > max
                        && n - nm > 8
                        && nm / (double)n < 0.1d
                        && (n + j >= seq3.length() || i + n >= seq5.length())) {

                    max = n - nm;
                    b3 = j;
                    b5 = i;
                    return tuple(b5, b3, max);
                }

            }
        }
        return tuple(b5, b3, max);
    }

    /**
     * Find the insertion
     * @param seq sequence
     * @param p position
     * @param ref map of reference bases
     * @param dir direction
     * @param chr chromosome name
     * @return Tuple of (BI (insert starting position), INS (insert sequence), BI2 ( = BI))
     */
    //perl version: 2206
    static Tuple.Tuple3<Integer, String, Integer> findInsertionBounds(String seq, int p, Map<Integer, Character> ref,
                                                                      final int dir, String chr) {
        final int dirExt = dir == -1 ? 1 : 0;
        int score = 0;
        int bi = 0;
        String ins = "";
        int bi2 = 0;

        for (int n = 6; n < seq.length(); n++) {
            if (p + 6 >= instance().chrLengths.get(chr)) {
                break;
            }
            int mm = 0;
            int i;
            Set<Character> m = new HashSet<>();
            for (i = 0; i + n < seq.length(); i++) {
                if (p + dir * i - dirExt < 1) {
                    break;
                }
                if (p + dir * i - dirExt > instance().chrLengths.get(chr)) {
                    break;
                }
                if (isNotEquals(seq.charAt(i + n), ref.get(p + dir * i - dirExt))) {
                    mm++;
                } else {
                    m.add(seq.charAt(i + n));
                }
                if (mm > MAX_MISMATCH_ALLOWED) {
                    break;
                }
            }
            int mnt = m.size();
            if (mnt < 2) { // at least three different NT for overhang sequences, weeding out low complexity seq
                continue;
            }
            if ((mnt >= 3 && i + n >= seq.length() - 1 && i >= 8 && mm / (double)i < 0.15)
                    || (mnt >= 2 && mm == 0 && i + n == seq.length() && n >= 20 && i >= 8)) {

                StringBuilder insert = new StringBuilder(substr(seq, 0, n));
                StringBuilder extra = new StringBuilder();
                int ept = 0;
                while (n + ept + 1 < seq.length() && (!Utils.isEquals(seq.charAt(n + ept), ref.get(p + ept * dir - dirExt))
                        || !Utils.isEquals(seq.charAt(n + ept + 1), ref.get(p + (ept + 1) * dir - dirExt)))) {
                    extra.append(seq.charAt(n + ept));
                    ept++;
                }
                if (dir == -1) {
                    insert.append(extra);
                    insert.reverse();
                    if (extra.length() > 0) {
                        insert.insert(insert.length() - extra.length(), "&");
                    }
                    if (mm == 0 && i + n == seq.length()) {
                        bi = p - 1 - extra.length();
                        ins = insert.toString();
                        bi2 = p - 1;
                        if (extra.length() == 0) {
                            Tuple.Tuple3<Integer, String, Integer> tpl = adjInsPos(bi, ins, ref);
                            bi = tpl._1;
                            ins = tpl._2;
                            bi2 = tpl._3;
                        }
                        return tuple(bi, ins, bi2);
                    } else if (i - mm > score) {
                        bi = p - 1 - extra.length();
                        ins = insert.toString();
                        bi2 = p - 1;
                        score = i - mm;
                    }
                } else {
                    int s = -1;
                    if (extra.length() > 0) {
                        insert.append("&").append(extra);
                    } else {
                        while (s >= -n && Utils.isEquals(charAt(insert, s), ref.get(p + s))) {
                            s--;
                        }
                        if (s < -1) {
                            String tins = substr(insert.toString(), s + 1, 1 - s);
                            insert.delete(insert.length() + s + 1, insert.length());
                            insert.insert(0, tins);
                        }

                    }
                    if (mm == 0 && i + n == seq.length()) {
                        bi = p + s;
                        ins = insert.toString();
                        bi2 = p + s + extra.length();
                        if (extra.length() == 0) {
                            Tuple.Tuple3<Integer, String, Integer> tpl = adjInsPos(bi, ins, ref);
                            bi = tpl._1;
                            ins = tpl._2;
                            bi2 = tpl._3;
                        }
                        return tuple(bi, ins, bi2);
                    } else if (i - mm > score) {
                        bi = p + s;
                        ins = insert.toString();
                        bi2 = p + s + extra.length();
                        score = i - mm;
                    }
                }
            }
        }
        if (bi2 == bi && ins.length() > 0 && bi != 0) {
            Tuple.Tuple3<Integer, String, Integer> tpl = adjInsPos(bi, ins, ref);
            bi = tpl._1;
            ins = tpl._2;
        }
        return tuple(bi, ins, bi2);
    }

    /**
     * Find breakpoint in a specified sequence
     * @param seq sequence
     * @param startPosition start position of sequence
     * @param dir direction
     * @return breakpoint position
     */
    //perl version 2298:  findbp()
    int findBreakpoint(String seq, int startPosition, int dir) {
        int bp = 0;
        int score = 0;
        int idx = instance().chrLengths.containsKey(chr) ? instance().chrLengths.get(chr) : 0;
        for (int n = 0; n < instance().conf.indelsize; n++) {
            int mm = 0;
            int i;
            Set<Character> m = new HashSet<>();
            for (i = 0; i < seq.length(); i++) {
                if (startPosition + dir * n + dir * i < 1) {
                    break;
                }
                if (startPosition + dir * n + dir * i > idx) {
                    break;
                }
                if (Utils.isEquals(seq.charAt(i), reference.get(startPosition + dir * n + dir * i))) {
                    m.add(seq.charAt(i));
                } else {
                    mm++;
                }
                if (mm > MAX_MISMATCH_ALLOWED - n / 100) {
                    break;
                }
            }
            if (m.size() < 3) {
                continue;
            }
            if (mm <= MAX_MISMATCH_ALLOWED - n / 100 && i >= seq.length() - 2 && i >= 8 + n / 10 && mm / (double)i < 0.12) {
                int lbp = startPosition + dir * n - (dir < 0 ? dir : 0);
                if (mm == 0 && i == seq.length()) {
                    if (instance().conf.y) {
                        System.err.printf("  Findbp: %s %s %s %s %s\n", seq, startPosition, lbp, mm, i);
                    }
                    return lbp;
                } else if (i - mm > score) {
                    bp = lbp;
                    score = i - mm;
                }
            }
        }
        if (instance().conf.y && bp != 0) {
            System.err.printf("  Findbp with mismatches: %s %s %s %s %s\n", seq, startPosition, bp, dir, score);
        }
        return bp;
    }

    private static int count(String str, char chr) {
        int cnt = 0;
        for (int i = 0; i < str.length(); i++) {
            if (str.charAt(i) == chr) {
                cnt++;
            }
        }
        return cnt;

    }

    /**
     * If any of the base is represented by 75%, it's a low complexity seq
     *
     * Low-complexity sequences are simple repeats such as ATATATATAT TODO it's wrong, isn't it?
     * or regions that are highly enriched for just one letter, e.g. AAACAAAAAAAAGAAAAAAC.
     * @param seq sequence
     * @return true if sequence has low complexity
     */
    //TODO we can count all nucleotides in one pass
    //perl version: 1754: islowcomplexseq
    static boolean isLowComplexSeq(String seq) {
        int len = seq.length();
        if (len == 0)
            return true;
        if (count(seq, 'A') / (double)len > 0.75)
            return true;
        if (count(seq, 'T') / (double)len > 0.75)
            return true;
        if (count(seq, 'G') / (double)len > 0.75)
            return true;
        return count(seq, 'C') / (double)len > 0.75;
    }

    /**
     * Forms the list of arrays with information about variations
     * @param changes map with position as a key and map with information about deletions of this position
     * and count of this position as a value
     * @return list of arrays [position, variation, count]
     */
    //perl version: 2653, 2020 (multiple usage)
    private static List<Object[]> makeListWithPositionToVariationCountInfo(Map<Integer, Map<String, Integer>> changes) {
        List<Object[]> tmp = new ArrayList<>();
        for (Map.Entry<Integer, Map<String, Integer>> positionToVariationCount : changes.entrySet()) {
            int position = positionToVariationCount.getKey();
            Map<String, Integer> v = positionToVariationCount.getValue();
            for (Map.Entry<String, Integer> entV : v.entrySet()) {
                String variation = entV.getKey();
                int count = entV.getValue();
                tmp.add(new Object[] { position, variation, count, /* ecnt */});
            }
        }
        tmp.sort(REALIGNDEL_COMPARATOR);
        return tmp;
    }

    /**
     * Returns the end position of the 3-end flanking sequence
     * @param position specified position in the reference
     * @param shiftLen shift to apply
     * @return end position of the 3-end flanking sequence
     */
    private int get3EndFlankingSeqEnd(Integer position, int shiftLen) {
        int sanend = position + shiftLen + 100;
        if (sanend > instance().chrLengths.get(chr)) {
            sanend = instance().chrLengths.get(chr);
        }
        return sanend;
    }

    /**
     * Returns the start position of the 5-end flanking sequence
     * @param position specified position in the reference
     * @param shiftLen shift to apply
     * @return start position of the 5-end flanking sequence
     */
    private int get5EndFlankingSeqStart(Integer position, int shiftLen) {
        int wustart = position - shiftLen - 100;
        if (wustart <= 1) {
            wustart = 1;
        }
        return wustart;
    }

    /**
     * Returns the variant from a specified position if it meets certain criteria, null otherwise
     * @param position specified position in the reference
     * @param endPosition end position of the reference softclip
     * @return variant
     */
    private Variation getVariant(Integer position, Integer endPosition) {
        Variation lref = null;
        if (endPosition > position && nonInsertionVariants.containsKey(position) &&
                reference.containsKey(position) && nonInsertionVariants.get(position).containsKey(reference.get(position).toString())) {
            lref = nonInsertionVariants.get(position).get(reference.get(position).toString());
        }
        return lref;
    }

    static class MismatchResult {
        private final List<Tuple.Tuple3<String, Integer, Integer>> mismatches;
        private final List<Integer> softClipPositions;
        private final int numberOfMismatches;
        private final int mismatchPosition;
        private final String mismatchVariant;

        public MismatchResult(List<Tuple.Tuple3<String, Integer, Integer>> mismatches, List<Integer> softClipPositions, int numberOfMismatches, int mismatchPosition, String mismatchVariant) {
            this.mismatches = mismatches;
            this.softClipPositions = softClipPositions;
            this.numberOfMismatches = numberOfMismatches;
            this.mismatchPosition = mismatchPosition;
            this.mismatchVariant = mismatchVariant;
        }

    }

}
