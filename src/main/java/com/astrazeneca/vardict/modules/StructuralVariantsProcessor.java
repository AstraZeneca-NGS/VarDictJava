package com.astrazeneca.vardict.modules;

import com.astrazeneca.vardict.Configuration;
import com.astrazeneca.vardict.data.ReferenceResource;
import com.astrazeneca.vardict.collection.DirectThreadExecutor;
import com.astrazeneca.vardict.data.Region;
import com.astrazeneca.vardict.collection.Tuple;
import com.astrazeneca.vardict.collection.VariationMap;
import com.astrazeneca.vardict.data.*;
import com.astrazeneca.vardict.data.scopedata.InitialData;
import com.astrazeneca.vardict.data.scopedata.RealignedVariationData;
import com.astrazeneca.vardict.data.scopedata.Scope;
import com.astrazeneca.vardict.data.scopedata.VariationData;
import com.astrazeneca.vardict.printers.VariantPrinter;
import com.astrazeneca.vardict.variations.*;
import htsjdk.samtools.util.SequenceUtil;

import java.time.LocalDateTime;
import java.util.*;

import static com.astrazeneca.vardict.data.scopedata.GlobalReadOnlyScope.getMode;
import static com.astrazeneca.vardict.data.scopedata.GlobalReadOnlyScope.instance;
import static com.astrazeneca.vardict.collection.VariationMap.getSV;
import static com.astrazeneca.vardict.data.ReferenceResource.isLoaded;
import static com.astrazeneca.vardict.collection.VariationMap.SV;
import static com.astrazeneca.vardict.modules.VariationRealigner.ismatchref;
import static com.astrazeneca.vardict.variations.VariationUtils.*;
import static com.astrazeneca.vardict.Utils.*;

import static java.lang.Math.abs;
import static java.util.Collections.reverseOrder;
import static java.util.Comparator.comparing;
import static java.util.stream.Collectors.toList;

/**
 * The main class for finding structural variants: INV, DEL and DUPs.
 */
public class StructuralVariantsProcessor implements Module<RealignedVariationData, RealignedVariationData>  {

    private Map<Integer, VariationMap<String, Variation>> nonInsertionVariants;
    private Map<Integer, VariationMap<String, Variation>> insertionVariants;
    private Map<Integer, List<Sclip>> SOFTP2SV;
    private Map<Integer, Integer> refCoverage;
    private Map<Integer, Sclip> softClips5End;
    private Map<Integer, Sclip> softClips3End;
    private Reference reference;
    private ReferenceResource referenceResource;
    private Region region;
    private Set<String> splice;
    private int maxReadLength;
    private String[] bams;
    private String bam;
    private SVStructures svStructures;
    private double duprate;
    private CurrentSegment CURSEG;
    private Scope<VariationData> previousScope;
    private VariantPrinter variantPrinter;

    private void initFromScope(Scope<RealignedVariationData> scope) {
        this.reference = scope.regionRef;
        this.CURSEG = scope.data.CURSEG;
        this.SOFTP2SV = scope.data.SOFTP2SV;
        this.referenceResource = scope.referenceResource;
        this.region = scope.region;
        this.nonInsertionVariants = scope.data.nonInsertionVariants;
        this.insertionVariants = scope.data.insertionVariants;
        this.refCoverage = scope.data.refCoverage;
        this.softClips5End = scope.data.softClips5End;
        this.softClips3End = scope.data.softClips3End;
        this.maxReadLength = scope.maxReadLength;
        this.bams = scope.bam != null ? scope.bam.split(":") : null;
        this.bam = scope.bam;
        this.splice = scope.splice;
        this.svStructures = scope.data.svStructures;
        this.duprate = scope.data.duprate;
        this.previousScope = scope.data.previousScope;
        this.variantPrinter = scope.out;
    }

    /**
     * Runs the structural variants processing on data got from Variation Realigner. Also print remaining soft-clipped
     * reads that haven't been used.
     * @param scope contains filled realigned variation data (maps of insertion and non-insertion variations)
     * @return updated realigned variation data contains structural variants
     */
    @Override
    public Scope<RealignedVariationData> process(Scope<RealignedVariationData> scope)  {
        initFromScope(scope);

        if (!instance().conf.disableSV) {
            findAllSVs();
        }
        adjSNV();
        if (instance().conf.y) {
            outputClipping();
            System.err.println("TIME: Finish realign:" + LocalDateTime.now());
        }
        return new Scope<>(
                scope,
                new RealignedVariationData(nonInsertionVariants, insertionVariants, softClips3End, softClips5End,
                        refCoverage, maxReadLength, duprate, svStructures, CURSEG, SOFTP2SV, previousScope));
    }

    /**
     * Find SVs for each structural variant structure
     */
    public void findAllSVs() {
        if (instance().conf.y) {
            System.err.println("Start Structural Variants: DEL\n");
        }
        findDEL();
        if (instance().conf.y) {
            System.err.println("Start Structural Variants: INV\n");
        }
        findINV();
        if (instance().conf.y) {
            System.err.println("Start Structural Variants\n");
        }
        findsv();
        if (instance().conf.y) {
            System.err.println("Start Structural Variants: DEL discordant pairs only\n");
        }
        findDELdisc() ;
        if (instance().conf.y) {
            System.err.println("Start Structural Variants: INV discordant pairs only\n");
        }
        findINVdisc();
        if (instance().conf.y) {
            System.err.println("Start Structural Variants: DUP discordant pairs only\n");
        }
        findDUPdisc();
    }


    /**
     * Find DEL SV
     */
    void findDEL() {
        int lastStart = 0;
        for (Sclip del : svStructures.svfdel) {
            try {
                lastStart = del.start;
                if (del.used) {
                    continue;
                }
                if (del.varsCount < instance().conf.minr) {
                    continue;
                }
                List<Tuple.Tuple2<Integer, Integer>> soft = del.soft.entrySet().stream()
                        .map(entry -> Tuple.tuple(entry.getKey(), entry.getValue()))
                        .sorted(comparing(tuple -> tuple._2, reverseOrder()))
                        .collect(toList());

                int softp = soft.isEmpty() ? 0 : soft.get(0)._1;
                if (instance().conf.y) {
                    System.err.printf("%n%nWorking DEL 5' %d mate cluster cnt: %d%n", softp, del.varsCount);
                }
                if (softp != 0) {
                    if (!softClips3End.containsKey(softp)) {
                        continue;
                    }
                    Sclip scv = softClips3End.get(softp);
                    if (scv.used) {
                        continue;
                    }
                    String seq = findconseq(scv, 0);
                    if (seq.isEmpty() || seq.length() < Configuration.SEED_2) {
                        continue;
                    }
                    if (!isLoaded(region.chr, del.mstart, del.mend, reference)) {
                        referenceResource.getReference(Region.newModifiedRegion(region, del.mstart, del.mend), 300, reference);
                        Region modifiedRegion = Region.newModifiedRegion(this.region, del.mstart - 200, del.mend + 200);
                        Scope<InitialData> currentScope = new Scope<>(bam, modifiedRegion, reference, referenceResource, maxReadLength, splice,
                                variantPrinter, new InitialData(nonInsertionVariants, insertionVariants, refCoverage, softClips3End, softClips5End));
                        getMode().partialPipeline(currentScope, new DirectThreadExecutor());
                    }
                    Match match = findMatch(seq, reference, softp, 1);
                    int bp = match.basePosition;
                    String extra = match.matchedSequence;
                    if (bp == 0) {
                        continue;
                    }
                    if (!(bp - softp > 30 && isOverlap(softp, bp, del.end, del.mstart, maxReadLength))) {
                        continue;
                    }
                    bp--;
                    int dellen = bp - softp + 1;
                    Map<Integer, Character> ref = reference.referenceSequences;
                    while (ref.containsKey(bp) && ref.containsKey(softp - 1) && ref.get(bp).equals(ref.get(softp - 1))) {
                        bp--;
                        if (bp != 0) {
                            softp--;
                        }
                    }

                    final Variation variation = getVariation(nonInsertionVariants, softp, "-" + dellen);
                    variation.varsCount = 0;

                    SV sv = getSV(nonInsertionVariants, softp);
                    sv.type = "DEL";
                    sv.pairs += del.varsCount;
                    sv.splits += scv.varsCount;
                    sv.clusters++;

                    if (!(refCoverage.containsKey(softp) && refCoverage.get(softp) > del.varsCount)) {
                        refCoverage.put(softp, del.varsCount);
                    }
                    if (refCoverage.containsKey(bp) && refCoverage.get(softp) < refCoverage.get(bp)) {
                        refCoverage.put(softp, refCoverage.get(bp));
                    }

                    adjCnt(variation, scv, nonInsertionVariants.containsKey(softp) && ref.containsKey(softp)
                            ? nonInsertionVariants.get(softp).get(ref.get(softp).toString()) : null);

                    int mcnt = del.varsCount;
                    Variation tv = new Variation();
                    tv.varsCount = mcnt;
                    tv.highQualityReadsCount = mcnt;
                    tv.varsCountOnForward = mcnt / 2;
                    tv.varsCountOnReverse = mcnt - mcnt / 2;
                    tv.meanQuality = del.meanQuality * mcnt / del.varsCount;
                    tv.meanPosition = del.meanPosition * mcnt / del.varsCount;
                    tv.meanMappingQuality = del.meanMappingQuality * mcnt / del.varsCount;
                    tv.numberOfMismatches = del.numberOfMismatches * mcnt / del.varsCount;
                    adjCnt(variation, tv);

                    del.used = true;
                    markSV(softp, bp, Arrays.asList(svStructures.svrdel), maxReadLength);
                    if (instance().conf.y) {
                        System.err.printf("    Found DEL SV from 5' softclip unhappy reads: %d -%d Cnt: %d AdjCnt: %d%n", bp, dellen, del.varsCount, mcnt);
                    }
                } else { // Look within a read length
                    if (instance().conf.y) {
                        System.err.printf("%n%nWorking DEL 5' no softp mate cluster cnt: %d%n", del.varsCount);
                    }
                    if (!isLoaded(region.chr, del.mstart, del.mend, reference)) {
                        referenceResource.getReference(Region.newModifiedRegion(region, del.mstart, del.mend), 300, reference);
                    }

                    for (Map.Entry<Integer, Sclip> entry : softClips3End.entrySet()) {
                        Integer i = entry.getKey();
                        Sclip scv = entry.getValue();

                        if (scv.used) {
                            continue;
                        }
                        if (!(i >= del.end - 3 && i - del.end < 3 * maxReadLength)) {
                            continue;
                        }
                        String seq = findconseq(scv, 0);
                        if (seq.isEmpty() || seq.length() < Configuration.SEED_2) {
                            continue;
                        }
                        softp = i;
                        Match match = findMatch(seq, reference, softp, 1);
                        int bp = match.basePosition;
                        String EXTRA = match.matchedSequence;
                        if (bp == 0) {
                            match = findMatch(seq, reference, softp, 1, Configuration.SEED_2, 0);
                            bp = match.basePosition;
                            EXTRA = match.matchedSequence;
                        }
                        if (bp == 0) {
                            continue;
                        }
                        if (!(bp - softp > 30 && isOverlap(softp, bp, del.end, del.mstart, maxReadLength)))
                            continue;
                        bp--;
                        int dellen = bp - softp + 1;

                        final Variation variation = getVariation(nonInsertionVariants, softp, "-" + dellen);
                        variation.varsCount = 0;

                        SV sv = getSV(nonInsertionVariants, softp);
                        sv.type = "DEL";
                        sv.pairs += del.varsCount;
                        sv.splits += scv.varsCount;
                        sv.clusters++;

                        if (!(refCoverage.containsKey(softp) && refCoverage.get(softp) > del.varsCount)) {
                            refCoverage.put(softp, del.varsCount);
                        }
                        if (refCoverage.containsKey(bp) && refCoverage.get(softp) < refCoverage.get(bp)) {
                            refCoverage.put(softp, refCoverage.get(bp));
                        }
                        adjCnt(variation, scv);
                        int mcnt = del.varsCount;
                        Variation tv = new Variation();
                        tv.varsCount = mcnt;
                        tv.highQualityReadsCount = mcnt;
                        tv.varsCountOnForward = mcnt / 2;
                        tv.varsCountOnReverse = mcnt - mcnt / 2;
                        tv.meanQuality = del.meanQuality * mcnt / del.varsCount;
                        tv.meanPosition = del.meanPosition * mcnt / del.varsCount;
                        tv.meanMappingQuality = del.meanMappingQuality * mcnt / del.varsCount;
                        tv.numberOfMismatches = del.numberOfMismatches * mcnt / del.varsCount;
                        adjCnt(variation, tv);
                        del.used = true;
                        markSV(softp, bp, Arrays.asList(svStructures.svrdel), maxReadLength);
                        if (instance().conf.y) {
                            System.err.printf("    Found DEL SV from 5' softclip happy reads: %d -%d Cnt: %d AdjCnt: %d%n", bp, dellen, del.varsCount, mcnt);
                        }
                        break;
                    }
                }
            } catch (Exception exception) {
                printExceptionAndContinue(exception, "variant", String.valueOf(lastStart), region);
            }
        }
        for (Sclip del : svStructures.svrdel) {
            try {
                lastStart = del.start;
                if (del.used) {
                    continue;
                }
                if (del.varsCount < instance().conf.minr) {
                    continue;
                }
                List<Tuple.Tuple2<Integer, Integer>> soft = del.soft.entrySet().stream()
                        .map(entry -> Tuple.tuple(entry.getKey(), entry.getValue()))
                        .sorted(comparing(tuple -> tuple._2, reverseOrder()))
                        .collect(toList());

                int softp = soft.isEmpty() ? 0 : soft.get(0)._1;
                if (softp != 0) {
                    if (instance().conf.y) {
                        System.err.printf("%n%nWorking DEL 3' %d mate cluster cnt: %s%n", softp, del.varsCount);
                    }
                    if (!softClips5End.containsKey(softp)) {
                        continue;
                    }
                    Sclip scv = softClips5End.get(softp);
                    if (scv.used) {
                        continue;
                    }
                    String seq = findconseq(scv, 0);
                    if (seq.isEmpty() || seq.length() < Configuration.SEED_2) {
                        continue;
                    }
                    if (!isLoaded(region.chr, del.mstart, del.mend, reference)) {
                        referenceResource.getReference(Region.newModifiedRegion(region, del.mstart, del.mend), 300, reference);
                        Region modifiedRegion = Region.newModifiedRegion(this.region, del.mstart - 200, del.mend + 200);
                        Scope<InitialData> currentScope = new Scope<>(bam, modifiedRegion, reference, referenceResource, maxReadLength, splice,
                                variantPrinter, new InitialData(nonInsertionVariants, insertionVariants, refCoverage, softClips3End, softClips5End));
                        getMode().partialPipeline(currentScope, new DirectThreadExecutor());
                    }
                    Match match = findMatch(seq, reference, softp, -1);
                    int bp = match.basePosition;
                    String EXTRA = match.matchedSequence;
                    if (bp == 0) {
                        match = findMatch(seq, reference, softp, -1, Configuration.SEED_2, 0);
                        bp = match.basePosition;
                        EXTRA = match.matchedSequence;
                    }
                    if (bp == 0) {
                        continue;
                    }
                    if (!(softp - bp > 30 && isOverlap(bp, softp, del.mend, del.start, maxReadLength)))
                        continue;
                    bp++;
                    softp--;
                    int dellen = softp - bp + 1;
                    final Variation variation = getVariation(nonInsertionVariants, bp, "-" + dellen);
                    variation.varsCount = 0;

                    SV sv = getSV(nonInsertionVariants, bp);
                    sv.type = "DEL";
                    sv.pairs += del.varsCount;
                    sv.splits += scv.varsCount;
                    sv.clusters++;

                    adjCnt(variation, scv);
                    if (!(refCoverage.containsKey(bp) && refCoverage.get(bp) > del.varsCount)) {
                        refCoverage.put(bp, del.varsCount);
                    }
                    if (refCoverage.containsKey(softp) && refCoverage.get(softp) > refCoverage.get(bp)) {
                        refCoverage.put(bp, refCoverage.get(softp));
                    }
                    int mcnt = del.varsCount;
                    Variation tv = new Variation();
                    tv.varsCount = mcnt;
                    tv.highQualityReadsCount = mcnt;
                    tv.varsCountOnForward = mcnt / 2;
                    tv.varsCountOnReverse = mcnt - mcnt / 2;
                    tv.meanQuality = del.meanQuality * mcnt / del.varsCount;
                    tv.meanPosition = del.meanPosition * mcnt / del.varsCount;
                    tv.meanMappingQuality = del.meanMappingQuality * mcnt / del.varsCount;
                    tv.numberOfMismatches = del.numberOfMismatches * mcnt / del.varsCount;
                    adjCnt(variation, tv);
                    del.used = true;
                    markSV(bp, softp, Arrays.asList(svStructures.svfdel), maxReadLength);
                    if (instance().conf.y) {
                        System.err.printf("    Found DEL SV from 3' softclip unhappy reads: %d -+%d Cnt: %d AdjCnt: %d%n", bp, dellen, del.varsCount, mcnt);
                    }
                } else {
                    if (instance().conf.y) {
                        System.err.printf("%n%nWorking DEL 3' no softp mate cluster %s %d %d cnt: %d%n", region.chr, del.mstart, del.mend, del.varsCount);
                    }
                    if (!isLoaded(region.chr, del.mstart, del.mend, reference)) {
                        referenceResource.getReference(Region.newModifiedRegion(region, del.mstart, del.mend), 300, reference);
                    }
                    for (Map.Entry<Integer, Sclip> entry : softClips5End.entrySet()) {
                        int i = entry.getKey();
                        Sclip scv = entry.getValue();
                        if (scv.used) {
                            continue;
                        }
                        if (!(i <= del.start + 3 && del.start - i < 3 * maxReadLength)) {
                            continue;
                        }
                        String seq = findconseq(scv, 0);
                        if (seq.isEmpty() || seq.length() < Configuration.SEED_2) {
                            continue;
                        }
                        softp = i;
                        Match match = findMatch(seq, reference, softp, -1);
                        int bp = match.basePosition;
                        String EXTRA = match.matchedSequence;
                        if (bp == 0) {
                            match = findMatch(seq, reference, softp, -1, Configuration.SEED_2, 0);
                            bp = match.basePosition;
                            EXTRA = match.matchedSequence;
                        }
                        if (bp == 0) {
                            continue;
                        }
                        if (!(softp - bp > 30 && isOverlap(bp, softp, del.mend, del.start, maxReadLength))) {
                            continue;
                        }
                        bp++;
                        softp--;
                        int dellen = softp - bp + 1;

                        final Variation variation = getVariation(nonInsertionVariants, bp, "-" + dellen);
                        variation.varsCount = 0;

                        SV sv = getSV(nonInsertionVariants, bp);
                        sv.type = "DEL";
                        sv.pairs += del.varsCount;
                        sv.splits += scv.varsCount;
                        sv.clusters++;

                        adjCnt(variation, scv);
                        if (!refCoverage.containsKey(bp)) {
                            refCoverage.put(bp, del.varsCount);
                        }
                        if (refCoverage.containsKey(softp) && refCoverage.get(softp) > refCoverage.get(bp)) {
                            refCoverage.put(bp, refCoverage.get(softp));
                        }
                        incCnt(refCoverage, bp, scv.varsCount);

                        int mcnt = del.varsCount;
                        Variation tv = new Variation();
                        tv.varsCount = mcnt;
                        tv.highQualityReadsCount = mcnt;
                        tv.varsCountOnForward = mcnt / 2;
                        tv.varsCountOnReverse = mcnt - mcnt / 2;
                        tv.meanQuality = del.meanQuality * mcnt / del.varsCount;
                        tv.meanPosition = del.meanPosition * mcnt / del.varsCount;
                        tv.meanMappingQuality = del.meanMappingQuality * mcnt / del.varsCount;
                        tv.numberOfMismatches = del.numberOfMismatches * mcnt / del.varsCount;
                        adjCnt(variation, tv);
                        del.used = true;
                        markSV(bp, softp, Arrays.asList(svStructures.svfdel), maxReadLength);
                        if (instance().conf.y) {
                            System.err.printf("    Found DEL SV from 3' softclip happy reads: %d -%d Cnt: %d AdjCnt: %d%n", bp, dellen, del.varsCount, mcnt);
                        }
                        break;
                    }
                }
            } catch (Exception exception) {
                printExceptionAndContinue(exception, "variant", String.valueOf(lastStart), region);
            }
        }
    }

    /**
     * Find INV SV for all structural variants structures
     */
    void findINV() {
        findINVsub(svStructures.svfinv5, 1 , Side._5);
        findINVsub(svStructures.svrinv5, -1 , Side._5);
        findINVsub(svStructures.svfinv3, 1 , Side._3);
        findINVsub(svStructures.svrinv3, -1 , Side._3);
    }

    /**
     * Find INV SV
     * Would only consider those with supports from both orientations
     * (svfinv5) --&gt; | &lt;-- (svrinv5) ..... (svfinv3) --&gt; | &lt;-- (svrinv3)
     * @param svref list of INV SVs
     * @param dir direction specified in findINV()
     * @param side 3 or 5 end
     * @return created variation
     */
    Variation findINVsub(Iterable<Sclip> svref, int dir, Side side) {
        int lastStart = 0;
        // dir = 1 means 3' soft clipping
        for (Sclip inv : svref) {
            try {
                lastStart = inv.start;
                if (inv.used) {
                    continue;
                }
                if (inv.varsCount < instance().conf.minr) {
                    continue;
                }
                List<Tuple.Tuple2<Integer, Integer>> soft = inv.soft.entrySet().stream()
                        .map(entry -> Tuple.tuple(entry.getKey(), entry.getValue()))
                        .sorted(comparing(tuple -> tuple._2, reverseOrder()))
                        .collect(toList());
                int softp = soft.isEmpty() ? 0 : soft.get(0)._1;
                Map<Integer, Sclip> sclip = dir == 1 ? softClips3End : softClips5End;

                if (instance().conf.y) {
                    System.err.printf("%n%nWorking INV %d %d %s pair_cnt: %d%n", softp, dir, side, inv.varsCount);
                }
                if (!isLoaded(region.chr, inv.mstart, inv.mend, reference)) {
                    referenceResource.getReference(Region.newModifiedRegion(region, inv.mstart, inv.mend), 500, reference);
                    Region modifiedRegion = Region.newModifiedRegion(this.region, inv.mstart - 200, inv.mend + 200);
                    Scope<InitialData> currentScope = new Scope<>(bam, modifiedRegion, reference, referenceResource, maxReadLength, splice,
                            variantPrinter, new InitialData(nonInsertionVariants, insertionVariants, refCoverage, softClips3End, softClips5End));
                    getMode().partialPipeline(currentScope, new DirectThreadExecutor());
                }
                int bp = 0;
                Sclip scv = new Sclip();
                String seq = "";
                String extra = "";
                if (softp != 0) {
                    if (!sclip.containsKey(softp)) {
                        continue;
                    }
                    scv = sclip.get(softp);
                    if (scv.used) {
                        continue;
                    }
                    seq = findconseq(scv, 0);
                    if (seq.isEmpty()) {
                        continue;
                    }
                    Match matchRev = findMatchRev(seq, reference, softp, dir);
                    bp = matchRev.basePosition;
                    extra = matchRev.matchedSequence;
                    if (bp == 0) {
                        matchRev = findMatchRev(seq, reference, softp, dir, Configuration.SEED_2, 0);
                        bp = matchRev.basePosition;
                        extra = matchRev.matchedSequence;
                    }
                    if (bp == 0) {
                        continue;
                    }
                } else {
                    // Look within 100bp to see whether a soft cliping can be found but not associated with discordant pairs
                    int sp = dir == 1 ? inv.end : inv.start; // starting position
                    for (int i = 1; i <= 2 * maxReadLength; i++) {
                        int cp = sp + i * dir;
                        if (!sclip.containsKey(cp)) {
                            continue;
                        }
                        scv = sclip.get(cp);
                        if (scv.used) {
                            continue;
                        }
                        seq = findconseq(scv, 0);
                        if (seq.isEmpty()) {
                            continue;
                        }
                        Match matchRev = findMatchRev(seq, reference, cp, dir);
                        bp = matchRev.basePosition;
                        extra = matchRev.matchedSequence;
                        if (bp == 0) {
                            matchRev = findMatchRev(seq, reference, cp, dir, Configuration.SEED_2, 0);
                            bp = matchRev.basePosition;
                            extra = matchRev.matchedSequence;
                        }
                        if (bp == 0) {
                            continue;
                        }
                        softp = cp;
                        if ((dir == 1 && abs(bp - inv.mend) < Configuration.MINSVCDIST * maxReadLength)
                                || (dir == -1 && abs(bp - inv.mstart) < Configuration.MINSVCDIST * maxReadLength))
                            break;
                    }
                    if (bp == 0) {
                        continue;
                    }
                }
                if (instance().conf.y) {
                    System.err.printf("    %d %d %d %s %s pair_cnt: %d soft_cnt: %d%n", softp, bp, dir, side, seq, inv.varsCount, scv.varsCount);
                }
                if (side == Side._5) {
                    if (dir == -1) {
                        bp--;
                    }
                } else {
                    if (dir == 1) {
                        bp++;
                        if (bp != 0) {
                            softp--;
                        }
                    } else {
                        softp--;
                    }
                }
                if (side == Side._3) {
                    int tmp = bp;
                    bp = softp;
                    softp = tmp;
                }
                Map<Integer, Character> ref = reference.referenceSequences;
                if ((dir == -1 && side == Side._5) || dir == 1 && side == Side._3) {
                    while (ref.containsKey(softp) && ref.containsKey(bp)
                            && ref.get(softp) == complement(ref.get(bp))) {
                        softp++;
                        if (softp != 0) {
                            bp--;
                        }
                    }
                }
                while (ref.containsKey(softp - 1) && ref.containsKey(bp + 1)
                        && ref.get(softp - 1) == complement(ref.get(bp + 1))) {
                    softp--;
                    if (softp != 0) {
                        bp++;
                    }
                }
                if (bp > softp && bp - softp > 150 && (bp - softp) / (double) abs(inv.mlen) < 1.5) {
                    int len = bp - softp + 1;
                    String ins5 = SequenceUtil.reverseComplement(joinRef(ref, bp - Configuration.SVFLANK + 1, bp));
                    String ins3 = SequenceUtil.reverseComplement(joinRef(ref, softp, softp + Configuration.SVFLANK - 1));
                    String ins = ins5 + "<inv" + (len - 2 * Configuration.SVFLANK) + ">" + ins3;
                    if (len - 2 * Configuration.SVFLANK <= 0) {
                        ins = SequenceUtil.reverseComplement(joinRef(ref, softp, bp));
                    }
                    if (dir == 1 && !extra.isEmpty()) {
                        extra = SequenceUtil.reverseComplement(extra);
                        ins = extra + ins;
                    } else if (dir == -1 && !extra.isEmpty()) {
                        ins = ins + extra;
                    }
                    String gt = "-" + len + "^" + ins;

                    final Variation vref = getVariation(nonInsertionVariants, softp, gt);
                    inv.used = true;
                    vref.pstd = true;
                    vref.qstd = true;

                    SV sv = getSV(nonInsertionVariants, softp);
                    sv.type = "INV";
                    sv.splits += scv.varsCount;
                    sv.pairs += inv.varsCount;
                    sv.clusters++;

                    Variation vrefSoftp = dir == -1
                            ? (nonInsertionVariants.containsKey(softp) && ref.containsKey(softp) ? nonInsertionVariants.get(softp).get(ref.get(softp).toString()) : null)
                            : null;
                    adjCnt(vref, scv, vrefSoftp);
                    Map<Integer, Map<String, Integer>> dels5 = new HashMap<>();
                    Map<String, Integer> map = new HashMap<>();
                    map.put(gt, inv.varsCount);
                    dels5.put(softp, map);
                    refCoverage.put(softp, refCoverage.containsKey(softp - 1) ? refCoverage.get(softp - 1) : inv.varsCount);
                    scv.used = true;

                    VariationRealigner variationRealigner = new VariationRealigner();
                    variationRealigner.initFromScope(previousScope);
                    variationRealigner.realigndel(bams, dels5);

                    if (instance().conf.y) {
                        System.err.printf(
                                "  Found INV SV: %s %d %s BP: %d cov: %d Cnt: %d EXTRA: %s %d %d %d cnt: %d %d\t DIR: %d Side: %s%n",
                                seq, softp, gt, bp, refCoverage.get(softp), inv.varsCount, extra, inv.mstart, inv.mend, inv.mlen, scv.varsCount, (bp - softp) / abs(inv.mlen), dir, side
                        );
                    }
                    return vref;
                }
            } catch (Exception exception) {
                printExceptionAndContinue(exception, "variant", String.valueOf(lastStart), region);
            }
        }
        return null;
    }

    /**
     * Find candidate SVs on 3' and 5' ends
     */
    void findsv() {
        Map<Integer, Character> ref = reference.referenceSequences;
        List<SortPositionSclip> tmp5 = fillAndSortTmpSV(softClips5End.entrySet());
        List<SortPositionSclip> tmp3 = fillAndSortTmpSV(softClips3End.entrySet());

        int lastPosition = 0;
        for (SortPositionSclip tuple5 : tmp5) {
            try {
                int p5 = tuple5.position;
                lastPosition = p5;
                Sclip sc5v = tuple5.softClip;
                int cnt5 = tuple5.count;
                if (cnt5 < instance().conf.minr) {
                    break;
                }
                if (sc5v.used) {
                    continue;
                }
                if (SOFTP2SV.containsKey(p5) && SOFTP2SV.get(p5).get(0).used) {
                    continue;
                }
                String seq = findconseq(sc5v, 0);
                if (seq.isEmpty() || seq.length() < Configuration.SEED_2) {
                    continue;
                }
                if (instance().conf.y) {
                    System.err.printf("  Finding SV 5': %s %s cnt: %s\n", seq, p5, cnt5);
                }
                Match match = findMatch(seq, reference, p5, -1);
                int bp = match.basePosition;
                String EXTRA = match.matchedSequence;

                if (bp != 0) {
                    // candidate deletion
                    if (bp < p5) {
                        PairsData pairsData =
                                checkPairs(region.chr, bp, p5, Arrays.asList(svStructures.svfdel, svStructures.svrdel), maxReadLength);
                        int pairs = pairsData.pairs;
                        double pmean = pairsData.pmean;
                        double qmean = pairsData.qmean;
                        double Qmean = pairsData.Qmean;
                        double nm = pairsData.nm;
                        if (pairs == 0) {
                            continue;
                        }
                        p5--;
                        bp++;
                        int dellen = p5 - bp + 1;

                        final Variation vref = getVariation(nonInsertionVariants, bp, "-" + dellen);
                        vref.varsCount = 0;

                        SV sv = getSV(nonInsertionVariants, bp);
                        sv.type = "DEL";
                        sv.pairs += pairs;
                        sv.splits += cnt5;
                        sv.clusters += pairs != 0 ? 1 : 0;

                        if (!refCoverage.containsKey(bp)) {
                            refCoverage.put(bp, pairs + sc5v.varsCount);
                        }
                        if (refCoverage.containsKey(p5 + 1) && refCoverage.get(bp) < refCoverage.get(p5 + 1)) {
                            refCoverage.put(bp, refCoverage.get(p5 + 1));
                        }
                        adjCnt(vref, sc5v);
                        Variation tmp = new Variation();
                        tmp.varsCount = pairs;
                        tmp.highQualityReadsCount = pairs;
                        tmp.varsCountOnForward = (pairs / 2);
                        tmp.varsCountOnReverse = pairs - (pairs / 2);
                        tmp.meanPosition = pmean;
                        tmp.meanQuality = qmean;
                        tmp.meanMappingQuality = Qmean;
                        tmp.numberOfMismatches = nm;
                        adjCnt(vref, tmp);
                        if (instance().conf.y) {
                            System.err.println("    Finding candidate deletion 5'");
                        }
                    } else { // candidate duplication
                    }
                } else { // candidate inversion
                    Match matchRev = findMatchRev(seq, reference, p5, -1);
                    bp = matchRev.basePosition;
                    EXTRA = matchRev.matchedSequence;
                    if (bp == 0) {
                        continue;
                    }
                    if (!(abs(bp - p5) > Configuration.SVFLANK)) {
                        continue;
                    }
                    if (bp > p5) { // bp at 3' side
                    } else { // bp at 5' side
                        int temp = bp;
                        bp = p5;
                        p5 = temp;
                    }
                    bp--;
                    while (ref.containsKey(bp + 1)
                            && isHasAndEquals(complement(ref.get(bp + 1)), ref, p5 - 1)) {
                        p5--;
                        if (p5 != 0) {
                            bp++;
                        }
                    }

                    String ins5 = SequenceUtil.reverseComplement(joinRef(ref, bp - Configuration.SVFLANK + 1, bp));
                    String ins3 = SequenceUtil.reverseComplement(joinRef(ref, p5, p5 + Configuration.SVFLANK - 1));
                    int mid = bp - p5 - ins5.length() - ins3.length() + 1;

                    String vn = "-" + (bp - p5 + 1) + "^" + ins5 + "<inv" + mid + ">" + ins3 + EXTRA;
                    if (mid <= 0) {
                        String tins = SequenceUtil.reverseComplement(joinRef(ref, p5, bp));
                        vn = "-" + (bp - p5 + 1) + "^" + tins + EXTRA;
                    }

                    final Variation vref = getVariation(nonInsertionVariants, p5, vn);

                    SV sv = getSV(nonInsertionVariants, p5);
                    sv.type = "INV";
                    sv.pairs += 0;
                    sv.splits += cnt5;
                    sv.clusters += 0;

                    adjCnt(vref, sc5v);
                    incCnt(refCoverage, p5, cnt5);
                    if (refCoverage.containsKey(bp) && refCoverage.get(p5) < refCoverage.get(bp)) {
                        refCoverage.put(p5, refCoverage.get(bp));
                    }
                    if (instance().conf.y) {
                        System.err.printf("    Found INV: %d %s Cnt:%d\n", p5, vn, cnt5);
                    }
                }
            } catch (Exception exception) {
                printExceptionAndContinue(exception, "variant", String.valueOf(lastPosition), region);
            }
        }
        for (SortPositionSclip tuple3 : tmp3) {
            try {
                int p3 = tuple3.position;
                lastPosition = p3;
                Sclip sc3v = tuple3.softClip;
                int cnt3 = tuple3.count;
                if (cnt3 < instance().conf.minr) {
                    break;
                }
                if (sc3v.used) {
                    continue;
                }
                if (SOFTP2SV.containsKey(p3) && SOFTP2SV.get(p3).get(0).used) {
                    continue;
                }
                String seq = findconseq(sc3v, 0);
                if (seq.isEmpty() || seq.length() < Configuration.SEED_2) {
                    continue;
                }
                if (instance().conf.y) {
                    System.err.printf("  Finding SV 3': %s %s cnt: %s\n", seq, p3, cnt3);
                }

                Match match = findMatch(seq, reference, p3, 1);
                int bp = match.basePosition;
                String EXTRA = match.matchedSequence;
                if (bp != 0) {
                    if (bp > p3) { // candidate deletion
                        PairsData pairsData =
                                checkPairs(region.chr, p3, bp, Arrays.asList(svStructures.svfdel, svStructures.svrdel), maxReadLength);
                        int pairs = pairsData.pairs;
                        double pmean = pairsData.pmean;
                        double qmean = pairsData.qmean;
                        double Qmean = pairsData.Qmean;
                        double nm = pairsData.nm;
                        if (pairs == 0) {
                            continue;
                        }
                        int dellen = bp - p3;
                        bp--;

                        while (isHasAndEquals(bp, ref, p3 - 1)) {
                            bp--;
                            if (bp != 0) {
                                p3--;
                            }
                        }

                        final Variation vref = getVariation(nonInsertionVariants, p3, "-" + dellen);
                        vref.varsCount = 0;

                        SV sv = getSV(nonInsertionVariants, p3);
                        sv.type = "DEL";
                        sv.pairs += pairs;
                        sv.splits += cnt3;
                        sv.clusters += pairs != 0 ? 1 : 0;

                        if (!refCoverage.containsKey(p3)) {
                            refCoverage.put(p3, pairs + sc3v.varsCount);
                        }
                        if (refCoverage.containsKey(bp) && refCoverage.get(bp) < refCoverage.get(p3)) {
                            refCoverage.put(bp, refCoverage.get(p3));
                        }

                        adjCnt(vref, sc3v);
                        Variation tmp = new Variation();
                        tmp.varsCount = pairs;
                        tmp.highQualityReadsCount = pairs;
                        tmp.varsCountOnForward = (int) (pairs / 2);
                        tmp.varsCountOnReverse = pairs - (int) (pairs / 2);
                        tmp.meanPosition = pmean;
                        tmp.meanQuality = qmean;
                        tmp.meanMappingQuality = Qmean;
                        tmp.numberOfMismatches = nm;
                        adjCnt(vref, tmp);
                        if (instance().conf.y) {
                            System.err.println("    Finding candidate deletion 3'");
                        }
                    } else { // candidate duplication
                    }
                } else { // candidate inversion
                    Match matchRev = findMatchRev(seq, reference, p3, 1);
                    bp = matchRev.basePosition;
                    EXTRA = matchRev.matchedSequence;

                    if (bp == 0) {
                        continue;
                    }
                    if (abs(bp - p3) <= Configuration.SVFLANK) {
                        continue;
                    }
                    if (bp < p3) { // bp at 5' side
                        int tmp = bp;
                        bp = p3;
                        p3 = tmp;
                        p3++;
                        bp--;
                    } else { // bp at 3' side
                    }

                    while (ref.containsKey(bp + 1)
                            && isHasAndEquals(complement(ref.get(bp + 1)), ref, p3 - 1)) {
                        p3--;
                        if (p3 != 0) {
                            bp++;
                        }
                    }
                    String ins5 = SequenceUtil.reverseComplement(joinRef(ref, bp - Configuration.SVFLANK + 1, bp));
                    String ins3 = SequenceUtil.reverseComplement(joinRef(ref, p3, p3 + Configuration.SVFLANK - 1));
                    int mid = bp - p3 - 2 * Configuration.SVFLANK + 1;

                    String vn = "-" + (bp - p3 + 1) + "^" + EXTRA + ins5 + "<inv" + mid + ">" + ins3;
                    if (mid <= 0) {
                        String tins = SequenceUtil.reverseComplement(joinRef(ref, p3, bp));
                        vn = "-" + (bp - p3 + 1) + "^" + EXTRA + tins;
                    }

                    final Variation vref = getVariation(nonInsertionVariants, p3, vn);

                    SV sv = getSV(nonInsertionVariants, p3);
                    sv.type = "INV";
                    sv.pairs += 0;
                    sv.splits += cnt3;
                    sv.clusters += 0;

                    adjCnt(vref, sc3v);
                    incCnt(refCoverage, p3, cnt3);

                    if (refCoverage.containsKey(bp) && refCoverage.get(p3) < refCoverage.get(bp)) {
                        refCoverage.put(p3, refCoverage.get(bp));
                    }
                    if (instance().conf.y) {
                        System.err.printf("    Found INV: %s BP: %s Cov: %s %s %s EXTRA: %s Cnt: %s\n",
                                p3, bp, refCoverage.get(p3), refCoverage.get(bp), vn, EXTRA, cnt3);
                    }
                }
            } catch (Exception exception) {
                printExceptionAndContinue(exception, "variant", String.valueOf(lastPosition), region);
            }
        }
    }

    private List<SortPositionSclip> fillAndSortTmpSV(Set<Map.Entry<Integer, Sclip>> entries) {
        List<SortPositionSclip> tmp = new ArrayList<>();
        for (Map.Entry<Integer, Sclip> entry : entries) {
            int position = entry.getKey();
            Sclip sclip = entry.getValue();
            if (sclip.used) {
                continue;
            }
            if (position < CURSEG.start || position > CURSEG.end) {
                continue;
            }
            tmp.add(new SortPositionSclip(position, sclip, sclip.varsCount));
        }
        tmp.sort(comparing((SortPositionSclip sortPositionSclip) -> sortPositionSclip.count).reversed());
        return tmp;
    }


    /**
     * Find DEL SV with discordant pairs only
     * (svfdel) --&gt; | ............ | &lt;-- (svrdel)
     */
    void findDELdisc() {
        int MINDIST = 8 * maxReadLength; // the minimum distance between two clusters
        int lastStart = 0;
        for (Sclip del : svStructures.svfdel) {
            try {
                lastStart = del.start;
                if (del.used) {
                    continue;
                }
                if (!splice.isEmpty() && abs(del.mlen) < 250000) {
                    continue; // more stringent for RNA-Seq
                }

                if (del.varsCount < instance().conf.minr + 5) {
                    continue;
                }
                if (del.mstart <= del.end + MINDIST) {
                    continue;
                }
                if (del.meanMappingQuality / del.varsCount <= Configuration.DISCPAIRQUAL) {
                    continue;
                }
                int mlen = del.mstart - del.end - maxReadLength / (del.varsCount + 1);
                if (!(mlen > 0 && mlen > MINDIST)) {
                    continue;
                }
                int bp = del.end + (maxReadLength / (del.varsCount + 1)) / 2;
                if (del.softp != 0) {
                    bp = del.softp;
                }

                if (!reference.referenceSequences.containsKey(bp)) {
                    referenceResource.getReference(Region.newModifiedRegion(region, bp - 150, bp + 150),
                            mlen < 1000 ? mlen : 1000, reference);
                }

                final Variation vref = getVariation(nonInsertionVariants, bp, "-" + mlen);
                vref.varsCount = 0;
                SV sv = getSV(nonInsertionVariants, bp);
                sv.type = "DEL";
                sv.splits += softClips3End.containsKey(del.end + 1) ? softClips3End.get(del.end + 1).varsCount : 0;
                sv.splits += softClips5End.containsKey(del.mstart) ? softClips5End.get(del.mstart).varsCount : 0;
                sv.pairs += del.varsCount;
                sv.clusters++;

                if (instance().conf.y) {
                    System.err.printf(
                            "  Found DEL with discordant pairs only: cnt: %d BP: %d Len: %d %d-%d<->%d-%d%n",
                            del.varsCount, bp, mlen, del.start, del.end, del.mstart, del.mend
                    );
                }
                Variation tv = new Variation();
                tv.varsCount = 2 * del.varsCount;
                tv.highQualityReadsCount = 2 * del.varsCount;
                tv.varsCountOnForward = del.varsCount;
                tv.varsCountOnReverse = del.varsCount;
                tv.meanQuality = 2 * del.meanQuality;
                tv.meanPosition = 2 * del.meanPosition;
                tv.meanMappingQuality = 2 * del.meanMappingQuality;
                tv.numberOfMismatches = 2 * del.numberOfMismatches;
                adjCnt(vref, tv);
                if (!refCoverage.containsKey(bp)) {
                    refCoverage.put(bp, 2 * del.varsCount);
                }
                del.used = true;
                markSV(del.end, del.mstart, Arrays.asList(svStructures.svrdel), maxReadLength);
            } catch (Exception exception) {
                printExceptionAndContinue(exception, "variant", String.valueOf(lastStart), region);
            }
        }
        for (Sclip del : svStructures.svrdel) {
            try {
                lastStart = del.start;
                if (del.used) {
                    continue;
                }
                if (!splice.isEmpty() && abs(del.mlen) < 250000) {
                    continue; // more stringent for RNA-Seq
                }
                if (del.varsCount < instance().conf.minr + 5) {
                    continue;
                }
                if (del.start <= del.mend + MINDIST) {
                    continue;
                }
                if (del.meanMappingQuality / del.varsCount <= Configuration.DISCPAIRQUAL) {
                    continue;
                }
                int mlen = del.start - del.mend - maxReadLength / (del.varsCount + 1);
                if (!(mlen > 0 && mlen > MINDIST)) {
                    continue;
                }
                int bp = del.mend + ((maxReadLength / (del.varsCount + 1)) / 2);

                if (!reference.referenceSequences.containsKey(bp)) {
                    referenceResource.getReference(Region.newModifiedRegion(region, bp - 150, bp + 150),
                            mlen < 1000 ? mlen : 1000, reference);
                }

                final Variation ref = getVariation(nonInsertionVariants, bp, "-" + mlen);
                ref.varsCount = 0;

                SV sv = getSV(nonInsertionVariants, bp);
                sv.type = "DEL";
                sv.splits += softClips3End.containsKey(del.mend + 1) ? softClips3End.get(del.mend + 1).varsCount : 0;
                sv.splits += softClips5End.containsKey(del.start) ? softClips5End.get(del.start).varsCount : 0;
                sv.pairs += del.varsCount;
                sv.clusters += 1;

                if (instance().conf.y) {
                    System.err.printf(
                            "  Found DEL with discordant pairs only (reverse): cnt: %d BP: %d Len: %d %d-%d<->%d-%d%n",
                            del.varsCount, bp, mlen, del.start, del.end, del.mstart, del.mend
                    );
                }
                if (del.softp != 0 && softClips5End.containsKey(del.softp)) {
                    softClips5End.get(del.softp).used = true;
                }
                Variation tv = new Variation();
                tv.varsCount = 2 * del.varsCount;
                tv.highQualityReadsCount = 2 * del.varsCount;
                tv.varsCountOnForward = del.varsCount;
                tv.varsCountOnReverse = del.varsCount;
                tv.meanQuality = 2 * del.meanQuality;
                tv.meanPosition = 2 * del.meanPosition;
                tv.meanMappingQuality = 2 * del.meanMappingQuality;
                tv.numberOfMismatches = 2 * del.numberOfMismatches;
                adjCnt(ref, tv);
                if (!refCoverage.containsKey(bp)) {
                    refCoverage.put(bp, 2 * del.varsCount);
                }
                if (refCoverage.containsKey(del.start) && refCoverage.get(bp) < refCoverage.get(del.start)) {
                    refCoverage.put(bp, refCoverage.get(del.start));
                }
                del.used = true;
                referenceResource.getReference(Region.newModifiedRegion(region, del.mstart - 100, del.mend + 100), 200, reference);
                markSV(del.mend, del.start, Arrays.asList(svStructures.svfdel), maxReadLength);
            } catch (Exception exception) {
                printExceptionAndContinue(exception, "variant", String.valueOf(lastStart), region);
            }
        }
    }

    /**
     * Find INV SV with discordant pairs only
     * Would only consider those with supports from both orientations
     * #  (svfinv5) --&gt; | &lt;-- (svrinv5) ..... (svfinv3) --&gt; | &lt;-- (svrinv3)
     */
    void findINVdisc () {
        Map<Integer, Character> ref = reference.referenceSequences;
        int lastStart = 0;
        for (Sclip invf5 : svStructures.svfinv5) {
            try {
                if (invf5.used) {
                    continue;
                }
                int cnt = invf5.varsCount;
                int me = invf5.mend;
                int ms = invf5.mstart;
                int end = invf5.end;
                int start = invf5.start;
                lastStart = start;
                double nm = invf5.numberOfMismatches;
                double pmean = invf5.meanPosition;
                double qmean = invf5.meanQuality;
                double Qmean = invf5.meanMappingQuality;
                if (!(Qmean / (double) cnt > Configuration.DISCPAIRQUAL)) {
                    continue;
                }
                for (Sclip invr5 : svStructures.svrinv5) {
                    try {
                        if (invr5.used) {
                            continue;
                        }
                        int rcnt = invr5.varsCount;
                        int rstart = invr5.start;
                        lastStart = rstart;
                        int rms = invr5.mstart;
                        double rnm = invr5.numberOfMismatches;
                        double rpmean = invr5.meanPosition;
                        double rqmean = invr5.meanQuality;
                        double rQmean = invr5.meanMappingQuality;
                        if (!(rQmean / (double) rcnt > Configuration.DISCPAIRQUAL)) {
                            continue;
                        }
                        if (!(cnt + rcnt > instance().conf.minr + 5)) {
                            continue;
                        }
                        if (isOverlap(end, me, rstart, rms, maxReadLength)) {
                            int bp = abs((end + rstart) / 2);
                            int pe = abs((me + rms) / 2);
                            if (!ref.containsKey(pe)) {
                                Region modifiedRegion = Region.newModifiedRegion(region, pe - 150, pe + 150);
                                referenceResource.getReference(modifiedRegion, 300, reference);
                            }
                            int len = pe - bp + 1;

                            String ins5 = SequenceUtil.reverseComplement(joinRef(ref, bp, bp + Configuration.SVFLANK - 1));
                            String ins3 = SequenceUtil.reverseComplement(joinRef(ref, pe - Configuration.SVFLANK + 1, pe));
                            String ins = ins3 + "<inv" + (len - 2 * Configuration.SVFLANK) + ">" + ins5;
                            if (len - 2 * Configuration.SVFLANK <= 0) {
                                ins = SequenceUtil.reverseComplement(joinRef(ref, bp, pe));
                            }
                            if (instance().conf.y) {
                                System.err.printf("  Found INV with discordant pairs only 5': cnt: %d Len: %d %d-%d<->%d-%d %s\n",
                                        cnt, len, end, rstart, me, rms, ins);
                            }
                            final Variation vref = getVariation(nonInsertionVariants, bp, "-" + len + "^" + ins);

                            invf5.used = true;
                            invr5.used = true;
                            vref.pstd = true;
                            vref.qstd = true;

                            Variation tmp = new Variation();
                            tmp.varsCount = cnt + rcnt;
                            tmp.highQualityReadsCount = cnt + rcnt;
                            tmp.varsCountOnForward = cnt;
                            tmp.varsCountOnReverse = rcnt;
                            tmp.meanQuality = qmean + rqmean;
                            tmp.meanPosition = pmean + rpmean;
                            tmp.meanMappingQuality = Qmean + rQmean;
                            tmp.numberOfMismatches = nm + rnm;
                            adjCnt(vref, tmp);

                            SV sv = getSV(nonInsertionVariants, bp);
                            sv.type = "INV";
                            sv.pairs += cnt;
                            sv.splits += softClips5End.containsKey(start) ? softClips5End.get(start).varsCount : 0;
                            sv.splits += softClips5End.containsKey(ms) ? softClips5End.get(ms).varsCount : 0;
                            sv.clusters++;

                            if (!refCoverage.containsKey(bp)) {
                                refCoverage.put(bp, 2 * cnt);
                            }
                            markSV(bp, pe, Arrays.asList(svStructures.svfinv3, svStructures.svrinv3), maxReadLength);
                        }
                    } catch (Exception exception) {
                        printExceptionAndContinue(exception, "variant", String.valueOf(lastStart), region);

                    }
                }
            } catch (Exception exception) {
                printExceptionAndContinue(exception, "variant", String.valueOf(lastStart), region);

            }
        }
        for (Sclip invf3 : svStructures.svfinv3) {
            try {
                if (invf3.used) {
                    continue;
                }
                int cnt = invf3.varsCount;
                int me = invf3.mend;
                int end = invf3.end;
                double nm = invf3.numberOfMismatches;
                double pmean = invf3.meanPosition;
                double qmean = invf3.meanQuality;
                double Qmean = invf3.meanMappingQuality;
                for (Sclip invr3 : svStructures.svrinv3) {
                    try {
                        if (invr3.used) {
                            continue;
                        }
                        int rcnt = invr3.varsCount;
                        int rstart = invr3.start;
                        lastStart = rstart;
                        int rms = invr3.mstart;
                        double rnm = invr3.numberOfMismatches;
                        double rpmean = invr3.meanPosition;
                        double rqmean = invr3.meanQuality;
                        double rQmean = invr3.meanMappingQuality;
                        if (!(rQmean / (double) rcnt > Configuration.DISCPAIRQUAL)) {
                            continue;
                        }
                        if (!(cnt + rcnt > instance().conf.minr + 5)) {
                            continue;
                        }
                        if (isOverlap(me, end, rms, rstart, maxReadLength)) {
                            int pe = abs((end + rstart) / 2);
                            int bp = abs((me + rms) / 2);
                            if (!ref.containsKey(bp)) {
                                Region modifiedRegion = Region.newModifiedRegion(region, bp - 150, bp + 150);
                                referenceResource.getReference(modifiedRegion, 300, reference);
                            }
                            int len = pe - bp + 1;
                            String ins5 = SequenceUtil.reverseComplement(joinRef(ref, bp, bp + Configuration.SVFLANK - 1));
                            String ins3 = SequenceUtil.reverseComplement(joinRef(ref, pe - Configuration.SVFLANK + 1, pe));

                            String ins = ins3 + "<inv" + (len - 2 * Configuration.SVFLANK) + ">" + ins5;
                            if (len - 2 * Configuration.SVFLANK <= 0) {
                                ins = SequenceUtil.reverseComplement(joinRef(ref, bp, pe));
                            }
                            if (instance().conf.y) {
                                System.err.printf("  Found INV with discordant pairs only 3': cnt: %d Len: %d %d-%d<->%d-%d %s\n",
                                        cnt, len, me, rms, end, rstart, ins);
                            }
                            final Variation vref = getVariation(nonInsertionVariants, bp, "-" + len + "^" + ins);

                            invf3.used = true;
                            invr3.used = true;
                            vref.pstd = true;
                            vref.qstd = true;

                            Variation tmp = new Variation();
                            tmp.varsCount = cnt + rcnt;
                            tmp.highQualityReadsCount = cnt + rcnt;
                            tmp.varsCountOnForward = cnt;
                            tmp.varsCountOnReverse = rcnt;
                            tmp.meanQuality = qmean + rqmean;
                            tmp.meanPosition = pmean + rpmean;
                            tmp.meanMappingQuality = Qmean + rQmean;
                            tmp.numberOfMismatches = nm + rnm;

                            adjCnt(vref, tmp);

                            SV sv = getSV(nonInsertionVariants, bp);
                            sv.type = "INV";
                            sv.pairs += cnt;
                            sv.splits += softClips3End.containsKey(end + 1) ? softClips3End.get(end + 1).varsCount : 0;
                            sv.splits += softClips3End.containsKey(me + 1) ? softClips3End.get(me + 1).varsCount : 0;
                            sv.clusters++;

                            if (!refCoverage.containsKey(bp)) {
                                refCoverage.put(bp, 2 * cnt);
                            }
                            markSV(bp, pe, Arrays.asList(svStructures.svfinv5, svStructures.svrinv5), maxReadLength);
                        }
                    } catch (Exception exception) {
                        printExceptionAndContinue(exception, "variant", String.valueOf(lastStart), region);
                    }
                }

            } catch (Exception exception) {
                printExceptionAndContinue(exception, "variant", String.valueOf(lastStart), region);
            }
        }
    }

    /**
     * Find DUP SVs with discordant pairs only
     * (svrdup) |&lt;--.........--&gt;| (svfdup)
     */
    void findDUPdisc (){
        Map<Integer, Character> ref = reference.referenceSequences;
        int lastStart = 0;
        for (Sclip dup : svStructures.svfdup) {
            try {
                if (dup.used) {
                    continue;
                }
                int ms = dup.mstart;
                int me = dup.mend;
                int cnt = dup.varsCount;
                int end = dup.end;
                int start = dup.start;
                lastStart = start;
                double pmean = dup.meanPosition;
                double qmean = dup.meanQuality;
                double Qmean = dup.meanMappingQuality;
                double nm = dup.numberOfMismatches;
                if (!(cnt >= instance().conf.minr + 5)) {
                    continue;
                }
                if (!(Qmean / cnt > Configuration.DISCPAIRQUAL)) {
                    continue;
                }
                int mlen = end - ms + maxReadLength / cnt;
                int bp = ms - (maxReadLength / cnt) / 2;
                int pe = end;

                if (!isLoaded(region.chr, ms, me, reference)) {
                    if (!ref.containsKey(bp)) {
                        referenceResource.getReference(Region.newModifiedRegion(region, bp - 150, bp + 150), 300, reference);
                    }
                    Region modifiedRegion = Region.newModifiedRegion(this.region, ms - 200, me + 200);
                    Scope<InitialData> currentScope = new Scope<>(bam, modifiedRegion, reference, referenceResource, maxReadLength, splice,
                            variantPrinter, new InitialData(nonInsertionVariants, insertionVariants, refCoverage, softClips3End, softClips5End));
                    getMode().partialPipeline(currentScope, new DirectThreadExecutor());
                }

                int cntf = cnt;
                int cntr = cnt;
                double qmeanf = qmean;
                double qmeanr = qmean;
                double Qmeanf = Qmean;
                double Qmeanr = Qmean;
                double pmeanf = pmean;
                double pmeanr = pmean;
                double nmf = nm;
                double nmr = nm;

                if (!dup.soft.isEmpty()) {
                    List<Map.Entry<Integer, Integer>> soft = new ArrayList<>(dup.soft.entrySet());
                    soft.sort(comparing((Map.Entry<Integer, Integer> entry) -> entry.getValue(),
                            Integer::compareTo).reversed());

                    if (soft.size() > 0) {
                        pe = soft.get(0).getKey();
                    }
                    if (!softClips3End.containsKey(pe)) {
                        continue;
                    }
                    if (softClips3End.get(pe).used) {
                        continue;
                    }
                    Sclip currentSclip3 = softClips3End.get(pe);
                    cntf = currentSclip3.varsCount;
                    qmeanf = currentSclip3.meanQuality;
                    Qmeanf = currentSclip3.meanMappingQuality;
                    pmeanf = currentSclip3.meanPosition;
                    nmf = currentSclip3.numberOfMismatches;

                    String seq = findconseq(currentSclip3, 0);
                    Match match = findMatch(seq, reference, bp, 1);
                    int tbp = match.basePosition;
                    String EXTRA = match.matchedSequence;

                    if (tbp != 0 && tbp < pe) {
                        currentSclip3.used = true;
                        while (isHasAndEquals(pe - 1, ref, tbp - 1)) {
                            tbp--;
                            if (tbp != 0) {
                                pe--;
                            }
                        }
                        mlen = pe - tbp;
                        bp = tbp;
                        pe--;
                        end = pe;
                        if (softClips5End.containsKey(bp)) {
                            Sclip currentSclip5 = softClips5End.get(bp);
                            cntr = currentSclip5.varsCount;
                            qmeanr = currentSclip5.meanQuality;
                            Qmeanr = currentSclip5.meanMappingQuality;
                            pmeanr = currentSclip5.meanPosition;
                            nmr = currentSclip5.numberOfMismatches;
                        }
                    }
                }

                String ins5 = joinRef(ref, bp, bp + Configuration.SVFLANK - 1);
                String ins3 = joinRef(ref, pe - Configuration.SVFLANK + 1, pe);
                String ins = ins5 + "<dup" + (mlen - 2 * Configuration.SVFLANK) + ">" + ins3;

                final Variation vref = getVariation(insertionVariants, bp, "+" + ins);
                vref.varsCount = 0;

                SV sv = getSV(nonInsertionVariants, bp);
                sv.type = "DUP";
                sv.pairs += cnt;
                sv.splits += dup.softp != 0 && softClips3End.get(dup.softp) != null
                        ? softClips3End.get(dup.softp).varsCount : 0;
                sv.clusters++;

                if (instance().conf.y) {
                    System.err.printf("  Found DUP with discordant pairs only (forward): cnt: %d BP: %d END: %d %s Len: %d %d-%d<->%d-%d\n",
                            cnt, bp, pe, ins, mlen, start, end, ms, me);
                }
                int tcnt = cntr + cntf;

                Variation tmp = new Variation();
                tmp.varsCount = tcnt;
                tmp.extracnt = tcnt;
                tmp.highQualityReadsCount = tcnt;
                tmp.varsCountOnForward = cntf;
                tmp.varsCountOnReverse = cntr;
                tmp.meanQuality = qmeanf + qmeanr;
                tmp.meanPosition = pmeanf + pmeanr;
                tmp.meanMappingQuality = Qmeanf + Qmeanr;
                tmp.numberOfMismatches = nmf + nmr;

                adjCnt(vref, tmp);
                dup.used = true;
                if (!refCoverage.containsKey(bp)) {
                    refCoverage.put(bp, tcnt);
                }
                if (refCoverage.containsKey(end) && refCoverage.get(bp) < refCoverage.get(end)) {
                    refCoverage.put(bp, refCoverage.get(end));
                }

                Tuple.Tuple2<Integer, Integer> tuple = markDUPSV(bp, pe, Collections.singletonList(svStructures.svrdup), maxReadLength);
                int clusters = tuple._1;
                sv.clusters += clusters;
            } catch (Exception exception) {
                printExceptionAndContinue(exception, "variant", String.valueOf(lastStart), region);
            }
        }

        for (Sclip dup : svStructures.svrdup) {
            try {
                if (dup.used) {
                    continue;
                }
                int ms = dup.mstart;
                int me = dup.mend;
                int cnt = dup.varsCount;
                int end = dup.end;
                int start = dup.start;
                lastStart = start;
                double pmean = dup.meanPosition;
                double qmean = dup.meanQuality;
                double Qmean = dup.meanMappingQuality;
                double nm = dup.numberOfMismatches;
                if (cnt < instance().conf.minr + 5) {
                    continue;
                }
                if (!(Qmean / cnt > Configuration.DISCPAIRQUAL)) {
                    continue;
                }
                int mlen = me - start + maxReadLength / cnt;
                int bp = start - (maxReadLength / cnt) / 2;
                int pe = mlen + bp - 1;
                int tpe = pe;
                if (!isLoaded(region.chr, ms, me, reference)) {
                    if (!ref.containsKey(pe)) {
                        referenceResource.getReference(Region.newModifiedRegion(region, pe - 150, pe + 150), 300, reference);
                    }
                    Region modifiedRegion = Region.newModifiedRegion(this.region, ms - 200, me + 200);
                    Scope<InitialData> currentScope = new Scope<>(bam, modifiedRegion, reference, referenceResource, maxReadLength, splice,
                            variantPrinter, new InitialData(nonInsertionVariants, insertionVariants, refCoverage, softClips3End, softClips5End));
                    getMode().partialPipeline(currentScope, new DirectThreadExecutor());
                }
                int cntf = cnt;
                int cntr = cnt;
                double qmeanf = qmean;
                double qmeanr = qmean;
                double Qmeanf = Qmean;
                double Qmeanr = Qmean;
                double pmeanf = pmean;
                double pmeanr = pmean;
                double nmf = nm;
                double nmr = nm;
                if (!dup.soft.isEmpty()) {
                    List<Map.Entry<Integer, Integer>> soft = new ArrayList<>(dup.soft.entrySet());
                    soft.sort(comparing((Map.Entry<Integer, Integer> entry) -> entry.getValue(),
                            Integer::compareTo).reversed());

                    if (soft.size() > 0) {
                        bp = soft.get(0).getKey();
                    }
                    if (!softClips5End.containsKey(bp)) {
                        continue;
                    }
                    Sclip currentSclip5 = softClips5End.get(bp);
                    if (currentSclip5.used) {
                        continue;
                    }
                    cntr = currentSclip5.varsCount;
                    qmeanr = currentSclip5.meanQuality;
                    Qmeanr = currentSclip5.meanMappingQuality;
                    pmeanr = currentSclip5.meanPosition;
                    nmr = currentSclip5.numberOfMismatches;
                    String seq = findconseq(currentSclip5, 0);
                    Match match = findMatch(seq, reference, pe, -1);
                    int tbp = match.basePosition;
                    String EXTRA = match.matchedSequence;
                    if (tbp != 0 && tbp > bp) {
                        currentSclip5.used = true;
                        pe = tbp;
                        mlen = pe - bp + 1;
                        tpe = pe + 1;
                        while (isHasAndEquals(tpe, ref, bp + (tpe - pe - 1))) {
                            tpe++;
                        }
                        if (softClips3End.containsKey(tpe)) {
                            Sclip currentSclip3 = softClips3End.get(tpe);
                            cntf = currentSclip3.varsCount;
                            qmeanf = currentSclip3.meanQuality;
                            Qmeanf = currentSclip3.meanMappingQuality;
                            pmeanf = currentSclip3.meanPosition;
                            nmf = currentSclip3.numberOfMismatches;
                        }
                    }
                }

                String ins5 = joinRef(ref, bp, bp + Configuration.SVFLANK - 1);
                String ins3 = joinRef(ref, pe - Configuration.SVFLANK + 1, pe);
                String ins = ins5 + "<dup" + (mlen - 2 * Configuration.SVFLANK) + ">" + ins3;

                final Variation vref = getVariation(insertionVariants, bp, "+" + ins);
                vref.varsCount = 0;

                SV sv = getSV(nonInsertionVariants, bp);
                sv.type = "DUP";
                sv.pairs += cnt;
                sv.splits += softClips5End.containsKey(bp) ? softClips5End.get(bp).varsCount : 0;
                sv.splits += softClips3End.containsKey(tpe) ? softClips3End.get(tpe).varsCount : 0;
                sv.clusters++;

                if (instance().conf.y) {
                    System.err.printf("  Found DUP with discordant pairs only (reverse): cnt: %d BP: %d Len: %d %d-%d<->%d-%d\n",
                            cnt, bp, mlen, start, end, ms, me);
                }
                int tcnt = cntr + cntf;
                Variation tmp = new Variation();
                tmp.varsCount = tcnt;
                tmp.extracnt = tcnt;
                tmp.highQualityReadsCount = tcnt;
                tmp.varsCountOnForward = cntf;
                tmp.varsCountOnReverse = cntr;
                tmp.meanQuality = qmeanf + qmeanr;
                tmp.meanPosition = pmeanf + pmeanr;
                tmp.meanMappingQuality = Qmeanf + Qmeanr;
                tmp.numberOfMismatches = nmf + nmr;
                adjCnt(vref, tmp);

                dup.used = true;
                if (!refCoverage.containsKey(bp)) {
                    refCoverage.put(bp, tcnt);
                }
                if (refCoverage.containsKey(me) && refCoverage.get(bp) < refCoverage.get(me)) {
                    refCoverage.put(bp, refCoverage.get(me));
                }
                Tuple.Tuple2<Integer, Integer> tuple = markDUPSV(bp, pe, Collections.singletonList(svStructures.svfdup), maxReadLength);
                int clusters = tuple._1;
                sv.clusters += clusters;
            } catch (Exception exception) {
                printExceptionAndContinue(exception, "variant", String.valueOf(lastStart), region);
            }
        }
    }

    /**
     * Mark SV clusters as used
     * @param start start of the region $s
     * @param end end of the region $e
     * @param structuralVariants_sv contains list of list of Structural variants $sv
     * @param rlen read length
     * @return tuple of coverage, count of SVs overlaping with mates, number of pairs
     */
    public static Tuple.Tuple3<Integer, Integer, Integer> markSV(int start,
                                                          int end,
                                                          List<List<Sclip>> structuralVariants_sv,
                                                          int rlen) {
        int cov = 0;
        int pairs = 0;
        int cnt = 0;

        for (List<Sclip> currentSclips_sr : structuralVariants_sv) {
            for (Sclip sv_r : currentSclips_sr) {
                int start2;
                int end2;
                if (sv_r.start < sv_r.mstart) {
                    start2 = sv_r.end;
                    end2 = sv_r.mstart;
                } else {
                    start2 = sv_r.mend;
                    end2 = sv_r.start;
                }
                if (instance().conf.y) {
                    System.err.printf("   Marking SV %s %s %s %s cnt: %s\n", start, end, start2, end2, sv_r.varsCount);
                }
                if (isOverlap(start, end, start2, end2, rlen) ) {
                    if (instance().conf.y) {
                        System.err.printf("       SV %s %s %s %s cnt: %s marked\n", start, end, start2, end2, sv_r.varsCount);
                    }
                    sv_r.used = true;
                    cnt++;
                    pairs += sv_r.varsCount;
                    if (sv_r.end != sv_r.start) {
                        cov += (int) ((sv_r.varsCount * rlen)/(sv_r.end - sv_r.start)) + 1;
                    }
                }
            }
        }
        return new Tuple.Tuple3<>(cov, cnt, pairs);
    }

    /**
     * Mark DUP clusters as used
     * @param start start of the region $s
     * @param end end of the region $e
     * @param structuralVariants_sv contains list of list of Structural variants $sv
     * @param rlen read length
     * @return tuple of count of SVs overlapping with mates, number of pairs
     */
    static Tuple.Tuple2<Integer, Integer> markDUPSV(int start,
                                                    int end,
                                                    List<List<Sclip>> structuralVariants_sv,
                                                    int rlen) {
        int cov = 0;
        int pairs = 0;
        int cnt = 0;

        for (List<Sclip> currentSclips_sr : structuralVariants_sv) {
            for (Sclip sv_r : currentSclips_sr) {
                int start2;
                int end2;
                if (sv_r.start < sv_r.mstart) {
                    start2 = sv_r.start;
                    end2 = sv_r.mend;
                } else {
                    start2 = sv_r.mstart;
                    end2 = sv_r.end;
                }
                if (instance().conf.y) {
                    System.err.printf("   Marking DUP SV %s %s %s %s cnt: %s\n", start, end, start2, end2, sv_r.varsCount);
                }
                if (isOverlap(start, end, start2, end2, rlen)) {
                    if (instance().conf.y) {
                        System.err.printf("       DUP SV %s %s %s %s cnt: %s marked\n", start, end, start2, end2, sv_r.varsCount);
                    }
                    sv_r.used = true;
                    cnt++;
                    pairs += sv_r.varsCount;
                    if (sv_r.end != sv_r.start) {
                        cov += (int) ((sv_r.varsCount * rlen)/(sv_r.end - sv_r.start)) + 1;
                    }
                }
            }
        }
        return new Tuple.Tuple2<>(cnt, pairs);
    }

    /**
     * Determine overlapping of the SVs
     * @param start1 start of the first SV $s1
     * @param end1 end of the first SV $e1
     * @param start2 start of second SV
     * @param end2 end of second SV
     * @param rlen read length
     * @return true if SV overlaps with mate
     */
    static boolean isOverlap(int start1,
                          int end1,
                          int start2,
                          int end2,
                          int rlen){
        if (start1 >= end2 || start2 >= end1) {
            return false;
        }
        List<Integer> positions = Arrays.asList(start1, end1, start2, end2);
        positions.sort(Integer::compareTo);

        int ins = positions.get(2) - positions.get(1);
        if ((end1 != start1) && (end2 != start2)
                && ins/(double)(end1 - start1) > 0.75 && ins/(double)(end2 - start2) > 0.75 ) {
            return true;
        }
        if (positions.get(1) - positions.get(0) + positions.get(3) - positions.get(2) < 3 * rlen) {
            return true;
        }
        return false;
    }

    /**
     * Given a candidate SV identified by clipped reads, check whether there're mates to support it
     * @param chr chromosome name
     * @param start start of the region
     * @param end end of the region
     * @param sv list of SVs in list of clusters
     * @param maxReadLength max read length $RLEN
     * @return PairsData of pairs, mean position, mean base quality, mean mapping quality and number of mismatches
     */
    PairsData checkPairs(String chr, int start, int end, List<List<Sclip>> sv, int maxReadLength) {
        int pairs = 0;
        double pmean = 0;
        double qmean = 0;
        double Qmean = 0;
        double nm = 0;

        PairsData pairsData = new PairsData(pairs, pmean, qmean, Qmean, nm);

        for (List<Sclip> svcluster: sv) {
            for(Sclip svr : svcluster) {
                if (svr.used) {
                    continue;
                }
                int s = (svr.start + svr.end) / 2;
                int e = (svr.mstart + svr.mend) / 2;
                if (s > e) {
                    int temp = s;
                    s = e;
                    e = temp;
                }
                if (!isOverlap(start, end, s, e, maxReadLength)) {
                    continue;
                }
                if (svr.varsCount > pairs) {
                    pairsData = new PairsData(svr.varsCount, svr.meanPosition, svr.meanQuality, svr.meanMappingQuality, svr.numberOfMismatches);
                    pairs = svr.varsCount;
                }
                svr.used = true;
                if (instance().conf.y) {
                    System.err.printf("      Pair [%s:%s-%s] overlapping [%s:%s-%s] found and marked.\n", chr, s, e, chr, start, end);
                }
            }
        }
        return pairsData;
    }


    /**
     * Find matching on reverse strand (without MM defined, by default MM = 3)
     * @param seq sequence to search match with reference
     * @param REF reference in a given region
     * @param position where match is searching
     * @param dir direction
     * @return tuple of base position with extra sequence that doesn't match. If there is no mismatch,
     * returns zero and empty string
     */
    Match findMatchRev(String seq, Reference REF, int position, int dir) {
        int MM = 3;
        return findMatchRev(seq, REF, position, dir, Configuration.SEED_1, MM);
    }

    /**
     * Find matching on reverse strand (with MM defined)
     * @param seq sequence to search match with reference
     * @param REF reference in a given region
     * @param position where match is searching
     * @param dir direction
     * @param SEED seed length from configuration
     * @param MM the minimum matches for a read to be considered
     * @return tuple of base position with extra sequence that doesn't match. If there is no mismatch,
     * returns zero and empty string
     */
    Match findMatchRev(String seq, Reference REF, int position, int dir, int SEED, int MM) {
        // dir = 1 means from 3' soft clip
        if (dir == 1) {
            seq = reverse(seq);
        }
        seq = complement(seq);
        if (instance().conf.y) {
            System.err.printf("    Working MatchRev %d %s %d%n", position, seq, dir);
        }
        String extra = "";
        for (int i = seq.length() - SEED; i >= 0; i--) {
            String seed = substr(seq, i, SEED);
            if (REF.seed.containsKey(seed)) {
                List<Integer> seeds = REF.seed.get(seed);
                if (seeds.size() == 1) {
                    Integer firstSeed = seeds.get(0);
                    int bp = dir == 1 ? firstSeed + seq.length() - i - 1 : firstSeed - i;
                    if (ismatchref(seq, REF.referenceSequences, bp, -1 * dir, MM)) {
                        if (instance().conf.y) {
                            System.err.printf(
                                    "      Found SV BP (reverse): %d BP: %d SEEDpos: %d %d %s %d %s \n",
                                    dir, bp, firstSeed, position, seed, i, seq
                            );
                        }
                        return new Match(bp, extra);
                    } else {
                        // for complex indels, allowing some mismatches at the end up to 15bp or 20% length
                        String sseq = seq;
                        int eqcnt = 0;
                        for (int j = 1; j <= 15; j++) {
                            bp -= dir;
                            sseq = dir == -1 ? substr(sseq, 1) : substr(sseq, 0, -1);
                            if (dir == -1) {
                                if (isHasAndNotEquals(charAt(sseq, 0), REF.referenceSequences, bp)) {
                                    continue;
                                }
                                eqcnt++;
                                if (isHasAndNotEquals(charAt(sseq, 1), REF.referenceSequences, bp+1)) {
                                    continue;
                                }
                                extra = substr(seq, 0, j);
                            } else {
                                if (isHasAndNotEquals(charAt(sseq, -1), REF.referenceSequences, bp)) {
                                    continue;
                                }
                                eqcnt++;
                                if (isHasAndNotEquals(charAt(sseq, -2), REF.referenceSequences, bp - 1 )) {
                                    continue;
                                }
                                extra = substr(seq, -j);
                            }
                            if (eqcnt >= 3 && eqcnt / (double)j > 0.5) {
                                break;
                            }
                            if (instance().conf.y) {
                                System.err.printf(
                                        "      FoundSEED SV BP (reverse): %d BP: %d SEEDpos%d %d %s %d %s EXTRA: %s%n",
                                        dir, bp, firstSeed, position, seed, i, seq, extra);
                            }
                            if (ismatchref(sseq, REF.referenceSequences, bp, -1 * dir, 1)) {
                                return new Match(bp, extra);
                            }
                        }
                    }
                }
            }
        }
        return new Match(0, "");
    }

    /**
     * Find matching on forward strand (without MM defined, by default MM = 3)
     * @param seq sequence to search match with reference
     * @param REF reference in a given region
     * @param position where match is searching
     * @param dir direction
     * @return tuple of base position with extra sequence that doesn't match. If there is no mismatch,
     * returns zero and empty string
     */
    public Match findMatch(String seq, Reference REF, int position, int dir) {
        int MM = 3;
        return findMatch(seq, REF, position, dir, Configuration.SEED_1, MM);
    }

    /**
     * Find matching on forward strand (with MM defined)
     * @param seq sequence to search match with reference
     * @param REF reference in a given region
     * @param position where match is searching
     * @param dir direction
     * @param SEED seed length from configuration
     * @param MM the minimum matches for a read to be considered
     * @return tuple of base position with extra sequence that doesn't match. If there is no mismatch,
     * returns zero and empty string
     */
     static Match findMatch(String seq, Reference REF, int position, int dir, int SEED, int MM) {
        if (dir == -1) {
            seq = reverse(seq); // dir==-1 means 5' clip
        }
        if (instance().conf.y) {
            System.err.printf("    Working Match %d %s %d SEED: %d\n", position, seq, dir, SEED);
        }
        String extra = "";
        for(int i = seq.length() - SEED; i >= 0; i--) {
            String seed = substr(seq, i, SEED);
            if (REF.seed.containsKey(seed)) {
                List<Integer> seeds = REF.seed.get(seed);
                if (seeds.size() == 1 ) {
                    Integer firstSeed = seeds.get(0);
                    int bp = dir == 1 ? firstSeed - i : firstSeed + seq.length() - i - 1;

                    if (ismatchref(seq, REF.referenceSequences, bp, dir, MM)) {
                        int mm = dir == -1 ? -1 : 0;
                        while(isHasAndNotEquals(charAt(seq, mm), REF.referenceSequences, bp)) {
                            extra += substr(seq, mm, 1);
                            bp += dir;
                            mm += dir;
                        }
                        if (!extra.isEmpty() && dir == -1 ) {
                            extra = reverse(extra);
                        }
                        if (instance().conf.y) {
                            System.err.printf("      Found SV BP: %d BP: %d SEEDpos %s %d %s %d %s extra: %s\n",
                                    dir, bp, firstSeed, position, seed, i, seq, extra);
                        }
                        return new Match(bp, extra);
                    } else { // for complex indels, allowing some mismatches at the end up to 15bp or 20% length
                        String sseq = seq;
                        int eqcnt = 0;
                        for(int ii = 1; ii <= 15; ii++) {
                            bp += dir;
                            sseq = dir == 1 ? substr(sseq, 1) : substr(sseq, 0, -1);
                            if(dir == 1) {
                                if (isHasAndNotEquals(charAt(sseq, 0), REF.referenceSequences, bp)) {
                                    continue;
                                }
                                eqcnt++;
                                if (isHasAndNotEquals(charAt(sseq, 1), REF.referenceSequences, bp + 1)) {
                                    continue;
                                }
                                extra = substr(seq, 0, ii);
                            } else {
                                if (isHasAndNotEquals(charAt(sseq, -1), REF.referenceSequences, bp)) {
                                    continue;
                                }
                                eqcnt++;
                                if (isHasAndNotEquals(charAt(sseq, -2), REF.referenceSequences, bp - 1)) {
                                    continue;
                                }
                                extra = substr(seq, -ii);
                            }
                            if  (eqcnt >= 3 && eqcnt/(double)ii > 0.5) {
                                break;
                            }
                            if (instance().conf.y) {
                                System.err.printf("      FoundSEED SV BP: %d BP: %d SEEDpos%s %d %s %d %s EXTRA: %s\n",
                                        dir, bp, firstSeed, position, seed, i, seq, extra);
                            }
                            if (ismatchref(sseq, REF.referenceSequences, bp, dir, 1) ) {
                                return new Match(bp, extra);
                            }
                        }
                    }
                }
            }
        }
        return new Match(0, "");
    }

    /**
     * Will rescue some SNVs that are sofly clipped due to at the end of the read
     */
    private void adjSNV() {
        for (Map.Entry<Integer, Sclip> entry: softClips5End.entrySet()) {
            int position = entry.getKey();
            Sclip sclip = entry.getValue();

            if (sclip.used) {
                continue;
            }
            String seq = findconseq(sclip, 0);
            if (seq.length() > 5) {
                continue;
            }
            String bp = substr(seq, 0, 1);

            int previousPosition = position - 1;
            if (nonInsertionVariants.containsKey(previousPosition)
                    && nonInsertionVariants.get(previousPosition).containsKey(bp)) {
                if (seq.length() > 1) {
                    if (isNotEquals(reference.referenceSequences.get(position - 2), seq.charAt(1))) {
                        continue;
                    }
                }
                adjCnt(nonInsertionVariants.get(previousPosition).get(bp), sclip);
                incCnt(refCoverage, previousPosition, sclip.varsCount);
            }
        }

        for (Map.Entry<Integer, Sclip> entry: softClips3End.entrySet()) {
            int position = entry.getKey();
            Sclip sclip = entry.getValue();

            if (sclip.used) {
                continue;
            }
            String seq = findconseq(sclip, 0);
            if (seq.length() > 5) {
                continue;
            }
            String bp = substr(seq, 0, 1);
            if (nonInsertionVariants.containsKey(position)
                    && nonInsertionVariants.get(position).containsKey(bp)) {
                if (seq.length() > 1) {
                    if (isNotEquals(reference.referenceSequences.get(position + 1), seq.charAt(1))) {
                        continue;
                    }
                }
                adjCnt(nonInsertionVariants.get(position).get(bp), sclip);
                incCnt(refCoverage, position, sclip.varsCount);
            }
        }
    }

    /**
     * Output remaining soft-clipped reads that haven't been used
     */
    public void outputClipping() {
        System.err.println("5' Remaining clipping reads");
        for (Map.Entry<Integer, Sclip> entry: softClips5End.entrySet()) {
            int position_p = entry.getKey();
            Sclip sclip_sc = entry.getValue();

            if (sclip_sc.used) {
                continue;
            }
            if (sclip_sc.varsCount < instance().conf.minr) {
                continue;
            }
            String seq = findconseq(sclip_sc, 0);
            if (!seq.isEmpty() && seq.length() > Configuration.SEED_2) {
                seq = new StringBuilder(seq).reverse().toString();
                System.err.printf("  P: %s Cnt: %s Seq: %s\n", position_p, sclip_sc.varsCount, seq);
            }
        }
        System.err.println("3' Remaining clipping reads");
        for (Map.Entry<Integer, Sclip> entry: softClips3End.entrySet()) {
            int position_p = entry.getKey();
            Sclip sclip_sc = entry.getValue();
            if (sclip_sc.used) {
                continue;
            }
            if (sclip_sc.varsCount < instance().conf.minr) {
                continue;
            }
            String seq = findconseq(sclip_sc, 0);
            if (!seq.isEmpty() && seq.length() > Configuration.SEED_2) {
                System.err.printf("  P: %s Cnt: %s Seq: %s\n", position_p, sclip_sc.varsCount, seq);
            }
        }
    }

    class PairsData {
        /**
         * Number of pairs
         */
        int pairs;
        /**
         * Mean position
         */
        double pmean;
        /**
         * Mean base quality
         */
        double qmean;
        /**
         * Mean mapping quality
         */
        double Qmean;
        /**
         * number of mismatches
         */
        double nm;

        public PairsData(int pairs, double pmean, double qmean, double Qmean, double nm) {
            this.pairs = pairs;
            this.pmean = pmean;
            this.qmean = qmean;
            this.Qmean = Qmean;
            this.nm = nm;
        }
    }
}
