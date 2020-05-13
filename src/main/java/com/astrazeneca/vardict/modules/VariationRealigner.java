package com.astrazeneca.vardict.modules;

import com.astrazeneca.vardict.Configuration;
import com.astrazeneca.vardict.data.ReferenceResource;
import com.astrazeneca.vardict.collection.DirectThreadExecutor;
import com.astrazeneca.vardict.data.*;
import com.astrazeneca.vardict.Utils;
import com.astrazeneca.vardict.collection.Tuple;
import com.astrazeneca.vardict.collection.VariationMap;
import com.astrazeneca.vardict.data.scopedata.InitialData;
import com.astrazeneca.vardict.data.scopedata.RealignedVariationData;
import com.astrazeneca.vardict.data.scopedata.Scope;
import com.astrazeneca.vardict.data.scopedata.VariationData;
import com.astrazeneca.vardict.printers.VariantPrinter;
import com.astrazeneca.vardict.variations.Cluster;
import com.astrazeneca.vardict.variations.Mate;
import com.astrazeneca.vardict.variations.Sclip;
import com.astrazeneca.vardict.variations.Variation;
import htsjdk.samtools.SAMRecord;

import java.time.LocalDateTime;
import java.util.*;
import java.util.regex.Matcher;

import static com.astrazeneca.vardict.data.scopedata.GlobalReadOnlyScope.getMode;
import static com.astrazeneca.vardict.data.scopedata.GlobalReadOnlyScope.instance;
import static com.astrazeneca.vardict.collection.VariationMap.getSV;
import static com.astrazeneca.vardict.collection.VariationMap.removeSV;
import static com.astrazeneca.vardict.data.Patterns.*;
import static com.astrazeneca.vardict.modules.CigarParser.getAlignedLength;
import static com.astrazeneca.vardict.modules.RecordPreprocessor.getChrName;
import static com.astrazeneca.vardict.data.ReferenceResource.isLoaded;
import static com.astrazeneca.vardict.modules.StructuralVariantsProcessor.*;
import static com.astrazeneca.vardict.variations.VariationUtils.*;
import static com.astrazeneca.vardict.Utils.*;
import static java.lang.String.format;
import static java.util.Collections.singletonMap;
import static java.util.Comparator.comparing;

/**
 * The step of pipeline which try to realign variations: softclips, indels, long insertions.
 */
public class VariationRealigner implements Module<VariationData, RealignedVariationData>  {
    private Map<Integer, List<Sclip>> SOFTP2SV = new HashMap<>();

    private Map<Integer, VariationMap<String, Variation>> nonInsertionVariants;
    private Map<Integer, VariationMap<String, Variation>> insertionVariants;
    private Map<Integer, Map<String, Integer>> positionToInsertionCount;
    private Map<Integer, Map<String, Integer>> positionToDeletionsCount;
    private Map<Integer, Integer> refCoverage;
    private Map<Integer, Sclip> softClips5End;
    private Map<Integer, Sclip> softClips3End;
    private ReferenceResource referenceResource;
    private Reference reference;
    private Region region;
    private Set<String> splice;
    private String chr;
    private int maxReadLength;
    private double duprate;
    private String[] bams;
    private String bam;
    private Map<Integer, Map<String, Integer>> mnp;
    private SVStructures svStructures;
    private VariantPrinter variantPrinter;

    /**
     * Starts the filtering of prepared SV structures, adjusting counts of MNPs and realining of indels and softclips.
     * @param scope variation data from the CigarParser that can be realigned. Contains maps of insertion
     *              and non-insertion variations.
     * @return realigned variations maps
     */
    @Override
    public Scope<RealignedVariationData> process(Scope<VariationData> scope) {
        initFromScope(scope);
        CurrentSegment CURSEG = new CurrentSegment(region.chr, region.start, region.end);

        if (!instance().conf.disableSV) {
            filterAllSVStructures();
        }

        adjustMNP();

        if (instance().conf.y) {
            System.err.println("TIME: Start realign: " + LocalDateTime.now());
        }

        if (instance().conf.performLocalRealignment) {
            realignIndels();
        }

        return new Scope<>(
                scope,
                new RealignedVariationData(nonInsertionVariants, insertionVariants, softClips3End, softClips5End,
                        refCoverage, maxReadLength, duprate, svStructures, CURSEG, SOFTP2SV, scope));
    }

    public void initFromScope(Scope<VariationData> scope) {
        this.region = scope.region;
        this.nonInsertionVariants = scope.data.nonInsertionVariants;
        this.insertionVariants = scope.data.insertionVariants;
        this.positionToInsertionCount = scope.data.positionToInsertionCount;
        this.positionToDeletionsCount = scope.data.positionToDeletionsCount;
        this.refCoverage = scope.data.refCoverage;
        this.softClips5End = scope.data.softClips5End;
        this.softClips3End = scope.data.softClips3End;
        this.reference = scope.regionRef;
        this.referenceResource = scope.referenceResource;
        this.chr = getChrName(scope.region);
        this.maxReadLength = scope.maxReadLength;
        this.bams = scope.bam != null ? scope.bam.split(":") : null;
        this.bam = scope.bam;
        this.mnp = scope.data.mnp;
        this.splice = scope.splice;
        this.svStructures = scope.data.svStructures;
        this.duprate = scope.data.duprate;
        this.variantPrinter = scope.out;
    }

    private static final Comparator<SortPositionSclip> COMP2 = new Comparator<SortPositionSclip>() {
        @Override
        public int compare(SortPositionSclip o1, SortPositionSclip o2) {
            int f = Integer.compare(o2.softClip.varsCount, o1.softClip.varsCount);
            if (f != 0)
                return f;
            return Integer.compare(o1.position, o2.position);
        }
    };

    private static final Comparator<SortPositionSclip> COMP3 = new Comparator<SortPositionSclip>() {
        @Override
        public int compare(SortPositionSclip o1, SortPositionSclip o2) {
            int f = Integer.compare(o2.count, o1.count);
            if (f != 0)
                return f;
            return Integer.compare(o1.position, o2.position);
        }
    };

    /**
     * Filter all possible structural variants structures.
     * All SVs clusters will be checked for forward and reverse direction.
     */
    public void filterAllSVStructures() {
        filterSV(svStructures.svfinv3);
        filterSV(svStructures.svrinv3);
        filterSV(svStructures.svfinv5);
        filterSV(svStructures.svrinv5);
        filterSV(svStructures.svfdel);
        filterSV(svStructures.svrdel);
        filterSV(svStructures.svfdup);
        filterSV(svStructures.svrdup);

        for (Map.Entry<String, List<Sclip>> svv : svStructures.svffus.entrySet()) {
            filterSV(svv.getValue());
        }
        for (Map.Entry<String, List<Sclip>> svv : svStructures.svrfus.entrySet()) {
            filterSV(svv.getValue());
        }
        for (Map.Entry<Integer, List<Sclip>> entry : SOFTP2SV.entrySet()) {
            Integer key = entry.getKey();
            List<Sclip> sclips = entry.getValue();
            sclips.sort(comparing((Sclip sclip) -> sclip.varsCount).reversed());
            SOFTP2SV.put(key, sclips);
        }
    }

    /**
     * Filter possible SVs by checking created clusters, Updates data of possible SV from cluster data.
     * @param svList_sva contains list of SVs to filter
     */
    void filterSV(List<Sclip> svList_sva) {
        for (Sclip sv: svList_sva) {
            try {
                Cluster cluster = checkCluster(sv.mates, maxReadLength);

                if (cluster.mateStart_ms != 0) {
                    sv.mstart = cluster.mateStart_ms;
                    sv.mend = cluster.mateEnd_me;
                    sv.varsCount = cluster.cnt;
                    sv.mlen = cluster.mateLength_mlen;
                    sv.start = cluster.start_s;
                    sv.end = cluster.end_e;
                    sv.meanPosition = cluster.pmean_rp;
                    sv.meanQuality = cluster.qmean_q;
                    sv.meanMappingQuality = cluster.Qmean_Q;
                    sv.numberOfMismatches = cluster.nm;
                } else {
                    sv.used = true;
                }
                // Too many unhappy mates are false positive
                if (sv.disc != 0 && sv.varsCount / (double) sv.disc < 0.5) {
                    if (!(sv.varsCount / (double) sv.disc >= 0.35 && sv.varsCount >= 5)) {
                        sv.used = true;
                    }
                }

                List<Map.Entry<Integer, Integer>> soft = new ArrayList<>(sv.soft.entrySet());
                soft.sort(comparing((Map.Entry<Integer, Integer> entry) -> entry.getValue(), Integer::compareTo).reversed());

                sv.softp = soft.size() > 0 ? soft.get(0).getKey() : 0;
                if (sv.softp != 0) {
                    List<Sclip> sclips = SOFTP2SV.getOrDefault(sv.softp, new ArrayList<>());
                    sclips.add(sv);
                    SOFTP2SV.put(sv.softp, sclips);
                }

                if (instance().conf.y) {
                    System.err.printf("SV cluster: %s %s %s %s Cnt: %s Discordant Cnt: %s Softp: %s Used: %s\n",
                            cluster.start_s, cluster.end_e,
                            cluster.mateStart_ms, cluster.mateEnd_me, sv.varsCount, sv.disc, sv.softp, sv.used);
                }
            } catch (Exception exception) {
                printExceptionAndContinue(exception, "cluster", String.valueOf(sv.start), region);
            }
        }
    }

    /**
     * Check cluster for start and end of SV mates
     * @param mates list of mates for current SV
     * @param rlen read length
     * @return Cluster with updated starts and ends
     */
    Cluster checkCluster(List<Mate> mates, int rlen) {
        mates.sort(comparing(mate -> mate.mateStart_ms, Integer::compareTo));

        List<Cluster> clusters = new ArrayList<>();
        Mate firstMate = mates.get(0);
        clusters.add(new Cluster(0, firstMate.mateStart_ms, firstMate.mateEnd_me, firstMate.start_s, firstMate.end_e));

        int cur = 0;
        for (Mate mate_m : mates) {
            Cluster currentCluster = clusters.get(cur);
            if (mate_m.mateStart_ms - currentCluster.mateEnd_me > Configuration.MINSVCDIST * rlen) {
                cur++;
                clusters.add(cur, new Cluster(0, mate_m.mateStart_ms, mate_m.mateEnd_me, mate_m.start_s, mate_m.end_e));
                currentCluster = clusters.get(cur);
            }

            currentCluster.cnt++;
            currentCluster.mateLength_mlen += mate_m.mateLength_mlen;

            if (mate_m.mateEnd_me > currentCluster.mateEnd_me) {
                currentCluster.mateEnd_me = mate_m.mateEnd_me;
            }
            if (mate_m.start_s < currentCluster.start_s) {
                currentCluster.start_s = mate_m.start_s;
            }
            if (mate_m.end_e > currentCluster.end_e) {
                currentCluster.end_e = mate_m.end_e;
            }
            currentCluster.pmean_rp += mate_m.pmean_rp;
            currentCluster.qmean_q += mate_m.qmean_q;
            currentCluster.Qmean_Q += mate_m.Qmean_Q;
            currentCluster.nm += mate_m.nm;
        }
        clusters.sort(comparing((Cluster cluster) -> cluster.cnt).reversed());

        if (instance().conf.y) {
            System.err.print("Clusters; ");
            clusters.forEach(cluster -> System.err.print(join("; ", cluster.cnt, cluster.start_s, cluster.end_e,
                    cluster.mateStart_ms, cluster.mateEnd_me, "")));
            System.err.println(join("; ", "; out of", mates.size()));
        }
        Cluster firstCluster = clusters.get(0);
        //TODO: refactor this to factory method?
        return firstCluster.cnt / (double) mates.size() >= 0.60
                ? new Cluster(firstCluster.mateStart_ms, firstCluster.mateEnd_me, firstCluster.cnt,
                firstCluster.mateLength_mlen/firstCluster.cnt, firstCluster.start_s, firstCluster.end_e,
                firstCluster.pmean_rp, firstCluster.qmean_q, firstCluster.Qmean_Q, firstCluster.nm)
                : new Cluster(0,0,0,0,0,0,0,0.0,0,0);
    }

    /**
     * Adjust MNP when there're breakpoints within MNP (multi-nucleotide polymorphism)
     */
    public void adjustMNP() {
        List<SortPositionDescription> tmp = fillAndSortTmp(mnp);
        for (SortPositionDescription tpl : tmp) {
            int lastPosition = 0;
            try {
                final Integer position = tpl.position;
                lastPosition = position;

                final String vn = tpl.descriptionString;
                final Map<String, Variation> varsOnPosition = nonInsertionVariants.get(position);
                if (varsOnPosition == null) {
                    continue;
                }
                final Variation vref = varsOnPosition.get(vn);
                if (vref == null) { // The variant is likely already been used by indel realignment
                    continue;
                }
                if (instance().conf.y) {
                    System.err.printf("  AdjMnt: %d %s %d\n", position, vn, vref.varsCount);
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
                        Variation tref = varsOnPosition.get(left);
                        if (tref != null) {
                            if (tref.varsCount <= 0) {
                                continue;
                            }
                            if (tref.varsCount < vref.varsCount && tref.meanPosition / tref.varsCount <= i + 1) {
                                if (instance().conf.y) {
                                    System.err.printf("    AdjMnt Left: %s %s Left: %s Cnt: %s\n", position, vn, left, tref.varsCount);
                                }
                                adjCnt(vref, tref);
                                varsOnPosition.remove(left);
                            }
                        }
                    }
                    if (nonInsertionVariants.containsKey(position + i + 1)) {
                        Variation tref = nonInsertionVariants.get(position + i + 1).get(right);
                        if (tref != null) {
                            if (tref.varsCount < 0) {
                                continue;
                            }
                            // #&& tref.pmean / tref.cnt <= mnt.length() - i - 1)
                            if (tref.varsCount < vref.varsCount) {
                                if (instance().conf.y) {
                                    System.err.printf("    AdjMnt Right: %s %s Right: %s Cnt: %s\n", position, vn, right, tref.varsCount);
                                }
                                adjCnt(vref, tref);
                                incCnt(refCoverage, position, tref.varsCount);
                                nonInsertionVariants.get(position + i + 1).remove(right);
                            }
                        }
                    }
                }
                if (softClips3End.containsKey(position)) {
                    final Sclip sc3v = softClips3End.get(position);
                    if (!sc3v.used) {
                        final String seq = findconseq(sc3v, 0);
                        if (seq.startsWith(mnt)) {
                            if (seq.length() == mnt.length()
                                    || ismatchref(seq.substring(mnt.length()), reference.referenceSequences, position + mnt.length(), 1)) {
                                adjCnt(nonInsertionVariants.get(position).get(vn), sc3v);
                                incCnt(refCoverage, position, sc3v.varsCount);
                                sc3v.used = true;
                            }
                        }
                    }
                }
                if (softClips5End.containsKey(position + mnt.length())) {
                    final Sclip sc5v = softClips5End.get(position + mnt.length());
                    if (!sc5v.used) {
                        String seq = findconseq(sc5v, 0);
                        if (!seq.isEmpty() && seq.length() >= mnt.length()) {
                            seq = new StringBuffer(seq).reverse().toString();
                            if (seq.endsWith(mnt)) {
                                if (seq.length() == mnt.length()
                                        || ismatchref(seq.substring(0, seq.length() - mnt.length()), reference.referenceSequences, position - 1, -1)) {
                                    adjCnt(nonInsertionVariants.get(position).get(vn), sc5v);
                                    incCnt(refCoverage, position, sc5v.varsCount);
                                    sc5v.used = true;
                                }
                            }
                        }

                    }
                }
            } catch (Exception exception) {
                printExceptionAndContinue(exception, "MNP", String.valueOf(lastPosition), region);
            }
        }
    }

    public void realignIndels()  {
        if (instance().conf.y)
            System.err.println("Start Realigndel");
        realigndel(bams, positionToDeletionsCount);
        if (instance().conf.y)
            System.err.println("Start Realignins");
        realignins(positionToInsertionCount);
        if (instance().conf.y)
            System.err.println("Start Realignlgdel");
        realignlgdel(svStructures.svfdel, svStructures.svrdel);
        if (instance().conf.y)
            System.err.println("Start Realignlgins30");
        realignlgins30();
        if (instance().conf.y)
            System.err.println("Start Realignlgins");
        realignlgins(svStructures.svfdup, svStructures.svrdup);
    }

    /**
     * Realign deletions if already present in alignment
     * @param positionToDeletionsCount deletion variants count on positions
     * @param bamsParameter BAM file list (can be null in few cases of running realigndel)
     */
    public void realigndel(String[] bamsParameter, Map<Integer, Map<String, Integer>> positionToDeletionsCount) {
        Map<Integer, Character> ref = reference.referenceSequences;
        String[] bams;
        if (bamsParameter == null) {
            bams = null;
        } else {
            bams = this.bams;
        }
        // In perl it doesn't commented, but it doesn't used
        // int longmm = 3; //Longest continued mismatches typical aligned at the end
        List<SortPositionDescription> tmp = fillAndSortTmp(positionToDeletionsCount);
        int lastPosition = 0;
        for (SortPositionDescription tpl : tmp) {
            try {
                Integer p = tpl.position;
                lastPosition = p;
                String vn = tpl.descriptionString;
                Integer dcnt = tpl.count;
                if (instance().conf.y) {
                    System.err.printf("  Realigndel for: %s %s %s cov: %s\n", p, vn, dcnt, refCoverage.get(p));
                }
                final Variation vref = getVariation(nonInsertionVariants, p, vn);
                int dellen = 0;
                Matcher mtch = BEGIN_MINUS_NUMBER.matcher(vn);
                if (mtch.find()) {
                    dellen = toInt(mtch.group(1));
                }
                mtch = UP_NUMBER_END.matcher(vn);
                if (mtch.find()) {
                    dellen += toInt(mtch.group(1));
                }
                String extrains = "";
                String extra = "";
                String inv5 = "";
                String inv3 = "";

                if ((mtch = MINUS_NUMBER_ATGNC_SV_ATGNC_END.matcher(vn)).find()) {
                    inv5 = mtch.group(1);
                    inv3 = mtch.group(2);
                } else if ((mtch = BEGIN_MINUS_NUMBER_ANY.matcher(vn)).find()) {
                    extra = mtch.group(1).replaceAll("\\^|&|#", "");
                    if ((mtch = CARET_ATGNC.matcher(vn)).find()) {
                        extrains = mtch.group(1);
                    }
                }

                int wustart = (p - 200) > 1 ? (p - 200) : 1;
                String wupseq = joinRef(ref, wustart, p - 1) + extra; // 5' flanking seq
                if (!inv3.isEmpty()) {
                    wupseq = inv3;
                }
                int sanend = (instance().chrLengths.get(chr) != null && p + 200 > instance().chrLengths.get(chr))
                        ? instance().chrLengths.get(chr)
                        : p + 200;

                // 3' flanking seq
                String sanpseq = extra + joinRef(ref, p + dellen + extra.length() - extrains.length(), sanend);
                if (!inv5.isEmpty()) {
                    sanpseq = inv5;
                }

                // mismatches, mismatch positions, 5 or 3 ends
                MismatchResult r3 = findMM3(ref, p, sanpseq);
                MismatchResult r5 = findMM5(ref, p + dellen + extra.length() - extrains.length() - 1, wupseq);

                List<Mismatch> mm3 = r3.getMismatches();
                List<Integer> sc3p = r3.getScp();
                int nm3 = r3.getNm();
                int misp3 = r3.getMisp();
                String misnt3 = r3.getMisnt();

                List<Mismatch> mm5 = r5.getMismatches();
                List<Integer> sc5p = r5.getScp();
                int nm5 = r5.getNm();
                int misp5 = r5.getMisp();
                String misnt5 = r5.getMisnt();
                if (instance().conf.y) {
                    System.err.printf("  Mismatches: misp3: %s-%s misp5: %s-%s sclip3: %s sclip5: %s\n",
                            misp3, misnt3, misp5, misnt5, Utils.toString(sc3p), Utils.toString(sc5p));
                }

                List<Mismatch> mmm = new ArrayList<>(mm3); //$mmm
                mmm.addAll(mm5);
                for (Mismatch mismatch : mmm) {
                    String mm = mismatch.mismatchSequence;
                    Integer mp = mismatch.mismatchPosition;
                    Integer me = mismatch.end;
                    if (mm.length() > 1) {
                        mm = mm.charAt(0) + "&" + mm.substring(1);
                    }
                    if (nonInsertionVariants.get(mp) == null) {
                        continue;
                    }
                    Variation tv = nonInsertionVariants.get(mp).get(mm);
                    if (tv == null) {
                        continue;
                    }
                    if (tv.varsCount == 0) {
                        continue;
                    }
                    if (tv.meanQuality / tv.varsCount < instance().conf.goodq) {
                        continue;
                    }

                    if (tv.meanPosition / tv.varsCount > (me == 3 ? nm3 + 4 : nm5 + 4)) {
                        continue;
                    }
                    if (tv.varsCount >= dcnt + dellen || tv.varsCount / dcnt >= 8) {
                        continue;
                    }
                    if (instance().conf.y) {
                        System.err.printf("  Realigndel Adj: %s %s %s %s %s %s %s %s cov: %s\n",
                                mm, mp, me, nm3, nm5, p, tv.varsCount, tv.meanQuality, refCoverage.get(p));
                    }
                    // Adjust ref cnt so that AF won't > 1
                    if (mp > p && me == 5) {
                        double f = tv.meanPosition != 0 ? (mp - p) / (tv.meanPosition / (double) tv.varsCount) : 1;
                        if (f > 1) {
                            f = 1;
                        }
                        incCnt(refCoverage, p, (int) (tv.varsCount * f));
                        adjRefCnt(tv, getVariationMaybe(nonInsertionVariants, p, ref.get(p)), dellen);
                    }
                    Variation lref = (mp > p && me == 3)
                            ? (nonInsertionVariants.containsKey(p) && nonInsertionVariants.get(p).containsKey(ref.get(p).toString())
                            ? nonInsertionVariants.get(p).get(ref.get(p).toString())
                            : null)
                            : null;
                    adjCnt(vref, tv, lref);
                    nonInsertionVariants.get(mp).remove(mm);
                    if (nonInsertionVariants.get(mp).isEmpty()) {
                        nonInsertionVariants.remove(mp);
                    }
                    if (instance().conf.y) {
                        System.err.printf("  Realigndel AdjA: %s %s %s %s %s %s %s %s cov: %s\n",
                                mm, mp, me, nm3, nm5, p, tv.varsCount, tv.meanQuality, refCoverage.get(p));
                    }
                }
                if (misp3 != 0 && mm3.size() == 1 && nonInsertionVariants.containsKey(misp3)
                        && nonInsertionVariants.get(misp3).containsKey(misnt3)
                        && nonInsertionVariants.get(misp3).get(misnt3).varsCount < dcnt) {
                    nonInsertionVariants.get(misp3).remove(misnt3);
                }
                if (misp5 != 0 && mm5.size() == 1 && nonInsertionVariants.containsKey(misp5)
                        && nonInsertionVariants.get(misp5).containsKey(misnt5)
                        && nonInsertionVariants.get(misp5).get(misnt5).varsCount < dcnt) {
                    nonInsertionVariants.get(misp5).remove(misnt5);
                }

                for (Integer sc5pp : sc5p) {
                    if (softClips5End.containsKey(sc5pp) && !softClips5End.get(sc5pp).used) {
                        Sclip tv = softClips5End.get(sc5pp);
                        String seq = findconseq(tv, 0);
                        //Make sure a couple of bogus mapping won't scoop up several fold soft-clip reads
                        if (dcnt <= 2 && tv.varsCount / dcnt > 5) {
                            continue;
                        }
                        if (instance().conf.y) {
                            System.err.printf("  Realigndel 5: %s %s seq: '%s' Wuseq: %s cnt: %s %s %s %s cov: %s\n",
                                    p, sc5pp, seq, new StringBuilder(wupseq).reverse(), tv.varsCount, dcnt, vn, p, refCoverage.get(p));
                        }
                        if (!seq.isEmpty() && ismatch(seq, wupseq, -1)) {
                            if (sc5pp > p) {
                                incCnt(refCoverage, p, tv.varsCount);
                            }
                            adjCnt(vref, tv);
                            softClips5End.get(sc5pp).used = true;
                            if (instance().conf.y) {
                                System.err.printf("  Realigndel 5: %s %s %s %s %s %s %s %s used cov: %s\n",
                                        p, sc5pp, seq, new StringBuilder(wupseq).reverse(), tv.varsCount, dcnt, vn, p, refCoverage.get(p));
                            }
                        }
                    }
                }

                for (Integer sc3pp : sc3p) {
                    if (softClips3End.containsKey(sc3pp) && !softClips3End.get(sc3pp).used) {
                        Sclip tv = softClips3End.get(sc3pp);
                        String seq = findconseq(tv, 0);
                        //Make sure a couple of bogus mapping won't scoop up several fold soft-clip reads
                        if (dcnt <= 2 && tv.varsCount / dcnt > 5) {
                            continue;
                        }
                        if (instance().conf.y) {
                            System.err.printf("  Realigndel 3: %s %s seq '%s' Sanseq: %s cnt: %s %s %s %s %s %s\n",
                                    p, sc3pp, seq, sanpseq, tv.varsCount, dcnt, vn, p, dellen, substr(sanpseq, sc3pp - p));
                        }
                        if (!seq.isEmpty() && ismatch(seq, substr(sanpseq, sc3pp - p), 1)) {
                            if (instance().conf.y) {
                                System.err.printf("  Realigndel 3: %s %s %s %s %s %s %s %s used\n",
                                        p, sc3pp, seq, sanpseq, tv.varsCount, dcnt, vn, p);
                            }
                            if (sc3pp <= p) {
                                incCnt(refCoverage, p, tv.varsCount);
                            }
                            Variation lref = sc3pp <= p ? null : getVariationMaybe(nonInsertionVariants, p, ref.get(p));
                            adjCnt(vref, tv, lref);
                            softClips3End.get(sc3pp).used = true;
                        }
                    }
                }
                // In perl it is commented too
                // int pe = p + dellen + extra.length() + compm.length();
                int pe = p + dellen + extra.length() - extrains.length();
                Variation h = getVariationMaybe(nonInsertionVariants, p, ref.get(p));
                // taking the size of gap into account
                if (bams != null && bams.length > 0
                        && pe - p >= 5
                        && pe - p < maxReadLength - 10
                        && h != null && h.varsCount != 0
                        && noPassingReads(chr, p, pe, bams)
                        && vref.varsCount > 2 * h.varsCount * (1 - (pe - p) / (double) maxReadLength)) {
                    adjCnt(vref, h, h);
                }
            } catch (Exception exception) {
                printExceptionAndContinue(exception, "variant", String.valueOf(lastPosition), region);
            }
        }

        for (int i = tmp.size() - 1; i > 0; i--) {
            try {
                SortPositionDescription tpl = tmp.get(i);
                int p = tpl.position;
                lastPosition = p;
                String vn = tpl.descriptionString;
                if (!nonInsertionVariants.containsKey(p)) {
                    continue;
                }
                Variation vref = nonInsertionVariants.get(p).get(vn);
                if (vref == null) {
                    continue;
                }
                Matcher matcher = MINUS_NUMBER_AMP_ATGCs_END.matcher(vn);
                if (matcher.find()) {
                    String tn = matcher.group(1);
                    Variation tref = nonInsertionVariants.get(p).get(tn);
                    if (tref != null) {
                        if (vref.varsCount < tref.varsCount) {
                            adjCnt(tref, vref);
                            nonInsertionVariants.get(p).remove(vn);
                        }
                    }
                }
            } catch (Exception exception) {
                printExceptionAndContinue(exception, "variant", String.valueOf(lastPosition), region);
            }
        }
    }

    /**
     * Realign insertions
     * @param positionToInsertionCount insertion variants count on positions
     * @return insertion sequence
     */
    public String realignins(Map<Integer, Map<String, Integer>> positionToInsertionCount) {
        Map<Integer, Character> ref = reference.referenceSequences;
        List<SortPositionDescription> tmp = fillAndSortTmp(positionToInsertionCount);
        String NEWINS = "";
        int lastPosition = 0;
        for (SortPositionDescription tpl : tmp) {
            try {
                Integer position = tpl.position;
                lastPosition = position;
                String vn = tpl.descriptionString;
                Integer insertionCount = tpl.count;
                if (instance().conf.y) {
                    System.err.println(format("  Realign Ins: %s %s %s", position, vn, insertionCount));
                }
                String insert;
                Matcher mtch = BEGIN_PLUS_ATGC.matcher(vn);
                if (mtch.find()) {
                    insert = mtch.group(1);
                } else {
                    continue;
                }
                String ins3 = "";
                int inslen = insert.length();

                mtch = DUP_NUM_ATGC.matcher(vn);
                if (mtch.find()) {
                    ins3 = mtch.group(2);
                    inslen += toInt(mtch.group(1)) + ins3.length();
                }
                String extra = "";
                mtch = AMP_ATGC.matcher(vn);
                if (mtch.find()) {
                    extra = mtch.group(1);
                }
                String compm = ""; // the match part for a complex variant
                mtch = HASH_ATGC.matcher(vn);
                if (mtch.find()) {
                    compm = mtch.group(1);
                }

                // In perl it doesn't commented, but not used
                String newins = ""; // the adjacent insertion
                mtch = CARET_ATGC_END.matcher(vn);
                if (mtch.find()) {
                    newins = mtch.group(1);
                }

                int newdel = 0; // the adjacent deletion
                mtch = UP_NUMBER_END.matcher(vn);
                if (mtch.find()) {
                    newdel = toInt(mtch.group(1));
                }
                String tn = vn.replaceFirst("^\\+", "")
                        .replaceFirst("&", "")
                        .replaceFirst("#", "")
                        .replaceFirst("\\^\\d+$", "")
                        .replaceFirst("\\^", "");

                int wustart = position - 150 > 1 ? (position - 150) : 1;
                String wupseq = joinRef(ref, wustart, position) + tn; // 5prime flanking seq
            /*
            5' flanking region is a region of DNA that is adjacent to the 5' end of the gene.
            The 5' flanking region contains the promoter, and may contain enhancers or other protein binding sites.
            It is the region of DNA that is not transcribed into RNA.
             */
                Integer tend = instance().chrLengths.get(chr);
                int sanend = position + vn.length() + 100;
                if (tend != null && tend < sanend) {
                    sanend = tend;
                }
                // 3prime flanking seq
                String sanpseq = "";
            /*
            3' flanking region is a region of DNA which is NOT copied into the mature mRNA, but which is present adjacent
            to 3' end of the gene. It was originally thought that the 3' flanking DNA was not transcribed at all,
            but it was discovered to be transcribed into RNA, but quickly removed during processing of the primary
            transcript to form the mature mRNA. The 3' flanking region often contains sequences which affect the
            formation of the 3' end of the message. It may also contain enhancers or other sites to which proteins may bind.
             */
                MismatchResult findMM3;

                if (!ins3.isEmpty()) {
                    int p3 = position + inslen - ins3.length() + Configuration.SVFLANK;
                    if (ins3.length() > Configuration.SVFLANK) {
                        sanpseq = substr(ins3, Configuration.SVFLANK - ins3.length());
                    }
                    sanpseq += joinRef(ref, position + 1, position + 101);
                    findMM3 = findMM3(ref, p3 + 1, sanpseq);
                } else {
                    sanpseq = tn + joinRef(ref, position + extra.length() + 1 + compm.length() + newdel, sanend);
                    findMM3 = findMM3(ref, position + 1, sanpseq);
                }

                // mismatches, mismatch positions, 5 or 3 ends
                MismatchResult findMM5 = findMM5(ref, position + extra.length() + compm.length() + newdel, wupseq);

                List<Mismatch> mm3 = findMM3.getMismatches();
                List<Integer> sc3p = findMM3.getScp();
                int nm3 = findMM3.getNm();
                int misp3 = findMM3.getMisp();
                String misnt3 = findMM3.getMisnt();

                List<Mismatch> mm5 = findMM5.getMismatches();
                List<Integer> sc5p = findMM5.getScp();
                int nm5 = findMM5.getNm();
                int misp5 = findMM5.getMisp();
                String misnt5 = findMM5.getMisnt();

                List<Mismatch> mmm = new ArrayList<>(mm3);
                mmm.addAll(mm5);
                Variation vref = getVariation(insertionVariants, position, vn);
                for (Mismatch mismatch : mmm) {
                    // $mm mismatch nucleotide
                    String mismatchBases = mismatch.mismatchSequence;
                    // $mp start position of clip that contains mm
                    Integer mismatchPosition = mismatch.mismatchPosition;
                    // $me end (3 or 5)
                    Integer mismatchEnd = mismatch.end;

                    if (mismatchBases.length() > 1) {
                        mismatchBases = mismatchBases.charAt(0) + "&" + mismatchBases.substring(1);
                    }
                    if (!nonInsertionVariants.containsKey(mismatchPosition)) {
                        continue;
                    }

                    Variation variation = nonInsertionVariants.get(mismatchPosition).get(mismatchBases);
                    if (variation == null) {
                        continue;
                    }
                    if (variation.varsCount == 0) {
                        continue;
                    }
                    if (variation.meanQuality / variation.varsCount < instance().conf.goodq) {
                        continue;
                    }
                    if (variation.meanPosition / variation.varsCount > (mismatchEnd == 3 ? nm3 + 4 : nm5 + 4)) { // opt_k;
                        continue;
                    }
                    if (variation.varsCount >= insertionCount + insert.length() || variation.varsCount / insertionCount >= 8) {
                        continue;
                    }
                    if (instance().conf.y) {
                        System.err.printf("    insMM: %s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
                                mismatchBases, mismatchPosition, mismatchEnd, nm3, nm5, vn, insertionCount, variation.varsCount, variation.meanQuality, variation.meanPosition, refCoverage.get(position));
                    }
                    // Adjust ref cnt so that AF won't > 1
                    if (mismatchPosition > position && mismatchEnd == 5) {
                        incCnt(refCoverage, position, variation.varsCount);
                    }

                    Variation lref = null;
                    if (mismatchPosition > position && mismatchEnd == 3 &&
                            nonInsertionVariants.containsKey(position) &&
                            ref.containsKey(position) &&
                            nonInsertionVariants.get(position).containsKey(ref.get(position).toString())) {

                        lref = nonInsertionVariants.get(position).get(ref.get(position).toString());
                    }
                    adjCnt(vref, variation, lref);
                    nonInsertionVariants.get(mismatchPosition).remove(mismatchBases);
                    if (nonInsertionVariants.get(mismatchPosition).isEmpty()) {
                        nonInsertionVariants.remove(mismatchPosition);
                    }
                }
                if (misp3 != 0 && mm3.size() == 1 && nonInsertionVariants.containsKey(misp3)
                        && nonInsertionVariants.get(misp3).containsKey(misnt3)
                        && nonInsertionVariants.get(misp3).get(misnt3).varsCount < insertionCount) {
                    nonInsertionVariants.get(misp3).remove(misnt3);
                }
                if (misp5 != 0 && mm5.size() == 1 && nonInsertionVariants.containsKey(misp5)
                        && nonInsertionVariants.get(misp5).containsKey(misnt5)
                        && nonInsertionVariants.get(misp5).get(misnt5).varsCount < insertionCount) {
                    nonInsertionVariants.get(misp5).remove(misnt5);
                }
                for (Integer sc5pp : sc5p) {
                    Sclip tv = softClips5End.get(sc5pp);
                    if (instance().conf.y) {
                        System.err.printf("    55: %s %s VN: '%s'  5' seq: ^%s^\n", position, sc5pp, vn, wupseq);
                    }
                    if (tv != null && !tv.used) {
                        String seq = findconseq(tv, 0);
                        if (instance().conf.y) {
                            System.err.printf("    ins5: %s %s %s %s VN: %s iCnt: %s vCnt: %s\n", position, sc5pp, seq, wupseq, vn, insertionCount, tv.varsCount);
                        }
                        if (!seq.isEmpty() && ismatch(seq, wupseq, -1)) {
                            if (instance().conf.y) {
                                System.err.printf("      ins5: %s %s %s %s VN: %s iCnt: %s cCnt: %s used\n", position, sc5pp, seq, wupseq, vn, insertionCount, tv.varsCount);
                            }
                            if (sc5pp > position) {
                                incCnt(refCoverage, position, tv.varsCount);
                            }
                            adjCnt(vref, tv);
                            tv.used = true;

                            //To find a case and implement later
                            if (insert.length() + 1 == vn.length() && sc5pp <= position) {
                                // Not commented in perl, created for special case
                                // System.err.printf(" %s %s\n", sc5pp, p);
                            }
                        }
                    }
                }
                for (Integer sc3pp : sc3p) {
                    Sclip tv = softClips3End.get(sc3pp);
                    if (instance().conf.y) {
                        System.err.printf("    33: %s %s VN: '%s'  3' seq: ^%s^\n", position, sc3pp, vn, sanpseq);
                    }
                    if (tv != null && !tv.used) {
                        String seq = findconseq(tv, 0);
                        if (instance().conf.y) {
                            System.err.printf("    ins3: %s %s %s %s VN: %s iCnt: %s vCnt: %s\n", position, sc3pp, seq, sanpseq, vn, insertionCount, tv.varsCount);
                        }
                        String mseq = !ins3.isEmpty() ? sanpseq : substr(sanpseq, sc3pp - position - 1);
                        if (!seq.isEmpty() && ismatch(seq, mseq, 1)) {
                            if (instance().conf.y) {
                                System.err.printf("      ins3: %s %s %s VN: %s iCnt: %s vCnt: %s used\n", position, sc3pp, seq, vn, insertionCount, tv.varsCount);
                            }
                            if (sc3pp <= position || insert.length() > tv.meanPosition / tv.varsCount) {
                                incCnt(refCoverage, position, tv.varsCount);
                            }
                            Variation lref = null;
                            if (sc3pp > position &&
                                    nonInsertionVariants.containsKey(position) &&
                                    ref.containsKey(position) &&
                                    nonInsertionVariants.get(position).containsKey(ref.get(position).toString())) {

                                lref = nonInsertionVariants.get(position).get(ref.get(position).toString());
                            }
                            if (insert.length() > tv.meanPosition / tv.varsCount) {
                                lref = null;
                            }
                            adjCnt(vref, tv, lref);
                            tv.used = true;
                            if (insert.length() + 1 == vn.length() && insert.length() > maxReadLength
                                    && sc3pp >= position + 1 + insert.length()) {
                                int flag = 0;
                                int offset = (sc3pp - position - 1) % insert.length();
                                String tvn = vn;
                                for (int seqi = 0; seqi < seq.length() && seqi + offset < insert.length(); seqi++) {
                                    if (!substr(seq, seqi, 1).equals(substr(insert, seqi + offset, 1))) {
                                        flag++;
                                        int shift = seqi + offset + 1;
                                        tvn = tvn.substring(0, shift) + substr(seq, seqi, 1) + tvn.substring(shift + 1);
                                    }
                                }
                                if (flag > 0) {
                                    Variation variation = insertionVariants.get(position).get(vn);
                                    insertionVariants.get(position).put(tvn, variation);
                                    insertionVariants.get(position).remove(vn);
                                    NEWINS = tvn;
                                }
                            }
                        }
                    }
                }
                int first3 = sc3p.get(0);
                int first5 = sc5p.get(0);
                if (!sc3p.isEmpty() && !sc5p.isEmpty()
                        && first3 > first5 + 3
                        && first3 - first5 < maxReadLength * 0.75) {
                    if (ref.containsKey(position) && nonInsertionVariants.containsKey(position)
                            && nonInsertionVariants.get(position).containsKey(ref.get(position).toString())) {
                        adjRefFactor(nonInsertionVariants.get(position).get(ref.get(position).toString()), (first3 - first5 - 1) / (double) maxReadLength);
                    }
                    adjRefFactor(vref, -(first3 - first5 - 1) / (double) maxReadLength);
                }
            } catch (Exception exception) {
                printExceptionAndContinue(exception, "variant", String.valueOf(lastPosition), region);
            }
        }

        for (int i = tmp.size() - 1; i > 0; i--) {
            try {
                SortPositionDescription tpl = tmp.get(i);
                Integer p = tpl.position;
                lastPosition = p;
                String vn = tpl.descriptionString;
                if (!insertionVariants.containsKey(p)) {
                    continue;
                }
                Variation vref = insertionVariants.get(p).get(vn);
                if (vref == null) {
                    continue;
                }
                Matcher mtch = ATGSs_AMP_ATGSs_END.matcher(vn);
                if (mtch.find()) {
                    String tn = mtch.group(1);
                    Variation tref = insertionVariants.get(p).get(tn);
                    if (tref != null) {
                        if (vref.varsCount < tref.varsCount) {
                            adjCnt(tref, vref, getVariationMaybe(nonInsertionVariants, p, ref.get(p)));
                            insertionVariants.get(p).remove(vn);
                        }
                    }
                }
            } catch (Exception exception) {
                printExceptionAndContinue(exception, "variant", String.valueOf(lastPosition), region);
            }
        }
        return NEWINS;
    }


    /**
     * Realign large deletions that are not present in alignment
     * @param svfdel list of DEL SVs in forward strand
     * @param svrdel list of DEL SVs in reverse strand
     */
    public void realignlgdel(List<Sclip> svfdel,
                             List<Sclip> svrdel) {
        final int longmm = 3;
        Map<Integer, Character> ref = reference.referenceSequences;

        List<SortPositionSclip> tmp = new ArrayList<>();
        for (Map.Entry<Integer, Sclip> ent5 : softClips5End.entrySet()) {
            int p = ent5.getKey();
            if (p < region.start - Configuration.EXTENSION
                    || p > region.end + Configuration.EXTENSION) {
                continue;
            }
            tmp.add(new SortPositionSclip(ent5.getKey(), ent5.getValue(), 0));
        }
        Collections.sort(tmp, COMP2);
        int svcov = 0;
        int clusters = 0;
        int pairs = 0;
        int lastPosition = 0;
        for (SortPositionSclip t : tmp) {
            try {
                int p = t.position;
                lastPosition = p;
                Sclip sc5v = t.softClip;
                int cnt = t.softClip.varsCount;
                if (cnt < instance().conf.minr) {
                    break;
                }
                //already been used in
                if (sc5v.used) {
                    continue;
                }
                String seq = findconseq(sc5v, 5);
                if (seq.isEmpty()) {
                    continue;
                }

                if (seq.length() < 7) {
                    continue;
                }

                if (instance().conf.y) {
                    System.err.printf("  Working Realignlgdel: 5' %s '%s' %s\n", p, seq, cnt);
                }

                int bp = findbp(seq, p - 5, ref, -1, chr);

                String extra = "";
                String EXTRA = "";

                if (bp == 0) {
                    //next if ($SOFTP2SV{ $p } && $SOFTP2SV{ $p }->[0]->{ used });
                    if (islowcomplexseq(seq)) {
                        continue;
                    }
                    Match match = findMatch(seq, reference, p, -1, Configuration.SEED_1, 1);
                    bp = match.basePosition;
                    EXTRA = match.matchedSequence;
                    if (!(bp != 0 && p - bp > 15 && p - bp < Configuration.SVMAXLEN)) {
                        continue;
                    }
                    bp++;
                    Tuple.Tuple3<Integer, Integer, Integer> svMark = markSV(bp, p, Arrays.asList(svfdel, svrdel), maxReadLength);
                    svcov = svMark._1;
                    clusters = svMark._2;
                    pairs = svMark._3;
                    if (svcov == 0) {
                        if (cnt <= instance().conf.minr) {
                            continue;
                        }
                    }

                    VariationMap.SV sv = getSV(nonInsertionVariants, bp);
                    sv.type = "DEL";
                    sv.pairs += pairs;
                    sv.splits += cnt;
                    sv.clusters += clusters;

                    if (bp < region.start) {
                        int tts = bp - maxReadLength;
                        int tte = bp + maxReadLength;
                        if (bp + maxReadLength >= region.start) {
                            tte = region.start - 1;
                        }
                        Region modifiedRegion = Region.newModifiedRegion(region, tts, tte);
                        if (!isLoaded(chr, tts, tte, reference)) {
                            referenceResource.getReference(modifiedRegion, maxReadLength, reference);
                        }
                        Scope<InitialData> currentScope = new Scope<>(bam, modifiedRegion, reference, referenceResource, maxReadLength, splice,
                                variantPrinter, new InitialData(nonInsertionVariants, insertionVariants, refCoverage, softClips3End, softClips5End));
                        getMode().partialPipeline(currentScope, new DirectThreadExecutor());
                    }
                }
                int dellen = p - bp;
                int en = 0;
                String gt = "-" + dellen;

                if (EXTRA.isEmpty()) {
                    while (en < seq.length() && isNotEquals(seq.charAt(en), ref.get(bp - en - 1))) {
                        extra += substr(seq, en, 1);
                        en++;
                    }
                    if (!extra.isEmpty()) {
                        extra = reverse(extra);
                        gt = "-" + dellen + "&" + extra;
                        bp -= extra.length();
                    }
                } else {
                    dellen -= EXTRA.length();
                    gt = dellen == 0 ? "-" + EXTRA.length() + "^" + EXTRA : "-" + dellen + "&" + EXTRA;
                }
                if (instance().conf.y) {
                    System.err.printf("  Found Realignlgdel: %s %s 5' %s %s %s\n", bp, gt, p, seq, cnt);
                }

                // Work on the softclipped read at 3'
                int n = 0;

                if (extra.isEmpty() && EXTRA.isEmpty()) {
                    while (ref.containsKey(bp + n)
                            && ref.containsKey(bp + dellen + n)
                            && isEquals(ref.get(bp + n), ref.get(bp + dellen + n))) {
                        n++;
                    }
                }
                int sc3p = bp + n;
                StringBuilder str = new StringBuilder();
                int mcnt = 0;
                while (mcnt <= longmm
                        && ref.containsKey(bp + n)
                        && ref.containsKey(bp + dellen + n)
                        && isNotEquals(ref.get(bp + n), ref.get(bp + dellen + n))) {
                    str.append(ref.get(bp + dellen + n));
                    n++;
                    mcnt++;
                }
                if (str.length() == 1) {
                    int nm = 0;
                    while (ref.containsKey(bp + n)
                            && ref.containsKey(bp + dellen + n)
                            && isEquals(ref.get(bp + n), ref.get(bp + dellen + n))) {
                        n++;
                        if (n != 0) {
                            nm++;
                        }
                    }
                    if (nm >= 3 && !softClips3End.containsKey(sc3p)) {
                        sc3p = bp + n;
                    }
                }

                // likely a false positive
                if (nonInsertionVariants.containsKey(bp) && nonInsertionVariants.get(bp).sv != null && !softClips3End.containsKey(sc3p)) {
                    if (svcov == 0 && cnt <= instance().conf.minr) {
                        removeSV(nonInsertionVariants, bp);
                        continue;
                    }
                }
                final Variation tv = getVariation(nonInsertionVariants, bp, gt);
                //$hash->{$bp}->{ $gt }->{ cnt } = 0 unless( $hash->{$bp}->{ $gt }->{ cnt } );
                tv.qstd = true; // more accurate implementation lat
                tv.pstd = true; // more accurate implementation lat

                adjCnt(tv, sc5v);
                // $sc5v->{ used } = $bp;
                sc5v.used = bp != 0;
                //$cov->{ $bp } = $cov->{ $p } unless( $cov->{ $bp } );
                if (!refCoverage.containsKey(bp) && refCoverage.containsKey(p)) {
                    refCoverage.put(bp, refCoverage.get(p));
                }
                if (dellen < instance().conf.indelsize) {
                    for (int tp = bp; tp < bp + dellen; tp++) {
                        incCnt(refCoverage, tp, sc5v.varsCount);
                    }
                }

                if (softClips3End.containsKey(sc3p) && !softClips3End.get(sc3p).used) {
                    Sclip sclip = softClips3End.get(sc3p);
                    if (sc3p > bp) {
                        adjCnt(tv, sclip, getVariationMaybe(nonInsertionVariants, bp, ref.get(bp)));
                    } else {
                        adjCnt(tv, sclip);
                    }

                    if (sc3p == bp) {
                        if (dellen < instance().conf.indelsize) {
                            for (int tp = bp; tp < bp + dellen; tp++) {
                                incCnt(refCoverage, tp, sclip.varsCount);
                            }
                        }
                    }
                    for (int ip = bp + 1; ip < sc3p; ip++) {
                        Variation vv = getVariation(nonInsertionVariants, ip, ref.get(dellen + ip).toString());
                        rmCnt(vv, sclip);
                        if (vv.varsCount == 0) {
                            nonInsertionVariants.get(ip).remove(ref.get(dellen + ip).toString());
                        }
                        if (nonInsertionVariants.get(ip).isEmpty()) {
                            nonInsertionVariants.remove(ip);
                        }
                    }
                    sclip.used = bp != 0;
                }
                Map<Integer, Map<String, Integer>> dels5 = singletonMap(bp, singletonMap(gt, tv.varsCount));
                realigndel(bams, dels5);

                if (nonInsertionVariants.containsKey(bp) && nonInsertionVariants.get(bp).sv != null) {
                    nonInsertionVariants.get(bp).sv.splits += tv.varsCount - dels5.get(bp).get(gt);
                }
                if (svcov > tv.varsCount) {
                    addVarFactor(tv, (svcov - tv.varsCount) / (double) tv.varsCount);
                }
                if (instance().conf.y) {
                    System.err.printf("  Found lgdel done: %s %s %s 5' %s %s\n\n", bp, gt, p, seq, tv.varsCount);
                }
            } catch (Exception exception) {
                printExceptionAndContinue(exception, "variant", String.valueOf(lastPosition), region);
            }
        }

        // Work on 3' clipped reads
        tmp = new ArrayList<>();
        for (Map.Entry<Integer, Sclip> ent3 : softClips3End.entrySet()) {
            int p = ent3.getKey();
            if (p < region.start - Configuration.EXTENSION
                    || p > region.end + Configuration.EXTENSION) {
                continue;
            }
            tmp.add(new SortPositionSclip(ent3.getKey(), ent3.getValue(), 0));
        }
        Collections.sort(tmp, COMP2);
        svcov = 0;
        clusters = 0;
        pairs = 0;
        for (SortPositionSclip t : tmp) {
            try {
                int p = t.position;
                lastPosition = p;
                final Sclip sc3v = t.softClip;
                final int cnt = sc3v.varsCount;
                if (cnt < instance().conf.minr) {
                    break;
                }
                if (sc3v.used) {
                    continue;
                }
                String seq = findconseq(sc3v, 3);
                if (seq.isEmpty()) {
                    continue;
                }
                if (seq.length() < 7) {
                    continue;
                }

                if (instance().conf.y) {
                    System.err.printf("  Working Realignlgdel: 3' %s '%s' %s\n", p, seq, cnt);
                }

                int bp = findbp(seq, p + 5, ref, 1, chr);
                String extra = "";
                String EXTRA = "";

                if (bp == 0) {
                    //#next if ($SOFTP2SV{ $p } && $SOFTP2SV{ $p }->[0]->{ used });
                    if (islowcomplexseq(seq)) {
                        continue;
                    }
                    Match match = findMatch(seq, reference, p, 1, Configuration.SEED_1, 1);
                    bp = match.basePosition;
                    EXTRA = match.matchedSequence;
                    if (!(bp != 0 && bp - p > 15 && p - bp < Configuration.SVMAXLEN)) {
                        continue;
                    }

                    Tuple.Tuple3<Integer, Integer, Integer> svMark = markSV(p, bp, Arrays.asList(svfdel, svrdel), maxReadLength);
                    svcov = svMark._1;
                    clusters = svMark._2;
                    pairs = svMark._3;
                    if (svcov == 0) {
                        if (cnt <= instance().conf.minr) { //a little more stringent
                            continue;
                        }
                    }

                    VariationMap.SV sv = getSV(nonInsertionVariants, p);
                    sv.type = "DEL";
                    sv.pairs += pairs;
                    sv.splits += cnt;
                    sv.clusters += clusters;

                    if (bp > region.end) {
                        int tts = bp - maxReadLength;
                        int tte = bp + maxReadLength;
                        if (bp - maxReadLength <= region.end) {
                            tts = region.end + 1;
                        }
                        Region modifiedRegion = Region.newModifiedRegion(region, tts, tte);
                        if (!isLoaded(chr, tts, tte, reference)) {
                            referenceResource.getReference(modifiedRegion, maxReadLength, reference);
                        }
                        Scope<InitialData> currentScope = new Scope<>(bam, modifiedRegion, reference, referenceResource, maxReadLength, splice,
                                variantPrinter, new InitialData(nonInsertionVariants, insertionVariants, refCoverage, softClips3End, softClips5End));
                        getMode().partialPipeline(currentScope, new DirectThreadExecutor());
                    }
                }

                int dellen = bp - p;
                int en = 0;
                if (!EXTRA.isEmpty()) {
                    dellen -= EXTRA.length();
                } else {
                    while (en < seq.length() && isNotEquals(seq.charAt(en), ref.get(bp + en))) {
                        extra += seq.charAt(en);
                        en++;
                    }
                }
                String gt = "-" + dellen;
                int sc5p = bp;
                bp = p; // Set it to 5'
                if (!extra.isEmpty()) {
                    gt = "-" + dellen + "&" + extra;
                    sc5p += extra.length();
                } else if (!EXTRA.isEmpty()) {
                    gt = "-" + dellen + "&" + EXTRA;
                } else {
                    // 5' adjustment
                    while (ref.containsKey(bp - 1) && ref.containsKey(bp + dellen - 1)
                            && isEquals(ref.get(bp - 1), ref.get(bp + dellen - 1))) {
                        bp--;
                        if (bp != 0) {
                            sc5p--;
                        }
                    }
                    if (bp != p && nonInsertionVariants.containsKey(p) && nonInsertionVariants.get(p).sv != null) {
                        VariationMap.SV sv = getSV(nonInsertionVariants, bp);
                        sv.clusters = nonInsertionVariants.get(p).sv.clusters;
                        sv.pairs = nonInsertionVariants.get(p).sv.pairs;
                        sv.splits = nonInsertionVariants.get(p).sv.splits;
                        sv.type = nonInsertionVariants.get(p).sv.type;
                        removeSV(nonInsertionVariants, p);
                    }
                }

                if (nonInsertionVariants.containsKey(bp) && nonInsertionVariants.get(bp).sv != null && !softClips5End.containsKey(sc5p)) {
                    if (svcov == 0 && cnt <= instance().conf.minr) {
                        removeSV(nonInsertionVariants, bp);
                        continue;
                    }
                }
                if (instance().conf.y) {
                    System.err.printf("  Found Realignlgdel: bp: %s %s 3' %s 5'clip: %s '%s' %s\n", bp, gt, p, sc5p, seq, cnt);
                }

                Variation tv = getVariation(nonInsertionVariants, bp, gt);
                tv.qstd = true; // more accurate implementation later
                tv.pstd = true; // more accurate implementation later
                if (dellen < instance().conf.indelsize) {
                    for (int tp = bp; tp < bp + dellen + extra.length() + EXTRA.length(); tp++) {
                        incCnt(refCoverage, tp, sc3v.varsCount);
                    }
                }
                //$cov->{$bp} = $cov->{ $p - 1 } ? $cov->{ $p - 1 } : $sc3v->{ cnt } unless( $cov->{ $bp } );
                if (!refCoverage.containsKey(bp)) {
                    if (refCoverage.containsKey(p - 1)) {
                        refCoverage.put(bp, refCoverage.get(p - 1));
                    } else refCoverage.put(bp, sc3v.varsCount);
                }
                sc3v.meanPosition += dellen * sc3v.varsCount;
                adjCnt(tv, sc3v);
                sc3v.used = p + dellen != 0;

                Map<Integer, Map<String, Integer>> dels5 = new HashMap<>();
                HashMap<String, Integer> map = new HashMap<>();
                map.put(gt, tv.varsCount);
                dels5.put(bp, map);
                realigndel(bams, dels5);

                if (nonInsertionVariants.containsKey(bp) && nonInsertionVariants.get(bp).sv != null) {
                    nonInsertionVariants.get(bp).sv.splits += tv.varsCount - dels5.get(bp).get(gt);
                }
                if (instance().conf.y) {
                    System.err.printf("  Found lgdel: %s %s %s 3' '%s' %s\n\n", bp, gt, p, seq, tv.varsCount);
                }
                if (svcov > tv.varsCount) {
                    addVarFactor(tv, (svcov - tv.varsCount) / (double) tv.varsCount);
                }
            } catch (Exception exception) {
                printExceptionAndContinue(exception, "variant", String.valueOf(lastPosition), region);
            }
        }
        if (instance().conf.y) {
            System.err.println("  Done: Realignlgdel\n");
        }
    }

    /**
     * This will try to realign large insertions (typically larger than 30bp)
     */
    void realignlgins30() {
        Map<Integer, Character> ref = reference.referenceSequences;

        List<SortPositionSclip> tmp5 = new ArrayList<>();
        for (Map.Entry<Integer, Sclip> ent5 : softClips5End.entrySet()) {
            //ent5: (position, soft clipping 5' variant)
            int p = ent5.getKey();
            if (p < region.start - Configuration.EXTENSION
                    || p > region.end + Configuration.EXTENSION) {
                continue;
            }
            tmp5.add(new SortPositionSclip(ent5.getKey(), ent5.getValue(), ent5.getValue().varsCount));
        }
        Collections.sort(tmp5, COMP3);

        List<SortPositionSclip> tmp3 = new ArrayList<>();
        for (Map.Entry<Integer, Sclip> ent3 : softClips3End.entrySet()) {
            int p = ent3.getKey();
            if (p < region.start - Configuration.EXTENSION
                    || p > region.end + Configuration.EXTENSION) {
                continue;
            }
            tmp3.add(new SortPositionSclip(ent3.getKey(), ent3.getValue(), ent3.getValue().varsCount));
        }
        Collections.sort(tmp3, COMP3);
        int lastPosition = 0;
        for (SortPositionSclip t5 : tmp5) {
            try {
                final int p5 = t5.position;
                lastPosition = p5;
                final Sclip sc5v = t5.softClip;
                final int cnt5 = t5.count;
                if (sc5v.used) {
                    continue;
                }

                for (SortPositionSclip t3 : tmp3) {
                    try {
                        final int p3 = t3.position;
                        lastPosition = p3;
                        final Sclip sc3v = t3.softClip;
                        final int cnt3 = t3.count;
                        if (sc5v.used) {
                            break;
                        }
                        if (sc3v.used) {
                            continue;
                        }
                        if (p5 - p3 > maxReadLength * 2.5) {
                            continue;
                        }
                        if (p3 - p5 > maxReadLength - 10) { // if they're too far away, don't even try
                            continue;
                        }
                        final String seq5 = findconseq(sc5v, 5);
                        final String seq3 = findconseq(sc3v, 3);
                        //next until at least one of consensus sequences has length > 10
                        if (seq5.length() <= 10 || seq3.length() <= 10) {
                            continue;
                        }
                        if (instance().conf.y) {
                            System.err.printf("  Working lgins30: %s %s 3: %s %s 5: %s %s\n",
                                    p3, p5, seq3, cnt3, new StringBuilder(seq5).reverse(), cnt5);
                        }
                        // Don't even try if there're extreme bias

                        if (!(cnt5 / (double) cnt3 >= 0.08 && cnt5 / (double) cnt3 <= 12)) {
                            continue;
                        }
                        Match35 match35 = find35match(seq5, seq3);
                        int bp5 = match35.matched5end;
                        int bp3 = match35.matched3End;
                        //length of match
                        int score = match35.maxMatchedLength;

                        if (score == 0) {
                            continue;
                        }

                        //to ensure higher quality of bases are used as read ends are usually poor quality
                        int smscore = (int) (score / 2);
                        //add matched part of seq3 (left clip)
                        String ins = bp3 + smscore > 1 ? substr(seq3, 0, -(bp3 + smscore) + 1) : seq3;
                        if (bp5 + smscore > 0) { //add not matched part of seq5 (right clip)
                            ins += reverse(substr(seq5, 0, bp5 + smscore));
                        }
                        if (islowcomplexseq(ins)) {
                            if (instance().conf.y) {
                                System.err.println("  Discard low complexity insertion found " + ins + ".");
                            }
                            continue;
                        }
                        int bi = 0;
                        Variation vref;
                        if (instance().conf.y) {
                            System.err.printf("  Found candidate lgins30: %s %s %s\n", p3, p5, ins);
                        }
                        if (p5 > p3) {
                            if (seq3.length() > ins.length()
                                    && !ismatch(substr(seq3, ins.length()), joinRef(ref, p5, p5 + seq3.length() - ins.length() + 2), 1)) {
                                continue;
                            }
                            if (seq5.length() > ins.length()
                                    && !ismatch(substr(seq5, ins.length()), joinRef(ref, p3 - (seq5.length() - ins.length()) - 2, p3 - 1), -1)) {
                                continue;
                            }
                            if (instance().conf.y) {
                                System.err.printf("  Found lgins30 complex: %s %s %s %s\n", p3, p5, ins.length(), ins);
                            }
                            String tmp = joinRef(ref, p3, p5 - 1);
                            if (tmp.length() > ins.length()) { // deletion is longer
                                ins = (p3 - p5) + "^" + ins;
                                bi = p3;
                                vref = getVariation(nonInsertionVariants, p3, ins);
                            } else if (tmp.length() < ins.length()) {
                                // In perl it isn't commented, but not used
//                        int p3s = p3 + tmp.length();
//                        int p3e = p3s + seq3.length() - ins.length() + 2;
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
                            if (seq3.length() > ins.length()
                                    && !ismatch(substr(seq3, ins.length()), joinRef(ref, p5, p5 + seq3.length() - ins.length() + 2), 1)) {
                                continue;
                            }
                            if (seq5.length() > ins.length()
                                    && !ismatch(substr(seq5, ins.length()), joinRef(ref, p3 - (seq5.length() - ins.length()) - 2, p3 - 1), -1)) {
                                continue;
                            }
                            String tmp = "";
                            if (ins.length() <= p3 - p5) { // Tandex duplication
                                int rpt = 2;
                                int tnr = 3;
                                while (((p3 - p5 + ins.length()) / (double) tnr) / (double) ins.length() > 1) {
                                    if ((p3 - p5 + ins.length()) % tnr == 0) {
                                        rpt++;
                                    }
                                    tnr++;
                                }
                                //TODO: maybe here in Perl we must cast position "to" to int?
                                // -1 because joinRef works to the bound exactly
                                tmp += joinRef(ref, p5, (p5 + (p3 - p5 + ins.length()) / (double) rpt - ins.length()));
                                ins = "+" + tmp + "" + ins;
                            } else {
                                // -1 because joinRef works to the bound exactly
                                tmp += joinRef(ref, p5, p3 - 1);
                                if ((ins.length() - tmp.length()) % 2 == 0) {
                                    int tex = (ins.length() - tmp.length()) / 2;
                                    ins = (tmp + substr(ins, 0, tex)).equals(substr(ins, tex))
                                            ? ("+" + substr(ins, tex))
                                            : "+" + tmp + "" + ins;
                                } else {
                                    ins = "+" + tmp + "" + ins;
                                }
                            }
                            if (instance().conf.y) {
                                System.err.printf("Found lgins30: %s %s %s %s + %s\n", p3, p5, ins.length(), tmp, ins);
                            }
                            bi = p5 - 1;
                            vref = getVariation(insertionVariants, bi, ins);
                        }
                        sc3v.used = true;
                        sc5v.used = true;
                        vref.pstd = true;
                        vref.qstd = true;
                        incCnt(refCoverage, bi, sc5v.varsCount);
                        if (instance().conf.y) {
                            System.err.printf(" lgins30 Found: '%s' %s %s %s\n", ins, bi, bp3, bp5);
                        }

                        if (ins.startsWith("+")) {
                            Variation mvref = getVariationMaybe(nonInsertionVariants, bi, ref.get(bi));
                            adjCnt(vref, sc3v, mvref);
                            adjCnt(vref, sc5v);
                            if (bams != null && bams.length > 0
                                    && p3 - p5 >= 5 && p3 - p5 < maxReadLength - 10
                                    && mvref != null && mvref.varsCount != 0
                                    && noPassingReads(chr, p5, p3, bams)
                                    && vref.varsCount > 2 * mvref.varsCount) {
                                adjCnt(vref, mvref, mvref);
                            }
                            Map<Integer, Map<String, Integer>> tins = new HashMap<>();
                            Map<String, Integer> map = new HashMap<>();
                            map.put(ins, vref.varsCount);
                            tins.put(bi, map);
                            realignins(tins);
                        } else if (ins.startsWith("-")) {
                            adjCnt(vref, sc3v, getVariationMaybe(nonInsertionVariants, bi, ref.get(bi)));
                            adjCnt(vref, sc5v);
                            Map<Integer, Map<String, Integer>> tdel = new HashMap<>();
                            Map<String, Integer> map = new HashMap<>();
                            map.put(ins, vref.varsCount);
                            tdel.put(bi, map);
                            realigndel(null, tdel);
                        } else {
                            adjCnt(vref, sc3v);
                            adjCnt(vref, sc5v);
                        }
                        break;
                    } catch (Exception exception) {
                        printExceptionAndContinue(exception, "variant", String.valueOf(lastPosition), region);
                    }
                }
            } catch (Exception exception) {
                printExceptionAndContinue(exception, "variant", String.valueOf(lastPosition), region);
            }
        }
        if (instance().conf.y) {
            System.err.println("Done: lgins30\n");
        }
    }

    /**
     * Realign large insertions that are not present in alignment
     * @param svfdup list of DUP SVs in forward strand
     * @param svrdup list of DUP SVs in reverse strand
     */
    public void realignlgins(List<Sclip> svfdup,
                             List<Sclip> svrdup)  {
        Map<Integer, Character> ref = reference.referenceSequences;

        List<SortPositionSclip> tmp = new ArrayList<>();
        for (Map.Entry<Integer, Sclip> ent5 : softClips5End.entrySet()) {
            int p = ent5.getKey();
            if (p < region.start - Configuration.EXTENSION
                    || p > region.end + Configuration.EXTENSION) {
                continue;
            }
            tmp.add(new SortPositionSclip(ent5.getKey(), ent5.getValue(), 0));
        }
        //sort by descending cnt
        Collections.sort(tmp, COMP2);
        int lastPosition = 0;
        for (SortPositionSclip t : tmp) {
            try {
                int p = t.position;
                lastPosition = p;
                Sclip sc5v = t.softClip;
                int cnt = t.softClip.varsCount;
                if (cnt < instance().conf.minr) {
                    break;
                }
                //already been used in
                if (sc5v.used) {
                    continue;
                }
                String seq = findconseq(sc5v, 0);
                if (seq.isEmpty()) {
                    continue;
                }
                if (instance().conf.y) {
                    System.err.println("  Working lgins: 5: " + p + " " + seq + " cnt: " + cnt);
                }
                if (seq.length() < 12) {
                    continue;
                }

                BaseInsertion tpl = findbi(seq, p, ref, -1, chr);
                int bi = tpl.baseInsert;
                String ins = tpl.insertionSequence;
                String EXTRA = "";

                if (bi == 0) {
                    //#next if ($SOFTP2SV{ $p } && $SOFTP2SV{ $p }->[0]->{ used });
                    if (islowcomplexseq(seq)) {
                        continue;
                    }
                    Match match = findMatch(seq, reference, p, -1, Configuration.SEED_1, 1);
                    bi = match.basePosition;
                    EXTRA = match.matchedSequence;
                    if (!(bi != 0 && bi - p > 15 && bi - p < Configuration.SVMAXLEN)) {
                        continue;
                    }

                    // For large insertions
                    if (bi > region.end) {
                        //my ($tts, $tte) = ($bi - $RLEN <= $END ? $END + 1 : $bi - $RLEN, $bi + $RLEN);
                        int tts = bi - maxReadLength;
                        int tte = bi + maxReadLength;
                        if (bi - maxReadLength <= region.end) {
                            tts = region.end + 1;
                        }
                        Region modifiedRegion = Region.newModifiedRegion(region, tts, tte);
                        if (!isLoaded(chr, tts, tte, reference)) {
                            referenceResource.getReference(modifiedRegion, maxReadLength, reference);
                        }
                        Scope<InitialData> currentScope = new Scope<>(bam, modifiedRegion, reference, referenceResource, maxReadLength, splice,
                                variantPrinter, new InitialData(nonInsertionVariants, insertionVariants, refCoverage, softClips3End, softClips5End));
                        getMode().partialPipeline(currentScope, new DirectThreadExecutor());
                    }
                    if (bi - p > instance().conf.SVMINLEN + 2 * Configuration.SVFLANK) {
                        ins = joinRef(ref, p, p + Configuration.SVFLANK - 1);
                        ins += "<dup" + (bi - p - 2 * Configuration.SVFLANK + 1) + ">";
                        ins += joinRefFor5Lgins(ref, bi - Configuration.SVFLANK + 1, bi, seq, EXTRA);
                    } else {
                        ins = joinRefFor5Lgins(ref, p, bi, seq, EXTRA);
                    }
                    ins += EXTRA;

                    Tuple.Tuple2<Integer, Integer> tp2 = markDUPSV(p, bi, Arrays.asList(svfdup, svrdup), maxReadLength);
                    int clusters = tp2._1;
                    int pairs = tp2._2;
                    if (!refCoverage.containsKey(p - 1) || (refCoverage.containsKey(bi) && refCoverage.containsKey(p - 1)
                            && refCoverage.get(p - 1) < refCoverage.get(bi))) {
                        if (refCoverage.containsKey(bi)) {
                            refCoverage.put(p - 1, refCoverage.get(bi));
                        } else {
                            refCoverage.put(p - 1, sc5v.varsCount);
                        }
                    } else {
                        if (sc5v.varsCount > refCoverage.get(p - 1)) {
                            incCnt(refCoverage, p - 1, sc5v.varsCount);
                        }
                    }

                    bi = p - 1;
                    VariationMap.SV sv = getSV(nonInsertionVariants, bi);
                    sv.type = "DUP";
                    sv.pairs += pairs;
                    sv.splits += cnt;
                    sv.clusters += clusters;
                }

                if (instance().conf.y) {
                    System.err.printf("  Found candidate lgins from 5: %s +%s %s %s\n", bi, ins, p, seq);
                }
                final Variation iref = getVariation(insertionVariants, bi, "+" + ins);
                iref.pstd = true;
                iref.qstd = true;
                adjCnt(iref, sc5v);
                boolean rpflag = true; // A flag to indicate whether an insertion is a repeat
                for (int i = 0; i < ins.length(); i++) {
                    if (!isEquals(ref.get(bi + 1 + i), ins.charAt(i))) {
                        rpflag = false;
                        break;
                    }
                }

                if (nonInsertionVariants.containsKey(bi) && nonInsertionVariants.get(bi).sv == null) {
                    incCnt(refCoverage, bi, sc5v.varsCount);
                }
                int len = ins.length();
                if (ins.indexOf('&') != -1) {
                    len--;
                }
                int seqLen = sc5v.seq.lastKey() + 1;
                for (int ii = len + 1; ii < seqLen; ii++) {
                    int pii = bi - ii + len;
                    if (!sc5v.seq.containsKey(ii)) {
                        continue;
                    }
                    for (Map.Entry<Character, Variation> ent : sc5v.seq.get(ii).entrySet()) {
                        Character tnt = ent.getKey();
                        Variation tv = ent.getValue();
                        Variation tvr = getVariation(nonInsertionVariants, pii, tnt.toString());
                        adjCnt(tvr, tv);
                        tvr.pstd = true;
                        tvr.qstd = true;
                        incCnt(refCoverage, pii, tv.varsCount);
                    }
                }
                sc5v.used = bi + len != 0;

                Map<Integer, Map<String, Integer>> tins = singletonMap(bi, singletonMap("+" + ins, iref.varsCount));
                String newins = realignins(tins);
                if (newins.isEmpty()) {
                    newins = "+" + ins;
                }
                Variation mref = getVariationMaybe(nonInsertionVariants, bi, ref.get(bi));
                Variation kref = getVariation(insertionVariants, bi, newins);

                if (nonInsertionVariants.containsKey(bi) && nonInsertionVariants.get(bi).sv != null) {
                    nonInsertionVariants.get(bi).sv.splits += kref.varsCount - tins.get(bi).get("+" + ins);
                }
                if (rpflag && bams.length > 0 && ins.length() >= 5
                        && ins.length() < maxReadLength - 10
                        && mref != null && mref.varsCount != 0
                        && noPassingReads(chr, bi, bi + ins.length(), bams)
                        && kref.varsCount > 2 * mref.varsCount) {
                    adjCnt(kref, mref, mref);
                }
            } catch (Exception exception) {
                printExceptionAndContinue(exception, "variant", String.valueOf(lastPosition), region);
            }
        }

        tmp = new ArrayList<>();
        for (Map.Entry<Integer, Sclip> ent3 : softClips3End.entrySet()) {
            int p = ent3.getKey();
            if (p < region.start - Configuration.EXTENSION
                    || p > region.end + Configuration.EXTENSION) {
                continue;
            }
            tmp.add(new SortPositionSclip(ent3.getKey(), ent3.getValue(), 0));
        }
        Collections.sort(tmp, COMP2);

        for (SortPositionSclip t : tmp) {
            try {
                int p = t.position;
                lastPosition = p;
                final Sclip sc3v = t.softClip;
                final int cnt = t.softClip.varsCount;
                if (cnt < instance().conf.minr) {
                    break;
                }
                if (sc3v.used) {
                    continue;
                }
                String seq = findconseq(sc3v, 0);
                if (seq.isEmpty()) {
                    continue;
                }
                if (instance().conf.y) {
                    System.err.println("  Working lgins 3: " + p + " " + seq + " cnt: " + cnt);
                }
                if (seq.length() < 12) {
                    continue;
                }

                BaseInsertion tpl = findbi(seq, p, ref, 1, chr);
                int bi = tpl.baseInsert;
                String ins = tpl.insertionSequence;
                String EXTRA = "";

                if (bi == 0) {
                    // #next if ($SOFTP2SV{ $p } && $SOFTP2SV{ $p }->[0]->{ used });
                    if (islowcomplexseq(seq)) {
                        continue;
                    }
                    Match match = findMatch(seq, reference, p, 1, Configuration.SEED_1, 1);
                    bi = match.basePosition;
                    EXTRA = match.matchedSequence;
                    if (!(bi != 0 && p - bi > 15 && p - bi < Configuration.SVMAXLEN)) {
                        continue;
                    }

                    // For large insertions
                    if (bi < region.start) {
                        //my ($tts, $tte) = ($bi - $RLEN, $bi + $RLEN >= $START ? $START - 1 : $bi + $RLEN);
                        //getREF($chr, $tts, $tte, $REF, $RLEN) unless( isLoaded( $chr, $tts, $tte, $REF ) );
                        //parseSAM($chr, $bi - $RLEN, $bi + $RLEN >= $START ? $START - 1 : $bi + $RLEN, $bams, $REF, $hash, $cov, $sclip5, $sclip3, 1);
                        int tts = bi - maxReadLength;
                        int tte = bi + maxReadLength;
                        if (bi + maxReadLength >= region.start) {
                            tte = region.start - 1;
                        }
                        Region modifiedRegion = Region.newModifiedRegion(region, tts, tte);
                        if (!isLoaded(chr, tts, tte, reference)) {
                            referenceResource.getReference(modifiedRegion, maxReadLength, reference);
                        }
                        Scope<InitialData> currentScope = new Scope<>(bam, modifiedRegion, reference, referenceResource, maxReadLength, splice,
                                variantPrinter, new InitialData(nonInsertionVariants, insertionVariants, refCoverage, softClips3End, softClips5End));
                        getMode().partialPipeline(currentScope, new DirectThreadExecutor());
                    }
                    int shift5 = 0;
                    while (ref.containsKey(p - 1) && ref.containsKey(bi - 1)
                            && ref.get(p - 1).equals(ref.get(bi - 1))) {
                        p--;
                        bi--;
                        shift5++;
                    }
                    if (p - bi > instance().conf.SVMINLEN + 2 * Configuration.SVFLANK) {
//                    join("", (map {$_-$bi >= $shift5 && $_-$bi-$shift5 < seq.length() - EXTRA.length()
//                            ? substr(seq, $_-bi-shift5+EXTRA.length(), 1)
//                            : $REF->{ $_ };} (bi .. (bi+Configuration.SVFLANK-1))));
                        ins = joinRefFor3Lgins(ref, bi, bi + Configuration.SVFLANK - 1, shift5, seq, EXTRA);
                        ins += "<dup" + (p - bi - 2 * Configuration.SVFLANK) + ">";
                        ins += joinRef(ref, p - Configuration.SVFLANK, p - 1);
                    } else {
                        ins = joinRefFor3Lgins(ref, bi, p - 1, shift5, seq, EXTRA);
                    }
                    ins += EXTRA;
                    Tuple.Tuple2<Integer, Integer> tp2 = markDUPSV(bi, p - 1, Arrays.asList(svfdup, svrdup), maxReadLength);
                    int clusters = tp2._1;
                    int pairs = tp2._2;
                    bi = bi - 1;

                    VariationMap.SV sv = getSV(nonInsertionVariants, bi);
                    sv.type = "DUP";
                    sv.pairs += pairs;
                    sv.splits += cnt;
                    sv.clusters += clusters;
                    if (!refCoverage.containsKey(bi) || (refCoverage.containsKey(p) && refCoverage.containsKey(bi)
                            && refCoverage.get(bi) < refCoverage.get(p))) {
                        if (refCoverage.containsKey(p)) {
                            refCoverage.put(bi, refCoverage.get(p));
                        } else {
                            refCoverage.put(bi, sc3v.varsCount);
                        }
                    } else {
                        if (sc3v.varsCount > refCoverage.get(bi)) {
                            incCnt(refCoverage, bi, sc3v.varsCount);
                        }
                    }
                }

                if (instance().conf.y) {
                    System.err.printf("  Found candidate lgins from 3: %s +%s %s %s\n", bi, ins, p, seq);
                }

                final Variation iref = getVariation(insertionVariants, bi, "+" + ins);
                iref.pstd = true;
                iref.qstd = true;
                Variation lref = getVariationMaybe(nonInsertionVariants, bi, ref.get(bi));
                if (p - bi > sc3v.meanPosition / cnt) {
                    lref = null;
                }
                adjCnt(iref, sc3v, lref);
                boolean rpflag = true;
                for (int i = 0; i < ins.length(); i++) {
                    if (!isEquals(ref.get(bi + 1 + i), ins.charAt(i))) {
                        rpflag = false;
                        break;
                    }
                }
                // In perl it doesn't commented, but doesn't used
                // final int offset = bi == be ? (p - bi - 1) : -(p + be - bi);
                int len = ins.length();
                if (ins.indexOf('&') != -1) {
                    len--;
                }
                int lenSeq = sc3v.seq.lastKey() + 1;
                for (int ii = len; ii < lenSeq; ii++) {
                    int pii = p + ii - len;
                    Map<Character, Variation> map = sc3v.seq.get(ii);
                    if (map == null) {
                        continue;
                    }
                    for (Map.Entry<Character, Variation> ent : map.entrySet()) {
                        Character tnt = ent.getKey();
                        Variation tv = ent.getValue();
                        Variation vref = getVariation(nonInsertionVariants, pii, tnt.toString());
                        adjCnt(vref, tv);
                        vref.pstd = true;
                        vref.qstd = true;
                        incCnt(refCoverage, pii, tv.varsCount);
                    }
                }
                sc3v.used = true;
                Map<Integer, Map<String, Integer>> tins = singletonMap(bi, singletonMap("+" + ins, iref.varsCount));
                realignins(tins);
                if (nonInsertionVariants.containsKey(bi) && nonInsertionVariants.get(bi).sv != null) {
                    nonInsertionVariants.get(bi).sv.splits += iref.varsCount - tins.get(bi).get("+" + ins);
                }
                Variation mref = getVariationMaybe(nonInsertionVariants, bi, ref.get(bi));
                if (rpflag && bams.length > 0 && ins.length() >= 5 && ins.length() < maxReadLength - 10
                        && mref != null && mref.varsCount != 0
                        && noPassingReads(chr, bi, bi + ins.length(), bams)
                        && iref.varsCount > 2 * mref.varsCount) {
                    adjCnt(iref, mref, mref);
                }
            } catch (Exception exception) {
                printExceptionAndContinue(exception, "variant", String.valueOf(lastPosition), region);
            }
        }
    }

    /**
     * Fill the temp tuple structure from hash tables of insertion or deletions positions to description string
     * and counts of variation.
     * @param changes initial hash map
     * @return tuple of (position, descriptions string and variation count).
     */
    List<SortPositionDescription> fillAndSortTmp(Map<Integer, Map<String, Integer>> changes) {
        //TODO: perl here have non-deterministic results because of hash, maybe we need to
        // make the sort in Perl more stringent (except simple sort b[2]<=>a[2]
        List<SortPositionDescription> tmp = new ArrayList<>();
        for (Map.Entry<Integer, Map<String, Integer>> entry : changes.entrySet()) {
            int position = entry.getKey();
            Map<String, Integer> v = entry.getValue();
            for (Map.Entry<String, Integer> entV : v.entrySet()) {
                String descriptionString = entV.getKey();
                int cnt = entV.getValue();
                // In perl it doesn't commented. but ecnt isn't used
                // int ecnt = 0;
                // Matcher mtch = ATGC_E.matcher(vn);
                // if (mtch.find()) {
                // ecnt = mtch.group(1).length();
                // }
                //, /* ecnt */
                tmp.add(new SortPositionDescription(position, descriptionString, cnt));
            }
        }
        Collections.sort(tmp, (o1, o2) -> {
            int x1 = o1.count;
            int x2 = o2.count;
            int f = Integer.compare(x2, x1);
            if (f != 0)
                return f;
            x1 = o1.position;
            x2 = o2.position;
            f =  Integer.compare(x1, x2);
            if (f != 0)
                return f;
            String s1 = o1.descriptionString;
            String s2 = o2.descriptionString;
            return s2.compareTo(s1);
        });
        return tmp;
    }

    /**
     * Temp structure from hash tables of insertion or deletions positions to description string
     * and counts of variation.
     */
    class SortPositionDescription {
        int position;
        /**
         * description string for variant
         */
        String descriptionString;
        /**
         * variation count
         */
        int count;

        public SortPositionDescription(int position, String descriptionString, int count) {
            this.position = position;
            this.descriptionString = descriptionString;
            this.count = count;
        }
    }

    /**
     * Test whether the two soft-clipped reads match
     * Returns the breakpoints in 5 and 3 prime soft-clipped reads
     * @param seq5 consensus sequence 5' strand
     * @param seq3 consensus sequence 3' strand
     * @return Match (start position of 5' strand, start position of 3' strand, max match length)
     */
    public Match35 find35match(String seq5, String seq3) {
        final int longMismatch = 2; //$longmm
        int maxMatchedLength = 0; //$max
        int b3 = 0;
        int b5 = 0;

        for (int i = 0; i < seq5.length() - 8; i++) {
            for (int j = 1; j < seq3.length() - 8; j++) {
                int numberOfMismatch = 0; //$nm
                int totalLength = 0; //$n
                while (totalLength + j <= seq3.length() && i + totalLength <= seq5.length()) {
                    if (!substr(seq3, -j - totalLength, 1).equals(substr(seq5, i + totalLength, 1))) {
                        numberOfMismatch++;
                    }
                    if (numberOfMismatch > longMismatch) {
                        break;
                    }
                    totalLength++;
                }
                if (totalLength - numberOfMismatch > maxMatchedLength
                        && totalLength - numberOfMismatch > 8
                        && numberOfMismatch / (double) totalLength < 0.1d
                        && (totalLength + j >= seq3.length() || i + totalLength >= seq5.length())) {

                    maxMatchedLength = totalLength - numberOfMismatch;
                    b3 = j;
                    b5 = i;
                    if (instance().conf.y) {
                        System.err.printf("      Found 35 Match, %s %s %s\n", seq5, seq3,
                                join("\t", totalLength + j, seq3.length(), i + totalLength, seq5.length(), totalLength, numberOfMismatch, i, j));
                    }
                    return new Match35(b5, b3, maxMatchedLength);
                }
            }
        }
        return new Match35(b5, b3, maxMatchedLength);
    }

    /**
     * Check whether there're reads supporting wild type in deletions
     * Only for indels that have micro-homology
     * @param chr chromosome name
     * @param start start position
     * @param end end position
     * @param bams BAM file list
     * @return true if any read was found in chr:s-e
     */
    boolean noPassingReads(String chr, int start, int end, String[] bams) {
        int cnt = 0;
        int midcnt = 0; // Reads end in the middle
        int dlen = end - start;
        String dlenqr = dlen + "D";
        Region region = new Region(chr, start, end, "");
        for (String bam : bams) {
            try (SamView reader = new SamView(bam, "0", region, instance().conf.validationStringency)) {
                SAMRecord record;
                while ((record = reader.read()) != null) {
                    if (record.getCigarString().contains(dlenqr)) {
                        continue;
                    }
                    int readStart = record.getAlignmentStart();
                    // The total aligned length, excluding soft-clipped bases and insertions
                    int readLengthIncludeMatchedAndDeleted = getAlignedLength(record.getCigar());
                    int readEnd = readStart + readLengthIncludeMatchedAndDeleted;
                    if (readEnd > end + 2 && readStart < start - 2) {
                        cnt++;
                    }
                    if (readStart < start - 2 && readEnd > start && readEnd < end) {
                        midcnt++;
                    }
                }

            }
        }
        if (instance().conf.y) {
            System.err.printf("    Passing Read CNT: %s %s %s %s %s\n", cnt, chr, start, end, midcnt);
        }
        return cnt <= 0 && midcnt + 1 > 0;
    }

    /**
     * Find if sequences match with no more than default number of mismatches (3).
     * @param seq1 first sequence
     * @param seq2 second sequence
     * @param dir direction of seq2 (1 or -1)
     * @return true if seq1 matches seq2 with no more than 3 mismatches
     */
    public boolean ismatch(String seq1,
                           String seq2,
                           int dir) {
        int MM = 3;
        return ismatch(seq1, seq2, dir, MM);
    }

    /**
     * $ismatch
     * Find if sequences match with no more than MM number of mismatches and total mismatches no more than 15% of
     * length of sequence
     * @param seq1 first sequence
     * @param seq2 second sequence
     * @param dir direction of seq2 (1 or -1)
     * @param MM length of mismatches for SV
     * @return true if seq1 matches seq2 with no more than MM number of mismatches
     */
    public boolean ismatch(String seq1,
                           String seq2,
                           int dir,
                           int MM) {
        if (instance().conf.y) {
            System.err.printf("    Matching two seqs %s %s %s %s\n", seq1, seq2, dir, MM);
        }
        seq2 = seq2.replaceAll("#|\\^", "");
        int mm = 0;
        for (int n = 0; n < seq1.length() && n < seq2.length(); n++) {
            if (seq1.charAt(n) != substr(seq2, dir * n - (dir == -1 ? 1 : 0), 1).charAt(0)) {
                mm++;
            }
        }
        return (mm <= MM && mm / (double)seq1.length() < 0.15);
    }

    /**
     * If any of the base is represented by 75%, it's a low complexity seq
     *
     * Low-complexity sequences are simple repeats such as ATATATATAT
     * or regions that are highly enriched for just one letter, e.g. AAACAAAAAAAAGAAAAAAC.
     * @param seq sequence
     * @return true if sequence has low complexity
     */
    public static boolean islowcomplexseq(String seq) {
        int len = seq.length();
        if (len == 0)
            return true;
        int ntcnt = 0;

        int a = count(seq, 'A');
        if (a > 0) ntcnt++;
        if (a / (double)len > 0.75)
            return true;

        int t = count(seq, 'T');
        if (t > 0) ntcnt++;
        if (t / (double)len > 0.75)
            return true;

        int g = count(seq, 'G');
        if (g > 0) ntcnt++;
        if (g / (double)len > 0.75)
            return true;

        int c = count(seq, 'C');
        if (c > 0) ntcnt++;
        if (c / (double)len > 0.75)
            return true;

        return ntcnt < 3;
    }

    /**
     * Counts the total count of base presence in the string
     * @param str string to count characters
     * @param chr character to seek in the string
     * @return number of characters in the string
     */
    static int count(String str, char chr) {
        int cnt = 0;
        for (int i = 0; i < str.length(); i++) {
            if (str.charAt(i) == chr) {
                cnt++;
            }
        }
        return cnt;
    }

    /**
     * Adjust the insertion position if necessary
     * @param bi starting position of insert
     * @param ins insert sequence
     * @param ref map of reference bases
     * @return Tuple of (int bi, String ins, int bi)
     */
    public static BaseInsertion adjInsPos(int bi, String ins, Map<Integer, Character> ref) {
        int n = 1;
        int len = ins.length();
        while (isEquals(ref.get(bi), ins.charAt(ins.length() - n))) {
            n++;
            if (n > len) {
                n = 1;
            }
            bi--;
        }
        if (n > 1) {
            ins = substr(ins, 1 - n) + substr(ins, 0, 1 - n);
        }
        return new BaseInsertion(bi, ins, bi);
    }

    /**
     * Find the insertion position
     * @param seq sequence
     * @param position start position of sequence
     * @param ref map of reference bases
     * @param dir direction
     * @param chr chromosome name
     * @return Tuple of (BI (insert starting position), INS (insert sequence), BI2 ( = BI))
     */
     BaseInsertion findbi(String seq,
                         int position,
                         Map<Integer, Character> ref,
                         final int dir,
                         String chr) {
        final int maxmm = 3; // maximum mismatches allowed
        final int dirExt = dir == -1 ? 1 : 0;
        int score = 0;
        int bi = 0;
        String ins = "";
        int bi2 = 0;

        for (int n = 6; n < seq.length(); n++) {
            if (position + 6 >= instance().chrLengths.get(chr)) {
                break;
            }
            int mm = 0;
            int i = 0;
            Set<Character> m = new HashSet<>();
            for (i = 0; i + n < seq.length(); i++) {
                if (position + dir * i - dirExt < 1) {
                    break;
                }
                if (position + dir * i - dirExt > instance().chrLengths.get(chr)) {
                    break;
                }
                if (isNotEquals(seq.charAt(i + n), ref.get(position + dir * i - dirExt))) {
                    mm++;
                } else {
                    m.add(seq.charAt(i + n));
                }
                if (mm > maxmm) {
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
                while (n + ept + 1 < seq.length() && (!isEquals(seq.charAt(n + ept), ref.get(position + ept * dir - dirExt))
                        || !isEquals(seq.charAt(n + ept + 1), ref.get(position + (ept + 1) * dir - dirExt)))) {
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
                        bi = position - 1 - extra.length();
                        ins = insert.toString();
                        bi2 = position - 1;
                        if (extra.length() == 0) {
                            BaseInsertion tpl = adjInsPos(bi, ins, ref);
                            bi = tpl.baseInsert;
                            ins = tpl.insertionSequence;
                            bi2 = tpl.baseInsert2;
                        }
                        return new BaseInsertion(bi, ins, bi2);
                    } else if (i - mm > score) {
                        bi = position - 1 - extra.length();
                        ins = insert.toString();
                        bi2 = position - 1;
                        score = i - mm;
                    }
                } else {
                    int s = -1;
                    if (extra.length() > 0) {
                        insert.append("&").append(extra);
                    } else {
                        while (s >= -n && isEquals(charAt(insert, s), ref.get(position + s))) {
                            s--;
                        }
                        if (s < -1) {
                            String tins = substr(insert.toString(), s + 1, 1 - s);
                            insert.delete(insert.length() + s + 1, insert.length());
                            insert.insert(0, tins);
                        }

                    }
                    if (mm == 0 && i + n == seq.length()) {
                        bi = position + s;
                        ins = insert.toString();
                        bi2 = position + s + extra.length();
                        if (extra.length() == 0) {
                            BaseInsertion tpl = adjInsPos(bi, ins, ref);
                            bi = tpl.baseInsert;
                            ins = tpl.insertionSequence;
                            bi2 = tpl.baseInsert2;
                        }
                        return new BaseInsertion(bi, ins, bi2);
                    } else if (i - mm > score) {
                        bi = position + s;
                        ins = insert.toString();
                        bi2 = position + s + extra.length();
                        score = i - mm;
                    }
                }
            }
        }
        if (bi2 == bi && ins.length() > 0 && bi != 0) {
            BaseInsertion tpl = adjInsPos(bi, ins, ref);
            bi = tpl.baseInsert;
            ins = tpl.insertionSequence;
        }
        return new BaseInsertion(bi, ins, bi2);
    }

    /**
     * Find breakpoint position in sequence
     * @param sequence sequence
     * @param startPosition $sp start position of sequence
     * @param ref map of reference bases
     * @param direction direction, reverse is -1, forward is 1.
     * @param chr chromosome name
     * @return breakpoint position
     */
     int findbp(String sequence,
                  int startPosition,
                  Map<Integer, Character> ref,
                  int direction,
                  String chr) {

        final int maxmm = 3; // maximum mismatches allowed
        int bp = 0;
        int score = 0;
        int idx = instance().chrLengths.containsKey(chr) ? instance().chrLengths.get(chr) : 0;
        for (int n = 0; n < instance().conf.indelsize; n++) {
            int mm = 0;
            int i = 0;
            Set<Character> m = new HashSet<>();
            for (i = 0; i < sequence.length(); i++) {
                if (startPosition + direction * n + direction * i < 1) {
                    break;
                }
                if (startPosition + direction * n + direction * i > idx) {
                    break;
                }
                if (isEquals(sequence.charAt(i), ref.get(startPosition + direction * n + direction * i))) {
                    m.add(sequence.charAt(i));
                } else {
                    mm++;
                }
                if (mm > maxmm - n / 100) {
                    break;
                }
            }
            if (m.size() < 3) {
                continue;
            }
            if (mm <= maxmm - n / 100 && i >= sequence.length() - 2 && i >= 8 + n / 10 && mm / (double)i < 0.12) {
                int lbp = startPosition + direction * n - (direction < 0 ? direction : 0);
                if (mm == 0 && i == sequence.length()) {
                    if (instance().conf.y) {
                        System.err.printf("  Findbp: %s %s %s %s %s\n", sequence, startPosition, lbp, mm, i);
                    }
                    return lbp;
                } else if (i - mm > score) {
                    bp = lbp;
                    score = i - mm;
                }
            }
        }
        if (instance().conf.y && bp != 0) {
            System.err.printf("  Findbp with mismatches: %s %s %s %s %s\n", sequence, startPosition, bp, direction, score);
        }
        return bp;
    }

    /**
     * Adjust the reference count.
     * @param tv non-reference variant
     * @param ref reference variant
     * @param len length for adhustment factor
     */
    void adjRefCnt(Variation tv,
                          Variation ref,
                          int len) {
        if (ref == null) {
            return;
        }

        if (instance().conf.y) {
            String refCnt = ref.varsCount != 0 ? String.valueOf(ref.varsCount) : "NA";
            System.err.printf("    AdjRefCnt: '+' %s %s %s %s Ref: %s\n",
                    ref.varsCount, tv.varsCount, ref.getDir(false), tv.getDir(false), refCnt);
            System.err.printf("    AdjRefCnt: '-' %s %s %s %s Ref: %s\n",
                    ref.varsCount, tv.varsCount, ref.getDir(true), tv.getDir(true), refCnt);
        }

        double f = tv.meanPosition != 0 ? (tv.meanPosition / (double)tv.varsCount - len + 1) / (tv.meanPosition / (double)tv.varsCount) : 0; // the adjustment factor
        if (f < 0) {
            return;
        }

        if (f > 1) {
            f = 1;
        }

        ref.varsCount -= (int) (f * tv.varsCount);
        ref.highQualityReadsCount -= (int) (f * tv.highQualityReadsCount);
        ref.lowQualityReadsCount -= (int) (f * tv.lowQualityReadsCount);
        ref.meanPosition -= f * tv.meanPosition;
        ref.meanQuality -= f * tv.meanQuality;
        ref.meanMappingQuality -= f * tv.meanMappingQuality;
        ref.numberOfMismatches -= f * tv.numberOfMismatches;
        ref.subDir(true, (int)(f * tv.getDir(true)));
        ref.subDir(false, (int)(f * tv.getDir(false)));
        correctCnt(ref);
    }

    /**
     * Adjust the reference by factor
     * @param ref reference variation
     * @param factor_f factor for adjustment counts of reference
     */
    void adjRefFactor(Variation ref, double factor_f) {
        if (ref == null) {
            return;
        }

        if (factor_f > 1) {
            factor_f = 1;
        }

        if (factor_f < -1) {
            return;
        }
        if (instance().conf.y) {
            System.err.printf("    AdjRefFactor: %s %s\n", ref.varsCount, factor_f);
        }
        int oldVarsCount = ref.varsCount;
        ref.varsCount -= (int) (factor_f * ref.varsCount);
        ref.highQualityReadsCount -= (int) (factor_f * ref.highQualityReadsCount);
        ref.lowQualityReadsCount -= (int) (factor_f * ref.lowQualityReadsCount);

        // Adjust mean mapping quality, mean quality and mean position only on number of changed counts
        double factorCnt = oldVarsCount != 0 ? Math.abs((ref.varsCount - oldVarsCount)) / (double) oldVarsCount : 1;
        // Factors must be the same sign
        if (factor_f < 0 && factorCnt > 0 || factor_f > 0 && factorCnt < 0){
            factorCnt = -factorCnt;
        }
        ref.meanPosition -= ref.meanPosition * factorCnt;
        ref.meanQuality -= ref.meanQuality * factorCnt;
        ref.meanMappingQuality -= ref.meanMappingQuality * factorCnt;

        ref.numberOfMismatches -= factor_f * ref.numberOfMismatches;
        ref.varsCountOnForward -= (int) (factor_f * ref.varsCountOnForward);
        ref.varsCountOnReverse -= (int) (factor_f * ref.varsCountOnReverse);

        correctCnt(ref);
    }

    /**
     * Add counts of variation by factor
     * @param vref variation  $ref
     * @param factor_f factor for adjustment counts of variation
     */
    void addVarFactor(Variation vref, double factor_f) {
        if (vref == null) {
            return;
        }

        if (factor_f < -1) {
            return;
        }

        vref.varsCount += (int) (factor_f * vref.varsCount);
        vref.highQualityReadsCount += (int) (factor_f * vref.highQualityReadsCount);
        vref.lowQualityReadsCount += (int) (factor_f * vref.lowQualityReadsCount);
        vref.meanPosition += factor_f * vref.meanPosition;
        vref.meanQuality += factor_f * vref.meanQuality;
        vref.meanMappingQuality += factor_f * vref.meanMappingQuality;
        vref.numberOfMismatches += factor_f * vref.numberOfMismatches;
        vref.varsCountOnForward += (int) (factor_f * vref.varsCountOnForward);
        vref.varsCountOnReverse += (int) (factor_f * vref.varsCountOnReverse);
    }


    /**
     * Given a variant sequence, find the mismatches and potential softclipping positions
     * @param ref map of reference bases
     * @param position position $p
     * @param wupseq sequence
     * @return MismatchResult contains mismatches lists and clipping positions
     */
     MismatchResult findMM5(Map<Integer, Character> ref,
                                  int position,
                                  String wupseq) {
        String seq = wupseq.replaceAll("#|\\^", "");
        int longmm = 3;
        List<Mismatch> mismatches = new ArrayList<>(); // $mm mismatches, mismatch positions, 5 or 3 ends
        int n = 0;
        int mn = 0;
        int mcnt = 0;
        StringBuilder str = new StringBuilder();
        List<Integer> sc5p = new ArrayList<>();
        while (isHasAndNotEquals(charAt(seq, -1 - n), ref, position - n) && mcnt < longmm) {
            str.insert(0, charAt(seq, -1 - n));
            mismatches.add(new Mismatch(str.toString(), position - n, 5));
            n++;
            mcnt++;
        }
        sc5p.add(position + 1);
        // Adjust clipping position if only one mismatch
        int misp = 0;
        Character misnt = null;
        if (str.length() == 1) {
            while (isHasAndEquals(charAt(seq, -1 - n), ref, position - n)) {
                n++;
                if (n != 0) {
                    mn++;
                }
            }
            if (mn > 1) {
                int n2 = 0;
                while (-1 - n - 1 - n2 >= 0
                        && isHasAndEquals(charAt(seq, -1 - n - 1 - n2), ref, position - n - 1 - n2)) {
                    n2++;
                }
                if (n2 > 2) {
                    sc5p.add(position - n - n2);
                    misp = position - n;
                    misnt = charAt(seq, -1 - n);
                    if (softClips5End.containsKey(position - n - n2)) {
                        softClips5End.get(position - n - n2).used = true;
                    }
                    mn += n2;
                } else {
                    sc5p.add(position - n);
                    if (softClips5End.containsKey(position - n)) {
                        softClips5End.get(position - n).used = true;
                    }
                }

            }
        }
        return new MismatchResult(mismatches, sc5p, mn, misp, misnt == null ? "" : misnt.toString());
    }

    /**
     * Given a variant sequence, find the mismatches and potential softclipping positions
     * @param ref map of reference bases
     * @param p position
     * @param sanpseq sequence
     * @return MismatchResult contains mismatches lists and clipping positions
     */
     MismatchResult findMM3(Map<Integer, Character> ref,
                                  int p,
                                  String sanpseq) {
        String seq = sanpseq.replaceAll("#|\\^", ""); // ~ s/#|\^//g;
        final int longmm = 3;
        // mismatches, mismatch positions, 5 or 3 ends
        List<Mismatch> mismatches = new ArrayList<>(); //$mm
        int n = 0;
        int mn = 0;
        int mcnt = 0;
        List<Integer> sc3p = new ArrayList<>();
        StringBuilder str = new StringBuilder();
        while (n < seq.length() && ref.containsKey(p + n) && isEquals(ref.get(p + n), seq.charAt(n))) {
            n++;
        }
        sc3p.add(p + n);
        int Tbp = p + n;
        while (mcnt <= longmm && n < seq.length() && isNotEquals(ref.get(p + n), seq.charAt(n))) {
            str.append(seq.charAt(n));
            mismatches.add(new Mismatch(str.toString(), Tbp, 3));
            n++;
            mcnt++;
        }
        // Adjust clipping position if only one mismatch
        int misp = 0;
        Character misnt = null;
        if (str.length() == 1) {
            while (n < seq.length() && isHasAndEquals(seq.charAt(n), ref, p + n)) {
                n++;
                if (n != 0) {
                    mn++;
                }
            }
            if (mn > 1) {
                int n2 = 0;
                while (n + n2 + 1 < seq.length() && isHasAndEquals(seq.charAt(n + n2 + 1), ref, p + n + 1 + n2)) {
                    n2++;
                }
                if (n2 > 2 && n + n2 + 1 < seq.length()) {
                    sc3p.add(p + n + n2);
                    misp = p + n;
                    misnt = seq.charAt(n);
                    if (softClips3End.containsKey(p + n + n2)) {
                        softClips3End.get(p + n + n2).used = true;
                    }
                    mn += n2;
                } else {
                    sc3p.add(p + n);
                    if (softClips3End.containsKey(p + n)) {
                        softClips3End.get(p + n).used = true;
                    }
                }
            }
        }
        return new MismatchResult(mismatches, sc3p, mn, misp, misnt == null ? "" : misnt.toString());
    }

    class MismatchResult {
        private final List<Mismatch> mismatches;
        private final List<Integer> scp;
        private final int nm;
        private final int misp;
        private final String misnt;

        public MismatchResult(List<Mismatch> mismatches, List<Integer> scp, int nm, int misp, String misnt) {
            this.mismatches = mismatches;
            this.scp = scp;
            this.nm = nm;
            this.misp = misp;
            this.misnt = misnt;
        }

        public List<Integer> getScp() {
            return scp;
        }

        public List<Mismatch> getMismatches() {
            return mismatches;
        }

        public int getNm() {
            return nm;
        }

        public int getMisp() {
            return misp;
        }

        public String getMisnt() {
            return misnt;
        }
    }

    class Mismatch {
        /**
         * mismatch string
         */
        String mismatchSequence; //$mm

        int mismatchPosition; //$mp
        /**
         * 3 or 5 end
         */
        int end; //# or 5

        public Mismatch(String mismatchSequence, int mismatchPosition, int end) {
            this.mismatchSequence = mismatchSequence;
            this.mismatchPosition = mismatchPosition;
            this.end = end;
        }
    }


    /**
     * Utility method for adjustMNP method with default number of mismatches = 3
     * @param sequence subsequence consensus sequence in soft-clipped reads  $seq
     * @param ref map of integer - characters (nucleotides) in reference sequence
     * @param dir direction (forward or reverse)
     * @param position key for MNP map     $p
     * @return true if sequence is matched to reference
     */
    public static boolean ismatchref(String sequence,
                                     Map<Integer, Character> ref,
                                     int position,
                                     int dir) {
        int MM = 3;
        return ismatchref(sequence, ref, position, dir, MM);
    }

    /**
     * Utility method for adjustMNP method
     * @param sequence subsequence consensus sequence in soft-clipped reads  $seq
     * @param ref map of integer - characters (nucleotides) in reference sequence
     * @param position key for MNP map     $p
     * @param dir direction (forward or reverse)
     * @param MM specific number of mismatches
     * @return true if sequence is matched to reference
     */
    public static boolean ismatchref(String sequence,
                                     Map<Integer, Character> ref,
                                     int position,
                                     int dir,
                                     int MM) {
        if (instance().conf.y) {
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

    /**
     * Subtract counts of one variation from another
     * @param vref variation where to subtract counts
     * @param tv variation which counts will be subtracted from the vref
     */
    void rmCnt(Variation vref,
               Variation tv) {
        vref.varsCount -= tv.varsCount;
        vref.highQualityReadsCount -= tv.highQualityReadsCount;
        vref.lowQualityReadsCount -= tv.lowQualityReadsCount;
        vref.meanPosition -= tv.meanPosition;
        vref.meanQuality -= tv.meanQuality;
        vref.meanMappingQuality -= tv.meanMappingQuality;
        vref.subDir(true, tv.getDir(true));
        vref.subDir(false, tv.getDir(false));
        correctCnt(vref);
    }
}
