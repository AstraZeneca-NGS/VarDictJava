package com.astrazeneca.vardict.modules;

import com.astrazeneca.vardict.Configuration;
import com.astrazeneca.vardict.data.Region;
import com.astrazeneca.vardict.Utils;
import com.astrazeneca.vardict.collection.Tuple;
import com.astrazeneca.vardict.collection.VariationMap;
import com.astrazeneca.vardict.data.Reference;
import com.astrazeneca.vardict.data.SamView;
import com.astrazeneca.vardict.variations.Sclip;
import com.astrazeneca.vardict.variations.Variation;
import htsjdk.samtools.SAMRecord;

import java.io.IOException;
import java.util.*;
import java.util.regex.Matcher;

import static com.astrazeneca.vardict.data.Patterns.*;
import static com.astrazeneca.vardict.modules.CigarUtils.getAlignedLength;
import static com.astrazeneca.vardict.modules.SAMFileParser.parseSAM;
import static com.astrazeneca.vardict.ReferenceResource.getREF;
import static com.astrazeneca.vardict.ReferenceResource.isLoaded;
import static com.astrazeneca.vardict.variations.VariationUtils.*;
import static com.astrazeneca.vardict.collection.Tuple.tuple;
import static com.astrazeneca.vardict.Utils.*;
import static java.lang.String.format;
import static java.util.Collections.singletonMap;

public class VariationRealigner {

    private static final Comparator<Tuple.Tuple2<Integer, Sclip>> COMP2 = new Comparator<Tuple.Tuple2<Integer, Sclip>>() {
        @Override
        public int compare(Tuple.Tuple2<Integer, Sclip> o1, Tuple.Tuple2<Integer, Sclip> o2) {
            int f = Integer.compare(o2._2.cnt, o1._2.cnt);
            if (f != 0)
                return f;
            return Integer.compare(o1._1, o2._1);
        }
    };

    private static final Comparator<Tuple.Tuple3<Integer, Sclip, Integer>> COMP3 = new Comparator<Tuple.Tuple3<Integer, Sclip, Integer>>() {
        @Override
        public int compare(Tuple.Tuple3<Integer, Sclip, Integer> o1, Tuple.Tuple3<Integer, Sclip, Integer> o2) {
            return Integer.compare(o2._3, o1._3);
        }
    };

    public static void realignIndels(Region region,
                                     Map<String, Integer> chrs,
                                     int rlen,
                                     Reference reference,
                                     Configuration conf,
                                     String bam,
                                     Map<Integer, VariationMap<String, Variation>> hash,
                                     Map<Integer, VariationMap<String, Variation>> iHash,
                                     Map<Integer, Integer> cov,
                                     Map<Integer, Sclip> sclip3,
                                     Map<Integer, Sclip> sclip5,
                                     String sample,
                                     Set<String> splice,
                                     String ampliconBasedCalling,
                                     Map<Integer, Map<String, Integer>> ins,
                                     Map<Integer, Map<String, Integer>> dels5) throws IOException {
        if (conf.y)
            System.err.println("Start Realigndel");
        realigndel(hash, dels5, cov, sclip5, sclip3, reference, region, chrs, rlen, bam, conf);
        if (conf.y)
            System.err.println("Start Realignins");
        realignins(hash, iHash, ins, cov, sclip5, sclip3, reference, region, chrs, rlen, conf);
        if (conf.y)
            System.err.println("Start Realignlgdel");
        realignlgdel(hash, cov, sclip5, sclip3, sample, splice, ampliconBasedCalling,
                reference, region, chrs, iHash, rlen, bam, conf);
        if (conf.y)
            System.err.println("Start Realignlgins");
        realignlgins(hash, iHash, cov, sclip5, sclip3, sample, splice, ampliconBasedCalling,
                reference, region, chrs, rlen, bam, conf);
        if (conf.y)
            System.err.println("Start Realignlgins30");
        realignlgins30(hash, iHash, cov, sclip5, sclip3, reference, region, chrs, rlen, bam, conf);
    }

    /**
     * Realign insertions
     * @param hash map of non-insertion variants produced by parseSAM method
     * @param iHash map of insertion variants produced by parseSAM method
     * @param ins insertion variants
     * @param cov coverage
     * @param sclip5 5' softclips
     * @param sclip3 3' softclips
     * @param reference object contains map of reference bases
     * @param region region contains chromosome name, start and end
     * @param chrs map of chromosome lengths
     * @param conf configuration
     */
    public static void realignins(Map<Integer, VariationMap<String, Variation>> hash,
                                  Map<Integer, VariationMap<String, Variation>> iHash,
                                  Map<Integer, Map<String, Integer>> ins,
                                  Map<Integer, Integer> cov,
                                  Map<Integer, Sclip> sclip5,
                                  Map<Integer, Sclip> sclip3,
                                  Reference reference,
                                  Region region,
                                  Map<String, Integer> chrs,
                                  int rlen,
                                  Configuration conf) {
        String chr = region.chr;
        Map<Integer, Character> ref = reference.referenceSequences;
        List<Object[]> tmp = fillTmp(ins);
        for (Object[] objects : tmp) {
            Integer p = (Integer)objects[0];
            String vn = (String)objects[1];
            Integer icnt = (Integer)objects[2];
            if (conf.y) {
                System.err.println(format("  Realign Ins: %s %s %s", p, vn, icnt));
            }
            String insert;
            Matcher mtch = BEGIN_PLUS_ATGC.matcher(vn);
            if (mtch.find()) {
                insert = mtch.group(1);
            } else {
                continue;
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
//            String newins = ""; // the adjacent insertion
//            mtch = CARET_ATGC_E.matcher(vn);
//            if (mtch.find()) {
//                newins = mtch.group(1);
//            }

            int newdel = 0; // the adjacent deletion
            try {
                newdel = toInt(vn);
            } catch (NumberFormatException ignore) {
            }
            String tn = vn.replaceFirst("^\\+", "")
                    .replaceFirst("&", "")
                    .replaceFirst("#", "")
                    .replaceFirst("\\^\\d+$", "")
                    .replaceFirst("\\^", "");

            int wustart = p - 100 - vn.length() + 1;
            if (wustart <= 1) {
                wustart = 1;
            }
            String wupseq = joinRef(ref, wustart, p) + tn; // 5prime flanking seq
            /*
            5' flanking region is a region of DNA that is adjacent to the 5' end of the gene.
            The 5' flanking region contains the promoter, and may contain enhancers or other protein binding sites.
            It is the region of DNA that is not transcribed into RNA.
             */
            Integer tend = chrs.get(chr);
            int sanend = p + vn.length() + 100;
            if (tend != null && tend < sanend) {
                sanend = tend;
            }
            // 3prime flanking seq
            String sanpseq = tn + joinRef(ref, p + extra.length() + 1 + compm.length() + newdel, sanend);
            /*
            3' flanking region is a region of DNA which is NOT copied into the mature mRNA, but which is present adjacent
            to 3' end of the gene. It was originally thought that the 3' flanking DNA was not transcribed at all,
            but it was discovered to be transcribed into RNA, but quickly removed during processing of the primary
            transcript to form the mature mRNA. The 3' flanking region often contains sequences which affect the
            formation of the 3' end of the message. It may also contain enhancers or other sites to which proteins may bind.
             */
            // mismatches, mismatch positions, 5 or 3 ends
            MismatchResult findMM3 = findMM3(ref, p + 1, sanpseq, insert.length() + compm.length(), sclip3);

            MismatchResult findMM5 = findMM5(ref, p + extra.length() + compm.length() + newdel, wupseq,
                    insert.length() + compm.length(), sclip5);

            List<Tuple.Tuple3<String, Integer, Integer>> mm3 = findMM3.getMm();
            List<Integer> sc3p = findMM3.getScp();
            int nm3 = findMM3.getNm();
            int misp3 = findMM3.getMisp();
            String misnt3 = findMM3.getMisnt();

            List<Tuple.Tuple3<String, Integer, Integer>> mm5 = findMM5.getMm();
            List<Integer> sc5p = findMM5.getScp();
            int nm5 = findMM5.getNm();
            int misp5 = findMM5.getMisp();
            String misnt5 = findMM5.getMisnt();

            List<Tuple.Tuple3<String, Integer, Integer>> mmm = new ArrayList<>(mm3);
            mmm.addAll(mm5);
            Variation vref = getVariation(iHash, p, vn);
            for (Tuple.Tuple3<String, Integer, Integer> tuple3 : mmm) {
                //mismatch nucleotide
                String mm = tuple3._1;
                //start position of clip that contains mm
                Integer mp = tuple3._2;
                //end (3 or 5)
                Integer me = tuple3._3;

                if (mm.length() > 1) {
                    mm = mm.charAt(0) + "&" + mm.substring(1);
                }
                if (!hash.containsKey(mp)) {
                    continue;
                }

                Variation tv = hash.get(mp).get(mm);
                if (tv == null) {
                    continue;
                }
                if (tv.cnt == 0) {
                    continue;
                }
                if (tv.qmean / tv.cnt < conf.goodq) {
                    continue;
                }
                if (tv.pmean / tv.cnt > (me == 3 ? nm3 + 4 : nm5 + 4)) { // opt_k;
                    continue;
                }
                if (tv.cnt >= icnt + insert.length() || tv.cnt / icnt >= 8) {
                    continue;
                }
                if (conf.y) {
                    System.err.printf("    insMM: %s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
                            mm, mp, me, nm3, nm5, vn, icnt, tv.cnt, tv.qmean, tv.pmean, cov.get(p));
                }
                // Adjust ref cnt so that AF won't > 1
                if (mp > p && me == 5) {
                    incCnt(cov, p, tv.cnt);
                }

                Variation lref = null;
                if (mp > p && me == 3 &&
                        hash.containsKey(p) &&
                        ref.containsKey(p) &&
                        hash.get(p).containsKey(ref.get(p).toString())) {

                    lref = hash.get(p).get(ref.get(p).toString());
                }
                adjCnt(vref, tv, lref, conf);
                hash.get(mp).remove(mm);
                if (hash.get(mp).isEmpty()) {
                    hash.remove(mp);
                }
            }
            if (misp3 != 0 && mm3.size() == 1 && hash.containsKey(misp3)) {
                hash.get(misp3).remove(misnt3);
            }
            if (misp5 != 0 && mm5.size() == 1 && hash.containsKey(misp5)) {
                hash.get(misp5).remove(misnt5);
            }
            for (Integer sc5pp : sc5p) {
                Sclip tv = sclip5.get(sc5pp);
                if (conf.y) {
                    System.err.printf("    55: %s %s VN: '%s'  5' seq: ^%s^\n", p, sc5pp, vn, wupseq);
                }
                if (tv != null && !tv.used) {
                    String seq = findconseq(tv, conf);
                    if (conf.y) {
                        System.err.printf("    ins5: %s %s %s %s %s %s\n", p, sc5pp, seq, wupseq, icnt, tv.cnt);
                    }
                    if (!seq.isEmpty() && ismatch(seq, wupseq, -1, conf.y)) {
                        if (conf.y) {
                            System.err.printf("    ins5: %s %s %s %s %s %s used\n", p, sc5pp, seq, wupseq, icnt, tv.cnt);
                        }
                        if (sc5pp > p) {
                            incCnt(cov, p, tv.cnt);
                        }
                        adjCnt(vref, tv, conf);
                        tv.used = true;
                    }
                }
            }
            for (Integer sc3pp : sc3p) {
                Sclip tv = sclip3.get(sc3pp);
                if (conf.y) {
                    System.err.printf("    33: %s %s VN: '%s'  3' seq: ^%s^\n", p, sc3pp, vn, sanpseq);
                }
                if (tv != null && !tv.used) {
                    String seq = findconseq(tv, conf);
                    if (conf.y) {
                        System.err.printf("    ins3: %s %s %s %s %s %s %s\n", p, sc3pp, seq, sanpseq, vn, icnt, tv.cnt);
                    }
                    if (!seq.isEmpty() && ismatch(seq, substr(sanpseq, sc3pp - p - 1), 1, conf.y)) {
                        if (conf.y) {
                            System.err.printf("    ins3: %s %s %s %s %s %s used\n", p, sc3pp, seq, vn, icnt, tv.cnt);
                        }
                        if (sc3pp <= p) {
                            incCnt(cov, p, tv.cnt);
                        }
                        Variation lref = null;
                        if (sc3pp > p &&
                                hash.containsKey(p) &&
                                ref.containsKey(p) &&
                                hash.get(p).containsKey(ref.get(p).toString())) {

                            lref = hash.get(p).get(ref.get(p).toString());
                        }
                        adjCnt(vref, tv, lref, conf);
                        tv.used = true;
                    }

                }
            }
            int first3 = sc3p.get(0);
            int first5 = sc5p.get(0);
            if (!sc3p.isEmpty() && !sc5p.isEmpty() && ref.containsKey(p)
                    && hash.containsKey(p) && hash.get(p).containsKey(ref.get(p).toString())
                    && first3 > first5 + 3
                    && first3 - first5 < rlen * 0.75) {
                adjRefFactor(hash.get(p).get(ref.get(p).toString()), (first3 - first5 - 1) / (double) rlen, conf.y);
                adjRefFactor(vref, -(first3 - first5 - 1)/ (double) rlen, conf.y);
            }
        }

        for (Object[] objects : tmp) {
            Integer p = (Integer)objects[0];
            String vn = (String)objects[1];
            if (!iHash.containsKey(p)) {
                continue;
            }
            Variation vref = iHash.get(p).get(vn);
            if (vref == null) {
                continue;
            }
            Matcher mtch = ATGSs_AMP_ATGSs_END.matcher(vn);
            if (mtch.find()) {
                String tn = mtch.group(1);
                Variation tref = iHash.get(p).get(tn);
                if (tref != null) {
                    if (vref.cnt < tref.cnt) {
                        adjCnt(tref, vref, getVariationMaybe(hash, p, ref.get(p)), conf);
                        iHash.get(p).remove(vn);
                    }
                }
            }
        }

    }

    /**
     * Realign deletions if already present in alignment
     * @param hash variant structure produced by parseSam method
     * @param dels5 deletion variants
     * @param cov coverage
     * @param sclip5 5' softclips
     * @param sclip3 3' softclips
     * @param reference object contains map of reference bases
     * @param region region contains chromosome name, start and end
     * @param chrs map of chromosome lengths
     * @param rlen read length
     * @param bam BAM file list
     * @param conf configuration
     * @throws IOException
     */
    public static void realigndel(Map<Integer, VariationMap<String, Variation>> hash,
                           Map<Integer, Map<String, Integer>> dels5,
                           Map<Integer, Integer> cov,
                           Map<Integer, Sclip> sclip5,
                           Map<Integer, Sclip> sclip3,
                           Reference reference,
                           Region region,
                           Map<String, Integer> chrs,
                           final int rlen,
                           String bam,
                           Configuration conf) throws IOException {
        String chr = region.chr;
        Map<Integer, Character> ref = reference.referenceSequences;
        String[] bams = bam.split(":");

        // In perl it doesn't commented, but it doesn't used
        // int longmm = 3; //Longest continued mismatches typical aligned at the end
        List<Object[]> tmp = fillTmp(dels5);

        for (Object[] objects : tmp) {
            Integer p = (Integer)objects[0];
            String vn = (String)objects[1];
            Integer dcnt = (Integer)objects[2];
            if (conf.y) {
                System.err.printf("  Realigndel for: %s %s %s cov: %s\n", p, vn, dcnt, cov.get(p));
            }
            final Variation vref = getVariation(hash, p, vn);
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
            mtch = CARET_ATGNC.matcher(vn);
            if (mtch.find()) {
                extrains = mtch.group(1);
            }
            String extra = "";
            mtch = BEGIN_MINUS_NUMBER_ANY.matcher(vn);
            if (mtch.find()) {
                extra = mtch.group(1).replaceAll("\\^|&|#", "");
            }

            int wustart = p - dellen - 100;
            if (wustart <= 1) {
                wustart = 1;
            }

            String wupseq = joinRef(ref, wustart, p - 1) + extra; // 5' flanking seq
            int sanend = p + 2 * dellen + 100;
            if (sanend > chrs.get(chr)) {
                sanend = chrs.get(chr);
            }
            // 3' flanking seq
            String sanpseq = extra + joinRef(ref, p + dellen + extra.length() - extrains.length(), sanend);
            // mismatches, mismatch positions, 5 or 3 ends
            MismatchResult r3 = findMM3(ref, p, sanpseq, dellen, sclip3);
            MismatchResult r5 = findMM5(ref, p + dellen + extra.length() - extrains.length() - 1, wupseq, dellen, sclip5);

            List<Tuple.Tuple3<String, Integer, Integer>> mm3 = r3.getMm();
            List<Integer> sc3p = r3.getScp();
            int nm3 = r3.getNm();
            int misp3 = r3.getMisp();
            String misnt3 = r3.getMisnt();

            List<Tuple.Tuple3<String, Integer, Integer>> mm5 = r5.getMm();
            List<Integer> sc5p = r5.getScp();
            int nm5 = r5.getNm();
            int misp5 = r5.getMisp();
            String misnt5 = r5.getMisnt();
            if (conf.y) {
                System.err.printf("  Mismatches: misp3: %s-%s misp5: %s-%s sclip3: %s sclip5: %s\n",
                        misp3, misnt3, misp5, misnt5, Utils.toString(sc3p), Utils.toString(sc5p));
            }

            List<Tuple.Tuple3<String, Integer, Integer>> mmm = new ArrayList<>(mm3);
            mmm.addAll(mm5);
            for (Tuple.Tuple3<String, Integer, Integer> tuple : mmm) {
                String mm = tuple._1;
                Integer mp = tuple._2;
                Integer me = tuple._3;
                if (mm.length() > 1) {
                    mm = mm.charAt(0) + "&" + mm.substring(1);
                }
                if (hash.get(mp) == null) {
                    continue;
                }
                Variation tv = hash.get(mp).get(mm);
                if (tv == null) {
                    continue;
                }
                if (tv.cnt == 0) {
                    continue;
                }
                if (tv.qmean / tv.cnt < conf.goodq) {
                    continue;
                }

                if (tv.pmean / tv.cnt > (me == 3 ? nm3 + 4 : nm5 + 4)) {
                    continue;
                }
                if (tv.cnt >= dcnt + dellen || tv.cnt / dcnt >= 8) {
                    continue;
                }
                if (conf.y) {
                    System.err.printf("  Realigndel Adj: %s %s %s %s %s %s %s %s cov: %s\n",
                            mm, mp, me, nm3, nm5, p, tv.cnt, tv.qmean, cov.get(p));
                }
                // Adjust ref cnt so that AF won't > 1
                if (mp > p && me == 5) {
                    double f = tv.pmean != 0 ? (mp - p) / (tv.pmean / (double)tv.cnt) : 1;
                    if (f > 1) {
                        f = 1;
                    }
                    incCnt(cov, p, (int)(tv.cnt * f));
                    adjRefCnt(tv, getVariationMaybe(hash, p, ref.get(p)), dellen);
                }
                Variation lref = (mp > p && me == 3)
                        ? (hash.containsKey(p) && hash.get(p).containsKey(ref.get(p).toString())
                        ? hash.get(p).get(ref.get(p).toString())
                        : null)
                        : null;
                adjCnt(vref, tv, lref, conf);
                hash.get(mp).remove(mm);
                if (hash.get(mp).isEmpty()) {
                    hash.remove(mp);
                }
                if (conf.y) {
                    System.err.printf("  Realigndel AdjA: %s %s %s %s %s %s %s %s cov: %s\n",
                            mm, mp, me, nm3, nm5, p, tv.cnt, tv.qmean, cov.get(p));
                }
            }
            if (misp3 != 0 && mm3.size() == 1 && hash.containsKey(misp3)) {
                hash.get(misp3).remove(misnt3);
            }
            if (misp5 != 0 && mm5.size() == 1 && hash.containsKey(misp5)) {
                hash.get(misp5).remove(misnt5);
            }

            for (Integer sc5pp : sc5p) {
                if (sclip5.containsKey(sc5pp) && !sclip5.get(sc5pp).used) {
                    Sclip tv = sclip5.get(sc5pp);
                    String seq = findconseq(tv, conf);
                    if (conf.y) {
                        System.err.printf("  Realigndel 5: %s %s Seq: '%s' %s %s %s %s %s cov: %s\n",
                                p, sc5pp, seq, new StringBuilder(wupseq).reverse(), tv.cnt, dcnt, vn, p, cov.get(p));
                    }
                    if (!seq.isEmpty() && ismatch(seq, wupseq, -1, conf.y)) {
                        if (sc5pp > p) {
                            incCnt(cov, p, tv.cnt);
                        }
                        adjCnt(vref, tv, conf);
                        sclip5.get(sc5pp).used = true;
                        if (conf.y) {
                            System.err.printf("  Realigndel 5: %s %s %s %s %s %s %s %s used cov: %s\n",
                                    p, sc5pp, seq, new StringBuilder(wupseq).reverse(), tv.cnt, dcnt, vn, p, cov.get(p));
                        }
                    }
                }
            }

            for (Integer sc3pp : sc3p) {
                if (sclip3.containsKey(sc3pp) && !sclip3.get(sc3pp).used) {
                    Sclip tv = sclip3.get(sc3pp);
                    String seq = findconseq(tv, conf);
                    if (conf.y) {
                        System.err.printf("  Realigndel 3: %s %s seq '%s' %s %s %s %s %s %s %s\n",
                                p, sc3pp, seq, sanpseq, tv.cnt, dcnt, vn, p, dellen, substr(sanpseq, sc3pp - p));
                    }
                    if (!seq.isEmpty() && ismatch(seq, substr(sanpseq, sc3pp - p), 1, conf.y)) {
                        if (conf.y) {
                            System.err.printf("  Realigndel 3: %s %s %s %s %s %s %s %s used\n",
                                    p, sc3pp, seq, sanpseq, tv.cnt, dcnt, vn, p);
                        }
                        if (sc3pp <= p) {
                            incCnt(cov, p, tv.cnt);
                        }
                        Variation lref = sc3pp <= p ? null : getVariationMaybe(hash, p, ref.get(p));
                        adjCnt(vref, tv, lref, conf);
                        sclip3.get(sc3pp).used = true;
                    }
                }
            }
            // In perl it is commented too
            // int pe = p + dellen + extra.length() + compm.length();
            int pe = p + dellen + extra.length() - extrains.length();
            Variation h = getVariationMaybe(hash, p, ref.get(p));
            if (bams != null && bams.length > 0
                    && pe - p >= 5
                    && pe - p < rlen - 10
                    && h != null && h.cnt != 0
                    && noPassingReads(chr, p, pe, bams, conf)
                    && vref.cnt > 2 * h.cnt * (1 - (pe - p) / (double)rlen)) {

                adjCnt(vref, h, h, conf);
            }
        }

        for (int i = tmp.size() - 1; i >= 0; i--) {
            Object[] os = tmp.get(i);
            int p = (Integer)os[0];
            String vn = (String)os[1];
            if (!hash.containsKey(p)) {
                continue;
            }
            Variation vref = hash.get(p).get(vn);
            if (vref == null) {
                continue;
            }
            Matcher matcher = MINUS_NUMBER_AMP_ATGCs_END.matcher(vn);
            if (matcher.find()) {
                String tn = matcher.group(1);
                Variation tref = hash.get(p).get(tn);
                if (tref != null) {
                    if (vref.cnt < tref.cnt) {
                        adjCnt(tref, vref, conf);
                        hash.get(p).remove(vn);
                    }
                }
            }
        }
    }

    /**
     * Realign large deletions that are not present in alignment
     * @param hash variant structure produced by parseSam method
     * @param cov coverage
     * @param sclip5 5' softclips
     * @param sclip3 3' softclips
     * @param reference object contains map of reference bases
     * @param region region contains chromosome name, start and end
     * @param chrs map of chromosome lengths
     * @param rlen read length
     * @param bam BAM file list
     * @param conf configuration
     * @throws IOException
     */
    public static void realignlgdel(Map<Integer, VariationMap<String, Variation>> hash,
                             Map<Integer, Integer> cov,
                             Map<Integer, Sclip> sclip5,
                             Map<Integer, Sclip> sclip3,
                             String sample, Set<String> splice,
                             String ampliconBasedCalling,
                             Reference reference,
                             Region region,
                             Map<String, Integer> chrs,
                             Map<Integer, VariationMap<String, Variation>> iHash,
                             final int rlen,
                             String bam,
                             Configuration conf) throws IOException {

        final int longmm = 3;
        String chr = region.chr;
        Map<Integer, Character> ref = reference.referenceSequences;

        List<Tuple.Tuple2<Integer, Sclip>> tmp = new ArrayList<>();
        for (Map.Entry<Integer, Sclip> ent5 : sclip5.entrySet()) {
            int p = ent5.getKey();
            if (p < region.start - Configuration.EXTENSION
                    || p > region.end + Configuration.EXTENSION) {
                continue;
            }
            tmp.add(tuple(ent5.getKey(), ent5.getValue()));
        }
        Collections.sort(tmp, COMP2);
        for (Tuple.Tuple2<Integer, Sclip> t : tmp) {
            int p = t._1;
            Sclip sc5v = t._2;
            int cnt = t._2.cnt;
            if (sc5v.used) {
                continue;
            }
            String seq = findconseq(sc5v, conf);
            if (seq.isEmpty()) {
                continue;
            }
            if (B_A7.matcher(seq).find() || B_T7.matcher(seq).find()) {
                continue;
            }
            if (seq.length() < 7) {
                continue;
            }
            if (islowcomplexseq(seq)) {
                continue;
            }
            int bp = findbp(seq, p - 5, ref, conf.indelsize, -1, chr, chrs, conf.y);
            final int dellen = p - bp;
            if (bp == 0) {
                continue;
            }
            // TODO: these parts in realignlgins and realignlgdel must be moved in appropriate parts
            // after implementing changes for realignment
            if (bp < region.start) {
                int tts = bp - rlen;
                int tte = bp + rlen;
                if (bp + rlen >= region.start) {
                    tte = region.start - 1;
                }
                Region modifiedRegion = Region.newModifiedRegion(region, tts, tte);
                if(!isLoaded(chr, tts, tte, reference)) {
                    getREF(modifiedRegion, chrs, conf, rlen, reference);
                }
                parseSAM(modifiedRegion,  bam, chrs, sample, splice, ampliconBasedCalling,
                        rlen, reference, conf, hash, iHash, cov, sclip5, sclip3, true);
            }

            StringBuilder extra = new StringBuilder();
            int en = 0;
            while (!Character.valueOf(seq.charAt(en)).equals(ref.get(bp - en - 1)) && en < seq.length()) {
                extra.append(seq.charAt(en));
                en++;
            }
            final String gt;
            if (extra.length() > 0) {
                extra = extra.reverse();
                gt = "-" + dellen + "&" + extra;
                bp -= extra.length();
            } else {
                gt = String.valueOf(-dellen);
            }
            if (conf.y) {
                System.err.printf("  Found Realignlgdel: %s %s 5' %s %s %s\n", bp, gt, p, seq, cnt);
            }
            final Variation tv = getVariation(hash, bp, gt);
            tv.qstd = true; // more accurate implementation lat
            tv.pstd = true; // more accurate implementation lat
            for (int tp = bp; tp < bp + dellen; tp++) {
                incCnt(cov, tp, sc5v.cnt);
            }
            adjCnt(tv, sc5v, conf);
            sc5v.used = true;

            // Work on the softclipped read at 3'
            int n = 0;
            while (ref.containsKey(bp + n)
                    && ref.containsKey(bp + dellen + n)
                    && isEquals(ref.get(bp + n), ref.get(bp + dellen + n))) {
                n++;
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
                while (ref.containsKey(bp + n)
                        && ref.containsKey(bp + dellen + n)
                        && isEquals(ref.get(bp + n), ref.get(bp + dellen + n))) {
                    n++;
                }
                sc3p = bp + n;
            }
            if (sclip3.containsKey(sc3p) && !sclip3.get(sc3p).used) {
                Sclip sclip = sclip3.get(sc3p);
                if (sc3p > bp) {
                    adjCnt(tv, sclip, getVariationMaybe(hash, bp, ref.get(bp)), conf);
                } else {
                    adjCnt(tv, sclip, conf);
                }

                if (sc3p == bp) {
                    for (int tp = bp; tp < bp + dellen; tp++) {
                        incCnt(cov, tp, sclip3.get(sc3p).cnt);
                    }
                }
                for (int ip = bp + 1; ip < sc3p; ip++) {
                    Variation vv = getVariation(hash, ip, ref.get(dellen + ip).toString());
                    rmCnt(vv, sclip);
                    if (vv.cnt == 0) {
                        hash.get(ip).remove(ref.get(dellen + ip).toString());
                    }
                    if (hash.get(ip).size() == 0) {
                        hash.remove(ip);
                    }
                }
                sclip.used = true;
            }
            Map<Integer, Map<String, Integer>> dels5 = singletonMap(bp, singletonMap(gt, tv.cnt));
            realigndel(hash, dels5, cov, sclip5, sclip3, reference, region, chrs, rlen, bam, conf);
            if (conf.y) {
                System.err.printf("  Found lgdel done: %s %s %s 5' %s %s\n\n", bp, gt, p, seq, tv.cnt);
            }
        }

        tmp = new ArrayList<>();
        for (Map.Entry<Integer, Sclip> ent3 : sclip3.entrySet()) {
            int p = ent3.getKey();
            if (p < region.start - Configuration.EXTENSION
                    || p > region.end + Configuration.EXTENSION) {
                continue;
            }
            tmp.add(tuple(ent3.getKey(), ent3.getValue()));
        }
        Collections.sort(tmp, COMP2);
        for (Tuple.Tuple2<Integer, Sclip> t : tmp) {
            int p = t._1;
            final Sclip sc3v = t._2;
            final int cnt = sc3v.cnt;
            if (sc3v.used) {
                continue;
            }
            String seq = findconseq(sc3v, conf);
            if (seq.isEmpty()) {
                continue;
            }
            if (B_A7.matcher(seq).find() || B_T7.matcher(seq).find()) {
                continue;
            }
            if (seq.length() < 7) {
                continue;
            }
            if (islowcomplexseq(seq)) {
                continue;
            }
            int bp = findbp(seq, p + 5, ref, conf.indelsize, 1, chr, chrs, conf.y);
            final int dellen = bp - p;
            if (bp == 0) {
                continue;
            }

            if (bp > region.end) {
//              my ($tts, $tte) = ($bp - $RLEN <= $END ? $END + 1 : $bp - $RLEN, $bp + $RLEN);
//		        getREF($chr, $tts, $tte, $REF, $RLEN) unless( isLoaded( $chr, $tts, $tte, $REF ) );
//		        parseSAM($chr, $bp - $RLEN <= $END ? $END + 1 : $bp - $RLEN, $bp + $RLEN, $bams, $REF, $hash, $cov, $sclip5, $sclip3, 1);
                int tts = bp - rlen;
                int tte = bp + rlen;
                if (bp - rlen <= region.end) {
                    tts = region.end + 1;
                }
                Region modifiedRegion = Region.newModifiedRegion(region, tts, tte);
                if(!isLoaded(chr, tts, tte, reference)) {
                    getREF(modifiedRegion, chrs, conf, rlen, reference);
                }
                parseSAM(modifiedRegion,  bam, chrs, sample, splice, ampliconBasedCalling,
                        rlen, reference, conf, hash, iHash, cov, sclip5, sclip3, true);
            }


            StringBuilder extra = new StringBuilder();
            int en = 0;
            while (en < seq.length() && isNotEquals(seq.charAt(en), ref.get(bp + en))) {
                extra.append(seq.charAt(en));
                en++;
            }
            bp = p; // Set it to 5'
            final String gt;
            if (extra.length() > 0) {
                gt = "-" + dellen + "&" + extra;
            } else {
                gt = "-" + dellen;
                while (isEquals(ref.get(bp - 1), ref.get(bp + dellen - 1))) {
                    bp--;
                }
            }
            if (conf.y) {
                System.err.printf("  Found Realignlgdel: %s %s 3' %s %s %s\n", bp, gt, p, seq, cnt);
            }
            Variation tv = getVariation(hash, bp, gt);
            tv.qstd = true; // more accurate implementation later
            tv.pstd = true; // more accurate implementation later
            for (int tp = bp; tp < bp + dellen; tp++) {
                incCnt(cov, tp, sc3v.cnt);
            }
            sc3v.pmean += dellen * sc3v.cnt;
            adjCnt(tv, sc3v, conf);
            sc3v.used = true;

            Map<Integer, Map<String, Integer>> dels5 = new HashMap<>();
            HashMap<String, Integer> map = new HashMap<>();
            map.put(gt, tv.cnt);
            dels5.put(bp, map);
            realigndel(hash, dels5, cov, sclip5, sclip3, reference, region, chrs, rlen, bam, conf);
            if (conf.y) {
                System.err.printf("  Found lgdel: %s %s $p 3' %s %s\n\n", bp, gt, p, tv.cnt);
            }
        }
        if (conf.y) {
            System.err.println("  Done: Realignlgdel\n");
        }
    }

    /**
     * Realign large insertions that are not present in alignment
     * @param hash map of non-insertion variants produced by parseSAM method
     * @param iHash map of insertion variants produced by parseSAM method
     * @param cov coverage
     * @param sclip5 5' softclips
     * @param sclip3 3' softclips
     * @param reference object contains map of reference bases
     * @param region region contains chromosome name, start and end
     * @param chrs map of chromosome lengths
     * @param rlen read length
     * @param bam BAM file list
     * @param conf configuration
     * @throws IOException
     */
    public static void realignlgins(Map<Integer, VariationMap<String, Variation>> hash,
                             Map<Integer, VariationMap<String, Variation>> iHash,
                             Map<Integer, Integer> cov,
                             Map<Integer, Sclip> sclip5,
                             Map<Integer, Sclip> sclip3,
                             String sample, Set<String> splice,
                             String ampliconBasedCalling,
                             Reference reference,
                             Region region,
                             Map<String, Integer> chrs,
                             int rlen,
                             String bam,
                             Configuration conf) throws IOException {
        String chr = region.chr;
        Map<Integer, Character> ref = reference.referenceSequences;
        String[] bams = bam.split(":");

        List<Tuple.Tuple2<Integer, Sclip>> tmp = new ArrayList<>();
        for (Map.Entry<Integer, Sclip> ent5 : sclip5.entrySet()) {
            int p = ent5.getKey();
            if (p < region.start - Configuration.EXTENSION
                    || p > region.end + Configuration.EXTENSION) {
                continue;
            }
            tmp.add(tuple(ent5.getKey(), ent5.getValue()));
        }
        //sort by descending cnt
        Collections.sort(tmp, COMP2);

        for (Tuple.Tuple2<Integer, Sclip> t : tmp) {
            int p = t._1;
            Sclip sc5v = t._2;
            //already been used in
            if (sc5v.used) {
                continue;
            }
            String seq = findconseq(sc5v, conf);
            if (seq.isEmpty()) {
                continue;
            }
            if (conf.y) {
                System.err.println("  Working lgins: 5: " + p + " " + seq);
            }
            if (B_A8.matcher(seq).find() || B_T8.matcher(seq).find()) {
                continue;
            }
            if (seq.length() < 12) {
                continue;
            }
            if (islowcomplexseq(seq)) {
                continue;
            }
            Tuple.Tuple3<Integer, String, Integer> tpl = findbi(seq, p, ref, -1, chr, chrs);
            final int bi = tpl._1;
            final String ins = tpl._2;
            if (bi == 0) {
                continue;
            }

            // For large insertions
            if (bi > region.end) {
                //my ($tts, $tte) = ($bi - $RLEN <= $END ? $END + 1 : $bi - $RLEN, $bi + $RLEN);
                int tts = region.end + 1;
                int tte = bi + rlen;
                if (bi - rlen <= region.end) {
                    tts = bi - rlen;
                }
                Region modifiedRegion = Region.newModifiedRegion(region, tts, tte);
                if(!isLoaded(chr, tts, tte, reference)) {
                    getREF(modifiedRegion, chrs, conf, rlen, reference);
                }
                parseSAM(modifiedRegion,  bam, chrs, sample, splice, ampliconBasedCalling,
                        rlen, reference, conf, hash, iHash, cov, sclip5, sclip3, true);
            }
            if (conf.y) {
                System.err.printf("  Found candidate lgins from 5: %s +%s %s %s\n", bi, ins, p, seq);
            }
            final Variation iref = getVariation(iHash, bi, "+" + ins);
            iref.pstd = true;
            iref.qstd = true;
            adjCnt(iref, sc5v, conf);
            boolean rpflag = true; // A flag to indicate whether an insertion is a repeat
            for (int i = 0; i < ins.length(); i++) {
                if (!isEquals(ref.get(bi + 1 + i), ins.charAt(i))) {
                    rpflag = false;
                    break;
                }
            }
            incCnt(cov, bi, sc5v.cnt);
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
                    Variation tvr = getVariation(hash, pii, tnt.toString());
                    adjCnt(tvr, tv, conf);
                    tvr.pstd = true;
                    tvr.qstd = true;
                    incCnt(cov, pii, tv.cnt);
                }
            }
            sc5v.used = true;
            Map<Integer, Map<String, Integer>> tins = singletonMap(bi, singletonMap("+" + ins, iref.cnt));
            realignins(hash, iHash, tins, cov, sclip5, sclip3, reference, region, chrs, rlen, conf);
            Variation mref = getVariationMaybe(hash, bi, ref.get(bi));
            if (rpflag && bams.length > 0 && ins.length() >= 5
                    && ins.length() < rlen - 10
                    && mref != null && mref.cnt != 0
                    && noPassingReads(chr, bi, bi + ins.length(), bams, conf)
                    && iref.cnt > 2 * mref.cnt) {
                adjCnt(iref, mref, mref, conf);
            }
        }

        tmp = new ArrayList<>();
        for (Map.Entry<Integer, Sclip> ent3 : sclip3.entrySet()) {
            int p = ent3.getKey();
            if (p < region.start - Configuration.EXTENSION
                    || p > region.end + Configuration.EXTENSION) {
                continue;
            }
            tmp.add(tuple(ent3.getKey(), ent3.getValue()));
        }
        Collections.sort(tmp, COMP2);
        for (Tuple.Tuple2<Integer, Sclip> t : tmp) {
            final int p = t._1;
            final Sclip sc3v = t._2;
            if (sc3v.used) {
                continue;
            }
            String seq = findconseq(sc3v, conf);
            if (seq.isEmpty()) {
                continue;
            }
            if (conf.y) {
                System.err.println("  Working lgins 3: " + seq);
            }
            if (B_A7.matcher(seq).find() || B_T7.matcher(seq).find()) {
                continue;
            }
            if (seq.length() < 12) {
                continue;
            }
            if (islowcomplexseq(seq)) {
                continue;
            }
            Tuple.Tuple3<Integer, String, Integer> tpl = findbi(seq, p, ref, 1, chr, chrs);
            final int bi = tpl._1;
            final String ins = tpl._2;

            if (bi == 0) {
                continue;
            }

            // For large insertions
            if (bi < region.start) {
                //my ($tts, $tte) = ($bi - $RLEN, $bi + $RLEN >= $START ? $START - 1 : $bi + $RLEN);
                //getREF($chr, $tts, $tte, $REF, $RLEN) unless( isLoaded( $chr, $tts, $tte, $REF ) );
                //parseSAM($chr, $bi - $RLEN, $bi + $RLEN >= $START ? $START - 1 : $bi + $RLEN, $bams, $REF, $hash, $cov, $sclip5, $sclip3, 1);
                int tts = bi - rlen;
                int tte = bi + rlen;
                if (bi + rlen >= region.start) {
                    tte = region.start - 1;
                }
                Region modifiedRegion = Region.newModifiedRegion(region, tts, tte);
                if(!isLoaded(chr, tts, tte, reference)) {
                    getREF(modifiedRegion, chrs, conf, rlen, reference);
                }
                parseSAM(modifiedRegion,  bam, chrs, sample, splice, ampliconBasedCalling,
                        rlen, reference, conf, hash, iHash, cov, sclip5, sclip3, true);
            }

            final Variation iref = getVariation(iHash, bi, "+" + ins);
            iref.pstd = true;
            iref.qstd = true;
            final Variation lref = getVariationMaybe(hash, bi, ref.get(bi));
            adjCnt(iref, sc3v, lref, conf);
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
                    Variation vref = getVariation(hash, pii, tnt.toString());
                    adjCnt(vref, tv, conf);
                    vref.pstd = true;
                    vref.qstd = true;
                    incCnt(cov, pii, tv.cnt);
                }
            }
            sc3v.used = true;
            Map<Integer, Map<String, Integer>> tins = singletonMap(bi, singletonMap("+" + ins, iref.cnt));
            realignins(hash, iHash, tins, cov, sclip5, sclip3, reference, region, chrs, rlen, conf);
            Variation mref = getVariationMaybe(hash, bi, ref.get(bi));
            if (rpflag && bams.length > 0 && ins.length() >= 5 && ins.length() < rlen - 10
                    && mref != null && mref.cnt != 0
                    && noPassingReads(chr, bi, bi + ins.length(), bams, conf)
                    && iref.cnt > 2 * mref.cnt) {

                adjCnt(iref, mref, mref, conf);
            }
        }
    }

    /**
     * this will try to realign large insertions (typically larger than 30bp)
     * @param hash map of non-insertion variants produced by parseSAM method
     * @param iHash map of insertion variants produced by parseSAM method
     * @param cov coverage
     * @param sclip5 5' softclips
     * @param sclip3 3' softclips
     * @param reference object contains map of reference bases
     * @param region region contains chromosome name, start and end
     * @param chrs map of chromosome lengths
     * @param rlen max read length
     * @param bam BAM file list
     * @param conf configuration
     * @throws IOException
     */
    static void realignlgins30(Map<Integer, VariationMap<String, Variation>> hash,
                               Map<Integer, VariationMap<String, Variation>> iHash,
                               Map<Integer, Integer> cov,
                               Map<Integer, Sclip> sclip5,
                               Map<Integer, Sclip> sclip3,
                               Reference reference,
                               Region region,
                               Map<String, Integer> chrs,
                               int rlen,
                               String bam,
                               Configuration conf) throws IOException {
        String chr = region.chr;
        Map<Integer, Character> ref = reference.referenceSequences;
        String[] bams = bam.split(":");

        List<Tuple.Tuple3<Integer, Sclip, Integer>> tmp5 = new ArrayList<>();
        for (Map.Entry<Integer, Sclip> ent5 : sclip5.entrySet()) {
            //ent5: (position, soft clipping 5' variant)
            int p = ent5.getKey();
            if (p < region.start - Configuration.EXTENSION
                    || p > region.end + Configuration.EXTENSION) {
                continue;
            }
            tmp5.add(tuple(ent5.getKey(), ent5.getValue(), ent5.getValue().cnt));
        }
        Collections.sort(tmp5, COMP3);

        List<Tuple.Tuple3<Integer, Sclip, Integer>> tmp3 = new ArrayList<>();
        for (Map.Entry<Integer, Sclip> ent3 : sclip3.entrySet()) {
            int p = ent3.getKey();
            if (p < region.start - Configuration.EXTENSION
                    || p > region.end + Configuration.EXTENSION) {
                continue;
            }
            tmp3.add(tuple(ent3.getKey(), ent3.getValue(), ent3.getValue().cnt));
        }
        Collections.sort(tmp3, COMP3);

        for (Tuple.Tuple3<Integer, Sclip, Integer> t5 : tmp5) {
            final int p5 = t5._1;
            final Sclip sc5v = t5._2;
            final int cnt5 = t5._3;
            if (sc5v.used) {
                continue;
            }

            for (Tuple.Tuple3<Integer, Sclip, Integer> t3 : tmp3) {
                final int p3 = t3._1;
                final Sclip sc3v = t3._2;
                final int cnt3 = t3._3;
                if (sc3v.used) {
                    continue;
                }
                if (p5 - p3 > rlen / 1.5) {
                    continue;
                }
                if (p3 - p5 > rlen - 10) { // if they're too far away, don't even try
                    continue;
                }
                final String seq5 = findconseq(sc5v, conf);
                final String seq3 = findconseq(sc3v, conf);
                //next until at least one of consensus sequences has length > 10
                if (seq5.length() <= 10 || seq3.length() <= 10) {
                    continue;
                }
                if (conf.y) {
                    System.err.printf("  Working lgins30: %s %s 3: %s %s 5: %s %s\n",
                            p3, p5, seq3, cnt3, new StringBuilder(seq5).reverse(), cnt5);
                }
                Tuple.Tuple3<Integer, Integer, Integer> tpl = find35match(seq5, seq3, p5, p3, ref);
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
                if (islowcomplexseq(ins)) {
                    if (conf.y) {
                        System.err.println("  Discard low complex insertion found " + ins + ".");
                    }
                    continue;
                }
                int bi = 0;
                Variation vref;
                if (conf.y) {
                    System.err.printf("  Found candidate lgins30: %s %s %s\n", p3, p5, ins);
                }
                if (p5 > p3) {
                    if (seq3.length() > ins.length()
                            && !ismatch(substr(seq3, ins.length()), joinRef(ref, p5, p5 + seq3.length() - ins.length() + 2), 1, conf.y)) {
                        continue;
                    }
                    if (seq5.length() > ins.length()
                            && !ismatch(substr(seq5, ins.length()), joinRef(ref, p3 - seq5.length() - ins.length() - 2, p3 - 1), -1, conf.y)) {
                        continue;
                    }
                    if (conf.y) {
                        System.err.printf("  Found lgins30 complex: %s %s %s %s\n", p3, p5, ins.length(), ins);
                    }
                    String tmp = joinRef(ref, p3, p5 - 1);
                    if (tmp.length() > ins.length()) { // deletion is longer
                        ins = (p3 - p5) + "^" + ins;
                        bi = p3;
                        vref = getVariation(hash, p3, ins);
                    } else if (tmp.length() < ins.length()) {
                        // In perl it isn't commented, but not used
//                        int p3s = p3 + tmp.length();
//                        int p3e = p3s + seq3.length() - ins.length() + 2;
                        ins = substr(ins, 0, ins.length() - tmp.length()) + "&" + substr(ins, p3 - p5);
                        ins = "+" + ins;
                        bi = p3 - 1;
                        vref = getVariation(iHash, bi, ins);
                    } else { // long MNP
                        ins = "-" + ins.length() + "^" + ins;
                        bi = p3;
                        vref = getVariation(hash, p3, ins);
                    }
                } else {
                    if (seq3.length() > ins.length()
                            && !ismatch(substr(seq3, ins.length()), joinRef(ref, p5, p5 + seq3.length() - ins.length() + 2), 1, conf.y)) {
                        continue;
                    }
                    if (seq5.length() > ins.length()
                            && !ismatch(substr(seq5, ins.length()), joinRef(ref, p3 - (seq5.length() - ins.length()) - 2, p3 - 1), -1, conf.y)) {
                        continue;
                    }
                    String tmp = ins.length() > p3 - p5
                            ? joinRef(ref, p5, p3)
                            : joinRef(ref, p5, p5 + (p3 - p5 - ins.length()) / 2); // Tandem duplication
                    if (conf.y) {
                        System.err.printf("Found lgins30: %s %s %s %s + %s\n", p3, p5, ins.length(), tmp, ins);
                    }
                    ins = "+" + tmp + ins;
                    bi = p5 - 1;
                    vref = getVariation(iHash, bi, ins);
                }
                sc3v.used = true;
                sc5v.used = true;
                vref.pstd = true;
                vref.qstd = true;
                incCnt(cov, bi, sc5v.cnt);
                if (conf.y) {
                    System.err.printf(" lgins30 Found: '%s' %s %s %s\n", ins, bi, bp3, bp5);
                }

                if (ins.startsWith("+")) {
                    Variation mvref = getVariationMaybe(hash, bi, ref.get(bi));
                    adjCnt(vref, sc3v, mvref, conf);
                    adjCnt(vref, sc5v, conf);
                    if (bams != null && bams.length > 0
                            && p3 - p5 >= 5 && p3 - p5 > rlen - 10
                            && mvref != null && mvref.cnt != 0
                            && vref.cnt > 2 * mvref.cnt
                            && noPassingReads(chr, p5, p3, bams, conf)) {
                        adjCnt(vref, mvref, mvref, conf);
                    }
                    Map<Integer, Map<String, Integer>> tins = new HashMap<>();
                    Map<String, Integer> map = new HashMap<>();
                    map.put(ins, vref.cnt);
                    tins.put(bi, map);
                    realignins(hash, iHash, tins, cov, sclip5, sclip3, reference, region, chrs, rlen, conf);
                } else if (ins.startsWith("-")) {
                    adjCnt(vref, sc3v, getVariationMaybe(hash, bi, ref.get(bi)), conf);
                    adjCnt(vref, sc5v, conf);
                    Map<Integer, Map<String, Integer>> tdel = new HashMap<>();
                    Map<String, Integer> map = new HashMap<>();
                    map.put(ins, vref.cnt);
                    tdel.put(bi, map);
                    realigndel(hash, tdel, cov, sclip5, sclip3, reference, region, chrs, rlen, bam, conf);
                } else {
                    adjCnt(vref, sc3v, conf);
                    adjCnt(vref, sc5v, conf);
                }
                break;
            }
        }
        if (conf.y) {
            System.err.println("Done: lgins30\n");
        }
    }

    static List<Object[]> fillTmp(Map<Integer, Map<String, Integer>> changes) {
        List<Object[]> tmp = new ArrayList<>();
        for (Map.Entry<Integer, Map<String, Integer>> ent : changes.entrySet()) {
            int p = ent.getKey();
            Map<String, Integer> v = ent.getValue();
            for (Map.Entry<String, Integer> entV : v.entrySet()) {
                String vn = entV.getKey();
                int cnt = entV.getValue();
                // In perl it doesn't commented. but ecnt isn't used
                // int ecnt = 0;
                // Matcher mtch = ATGC_E.matcher(vn);
                // if (mtch.find()) {
                // ecnt = mtch.group(1).length();
                // }
                tmp.add(new Object[] { p, vn, cnt, /* ecnt */});
            }
        }
        Collections.sort(tmp, REALIGNDEL_COMPARATOR);
        return tmp;
    }

    private static final Comparator<Object[]> REALIGNDEL_COMPARATOR = new Comparator<Object[]>() {
        @Override
        public int compare(Object[] o1, Object[] o2) {
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

        }
    };

    /**
     * Test whether the two soft-clipped reads match
     * Returns the breakpoints in 5 and 3 prime soft-clipped reads
     * @param seq5 consensus sequence 5' strand
     * @param seq3 consensus sequence 3' strand
     * @param p5 position 5'
     * @param p3 position 3'
     * @param ref map of reference bases (not used)
     * @return Tuple of (start position of 5' strand, start position of 3' strand, max match length)
     */
    static Tuple.Tuple3<Integer, Integer, Integer> find35match(String seq5,
                                                               String seq3,
                                                               int p5,
                                                               int p3,
                                                               Map<Integer, Character> ref) {
        final int longmm = 3;
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
                    if (nm > longmm) {
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
     * Check whether there're reads supporting wild type in deletions
     * Only for indels that have micro-homology
     * @param chr chromosome name
     * @param s start position
     * @param e end position
     * @param bams BAM file list
     * @param conf configuration
     * @return true if any read was found in chr:s-e
     * @throws IOException
     */
    static boolean noPassingReads(String chr,
                                  int s,
                                  int e,
                                  String[] bams,
                                  Configuration conf) throws IOException {
        int cnt = 0;
        int midcnt = 0; // Reads end in the middle
        int dlen = e - s;
        String dlenqr = dlen + "D";
        Region region = new Region(chr, s, e, "");
        for (String bam : bams) {
            try (SamView reader = new SamView(bam, "", region, conf.validationStringency)) {
                SAMRecord record;
                while ((record = reader.read()) != null) {
                    if (record.getCigarString().contains(dlenqr)) {
                        continue;
                    }
                    int rs = record.getAlignmentStart();
                    // The total aligned length, excluding soft-clipped bases and insertions
                    int rlen = getAlignedLength(record.getCigar());
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
        if (conf.y) {
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
     * @return true if seq1 matches seq2 with no more than 3 mismatches
     */
    static boolean ismatch(String seq1,
                           String seq2,
                           int dir,
                           boolean debugLog) {
        int MM = 3;
        return ismatch(seq1, seq2, dir, debugLog, MM);
    }

    /**
     * $ismatch
     * Find if sequences match
     * @param seq1 first sequence
     * @param seq2 second sequence
     * @param dir direction of seq2 (1 or -1)
     * @param debugLog print debug message if true
     * @param MM length of mismatches for SV
     * @return true if seq1 matches seq2 with no more than MM number of mismatches
     */
    static boolean ismatch(String seq1,
                           String seq2,
                           int dir,
                           boolean debugLog,
                           int MM) {
        if (debugLog) {
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
    public static Tuple.Tuple3<Integer, String, Integer> adjInsPos(int bi,
                                                                   String ins,
                                                                   Map<Integer, Character> ref) {
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
        return tuple(bi, ins, bi);
    }

    /**
     * Find the insertion
     * @param seq sequence
     * @param p position
     * @param ref map of reference bases
     * @param dir direction
     * @param chr chromosome name
     * @param chrs map of chromosome lengths
     * @return Tuple of (BI (insert starting position), INS (insert sequence), BI2 ( = BI))
     */
    static Tuple.Tuple3<Integer, String, Integer> findbi(String seq,
                                                         int p,
                                                         Map<Integer, Character> ref,
                                                         final int dir,
                                                         String chr,
                                                         Map<String, Integer> chrs) {
        final int maxmm = 3; // maximum mismatches allowed
        final int dirExt = dir == -1 ? 1 : 0;
        int score = 0;
        int bi = 0;
        String ins = "";
        int bi2 = 0;

        for (int n = 6; n < seq.length(); n++) {
            if (p + 6 >= chrs.get(chr)) {
                break;
            }
            int mm = 0;
            int i = 0;
            Set<Character> m = new HashSet<>();
            for (i = 0; i + n < seq.length(); i++) {
                if (p + dir * i - dirExt < 1) {
                    break;
                }
                if (p + dir * i - dirExt > chrs.get(chr)) {
                    break;
                }
                if (isNotEquals(seq.charAt(i + n), ref.get(p + dir * i - dirExt))) {
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
                while (n + ept + 1 < seq.length() && (!isEquals(seq.charAt(n + ept), ref.get(p + ept * dir - dirExt))
                        || !isEquals(seq.charAt(n + ept + 1), ref.get(p + (ept + 1) * dir - dirExt)))) {
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
                        while (s >= -n && isEquals(charAt(insert, s), ref.get(p + s))) {
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
     * Find breakpoint
     * @param seq sequence
     * @param sp start position of sequence
     * @param ref map of reference bases
     * @param dis distance (always = $INDELSIZE)
     * @param dir direction
     * @param chr chromosome name
     * @param chrs map of chromosome lengths
     * @param debugLog print debug message if true
     * @return breakpoint position
     */
    static int findbp(String seq,
                      int sp,
                      Map<Integer, Character> ref,
                      int dis,
                      int dir,
                      String chr,
                      Map<String, Integer> chrs,
                      boolean debugLog) {

        final int maxmm = 3; // maximum mismatches allowed
        int bp = 0;
        int score = 0;
        int idx = chrs.containsKey(chr) ? chrs.get(chr) : 0;
        for (int n = 0; n < dis; n++) {
            int mm = 0;
            int i = 0;
            Set<Character> m = new HashSet<>();
            for (i = 0; i < seq.length(); i++) {
                if (sp + dir * n + dir * i < 1) {
                    break;
                }
                if (sp + dir * n + dir * i > idx) {
                    break;
                }
                if (isEquals(seq.charAt(i), ref.get(sp + dir * n + dir * i))) {
                    m.add(seq.charAt(i));
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
            if (mm <= maxmm - n / 100 && i >= seq.length() - 2 && i >= 8 + n / 10 && mm / (double)i < 0.12) {
                int lbp = sp + dir * n - (dir < 0 ? dir : 0);
                if (mm == 0 && i == seq.length()) {
                    if (debugLog) {
                        System.err.printf("  Findbp: %s %s %s %s %s\n", seq, sp, lbp, mm, i);
                    }
                    return lbp;
                } else if (i - mm > score) {
                    bp = lbp;
                    score = i - mm;
                }
            }
        }
        if (debugLog && bp != 0) {
            System.err.printf("  Findbp with mismatches: %s %s %s %s %s\n", seq, sp, bp, dir, score);
        }
        return bp;
    }

    static void rmCnt(Variation vref,
                      Variation tv) {
        vref.cnt -= tv.cnt;
        vref.hicnt -= tv.hicnt;
        vref.locnt -= tv.locnt;
        vref.pmean -= tv.pmean;
        vref.qmean -= tv.qmean;
        vref.Qmean -= tv.Qmean;
        vref.subDir(true, tv.getDir(true));
        vref.subDir(false, tv.getDir(false));
        correctCnt(vref);
    }

    /**
     * Adjust the reference count.
     * @param tv non-reference variant
     * @param ref reference variant
     * @param len length for adhustment factor
     */
    static void adjRefCnt(Variation tv,
                          Variation ref,
                          int len) {
        if (ref == null) {
            return;
        }
        double f = tv.pmean != 0 ? (tv.pmean / (double)tv.cnt - len + 1) / (tv.pmean / (double)tv.cnt) : 0; // the adjustment factor
        if (f < 0) {
            return;
        }

        if (f > 1) {
            f = 1;
        }

        ref.cnt -= (int) (f * tv.cnt);
        ref.hicnt -= (int) (f * tv.hicnt);
        ref.locnt -= (int) (f * tv.locnt);
        ref.pmean -= f * tv.pmean;
        ref.qmean -= f * tv.qmean;
        ref.Qmean -= f * tv.Qmean;
        ref.nm -= f * tv.nm;
        ref.subDir(true, (int)(f * tv.getDir(true)));
        ref.subDir(false, (int)(f * tv.getDir(false)));
        correctCnt(ref);
    }

    /**
     * Adjust the reference by factor
     * @param ref
     * @param factor_f
     * @param debugLog boolean parameter of 'y' option of Configuartion (print debug message if true)
     */
    static void adjRefFactor(Variation ref, double factor_f, boolean debugLog) {
        if (ref == null) {
            return;
        }

        if (factor_f > 1) {
            factor_f = 1;
        }

        if (factor_f < -1) {
            return;
        }
        if (debugLog) {
            System.err.printf("adjRefFactor(%s, %s);\n", ref.cnt, factor_f);
        }

        ref.cnt -= (int) (factor_f * ref.cnt);
        ref.hicnt -= (int) (factor_f * ref.hicnt);
        ref.locnt -= (int) (factor_f * ref.locnt);
        ref.pmean -= factor_f * ref.pmean;
        ref.qmean -= factor_f * ref.qmean;
        ref.Qmean -= factor_f * ref.Qmean;
        ref.nm -= factor_f * ref.nm;
        ref.dirPlus -= (int) (factor_f * ref.dirPlus);
        ref.dirMinus -= (int) (factor_f * ref.dirMinus);

        correctCnt(ref);
    }

    /**
     * Add variant by factor
     * @param vref $ref
     * @param factor_f
     */
    static void addVarFactor(Variation vref, double factor_f) {
        if (vref == null) {
            return;
        }

        if (factor_f < -1) {
            return;
        }

        vref.cnt += (int) (factor_f * vref.cnt);
        vref.hicnt += (int) (factor_f * vref.hicnt);
        vref.locnt += (int) (factor_f * vref.locnt);
        vref.pmean += factor_f * vref.pmean;
        vref.qmean += factor_f * vref.qmean;
        vref.Qmean += factor_f * vref.Qmean;
        vref.nm += factor_f * vref.nm;
        vref.dirPlus += (int) (factor_f * vref.dirPlus);
        vref.dirMinus += (int) (factor_f * vref.dirMinus);
    }


    /**
     * Given a variant sequence, find the mismatches and potential softclipping positions
     * @param ref map of reference bases
     * @param p position
     * @param wupseq sequence
     * @param len not used
     * @param sclip5 map of 5' softclips
     * @return MismatchResult contains mismatches lists and clipping positions
     */
    static MismatchResult findMM5(Map<Integer, Character> ref,
                                  int p,
                                  String wupseq,
                                  int len,
                                  Map<Integer, Sclip> sclip5) {
        String seq = wupseq.replaceAll("#|\\^", "");
        int longmm = 3;
        List<Tuple.Tuple3<String, Integer, Integer>> mm = new ArrayList<>(); // mismatches, mismatch positions, 5 or 3 ends
        int n = 0;
        int mn = 0;
        int mcnt = 0;
        StringBuilder str = new StringBuilder();
        List<Integer> sc5p = new ArrayList<>();
        while (isHasAndNotEquals(charAt(seq, -1 - n), ref, p - n) && mcnt < longmm) {
            str.insert(0, charAt(seq, -1 - n));
            mm.add(tuple(str.toString(), p - n, 5));
            n++;
            mcnt++;
        }
        sc5p.add(p + 1);
        // Adjust clipping position if only one mismatch
        int misp = 0;
        Character misnt = null;
        if (str.length() == 1) {
            while (isHasAndEquals(charAt(seq, -1 - n), ref, p - n)) {
                n++;
                mn++;
            }
            if (mn > 1) {
                int n2 = 0;
                while (-1 - n - 1 - n2 >= 0
                        && isHasAndEquals(charAt(seq, -1 - n - 1 - n2), ref, p - n - 1 - n2)) {
                    n2++;
                }
                if (n2 > 2) {
                    sc5p.add(p - n - n2);
                    misp = p - n;
                    misnt = charAt(seq, -1 - n);
                    if (sclip5.containsKey(p - n - n2)) {
                        sclip5.get(p - n - n2).used = true;
                    }
                    mn += n2;
                } else {
                    sc5p.add(p - n);
                    if (sclip5.containsKey(p - n)) {
                        sclip5.get(p - n).used = true;
                    }
                }

            }
        }
        return new MismatchResult(mm, sc5p, mn, misp, misnt == null ? "" : misnt.toString());
    }

    /**
     * Given a variant sequence, find the mismatches and potential softclipping positions
     * @param ref map of reference bases
     * @param p position
     * @param sanpseq sequence
     * @param len not used
     * @param sclip3 map of 3' softclips
     * @return MismatchResult contains mismatches lists and clipping positions
     */
    static MismatchResult findMM3(Map<Integer, Character> ref,
                                  int p,
                                  String sanpseq,
                                  int len,
                                  Map<Integer, Sclip> sclip3) {
        String seq = sanpseq.replaceAll("#|\\^", ""); // ~ s/#|\^//g;
        final int longmm = 3;
        // mismatches, mismatch positions, 5 or 3 ends
        List<Tuple.Tuple3<String, Integer, Integer>> mm = new ArrayList<>();
        int n = 0;
        int mn = 0;
        int mcnt = 0;
        List<Integer> sc3p = new ArrayList<>();
        StringBuilder str = new StringBuilder();
        while (n < seq.length() && isEquals(ref.get(p + n), seq.charAt(n))) {
            n++;
        }
        sc3p.add(p + n);
        int Tbp = p + n;
        while (mcnt <= longmm && n < seq.length() && isNotEquals(ref.get(p + n), seq.charAt(n))) {
            str.append(seq.charAt(n));
            mm.add(tuple(str.toString(), Tbp, 3));
            n++;
            mcnt++;
        }
        // Adjust clipping position if only one mismatch
        int misp = 0;
        Character misnt = null;
        if (str.length() == 1) {
            while (n < seq.length() && isHasAndEquals(seq.charAt(n), ref, p + n)) {
                n++;
                mn++;
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
                    if (sclip3.containsKey(p + n + n2)) {
                        sclip3.get(p + n + n2).used = true;
                    }
                    mn += n2;
                } else {
                    sc3p.add(p + n);
                    if (sclip3.containsKey(p + n)) {
                        sclip3.get(p + n).used = true;
                    }
                }
            }
        }
        return new MismatchResult(mm, sc3p, mn, misp, misnt == null ? "" : misnt.toString());
    }

    static class MismatchResult {
        private final List<Tuple.Tuple3<String, Integer, Integer>> mm;
        private final List<Integer> scp;
        private final int nm;
        private final int misp;
        private final String misnt;

        public MismatchResult(List<Tuple.Tuple3<String, Integer, Integer>> mm, List<Integer> scp, int nm, int misp, String misnt) {
            this.mm = mm;
            this.scp = scp;
            this.nm = nm;
            this.misp = misp;
            this.misnt = misnt;
        }

        public List<Integer> getScp() {
            return scp;
        }

        public List<Tuple.Tuple3<String, Integer, Integer>> getMm() {
            return mm;
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

}
