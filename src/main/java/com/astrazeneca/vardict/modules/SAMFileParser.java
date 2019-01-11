package com.astrazeneca.vardict.modules;

import com.astrazeneca.vardict.collection.VariationMap;
import com.astrazeneca.vardict.data.Region;
import com.astrazeneca.vardict.data.scopedata.InitialData;
import com.astrazeneca.vardict.data.scopedata.Scope;
import com.astrazeneca.vardict.variations.Sclip;
import com.astrazeneca.vardict.variations.Variation;

import java.util.*;

import static com.astrazeneca.vardict.data.scopedata.GlobalReadOnlyScope.instance;
import static com.astrazeneca.vardict.Utils.*;
import static com.astrazeneca.vardict.variations.VariationUtils.*;
import static java.lang.String.format;

public class SAMFileParser implements Module<InitialData, RecordPreprocessor> {
    @Override
    public Scope<RecordPreprocessor> process(Scope<InitialData> scope) {
        return new Scope<>(
                scope,
                new RecordPreprocessor(scope.bam.split(":"), scope.region, scope.data)
        );
    }

    /**
     * Adjust MNP when there're breakpoints within MNP (multi-nucleotide polymorphism)
     * @param hash map of variants produced by parseSAM method
     * @param mnp map of MNP variants
     * @param cov coverage
     * @param ref map of reference bases
     * @param sclip5 map of 5' softclips
     * @param sclip3 map of 3' softclips
     */
    public static void adjustMNP(Map<Integer, VariationMap<String, Variation>> hash,
                                 Map<Integer, Map<String, Integer>> mnp,
                                 Map<Integer, Integer> cov, Map<Integer, Character> ref,
                                 Map<Integer, Sclip> sclip3, Map<Integer, Sclip> sclip5,
                                 Region region) {

        for (Map.Entry<Integer, Map<String, Integer>> entry : mnp.entrySet()) {
            int lastPosition = 0;

            try {
                final Integer p = entry.getKey();
                lastPosition = p;
                Map<String, Integer> v = entry.getValue();

                for (Map.Entry<String, Integer> en : v.entrySet()) {
                    final String vn = en.getKey();
                    final Map<String, Variation> hashP = hash.get(p);
                    if (hashP == null) {
                        continue;
                    }
                    final Variation vref = hashP.get(vn);
                    if (vref == null) { // The variant is likely already been used by indel realignment
                        continue;
                    }
                    if (instance().conf.y) {
                        System.err.printf("  AdjMnt: %d %s %d\n", p, vn, vref.varsCount);
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
                                if (tref.varsCount <= 0) {
                                    continue;
                                }
                                if (tref.varsCount < vref.varsCount && tref.meanPosition / tref.varsCount <= i + 1) {
                                    if (instance().conf.y) {
                                        System.err.printf("    AdjMnt Left: %s %s Left: %s Cnt: %s\n", p, vn, left, tref.varsCount);
                                    }
                                    adjCnt(vref, tref);
                                    hashP.remove(left);
                                }
                            }
                        }
                        if (hash.containsKey(p + i + 1)) {
                            Variation tref = hash.get(p + i + 1).get(right);
                            if (tref != null) {
                                if (tref.varsCount < 0) {
                                    continue;
                                }
                                // #&& tref.pmean / tref.cnt <= mnt.length() - i - 1)
                                if (tref.varsCount < vref.varsCount) {
                                    if (instance().conf.y) {
                                        System.err.printf("    AdjMnt Right: %s %s Right: %s Cnt: %s\n", p, vn, right, tref.varsCount);
                                    }
                                    adjCnt(vref, tref);
                                    incCnt(cov, p, tref.varsCount);
                                    hash.get(p + i + 1).remove(right);
                                }
                            }
                        }
                    }
                    if (sclip3.containsKey(p)) {
                        final Sclip sc3v = sclip3.get(p);
                        if (!sc3v.used) {
                            final String seq = findconseq(sc3v, 0);
                            if (seq.startsWith(mnt)) {
                                if (seq.length() == mnt.length()
                                        || ismatchref(seq.substring(mnt.length()), ref, p + mnt.length(), 1)) {
                                    adjCnt(hash.get(p).get(vn), sc3v);
                                    incCnt(cov, p, sc3v.varsCount);
                                    sc3v.used = true;
                                }
                            }
                        }
                    }
                    if (sclip5.containsKey(p + mnt.length())) {
                        final Sclip sc5v = sclip5.get(p + mnt.length());
                        if (!sc5v.used) {
                            String seq = findconseq(sc5v, 0);
                            if (!seq.isEmpty() && seq.length() >= mnt.length()) {
                                seq = new StringBuffer(seq).reverse().toString();
                                if (seq.endsWith(mnt)) {
                                    if (seq.length() == mnt.length()
                                            || ismatchref(seq.substring(0, seq.length() - mnt.length()), ref, p - 1, -1)) {
                                        adjCnt(hash.get(p).get(vn), sc5v);
                                        incCnt(cov, p, sc5v.varsCount);
                                        sc5v.used = true;
                                    }
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

    /**
     * Utility method for adjustMNP method with default number of mismatches = 3
     * @param sequence subsequence consensus sequence in soft-clipped reads  $seq
     * @param ref map of integer - characters (nucleotides) in reference sequence
     * @param dir direction (forward or reverse)
     * @param position key for MNP map     $p
     * @return true if sequence is matched to reference
     */
    public static boolean ismatchref(String sequence, Map<Integer, Character> ref, int position, int dir) {
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
}
