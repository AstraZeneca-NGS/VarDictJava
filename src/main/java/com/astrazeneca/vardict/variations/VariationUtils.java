package com.astrazeneca.vardict.variations;

import com.astrazeneca.vardict.Configuration;
import com.astrazeneca.vardict.collection.VariationMap;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;

import static com.astrazeneca.vardict.data.Patterns.B_A7;
import static com.astrazeneca.vardict.data.Patterns.B_T7;
import static com.astrazeneca.vardict.modules.VariationRealigner.islowcomplexseq;

public class VariationUtils {
    /**
     * VARN - array of variants, REF - reference variant, VAR - map of (variant description string =&gt; variant)
     */
    public enum VarsType {
        varn, ref, var,
    }

    /**
     * Increase count for given key
     * @param counts map of counts $cnts
     * @param key key to add count
     * @param add amount to add
     */
    @SuppressWarnings({ "unchecked", "rawtypes" })
    public static void incCnt(Map counts,
                              Object key,
                              int add) {
        Integer integer = (Integer)counts.get(key);
        if (integer == null) {
            counts.put(key, add);
        } else {
            counts.put(key, integer + add);
        }
    }

    /**
     * Calculate strand bias flag
     * @param forwardCount Variant count for forward strand  $fwd
     * @param reverseCount Variant count for reverse strand  $rev
     * @param bias the cutoff from configuration
     * @param minBiasReads minimum reads for bias   $minb
     * @return 0 - small total count, only one of strands, 1 - strand bias, 2 - no strand bias
     */
    public static int strandBias(int forwardCount,
                                 int reverseCount,
                                 double bias,
                                 int minBiasReads) {
        // using p=0.01, because prop.test(1,12) = 0.01
        if (forwardCount + reverseCount <= 12) {
            return forwardCount * reverseCount > 0 ? 2 : 0;
        }

        return (forwardCount / (double)(forwardCount + reverseCount) >= bias
                && reverseCount / (double)(forwardCount + reverseCount) >= bias
                && forwardCount >= minBiasReads
                && reverseCount >= minBiasReads) ? 2 : 1;
    }

    /**
     * Find the consensus sequence in soft-clipped reads. Consensus is called if
     * the matched nucleotides are &gt;90% of all softly clipped nucleotides.
     * @param softClip soft-clipped sequences    $scv
     * @param conf Configuration for debug and adseed options
     * @return consensus sequence
     */
    public static String findconseq(Sclip softClip,
                                    Configuration conf) {
        if (softClip.sequence != null) {
            return softClip.sequence;
        }

        int total = 0;
        int match = 0;
        StringBuilder seq = new StringBuilder();
        boolean flag = false;
        for (Map.Entry<Integer, Map<Character, Integer>> nve : softClip.nt.entrySet()) {
            Integer i = nve.getKey();
            Map<Character, Integer> nv = nve.getValue();
            int max = 0;
            double maxq = 0;
            Character mnt = null;
            int tt = 0;
            for (Map.Entry<Character, Integer> ent : nv.entrySet()) {
                Character nt = ent.getKey();
                int ncnt = ent.getValue();
                tt += ncnt;
                if (ncnt > max && softClip.seq.containsKey(i) && softClip.seq.get(i).containsKey(nt) && softClip.seq.get(i).get(nt).qmean > maxq) {
                    max = ncnt;
                    mnt = nt;
                    maxq = softClip.seq.get(i).get(nt).qmean;
                }
            }
            if (i == 3 && softClip.nt.size() >= 6 && tt/(double)softClip.cnt < 0.2 && tt <= 2) {
                break;
            }
            if ((tt - max > 2 || max <= tt - max) && max / (double)tt < 0.8) {
                if (flag) {
                    break;
                }
                flag = true;
            }
            total += tt;
            match += max;
            if (mnt != null) {
                seq.append(mnt);
            }
        }

        Integer ntSize = softClip.nt.lastKey();
        if (total != 0
                && match / (double)total > 0.9
                && seq.length() / 1.5 > ntSize - seq.length()
                && (seq.length() / (double)ntSize > 0.8
                || ntSize - seq.length() < 10
                || seq.length() > 25)) {
            softClip.sequence = seq.toString();
        } else {
            softClip.sequence = "";
        }

        if (!softClip.sequence.isEmpty() && seq.length() > Configuration.SEED_2) {
            Matcher mm1  =  B_A7.matcher(seq);
            Matcher mm2  =  B_T7.matcher(seq);
            if (mm1.find() || mm2.find()) {
                softClip.used = true;
            }
            if (islowcomplexseq(seq.toString())) {
                softClip.used = true;
            }
        }
        if (conf.y) {
            System.err.printf("  Candidate consensus: %s Reads: %s M: %s T: %s Final: %s\n", seq, softClip.cnt, match, total, softClip.sequence);
        }
        return softClip.sequence;
    }

    // TODO: Need to implement in adaptor task
//    static String findconseq(Sclip scv, Configuration conf, int dir) {
//        String sequence = findconseq(scv, conf);
//
//        if (!sequence.isEmpty() && sequence.length() >= conf.adseed) {
//            if ( dir == 3 ) { // 3'
//                if (adaptor.containsKey(substr(sequence, 0, conf.adseed))) {
//                    scv.sequence = "";
//                }
//            } else if ( dir == 5 ) { // 5'
//                if (adaptor_rev.containsKey(new StringBuilder(substr(sequence, 0, conf.adseed)).reverse().toString())) {
//                    scv.sequence = "";
//                }
//            }
//        }
//        return scv.sequence;
//    }

    /**
     * Construct reference sequence
     * @param baseToPosition map of reference bases    $ref
     * @param from start position
     * @param to end position
     * @return reference sequence starting at from and ending at to
     */
    public static String joinRef(Map<Integer, Character> baseToPosition,
                          int from,
                          int to) {
        StringBuilder sb = new StringBuilder();
        for (int i = from; i <= to; i++) {
            Character ch = baseToPosition.get(i);
            if (ch != null) {
                sb.append(ch);
            }
        }
        return sb.toString();
    }
    /**
     * Adjust count. Don't adjust reference
     * @param varToAdd variant to add counts    $vref
     * @param variation variant     $tv
     * @param conf configuration
     */
    public static void adjCnt(Variation varToAdd,
                              Variation variation,
                              Configuration conf) {
        adjCnt(varToAdd, variation, null, conf);
    }

    /**
     * Adjust the count,  If ref is not null, the count for reference is also adjusted.
     * Variant counts of tv are added to vref and removed from ref
     * @param varToAdd variant to add counts    $vref
     * @param variant variant    $tv
     * @param referenceVar reference variant  $ref
     * @param conf configuration
     */
    public static void adjCnt(Variation varToAdd,
                              Variation variant,
                              Variation referenceVar, Configuration conf) {
        varToAdd.cnt += variant.cnt;
        varToAdd.extracnt += variant.cnt;
        varToAdd.hicnt += variant.hicnt;
        varToAdd.locnt += variant.locnt;
        varToAdd.pmean += variant.pmean;
        varToAdd.qmean += variant.qmean;
        varToAdd.Qmean += variant.Qmean;
        varToAdd.nm += variant.nm;
        varToAdd.pstd = true;
        varToAdd.qstd = true;
        varToAdd.addDir(true, variant.getDir(true));
        varToAdd.addDir(false, variant.getDir(false));

        if (conf.y) {
            String refCnt = referenceVar != null ? String.valueOf(referenceVar.cnt) : "NA";
            System.err.printf("    AdjCnt: '+' %s %s %s %s Ref: %s\n",
                    varToAdd.cnt, variant.cnt, varToAdd.getDir(false), variant.getDir(false), refCnt);
            System.err.printf("    AdjCnt: '-' %s %s %s %s Ref: %s\n",
                    varToAdd.cnt, variant.cnt, varToAdd.getDir(true), variant.getDir(true), refCnt);
        }

        if (referenceVar == null)
            return;

        referenceVar.cnt -= variant.cnt;
        referenceVar.hicnt -= variant.hicnt;
        referenceVar.locnt -= variant.locnt;
        referenceVar.pmean -= variant.pmean;
        referenceVar.qmean -= variant.qmean;
        referenceVar.Qmean -= variant.Qmean;
        referenceVar.nm -= variant.nm;
        referenceVar.subDir(true, variant.getDir(true));
        referenceVar.subDir(false, variant.getDir(false));
        correctCnt(referenceVar);
    }

    /**
     * Correct counts for negative values
     * @param varToCorrect variant   $ref
     */
    public static void correctCnt(Variation varToCorrect) {
        if (varToCorrect.cnt < 0)
            varToCorrect.cnt = 0;
        if (varToCorrect.hicnt < 0)
            varToCorrect.hicnt = 0;
        if (varToCorrect.locnt < 0)
            varToCorrect.locnt = 0;
        if (varToCorrect.pmean < 0)
            varToCorrect.pmean = 0;
        if (varToCorrect.qmean < 0)
            varToCorrect.qmean = 0;
        if (varToCorrect.Qmean < 0)
            varToCorrect.Qmean = 0;
        if (varToCorrect.getDir(true) < 0)
            varToCorrect.addDir(true, -varToCorrect.getDir(true));
        if (varToCorrect.getDir(false) < 0)
            varToCorrect.addDir(false, -varToCorrect.getDir(false));
    }

    /**
     * Returns <tt>true</tt> whether a variant meet specified criteria
     * @param vref variant
     * @param referenceVar reference variant    $rref
     * @param type Type of variant
     * @param splice set of strings representing introns in splice
     * @param conf Configuration (contains preferences for min, freq, filter and etc)
     * @return <tt>true</tt> if variant meet specified criteria
     */
    public static boolean isGoodVar(Variant vref, Variant referenceVar, String type,
                                    Set<String> splice,
                                    Configuration conf) {
        if (vref == null || vref.refallele == null || vref.refallele.isEmpty())
            return false;

        if (type == null || type.isEmpty()) {
            type = vref.varType();
        }
        if (vref.freq < conf.freq
                || vref.hicnt < conf.minr
                || vref.pmean < conf.readPosFilter
                || vref.qual < conf.goodq) {
            return false;
        }

        if (referenceVar != null && referenceVar.hicnt > conf.minr && vref.freq < 0.25d) {
            //The reference allele has much better mean mapq than var allele, thus likely false positives
            double d = vref.mapq + vref.refallele.length() + vref.varallele.length();
            if ((d - 2 < 5 && referenceVar.mapq > 20)
                    || (1 + d) / (referenceVar.mapq + 1) < 0.25d) {
                return false;
            }
        }

        if (type.equals("Deletion") && splice.contains(vref.sp + "-" + vref.ep)) {
            return false;
        }
        if (vref.qratio < conf.qratio) {
            return false;
        }
        if (vref.freq > 0.30d) {
            return true;
        }
        if (vref.mapq < conf.mapq) {
            return false;
        }
        if (vref.msi >= 15 && vref.freq <= 0.25d && vref.msint == 1) {
            return false;
        }
        if (vref.msi >= 12 && vref.freq <= 0.1d && vref.msint > 1) {
            return false;
        }
        if (vref.bias.equals("2;1") && vref.freq < 0.20d) {
            if (type == null || type.equals("SNV") || (vref.refallele.length() < 3 && vref.varallele.length() < 3)) {
                return false;
            }
        }
        return true;
    }


    /**
     * Returns var from map of vars or null if map doesn't contain such key.
     * @param alignedVariants alignedVariants $vars
     */
    public static Variant getVarMaybe(Map<Integer, Vars> alignedVariants,
                               int key,
                               VarsType type,
                               Object... keys) {
        if (!alignedVariants.containsKey(key)) {
            return null;
        }
        Vars vars = alignedVariants.get(key);
        return getVarMaybe(vars, type, keys);
    }

    public static Variant getVarMaybe(Vars vars,
                               VarsType type,
                               Object... keys) {
        switch (type) {
            case var:
                if (vars.var.size() > (Integer)keys[0]) {
                    return vars.var.get((Integer)keys[0]);
                }
            case varn:
                return vars.varn.get(keys[0]);
            case ref:
                return vars.ref;
        }
        return null;
    }

    public static Vars getOrPutVars(Map<Integer, Vars> map,
                                     int position) {
        Vars vars = map.get(position);
        if (vars == null) {
            vars = new Vars();
            map.put(position, vars);
        }
        return vars;
    }

    /**
     * Get {@link Variation} from {@link Sclip#seq} field
     * @param softClip variation with soft clipped reads    $sclip
     * @param idx index to get variation
     * @param ch key to get variation from {@link Sclip#seq}
     * @return variation
     */
    public static Variation getVariationFromSeq(Sclip softClip,
                                                int idx,
                                                Character ch) {
        Map<Character, Variation> map = softClip.seq.get(idx);
        if (map == null) {
            map = new HashMap<>();
            softClip.seq.put(idx, map);
        }
        Variation variation = map.get(ch);
        if (variation == null) {
            variation = new Variation();
            map.put(ch, variation);
        }
        return variation;
    }

    public static Variation getVariation(Map<Integer, VariationMap<String, Variation>> hash,
                                         int start,
                                         String ref) {
        VariationMap<String, Variation> map = hash.get(start);
        if (map == null) {
            map = new VariationMap<>();
            hash.put(start, map);
        }
        Variation variation = map.get(ref);
        if (variation == null) {
            variation = new Variation();
            map.put(ref, variation);
        }
        return variation;
    }

    public static Variation getVariationMaybe(Map<Integer, VariationMap<String, Variation>> hash,
                                              int start,
                                              Character ref) {
        if (ref == null)
            return null;

        Map<String, Variation> map = hash.get(start);
        if (map == null) {
            return null;
        }

        return map.get(ref.toString());
    }

    public static boolean isHasAndEquals(char ch1, Map<Integer, Character> ref, int index) {
        Character refc = ref.get(index);
        if (refc == null)
            return false;
        return refc.equals(ch1);
    }

    public static boolean isHasAndEquals(int index, Map<Integer, Character> ref, int index2) {
        Character refc = ref.get(index);
        if (refc == null)
            return false;
        Character refc2 = ref.get(index2);
        if (refc2 == null)
            return false;
        return refc.equals(refc2);
    }

    public static boolean isHasAndEquals(Map<Integer, Character> ref, int index1, String str, int index2) {
        Character refc = ref.get(index1);
        if (refc == null)
            return false;
        //if (index2 < 0) index2 = index2 + str.length();
        return refc.equals(str.charAt(index2));
    }

    public static boolean isHasAndNotEquals(Character ch1, Map<Integer, Character> ref, int index) {
        Character refc = ref.get(index);
        if (refc == null)
            return false;
        return !refc.equals(ch1);
    }

    public static boolean isHasAndNotEquals(Map<Integer, Character> ref, int index1, String str, int index2) {
        Character refc = ref.get(index1);
        if (refc == null)
            return false;
        //if (index2 < 0) index2 = index2 + str.length();
        return !refc.equals(str.charAt(index2));
    }

    public static boolean isEquals(Character ch1, Character ch2) {
        if (ch1 == null && ch2 == null)
            return true;
        if (ch1 == null || ch2 == null)
            return false;
        return ch1.equals(ch2);
    }

    public static boolean isNotEquals(Character ch1, Character ch2) {
        return !isEquals(ch1, ch2);
    }

}
