package com.astrazeneca.vardict.variations;

import com.astrazeneca.utils.Tuple;
import com.astrazeneca.vardict.Configuration;

import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Set;

import static com.astrazeneca.GlobalReadOnlyScope.instance;
import static com.astrazeneca.utils.Tuple.tuple;
import static com.astrazeneca.utils.Utils.substr;
import static com.astrazeneca.utils.Utils.isEquals;

public class VariationUtils {
    /**
     * VARN - array of variants, REF - reference variant, VAR - map of (variant description string =&gt; variant)
     */
    public enum VarsType {
        varn, ref, var,
    }

    /**
     * Increase count for given key
     * @param counts map of counts
     * @param key key to add count
     * @param add amount to add
     */
    //perl version: 2447
    @SuppressWarnings({ "unchecked", "rawtypes" })
    public static void incCnt(Map counts, Object key, int add) {
        Integer integer = (Integer)counts.get(key);
        if (integer == null) {
            counts.put(key, add);
        } else {
            counts.put(key, integer + add);
        }
    }

    /**
     * Calculate strand bias flag
     * @param forwardCount Variant count for forward strand
     * @param reverseCount Variant count for reverse strand
     * @param bias
     * @param minBiasReads minimum reads for bias
     * @return 0 - small total count, only one of strands, 1 - strand bias, 2 - no strand bias
     */
    //perl version: 2775
    public static int strandBias(int forwardCount, int reverseCount, double bias, int minBiasReads) {

        if (forwardCount + reverseCount <= 12) { // using p=0.01, because prop.test(1,12) = 0.01
            return forwardCount * reverseCount > 0 ? 2 : 0;
        }

        return (forwardCount / (double)(forwardCount + reverseCount) >= bias && reverseCount / (double)(forwardCount + reverseCount) >= bias && forwardCount >= minBiasReads && reverseCount >= minBiasReads) ? 2 : 1;
    }

    /**
     * Find the consensus sequence in soft-clipped reads. Consensus is called if
     * the matched nucleotides are &gt;90% of all softly clipped nucleotides.
     * @param softClip soft-clipped sequences
     * @return consensus sequence
     */
    //perl version: 2169
    public static String findconseq(SoftClip softClip) {
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
                if (softClip.seq.containsKey(i) && softClip.seq.get(i).containsKey(nt) && softClip.seq.get(i).get(nt).meanQuality > maxq) {
                    max = ncnt;
                    mnt = nt;
                    maxq = softClip.seq.get(i).get(nt).meanQuality;
                }
            }
            if ((tt - max > 2 || max <= tt - max) && max / (double)tt < 0.8) {
                if (flag)
                    break;
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
        if (instance().conf.y) {
            System.err.printf("  candidate consensus: %s M: %s T: %s Final: %s\n", seq, match, total, softClip.sequence);
        }
        return softClip.sequence;

    }

    /**
     * Construct reference sequence
     * @param baseToPosition map of reference bases
     * @param from start position
     * @param to end position
     * @return reference sequence starting at from and ending at to
     */
    //perl version: 2674
    public static String joinRef(Map<Integer, Character> baseToPosition, int from, int to) {
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
     * @param varToAdd variant to add counts
     * @param variant variant
     */
    public static void adjCnt(Variation varToAdd, Variation variant) {
        adjCnt(varToAdd, variant, null);
    }

    /**
     * Adjust the count,  If <code>referenceVar</code> is not null, the count for reference is also adjusted.
     * Variant counts of <code>variant</code> are added to <code>varToAdd</code> and removed from <code>referenceVar</code>
     * @param varToAdd variant to add counts
     * @param variant variant
     * @param referenceVar reference variant
     */
    //perl version: 2519
    public static void adjCnt(Variation varToAdd, Variation variant, Variation referenceVar) {
        varToAdd.varsCount += variant.varsCount;
        varToAdd.extraCount += variant.varsCount;
        varToAdd.highQualityReadsCount += variant.highQualityReadsCount;
        varToAdd.lowQualityReadsCount += variant.lowQualityReadsCount;
        varToAdd.meanPosition += variant.meanPosition;
        varToAdd.meanQuality += variant.meanQuality;
        varToAdd.meanMappingQuality += variant.meanMappingQuality;
        varToAdd.numberOfMismatches += variant.numberOfMismatches;
        varToAdd.isAtLeastAt2Positions = true;
        varToAdd.hasAtLeast2DiffQualities = true;
        varToAdd.addCountForDir(true, variant.getCountForDir(true));
        varToAdd.addCountForDir(false, variant.getCountForDir(false));

        if (instance().conf.y) {
            String refCnt = referenceVar != null ? String.valueOf(referenceVar.varsCount) : "NA";
            System.err.printf("    AdjCnt: '+' %s %s %s %s Ref: %s\n", varToAdd.varsCount, variant.varsCount, varToAdd.getCountForDir(false), variant.getCountForDir(false), refCnt);
            System.err.printf("    AdjCnt: '-' %s %s %s %s Ref: %s\n", varToAdd.varsCount, variant.varsCount, varToAdd.getCountForDir(true), variant.getCountForDir(true), refCnt);
        }

        if (referenceVar == null)
            return;

        referenceVar.numberOfMismatches -= variant.numberOfMismatches;
        referenceVar.adjVariant(variant, 1);
    }

    public static Variant getVarMaybe(Map<Integer, Vars> alignedVariants, int key, VarsType type, Object... keys) {
        if (!alignedVariants.containsKey(key)) {
            return null;
        }
        Vars vars = alignedVariants.get(key);
        return getVarMaybe(vars, type, keys);
    }

    public static Variant getVarMaybe(Vars vars, VarsType type, Object... keys) {
        switch (type) {
            case var:
                if (vars.variants.size() > (Integer)keys[0]) {
                    return vars.variants.get((Integer)keys[0]);
                }
            case varn:
                return vars.varDescriptionStringToVariants.get(keys[0]);
            case ref:
                return vars.referenceVariant;
        }
        return null;
    }

    /**
     * Returns <tt>true</tt> whether a variant meet specified criteria
     * @param vref variant
     * @param referenceVar reference variant
     * @param type Type of variant
     * @param splice set of strings representing introns in splice
     * @return <tt>true</tt> if variant meet specified criteria
     */
    //perl version: 523
    public static boolean isGoodVar(Variant vref, Variant referenceVar, Variant.Type type,
                                    Set<String> splice) {
        Configuration config = instance().conf;
        if (vref == null || vref.refAllele == null || vref.refAllele.isEmpty()) {
            return false;
        }

        if (type == null || type == Variant.Type.noInfo) {
            type = vref.getType();
        }
        if (vref.frequency < config.freq
                || vref.highQualityReadsCount < config.minr
                || vref.meanPosition < config.readPosFilter
                || vref.meanQuality < config.goodq) {
            return false;
        }

        if (referenceVar != null && referenceVar.highQualityReadsCount > config.minr && vref.frequency < 0.25d) {
            //The referenceVar allele has much better mean meanMappingQuality than variants allele, thus likely false positives
            double d = vref.meanMappingQuality + vref.refAllele.length() + vref.varAllele.length();
            if ((d - 2 < 5 && referenceVar.meanMappingQuality > 20)
                    || (1 + d) / (referenceVar.meanMappingQuality + 1) < 0.25d) {
                return false;
            }
        }

        if (type == Variant.Type.deletion && splice.contains(vref.startPosition + "-" + vref.endPosition)) {
            return false;
        }
        if (vref.highQualityToLowQualityRatio < config.qratio) {
            return false;
        }
        if (vref.frequency > 0.35d) {
            return true;
        }
        if (vref.meanMappingQuality < config.mapq) {
            return false;
        }
        if (vref.msi >= 13 && vref.frequency <= 0.275d && vref.msint == 1) {
            return false;
        }
        if (vref.msi >= 8 && vref.frequency <= 0.2d && vref.msint > 1) {
            return false;
        }
        if (vref.strandBiasFlag.equals("2;1") && vref.frequency < 0.25d) {
            if (type == Variant.Type.SNV || (vref.refAllele.length() <= 3 && vref.varAllele.length() <= 3)) {
                return false;
            }
        }

        return true;
    }

    public static Variation getVariation(Map<Integer, Map<String, Variation>> hash, int start, String key) {
        Map<String, Variation> map = hash.get(start);
        if (map == null) {
            map = new LinkedHashMap<>();
            hash.put(start, map);
        }
        Variation variation = map.get(key);
        if (variation == null) {
            variation = new Variation();
            map.put(key, variation);
        }
        return variation;
    }

    //perl version: 1039
    public static Variation getVariationMaybe(Map<Integer, Map<String, Variation>> hash, int start, Character key) {
        if (key == null)
            return null;

        Map<String, Variation> map = hash.get(start);
        if (map == null) {
            return null;
        }

        return map.get(key.toString());
    }

    /**
     * Adjust the reference count.
     * @param variant non-reference variant
     * @param referenceVar reference variant
     * @param length length for adjustment factor
     */
    //perl version: 2551
    public static void adjRefCnt(Variation variant, Variation referenceVar, int length) {
        if (referenceVar == null) {
            return;
        }
        double adjustmentFactor = variant.meanPosition != 0
                ? (variant.meanPosition / (double)variant.varsCount - length + 1) / (variant.meanPosition / (double)variant.varsCount)
                : 0;

        if (adjustmentFactor < 0) {
            return;
        }

        if (adjustmentFactor > 1) {
            adjustmentFactor = 1;
        }

        referenceVar.numberOfMismatches -= adjustmentFactor * variant.numberOfMismatches;
        referenceVar.adjVariant(variant, adjustmentFactor);
    }

    /**
     * Adjust the insertion position if necessary
     * @param startPos starting position of insert
     * @param insertSeq insert sequence
     * @param baseToPosition map of reference bases
     * @return Tuple of (startPos, insertSeq, startPos)
     */
    //perl version: 2284
    public static Tuple.Tuple3<Integer, String, Integer> adjInsPos(int startPos, String insertSeq, Map<Integer, Character> baseToPosition) {
        int n = 1;
        int len = insertSeq.length();
        while (isEquals(baseToPosition.get(startPos), insertSeq.charAt(insertSeq.length() - n))) {
            n++;
            if (n > len) {
                n = 1;
            }
            startPos--;
        }
        if (n > 1) {
            insertSeq = substr(insertSeq, 1 - n) + substr(insertSeq, 0, 1 - n);
        }
        return tuple(startPos, insertSeq, startPos);
    }

    //perl version: 1040
    public static boolean isHasAndEquals(Character ch1, Map<Integer, Character> ref, int index) {
        Character refc = ref.get(index);
        if (refc == null)
            return false;
        return refc.equals(ch1);
    }

    public static boolean isHasAndEquals(Map<Integer, Character> ref, int index1, String str, int index2) {
        Character refc = ref.get(index1);
        if (refc == null)
            return false;
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
        return !refc.equals(str.charAt(index2));
    }

    public static boolean isNotEquals(Character ch1, Character ch2) {
        return !isEquals(ch1, ch2);
    }
}
