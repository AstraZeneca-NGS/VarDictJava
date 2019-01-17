package com.astrazeneca.vardict.variations;

import com.astrazeneca.vardict.Configuration;
import com.astrazeneca.vardict.collection.VariationMap;

import java.util.*;
import java.util.regex.Matcher;

import static com.astrazeneca.vardict.Utils.reverse;
import static com.astrazeneca.vardict.Utils.substr;
import static com.astrazeneca.vardict.data.scopedata.GlobalReadOnlyScope.instance;
import static com.astrazeneca.vardict.Utils.charAt;
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
     * @return 0 - small total count, only one of strands, 1 - strand bias, 2 - no strand bias
     */
    public static int strandBias(int forwardCount,
                                 int reverseCount) {
        // using p=0.01, because prop.test(1,12) = 0.01
        if (forwardCount + reverseCount <= 12) {
            return forwardCount * reverseCount > 0 ? 2 : 0;
        }

        return (forwardCount / (double)(forwardCount + reverseCount) >= instance().conf.bias
                && reverseCount / (double)(forwardCount + reverseCount) >= instance().conf.bias
                && forwardCount >= instance().conf.minBiasReads
                && reverseCount >= instance().conf.minBiasReads) ? 2 : 1;
    }

    /**
     * Find the consensus sequence in soft-clipped reads. Consensus is called if
     * the matched nucleotides are &gt;90% of all softly clipped nucleotides.
     * @param softClip soft-clipped sequences    $scv
     * @param dir not used now, will be used when adaptor option will be added
     * @return consensus sequence
     */
    public static String findconseq(Sclip softClip,
                                    int dir) {
        if (softClip.sequence != null) {
            return softClip.sequence;
        }

        int total = 0;
        int match = 0;
        StringBuilder seq = new StringBuilder();
        boolean flag = false;
        for (Map.Entry<Integer, TreeMap<Character, Integer>> nve : softClip.nt.entrySet()) {
            Integer positionInSclip = nve.getKey();
            Map<Character, Integer> nv = nve.getValue();
            int maxCount = 0; //$max
            double maxQuality = 0; //$maxq
            Character chosenBase = null; //$mnt
            int totalCount = 0; //$tt
            for (Map.Entry<Character, Integer> ent : nv.entrySet()) {
                Character currentBase = ent.getKey(); //$nt
                int currentCount = ent.getValue(); //$ncnt
                totalCount += currentCount;
                if (currentCount > maxCount || (softClip.seq.containsKey(positionInSclip) && softClip.seq.get(positionInSclip).containsKey(currentBase)
                        && softClip.seq.get(positionInSclip).get(currentBase).meanQuality > maxQuality)) {
                    maxCount = currentCount;
                    chosenBase = currentBase;
                    maxQuality = softClip.seq.get(positionInSclip).get(currentBase).meanQuality;
                }
            }
            if (positionInSclip == 3 && softClip.nt.size() >= 6 && totalCount/(double)softClip.varsCount < 0.2 && totalCount <= 2) {
                break;
            }
            if ((totalCount - maxCount > 2 || maxCount <= totalCount - maxCount) && maxCount / (double)totalCount < 0.8) {
                if (flag) {
                    break;
                }
                flag = true;
            }
            total += totalCount;
            match += maxCount;
            if (chosenBase != null) {
                seq.append(chosenBase);
            }
        }
        String SEQ;
        Integer ntSize = softClip.nt.size();
        if (total != 0
                && match / (double)total > 0.9
                && seq.length() / 1.5 > ntSize - seq.length()
                && (seq.length() / (double)ntSize > 0.8
                || ntSize - seq.length() < 10
                || seq.length() > 25)) {
            SEQ = seq.toString();
        } else {
            SEQ = "";
        }

        if (!SEQ.isEmpty() && SEQ.length() > Configuration.SEED_2) {
            Matcher mm1  =  B_A7.matcher(SEQ);
            Matcher mm2  =  B_T7.matcher(SEQ);
            if (mm1.find() || mm2.find()) {
                softClip.used = true;
            }
            if (islowcomplexseq(SEQ)) {
                softClip.used = true;
            }
        }

        if (!SEQ.isEmpty() && SEQ.length() >= Configuration.ADSEED) {
            if ( dir == 3 ) { // 3'
                if (instance().adaptorForward.containsKey(substr(SEQ, 0, Configuration.ADSEED))) {
                    SEQ = "";
                }
            } else if ( dir == 5 ) { // 5'
                if (instance().adaptorReverse.containsKey(reverse(substr(SEQ, 0, Configuration.ADSEED)))) {
                    SEQ = "";
                }
            }
        }

        softClip.sequence = SEQ;

        if (instance().conf.y) {
            System.err.printf("  Candidate consensus: %s Reads: %s M: %s T: %s Final: %s\n", seq, softClip.varsCount, match, total, SEQ);
        }
        return SEQ;
    }

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
     * Construct reference sequence
     * @param baseToPosition map of reference bases    $ref
     * @param from start position
     * @param to end position in double (for end with deletion)
     * @return reference sequence starting at from and ending at to
     */
    public static String joinRef(Map<Integer, Character> baseToPosition,
                                 int from,
                                 double to) {
        StringBuilder sb = new StringBuilder();

        for (int i = from; i < to; i++) {
            Character ch = baseToPosition.get(i);
            if (ch != null) {
                sb.append(ch);
            }
        }
        return sb.toString();
    }

    /**
     * Construct reference sequence
     * @param baseToPosition map of reference bases    $ref
     * @param from start position
     * @param to end position
     * @param seq sequence to go through found in findconseq()
     * @param EXTRA extra subsequence from findMatch()
     * @return reference sequence starting at from and ending at to
     */
    public static String joinRefFor5Lgins(Map<Integer, Character> baseToPosition,
                                          int from,
                                          int to,
                                          String seq,
                                          String EXTRA) {
        StringBuilder sb = new StringBuilder();
        for (int i = from; i <= to; i++) {
            if (to - i < seq.length() - EXTRA.length()) {
                sb.append(charAt(seq, to - i + EXTRA.length()));
            } else {
                Character ch = baseToPosition.get(i);
                if (ch != null) {
                    sb.append(ch);
                }
            }
        }
        return sb.toString();
    }

    /**
     * Construct reference sequence
     * @param baseToPosition map of reference bases    $ref
     * @param from start position
     * @param to end position
     * @param shift5 shift positions
     * @param seq sequence to go through found in findconseq()
     * @param EXTRA extra subsequence from findMatch()
     * @return reference sequence starting at from and ending at to
     */
    public static String joinRefFor3Lgins(Map<Integer, Character> baseToPosition,
                                          int from,
                                          int to,
                                          int shift5,
                                          String seq,
                                          String EXTRA) {
        StringBuilder sb = new StringBuilder();

        for (int i = from; i <= to; i++) {
            if (i - from >= shift5 && i - from - shift5 < seq.length() - EXTRA.length()) {
                sb.append(charAt(seq, i - from - shift5 + EXTRA.length()));
            } else {
                Character ch = baseToPosition.get(i);
                if (ch != null) {
                    sb.append(ch);
                }
            }
        }
        return sb.toString();
    }

    /**
     * Adjust count. Don't adjust reference
     * @param varToAdd variant to add counts    $vref
     * @param variation variant     $tv
     */
    public static void adjCnt(Variation varToAdd,
                              Variation variation) {
        adjCnt(varToAdd, variation, null);
    }

    /**
     * Adjust the count,  If ref is not null, the count for reference is also adjusted.
     * Variant counts of tv are added to vref and removed from ref
     * @param varToAdd variant to add counts    $vref
     * @param variant variant    $tv
     * @param referenceVar reference variant  $ref
     */
    public static void adjCnt(Variation varToAdd,
                              Variation variant,
                              Variation referenceVar) {
        varToAdd.varsCount += variant.varsCount;
        varToAdd.extracnt += variant.varsCount;
        varToAdd.highQualityReadsCount += variant.highQualityReadsCount;
        varToAdd.lowQualityReadsCount += variant.lowQualityReadsCount;
        varToAdd.meanPosition += variant.meanPosition;
        varToAdd.meanQuality += variant.meanQuality;
        varToAdd.meanMappingQuality += variant.meanMappingQuality;
        varToAdd.numberOfMismatches += variant.numberOfMismatches;
        varToAdd.pstd = true;
        varToAdd.qstd = true;
        varToAdd.addDir(true, variant.getDir(true));
        varToAdd.addDir(false, variant.getDir(false));

        if (instance().conf.y) {
            String refCnt = (referenceVar != null && referenceVar.varsCount != 0) ? String.valueOf(referenceVar.varsCount) : "NA";
            System.err.printf("    AdjCnt: '+' %s %s %s %s Ref: %s\n",
                    varToAdd.varsCount, variant.varsCount, varToAdd.getDir(false), variant.getDir(false), refCnt);
            System.err.printf("    AdjCnt: '-' %s %s %s %s Ref: %s\n",
                    varToAdd.varsCount, variant.varsCount, varToAdd.getDir(true), variant.getDir(true), refCnt);
        }

        if (referenceVar == null)
            return;

        referenceVar.varsCount -= variant.varsCount;
        referenceVar.highQualityReadsCount -= variant.highQualityReadsCount;
        referenceVar.lowQualityReadsCount -= variant.lowQualityReadsCount;
        referenceVar.meanPosition -= variant.meanPosition;
        referenceVar.meanQuality -= variant.meanQuality;
        referenceVar.meanMappingQuality -= variant.meanMappingQuality;
        referenceVar.numberOfMismatches -= variant.numberOfMismatches;
        referenceVar.subDir(true, variant.getDir(true));
        referenceVar.subDir(false, variant.getDir(false));
        correctCnt(referenceVar);
    }

    /**
     * Correct counts for negative values
     * @param varToCorrect variant   $ref
     */
    public static void correctCnt(Variation varToCorrect) {
        if (varToCorrect.varsCount < 0)
            varToCorrect.varsCount = 0;
        if (varToCorrect.highQualityReadsCount < 0)
            varToCorrect.highQualityReadsCount = 0;
        if (varToCorrect.lowQualityReadsCount < 0)
            varToCorrect.lowQualityReadsCount = 0;
        if (varToCorrect.meanPosition < 0)
            varToCorrect.meanPosition = 0;
        if (varToCorrect.meanQuality < 0)
            varToCorrect.meanQuality = 0;
        if (varToCorrect.meanMappingQuality < 0)
            varToCorrect.meanMappingQuality = 0;
        if (varToCorrect.getDir(true) < 0)
            varToCorrect.addDir(true, -varToCorrect.getDir(true));
        if (varToCorrect.getDir(false) < 0)
            varToCorrect.addDir(false, -varToCorrect.getDir(false));
    }

    /**
     * Returns var from map of vars or null if map doesn't contain such key.
     * @param alignedVariants alignedVariants $vars
     * @param key position on which vars are searching
     * @param type type of VarsType
     * @param keys often the string that variant must contain
     * @return variant from map of vars of specified VarsType or null
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

    /**
     * Returns variant from Vars object or null if map doesn't contain such key.
     * @param vars alignedVariants $vars
     * @param type type of VarsType
     * @param keys often the string that variant must contain
     * @return variant from vars of specified VarsType or null
     */
    public static Variant getVarMaybe(Vars vars,
                               VarsType type,
                               Object... keys) {
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
     * If map contains Vars on position, then returns it, else create Vars, put it on position and returns.
     * @param map map of positions and vars
     * @param position start position of vars
     * @return Vars object on position
     */
    public static Vars getOrPutVars(Map<Integer, Vars> map, int position) {
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

    /**
     * Get {@link Variation} from map. If map contains variation, return it. If not, add variation with description string
     * at specific start position.
     * @param hash map contain position and description string on variations.
     * @param start start position of Variation
     * @param descriptionString string contains information about variation (length and type of variation)
     * @return variation
     */
    public static Variation getVariation(Map<Integer, VariationMap<String, Variation>> hash,
                                         int start,
                                         String descriptionString) {
        VariationMap<String, Variation> map = hash.get(start);
        if (map == null) {
            map = new VariationMap<>();
            hash.put(start, map);
        }
        Variation variation = map.get(descriptionString);
        if (variation == null) {
            variation = new Variation();
            map.put(descriptionString, variation);
        }
        return variation;
    }

    /**
     * Get {@link Variation} from map. If map contains variation, return it. If not, return null.
     * @param hash map contain position and description string on variations.
     * @param start start position of Variation
     * @param refBase contains nucleotide from the reference
     * @return variation
     */
    public static Variation getVariationMaybe(Map<Integer, VariationMap<String, Variation>> hash,
                                              int start,
                                              Character refBase) {
        if (refBase == null)
            return null;

        Map<String, Variation> map = hash.get(start);
        if (map == null) {
            return null;
        }

        return map.get(refBase.toString());
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
