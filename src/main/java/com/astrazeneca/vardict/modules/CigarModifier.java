package com.astrazeneca.vardict.modules;

import com.astrazeneca.vardict.Configuration;
import com.astrazeneca.vardict.collection.Tuple;
import com.astrazeneca.vardict.data.ModifiedCigar;
import com.astrazeneca.vardict.data.Reference;
import com.astrazeneca.vardict.data.Region;
import jregex.Replacer;

import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;

import static com.astrazeneca.vardict.Utils.*;
import static com.astrazeneca.vardict.collection.Tuple.tuple;
import static com.astrazeneca.vardict.data.Patterns.*;
import static com.astrazeneca.vardict.data.scopedata.GlobalReadOnlyScope.instance;
import static com.astrazeneca.vardict.variations.VariationUtils.isHasAndEquals;
import static com.astrazeneca.vardict.variations.VariationUtils.isHasAndNotEquals;

/**
 * Class for process CIGAR string of record in SAM/BAM file and modify it if it fits to matchers.
 */
public class CigarModifier {
    private int position;
    private String cigarStr;
    private final String originalCigar;
    private String querySequence;
    private String queryQuality;
    private Map<Integer, Character> reference;
    private Reference ref;
    private int indel;
    private int maxReadLength;
    private Region region;

    CigarModifier(int position, String cigarStr, String querySequence,
                  String queryQuality, Reference ref, int indel, Region region, int maxReadLength) {
        this.position = position;
        this.cigarStr = cigarStr;
        this.originalCigar = cigarStr;
        this.querySequence = querySequence;
        this.queryQuality = queryQuality;
        this.ref = ref;
        this.reference = ref.referenceSequences;
        this.indel = indel;
        this.region = region;
        this.maxReadLength = maxReadLength;
    }

    /**
     * Method modify cigar by applying different matchers.
     * @return modified position and cigar string
     */
    public ModifiedCigar modifyCigar() {
        //flag is set to true if CIGAR string is modified and should be looked at again
        boolean flag = true;
        try {
            // if CIGAR starts with deletion cut it off
            jregex.Matcher mm = BEGIN_NUMBER_D.matcher(cigarStr);
            if (mm.find()) {
                position += toInt(mm.group(1));
                Replacer r = BEGIN_NUMBER_D.replacer("");
                cigarStr = r.replace(cigarStr);
            }
            mm = END_NUMBER_D.matcher(cigarStr);
            if (mm.find()) {
                Replacer r = END_NUMBER_D.replacer("");
                cigarStr = r.replace(cigarStr);
            }
            // replace insertion at the beginning and end with soft clipping
            mm = BEGIN_NUMBER_I.matcher(cigarStr);
            if (mm.find()) {
                Replacer r = BEGIN_NUMBER_I.replacer(toInt(mm.group(1)) + "S");
                cigarStr = r.replace(cigarStr);
            }
            mm = END_NUMBER_I.matcher(cigarStr);
            if (mm.find()) {
                Replacer r = END_NUMBER_I.replacer(toInt(mm.group(1)) + "S");
                cigarStr = r.replace(cigarStr);
            }
            Matcher sc5_mm = SA_CIGAR_D_S_5clip_GROUP.matcher(cigarStr);
            Matcher sc3_mm = SA_CIGAR_D_S_3clip_GROUP.matcher(cigarStr);
            Map<String, List<Integer>> referenceSeedMap = ref.seed;
            if (sc5_mm.find()) {
               int cigarElement = toInt(sc5_mm.group(1));
               if (!instance().conf.chimeric && cigarElement >= Configuration.SEED_2 ) {
                   String sseq = substr(querySequence, 0, cigarElement);
                   String sequence = complement(reverse(sseq));
                   String reverseComplementedSeed = sequence.substring(0, Configuration.SEED_2);

                   if (referenceSeedMap.containsKey(reverseComplementedSeed)) {
                       List<Integer> positions = referenceSeedMap.get(reverseComplementedSeed);

                       if (positions.size() == 1 && Math.abs(position - positions.get(0)) < 2 * maxReadLength) {
                           Replacer r = SA_CIGAR_D_S_5clip_GROUP_Repl.replacer(""); // $a[5] = ~s / ^\d + S//;
                           cigarStr = r.replace(cigarStr);
                           querySequence = substr(querySequence, cigarElement);
                           queryQuality = substr(queryQuality, cigarElement);
                           if (instance().conf.y) {
                               System.err.println(sequence + " at 5' is a chimeric at "
                                       + position + " by SEED " + Configuration.SEED_2);
                           }
                       }
                   }
               }
           } else if (sc3_mm.find()) {
                int cigarElement = toInt(sc3_mm.group(1));
                if (!instance().conf.chimeric && cigarElement >= Configuration.SEED_2 ) {
                    String sseq = substr(querySequence, -cigarElement, cigarElement);
                    String sequence = complement(reverse(sseq));
                    String reverseComplementedSeed = substr(sequence, -Configuration.SEED_2, Configuration.SEED_2);

                    if (referenceSeedMap.containsKey(reverseComplementedSeed)) {
                        List<Integer> positions = referenceSeedMap.get(reverseComplementedSeed);

                        if (positions.size() == 1 && Math.abs(position - positions.get(0)) < 2* maxReadLength ) {
                            Replacer r = SA_CIGAR_D_S_3clip_GROUP_Repl.replacer("");  // $a[5] = ~s /\d\d + S$//;
                            cigarStr = r.replace(cigarStr);
                            querySequence = substr(querySequence, 0, querySequence.length() - cigarElement);
                            queryQuality = substr(queryQuality, 0, queryQuality.length() - cigarElement);
                            if (instance().conf.y) {
                                System.err.println(sequence + " at 3' is a chimeric at "
                                        + position + " by SEED " + Configuration.SEED_2);
                            }
                        }
                   }
               }
           }

            while (flag && indel > 0) {
                flag = false;
                mm = BEGIN_NUMBER_S_NUMBER_IorD.matcher(cigarStr);
                if (mm.find()) { // If CIGAR starts with soft-clipping followed by insertion or deletion
                    /*
                    If insertion follows soft-clipping, add the inserted sequence to soft-clipped start
                    Otherwise increase position by number of deleted bases
                     */
                    String tslen = toInt(mm.group(1)) + (mm.group(3).equals("I") ? toInt(mm.group(2)) : 0) + "S";
                    position += mm.group(3).equals("D") ? toInt(mm.group(2)) : 0;
                    //Regexp replaces found CIGAR sequence with tslen (number + S)
                    Replacer r = BEGIN_NUMBER_S_NUMBER_IorD.replacer(tslen);
                    cigarStr = r.replace(cigarStr);
                    flag = true;
                }
                mm = NUMBER_IorD_NUMBER_S_END.matcher(cigarStr);
                if (mm.find()) { // If CIGAR ends with insertion or deletion followed by soft-clipping
                    //Replace insertion or deletion with soft-clipping
                    String tslen = toInt(mm.group(3)) + (mm.group(2).equals("I") ? toInt(mm.group(1)) : 0) + "S";
                    //Regexp replaces found CIGAR sequence with $tslen (number + S)
                    Replacer r = NUMBER_IorD_NUMBER_S_END.replacer(tslen);
                    cigarStr = r.replace(cigarStr);
                    flag = true;
                }
                mm = BEGIN_NUMBER_S_NUMBER_M_NUMBER_IorD.matcher(cigarStr);
                if (mm.find()) { // If CIGAR starts with soft-clipping followed by matched sequence and insertion or deletion
                    int tmid = toInt(mm.group(2));
                    if (tmid <= 10) { // If matched sequence length is no more than 10, replace everything with soft-clipping
                        String tslen = toInt(mm.group(1)) + tmid + (mm.group(4).equals("I") ? toInt(mm.group(3)) : 0) + "S";
                        position += tmid + (mm.group(4).equals("D") ? toInt(mm.group(3)) : 0);
                        Replacer r = BEGIN_NUMBER_S_NUMBER_M_NUMBER_IorD.replacer(tslen);
                        cigarStr = r.replace(cigarStr);
                        flag = true;
                    }
                }
                mm = NUMBER_IorD_NUMBER_M_NUMBER_S_END.matcher(cigarStr);
                if (mm.find()) { // If CIGAR ends with insertion or deletion, matched sequence and soft-clipping
                    int tmid = toInt(mm.group(3));
                    if (tmid <= 10) { // If matched sequence length is no more than 10, replace everything with soft-clipping
                        String tslen = toInt(mm.group(4)) + tmid + (mm.group(2).equals("I") ? toInt(mm.group(1)) : 0) + "S";
                        Replacer r = NUMBER_IorD_NUMBER_M_NUMBER_S_END.replacer(tslen);
                        cigarStr = r.replace(cigarStr);
                        flag = true;
                    }
                }

                // The following two clauses to make indels at the end of reads as softly
                // clipped reads and let VarDict's algorithm identify indels

                //If CIGAR starts with 1-9 bases long matched sequence, insertion or deletion and matched sequence
                mm = BEGIN_DIGIT_M_NUMBER_IorD_NUMBER_M.matcher(cigarStr);
                if (mm.find()) {
                    flag = beginDigitMNumberIorDNumberM(mm);
                }
                mm = NUMBER_IorD_DIGIT_M_END.matcher(cigarStr);
                if (mm.find()) { //If CIGAR ends with insertion or deletion and 1-9 bases long matched sequence
                    int tmid = toInt(mm.group(3));
                    String tslen = tmid + (mm.group(2).equals("I") ? toInt(mm.group(1)) : 0) + "S";
                    Replacer r = NUMBER_IorD_NUMBER_M_END.replacer(tslen);
                    cigarStr = r.replace(cigarStr);
                    flag = true;
                }

                // Combine two deletions and insertion into one complex if they are close
                mm = D_M_D_DD_M_D_I_D_M_D_DD.matcher(cigarStr);
                jregex.Matcher threeIndelsMatcher = threeIndelsPattern.matcher(cigarStr);
                jregex.Matcher threeDeletionsMatcher = threeDeletionsPattern.matcher(cigarStr);
                if (mm.find()) {
                    flag = twoDeletionsInsertionToComplex(mm, flag);
                } else if (threeDeletionsMatcher.find()) {
                    flag = threeDeletions(threeDeletionsMatcher, flag);
                } else if (threeIndelsMatcher.find()) {
                    flag = threeIndels(threeIndelsMatcher, flag);
                }

                Matcher cm = DIG_D_DIG_M_DIG_DI_DIGI.matcher(cigarStr);
                if (cm.find()) {
                    flag  = combineToCloseToCorrect(cm, flag);
                }

                cm = NOTDIG_DIG_I_DIG_M_DIG_DI_DIGI.matcher(cigarStr);
                if (cm.find() && !"D".equals(cm.group(1)) && !"H".equals(cm.group(1))) {
                    flag = combineToCloseToOne(cm, flag);
                }
                cm = DIG_D_DIG_D.matcher(cigarStr);
                if (cm.find()) {
                    int dlen = toInt(cm.group(1)) + toInt(cm.group(2));
                    cigarStr = cm.replaceFirst(dlen + "D");
                    flag = true;
                }
                cm = DIG_I_DIG_I.matcher(cigarStr);
                if (cm.find()) {
                    int ilen = toInt(cm.group(1)) + toInt(cm.group(2));
                    cigarStr = cm.replaceFirst(ilen + "I");
                    flag = true;
                }
            }

            Matcher mtch = ANY_NUMBER_M_NUMBER_S_END.matcher(cigarStr);
            if (mtch.find()) {
                captureMisSoftlyMS(mtch);
            } else if ((mtch = BEGIN_ANY_DIG_M_END.matcher(cigarStr)).find()) {
                captureMisSoftly3Mismatches(mtch);
            }

            mtch = DIG_S_DIG_M.matcher(cigarStr);
            if (mtch.find()) {
                combineDigSDigM(mtch);
            } else if ((mtch = BEGIN_DIG_M.matcher(cigarStr)).find()) {
                combineBeginDigM(mtch);
            }
            return new ModifiedCigar(position, cigarStr, querySequence, queryQuality);
        } catch (Exception exception) {
            printExceptionAndContinue(exception, "cigar", String.valueOf(position) + " " + originalCigar, region);
        }
        return new ModifiedCigar(position, cigarStr, querySequence, queryQuality);
    }

    /**
     * Combine matched sequence to soft clip and matched and realign position. Make >=3 mismatches in the end as soft clipping
     * @param matcher Regexp find matched sequence at the start of the CIGAR
     */
    private void combineBeginDigM(Matcher matcher) {
        int mch = toInt(matcher.group(1));
        int rn = 0;
        int rrn = 0;
        int rmch = 0;
        while (rrn < mch && rn < mch) {
            if (!reference.containsKey(position + rrn)) {
                break;
            }
            if (isHasAndNotEquals(reference, position + rrn, querySequence, rrn)) {
                rn = rrn + 1;
                rmch = 0;
            } else if (isHasAndEquals(reference, position + rrn, querySequence, rrn)) {
                rmch++;
            }
            rrn++;
            if (rmch >= 3) {
                break;
            }
        }
        if (rn > 0 && rn <= 3) {
            mch -= rn;
            cigarStr = matcher.replaceFirst(rn + "S" + mch + "M");
            position += rn;
        }
    }

    /**
     * Combine soft clip and matched sequence and realign position
     * @param matcher Regexp find softclip and matched sequence at the start of the CIGAR
     */
    private void combineDigSDigM(Matcher matcher) {
        //length of matched sequence
        int mch = toInt(matcher.group(2));
        //length of soft-clipping
        int soft = toInt(matcher.group(1));
        //number of bases before matched sequence that match in reference and read sequences
        int rn = 0;
        Set<Character> RN = new HashSet<>();
        while (rn < soft && isHasAndEquals(reference, position - rn - 1, querySequence, soft - rn - 1)
                && queryQuality.charAt(soft - rn - 1) - 33 > Configuration.LOWQUAL) {
            rn++;
        }

        if (rn > 0) {
            mch += rn;
            soft -= rn;
            if (soft > 0) {
                cigarStr = matcher.replaceFirst(soft + "S" + mch + "M");
            } else {
                cigarStr = matcher.replaceFirst(mch + "M");
            }
            position -= rn;
            rn = 0;
        }
        if (soft > 0) {
            while (rn + 1 < soft && isHasAndEquals(reference, position - rn - 2, querySequence, soft - rn - 2) &&
                    queryQuality.charAt(soft - rn - 2) - 33 > Configuration.LOWQUAL) {
                rn++;
                RN.add(reference.get(position - rn - 2));
            }
            int rn_nt = RN.size();  // Don't adjust for homopolymers
            if ((rn > 4 && rn_nt > 1) || isHasAndEquals(reference, position - 1, querySequence, soft - 1)) {
                mch += rn + 1;
                soft -= rn + 1;
                if (soft > 0) {
                    cigarStr = DIG_S_DIG_M.matcher(cigarStr).replaceFirst(soft + "S" + mch + "M");
                } else {
                    cigarStr = DIG_S_DIG_M.matcher(cigarStr).replaceFirst(mch + "M");
                }
                position -= rn + 1;
            }
            if (rn == 0) {
                int rrn = 0;
                int rmch = 0;
                while (rrn < mch && rn < mch) {
                    if (!reference.containsKey(position + rrn)) {
                        break;
                    }
                    if (isHasAndNotEquals(reference, position + rrn, querySequence, soft + rrn)) {
                        rn = rrn + 1;
                        rmch = 0;
                    } else if (isHasAndEquals(reference, position + rrn, querySequence, soft + rrn)) {
                        rmch++;
                    }
                    rrn++;
                    if (rmch >= 3) {
                        break;
                    }
                }
                if (rn > 0 && rn < mch) {
                    soft += rn;
                    mch -= rn;
                    cigarStr = DIG_S_DIG_M.matcher(cigarStr).replaceFirst(soft + "S" + mch + "M");
                    position += rn;
                }
            }
        }
    }

    /**
     * Trying to capture sometimes mis-softly clipped reads by aligner. Make >=3 mismatches in the end as soft clipping
     * @param matcher Regexp find the matched sequence only
     */
    private void captureMisSoftly3Mismatches(Matcher matcher) {
        String ov5 = matcher.group(1);
        int mch = toInt(matcher.group(2));
        int refoff = position + toInt(matcher.group(2));
        int rdoff = toInt(matcher.group(2));
        if (!ov5.isEmpty()) { //If prefix is present
            //Add all matched and deletion lengths to reference position
            refoff += sum(globalFind(ALIGNED_LENGTH_MND, ov5)); // reference position
            // Add all matched, insertion and soft-clipped lengths to read position
            rdoff += sum(globalFind(SOFT_CLIPPED, ov5)); // read position
        }
        int rn = 0;
        int rrn = 0;
        int rmch = 0;
        while (rrn < mch && rn < mch) {
            if (!reference.containsKey(refoff - rrn - 1)) {
                break;
            }
            if (rrn < rdoff && isHasAndNotEquals(reference, refoff - rrn - 1, querySequence, rdoff - rrn - 1)) {
                rn = rrn + 1;
                rmch = 0;
            } else if (rrn < rdoff && isHasAndEquals(reference, refoff - rrn - 1, querySequence, rdoff - rrn - 1)) {
                rmch++;
            }
            rrn++;
            // Stop at three consecure matches
            if (rmch >= 3) {
                break;
            }
        }
        mch -= rn;
        if (rn > 0 && rn <= 3) {
            cigarStr = DIG_M_END.matcher(cigarStr).replaceFirst(mch + "M" + rn + "S");
        }
    }

    /**
     * Trying to capture sometimes mis-softly clipped reads by aligner
     * @param matcher Regexp finds number-M-number-S at end of CIGAR string ^(.*?) captures everything before M-S complex
     */
    private void captureMisSoftlyMS(Matcher matcher) {
        //prefix of CIGAR string before last matched sequence
        String ov5 = matcher.group(1);
        //length of matched sequence
        int mch = toInt(matcher.group(2));
        //length of soft-clipping
        int soft = toInt(matcher.group(3));
        //offset of soft-clipped sequence in the reference string (position + length of matched)
        int refoff = position + toInt(matcher.group(2));
        //offset of soft-clipped sequence in the read
        int rdoff = toInt(matcher.group(2));
        if (!ov5.isEmpty()) { //If prefix is present
            //Add all matched and deletion lengths to reference position
            refoff += sum(globalFind(ALIGNED_LENGTH_MND, ov5)); // reference position
            //Add all matched, insertion and soft-clipped lengths to read position
            rdoff += sum(globalFind(SOFT_CLIPPED, ov5)); // read position
        }
        //number of bases after refoff/rdoff that match in reference and read sequences
        int rn = 0;
        Set<Character> RN = new HashSet<>();
        while (rn < soft && isHasAndEquals(reference, refoff + rn, querySequence, rdoff + rn)
                && queryQuality.charAt(rdoff + rn) - 33 > Configuration.LOWQUAL) {
            rn++;
        }
        if (rn > 0) {
            mch += rn;
            soft -= rn;
            if (soft > 0) {
                cigarStr = DIG_M_DIG_S_END.matcher(cigarStr).replaceFirst(mch + "M" + soft + "S");
            } else {
                cigarStr = DIG_M_DIG_S_END.matcher(cigarStr).replaceFirst(mch + "M");
            }
            rn = 0;
        }

        if (soft > 0) {
            while (rn + 1 < soft && isHasAndEquals(reference, refoff + rn + 1, querySequence, rdoff + rn + 1) &&
                    queryQuality.charAt(rdoff + rn + 1) - 33 > Configuration.LOWQUAL) {
                rn++;
                RN.add(reference.get(refoff + rn + 1)); // to my mind in perl value increasing here ? not adding el
            }
            int rn_nt = RN.size(); // don't adjust if homopolymer
            if (rn > 4 && rn_nt > 1) {
                //If more than 3 bases
                // match after refoff/rdoff or base at refoff/rdoff match
                mch += rn + 1;
                soft -= rn + 1;
                if (soft > 0) {
                    cigarStr = DIG_M_DIG_S_END.matcher(cigarStr).replaceFirst(mch + "M" + soft + "S");
                } else {
                    cigarStr = DIG_M_DIG_S_END.matcher(cigarStr).replaceFirst(mch + "M");
                }
            }
            if (rn == 0) {
                int rrn = 0;
                int rmch = 0;
                while (rrn < mch && rn < mch) {
                    if (!reference.containsKey(refoff - rrn - 1)) {
                        break;
                    }
                    if (rrn < rdoff && isHasAndNotEquals(reference, refoff - rrn - 1, querySequence, rdoff - rrn - 1)) {
                        rn = rrn + 1;
                        rmch = 0;
                    } else if (rrn < rdoff && isHasAndEquals(reference, refoff - rrn - 1, querySequence, rdoff - rrn - 1)) {
                        rmch++;
                    }
                    rrn++;
                    // Stop at three consecure matches
                    if (rmch >= 3) {
                        break;
                    }
                }

                if (rn > 0 && rn < mch) {
                    soft += rn;
                    mch -= rn;
                    cigarStr = DIG_M_DIG_S_END.matcher(cigarStr).replaceFirst(mch + "M" + soft + "S");
                }
            }
        }
    }

    /**
     * Combine two close indels (<15bp) into one
     * @param matcher Regexp finds not_digit-number-I-number-M-number-D and optionally number-I
     * @param flag flag to determine if read was processed before
     * @return true if length of internal matched sequences is no more than 15
     */
    private boolean combineToCloseToOne(Matcher matcher, boolean flag) {
        //If CIGAR string contains any letter except D, matched sequence, deletion and possibly insertion
        String op = matcher.group(5);
        int g2 = toInt(matcher.group(2));
        int g3 = toInt(matcher.group(3));
        int g4 = toInt(matcher.group(4));

        if (g3 <= 15) {
            //length of matched sequence and deletion
            int dlen = g3;
            //length of first insertion and matched sequence
            int ilen = g2 + g3;
            if (op.equals("I")) {
                ilen += g4;
            } else if (op.equals("D")) {
                dlen += g4;
                //last insertion string
                if (matcher.groupCount() > 5) {
                    String istr = matcher.group(6);
                    if (istr != null) { //If digit and after it Insertion is present after deletion,
                        // add its length to ilen
                        ilen += toInt(istr.substring(0, istr.length() - 1));
                    }
                }
            }
            matcher = DIG_I_DIG_M_DIG_DI_DIGI.matcher(cigarStr);
            //Replace I-M-D-I? complex with deletion and insertion
            cigarStr = matcher.replaceFirst(dlen + "D" + ilen + "I");
            flag = true;
        }
        return flag;
    }

    /**
     * Combine two close deletions (<15bp) into one correct from <10 to <15.
     * If CIGAR string contains deletion, short (<= 9 bases) matched sequence, deletion and possibly insertion,
     * replace with deletion and insertion.
     * @param matcher Regexp finds number-D-digit-M-number-D and optionally number-I
     * @param flag flag to determine if read was processed before
     * @return true if length of internal matched sequences is no more than 15
     */
    private boolean combineToCloseToCorrect(Matcher matcher, boolean flag) {
        int g1 = toInt(matcher.group(1));
        int g2 = toInt(matcher.group(2));
        int g3 = toInt(matcher.group(3));
        if (g2 <= 15) {
            String op = matcher.group(4);
            //length of both deletions and matched sequence
            int dlen = g1 + g2;
            //matched sequence length
            int ilen = g2;
            if (op.equals("I")) {
                ilen += g3;
            } else if (op.equals("D")) {
                dlen += g3;
                if (matcher.groupCount() > 4) {
                    //insertion string
                    String istr = matcher.group(5);
                    if (istr != null) { //If insertion is present after 2nd deletion, add its length to $ilen
                        ilen += toInt(istr.substring(0, istr.length() - 1));
                    }
                }
            }
            cigarStr = matcher.replaceFirst(dlen + "D" + ilen + "I");
            flag = true;
        }
        return flag;
    }

    /**
     * If CIGAR string contains of 3 deletions/insertions with matched sequences, replace with deletion and matched,
     * or deletion, insertion and matched if matched and insertion more then zero
     * @param matcher regexp found 3 deletions/insertions with matched sequences between them
     * @param flag to determine if read was processed before
     * @return true if length of internal matched sequences is no more than 15
     */
    private boolean threeIndels(jregex.Matcher matcher, boolean flag) {
        int tslen = toInt(matcher.group(5)) + toInt(matcher.group(8));
        if ("I".equals(matcher.group(4))) {
            tslen += toInt(matcher.group(3));
        }
        if ("I".equals(matcher.group(7))) {
            tslen += toInt(matcher.group(6));
        }
        if ("I".equals(matcher.group(10))) {
            tslen += toInt(matcher.group(9));
        }

        int dlen = toInt(matcher.group(5)) + toInt(matcher.group(8));
        if ("D".equals(matcher.group(4))) {
            dlen += toInt(matcher.group(3));
        }
        if ("D".equals(matcher.group(7))) {
            dlen += toInt(matcher.group(6));
        }
        if ("D".equals(matcher.group(10))) {
            dlen += toInt(matcher.group(9));
        }

        int mid = toInt(matcher.group(5)) + toInt(matcher.group(8));
        String ov5 = matcher.group(1);
        int refoff = position + toInt(matcher.group(2));
        int rdoff = toInt(matcher.group(2));
        int RDOFF = toInt(matcher.group(2));
        int rm = toInt(matcher.group(11));

        if (!ov5.isEmpty()) { //If the complex is not at start of CIGAR string
            refoff += sum(globalFind(ALIGNED_LENGTH_MND, ov5)); // reference position
            rdoff += sum(globalFind(SOFT_CLIPPED, ov5)); // read position
        }

        int rn = 0;
        while (rdoff + rn < querySequence.length() && isHasAndEquals(querySequence.charAt(rdoff + rn), reference, refoff + rn)) {
            rn++;
        }
        RDOFF += rn;
        dlen -= rn;
        tslen -= rn;

        String newCigarStr = RDOFF + "M";
        if (tslen <= 0) {
            dlen -= tslen;
            rm += tslen;
            if (dlen == 0) {
                RDOFF = RDOFF + rm;
                newCigarStr = RDOFF + "M";
            } else if (dlen < 0) {
                tslen = -dlen;
                rm += dlen;
                if (rm < 0) {
                    RDOFF = RDOFF + rm;
                    newCigarStr = RDOFF + "M" + tslen + "I";
                } else {
                    newCigarStr += tslen + "I" + rm + "M";
                }
            } else {
                newCigarStr += dlen + "D" + rm + "M";
            }
        } else {
            if (dlen == 0) {
                newCigarStr += tslen + "I" + rm + "M";
            } else if (dlen < 0) {
                rm += dlen;
                newCigarStr += tslen + "I" + rm + "M";
            } else {
                newCigarStr += dlen + "D" + tslen + "I" + rm + "M";
            }
        }
        if (mid <= 15) {
            cigarStr = DIGM_D_DI_DIGM_D_DI_DIGM_DI_DIGM.matcher(cigarStr).replaceFirst(newCigarStr);
            flag = true;
        }
        return flag;
    }

    /**
     * If CIGAR string contains of 3 deletions with matched sequences, replace with deletion and matched
     * or deletion, insertion and matched if matched and insertion more then zero
     * @param matcher regexp found 3 deletions with matched sequences between them
     * @param flag to determine if read was processed before
     * @return true if length of internal matched sequences is no more than 15
     */
    private boolean threeDeletions(jregex.Matcher matcher, boolean flag) {
        //length of both matched sequences and insertion
        int tslen = toInt(matcher.group(4)) + toInt(matcher.group(6));
        //length of deletions and internal matched sequences
        int dlen = toInt(matcher.group(3)) + toInt(matcher.group(4)) +
                toInt(matcher.group(5)) + toInt(matcher.group(6)) +
                toInt(matcher.group(7));

        //length of internal matched sequences
        int mid = toInt(matcher.group(4)) + toInt(matcher.group(6));

        //prefix of CIGAR string before the M-D-M-I-M-D complex
        String ov5 = matcher.group(1);

        //offset of first deletion in the reference sequence
        int refoff = position + toInt(matcher.group(2));

        //offset of first deletion in the read
        int rdoff = toInt(matcher.group(2));

        //offset of first deletion in the read corrected by possibly matching bases
        int RDOFF = toInt(matcher.group(2));

        int rm = toInt(matcher.group(8));

        if (!ov5.isEmpty()) { //If the complex is not at start of CIGAR string
            refoff += sum(globalFind(ALIGNED_LENGTH_MND, ov5)); // reference position // - 937 -
            rdoff += sum(globalFind(SOFT_CLIPPED, ov5)); // read position
            // MND changes from 27.04.2018
        }

        //number of bases after refoff/rdoff that match in reference and read
        int rn = 0;
        while (rdoff + rn < querySequence.length() && isHasAndEquals(querySequence.charAt(rdoff + rn), reference, refoff + rn)) {
            rn++;
        }
        RDOFF += rn;
        dlen -= rn;
        tslen -= rn;
        String newCigarStr = RDOFF + "M";
        if (tslen <= 0) {
            dlen -= tslen;
            rm += tslen;
            newCigarStr += dlen + "D" + rm + "M";
        } else {
            newCigarStr += dlen + "D" + tslen + "I" + rm + "M";
        }
        if (mid <= 15) {
            cigarStr = DM_DD_DM_DD_DM_DD_DM.matcher(cigarStr).replaceFirst(newCigarStr);
            flag = true;
        }
        return flag;
    }

    /**
     * If CIGAR string contains matched sequence, deletion, short (<=9 bases) matched
       sequence, insertion, short (<=9 bases) matched sequence and deletion, replace with deletion and matched
     * @param matcher Regexp finds number-M-number-D-digit-M-number-I-digit-M-number-D. ^(.*?) captures everything before the
     * M-D-M-I-M-D complex
     * @param flag to determine if read was processed before
     * @return true if length of internal matched sequences is no more than 15
     */
    private boolean twoDeletionsInsertionToComplex(jregex.Matcher matcher, boolean flag) {
        //length of both matched sequences and insertion
        int tslen = toInt(matcher.group(4)) + toInt(matcher.group(5)) + toInt(matcher.group(6));

        //length of deletions and internal matched sequences
        int dlen = toInt(matcher.group(3)) + toInt(matcher.group(4)) + toInt(matcher.group(6)) + toInt(matcher.group(7));

        //length of internal matched sequences
        int mid = toInt(matcher.group(4)) + toInt(matcher.group(6));

        //prefix of CIGAR string before the M-D-M-I-M-D complex
        String ov5 = matcher.group(1);

        //offset of first deletion in the reference sequence
        int refoff = position + toInt(matcher.group(2));

        //offset of first deletion in the read
        int rdoff = toInt(matcher.group(2));

        //offset of first deletion in the read corrected by possibly matching bases
        int RDOFF = toInt(matcher.group(2));

        int rm = toInt(matcher.group(8));

        if (!ov5.isEmpty()) { //If the complex is not at start of CIGAR string
            refoff += sum(globalFind(ALIGNED_LENGTH_MND, ov5)); // reference position
            rdoff += sum(globalFind(SOFT_CLIPPED, ov5)); // read position
        }

        //number of bases after refoff/rdoff that match in reference and read
        int rn = 0;
        while (rdoff + rn < querySequence.length() && isHasAndEquals(querySequence.charAt(rdoff + rn), reference, refoff + rn)) {
            rn++;
        }
        RDOFF += rn;
        dlen -= rn;
        tslen -= rn;

        String newCigarStr = RDOFF + "M";
        if (tslen <= 0) {
            dlen -= tslen;
            rm += tslen;
            newCigarStr += dlen + "D" + rm + "M";
        } else {
            newCigarStr += dlen + "D" + tslen + "I" + rm + "M";
        }
        if (mid <= 15) {
            //If length of internal matched sequences is no more than 15, replace M-D-M-I-M-D complex with M-D-I
            cigarStr = D_M_D_DD_M_D_I_D_M_D_DD_prim.matcher(cigarStr).replaceFirst(newCigarStr);
            flag = true;
        }
        return flag;
    }

    /**
     * If CIGAR starts with 1-9 bases long matched sequence, insertion or deletion and matched sequence,
     * sequence and insertion/deletion are replaced with soft-clipping up to 1st matching base in last matched sequence.
     * For deletion position is adjusted by deletion length.
     * @param matcher Regexp finds digit-M-number-(I or D)-number-M
     * @return true (cigar was processed)
     */
    private boolean beginDigitMNumberIorDNumberM(jregex.Matcher matcher) {
        int tmid = toInt(matcher.group(1));
        int mlen = toInt(matcher.group(4));
        int tn = 0;

        int tslen = tmid + (matcher.group(3).equals("I") ? toInt(matcher.group(2)) : 0);
        position += tmid + (matcher.group(3).equals("D") ? toInt(matcher.group(2)) : 0);
        while ((tn < mlen) && isHasAndNotEquals(querySequence.charAt(tslen + tn), reference, position + tn)) {
            tn++;
        }
        tslen += tn;
        mlen -= tn;
        position += tn;
        cigarStr = BEGIN_DIGIT_M_NUMBER_IorD_NUMBER_M_.matcher(cigarStr).replaceFirst(tslen + "S" + mlen + "M");
        return true;
    }
}
