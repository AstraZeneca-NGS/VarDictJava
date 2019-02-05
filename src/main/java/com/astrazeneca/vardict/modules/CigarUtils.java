package com.astrazeneca.vardict.modules;

import com.astrazeneca.vardict.Configuration;
import com.astrazeneca.vardict.collection.Tuple;
import com.astrazeneca.vardict.data.Region;
import jregex.Replacer;

import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;

import static com.astrazeneca.vardict.Utils.*;
import static com.astrazeneca.vardict.collection.Tuple.tuple;
import static com.astrazeneca.vardict.data.Patterns.*;
import static com.astrazeneca.vardict.variations.VariationUtils.isHasAndEquals;
import static com.astrazeneca.vardict.variations.VariationUtils.isHasAndNotEquals;

public class CigarUtils {
    private int position;
    private String cigarStr;
    private final String originalCigar;
    private String querySequence;
    private String queryQuality;
    private Map<Integer, Character> reference;
    private int indel;
    private Region region;

    CigarUtils(int position, String cigarStr, String querySequence,
               String queryQuality, Map<Integer, Character> reference, int indel, Region region) {
        this.position = position;
        this.cigarStr = cigarStr;
        this.originalCigar = cigarStr;
        this.querySequence = querySequence;
        this.queryQuality = queryQuality;
        this.reference = reference;
        this.indel = indel;
        this.region = region;
    }

    public Tuple.Tuple2<Integer, String> modifyCigar() {
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
                if (mm.find()) { //If CIGAR string contains matched sequence, deletion, short (<=9 bases) matched
                    // sequence, insertion, short (<=9 bases) matched sequence and deletion
                    flag = twoDeletionsInsertionToComplex(mm, flag);
                } else if (threeDeletionsMatcher.find()) {
                    flag = threeDeletions(threeDeletionsMatcher, flag);
                } else if (threeIndelsMatcher.find()) {
                    flag = threeIndels(threeIndelsMatcher, flag);
                }
                // Combine two close deletions (<15bp) into one correct from <10 to <15
                Matcher cm = DIG_D_DIG_M_DIG_DI_DIGI.matcher(cigarStr);
                if (cm.find()) { //If CIGAR string contains deletion, short (<= 9 bases) matched sequence, deletion and
                    // possibly insertion
                    flag  = combineToCloseToCorrect(cm, flag);
                }

                // Combine two close indels (<15bp) into one
                cm = NOTDIG_DIG_I_DIG_M_DIG_DI_DIGI.matcher(cigarStr);
                if (cm.find() && !"D".equals(cm.group(1))) {
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
                    cigarStr = cm.replaceFirst(ilen + "D");
                    flag = true;
                }
            }

            //The following two clauses to capture sometimes mis-softly clipped reads by aligner
            Matcher mtch = ANY_NUMBER_M_NUMBER_S_END.matcher(cigarStr);
            if (mtch.find()) {
                captureMisSoftlyMS(mtch);
            } else if ((mtch = BEGIN_ANY_DIG_M_END.matcher(cigarStr)).find()) { // Make >=3 mismatches in the end as soft
                // clipping
                captureMisSoftly3Mismatches(mtch);
            }

            mtch = DIG_S_DIG_M.matcher(cigarStr);
            if (mtch.find()) {
                combineDigSDigM(mtch);
            } else if ((mtch = BEGIN_DIG_M.matcher(cigarStr)).find()) {// Make >=3 mismatches in the end as soft clipping
                combineBeginDigM(mtch);
            }
            return tuple(position, cigarStr);
        } catch (Exception exception) {
            printExceptionAndContinue(exception, "cigar", String.valueOf(position) + " " + originalCigar, region);
        }
        return tuple(position, cigarStr);
    }

    private void combineBeginDigM(Matcher mtch) {
        int mch = toInt(mtch.group(1));
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
        if (rn >= 3) {
            mch -= rn;
            cigarStr = mtch.replaceFirst(rn + "S" + mch + "M");
            position += rn;
        }
    }

    private void combineDigSDigM(Matcher mtch) {
        //length of matched sequence
        int mch = toInt(mtch.group(2));
        //length of soft-clipping
        int soft = toInt(mtch.group(1));
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
                cigarStr = mtch.replaceFirst(soft + "S" + mch + "M");
            } else {
                cigarStr = mtch.replaceFirst(mch + "M");
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

    private void captureMisSoftly3Mismatches(Matcher mtch) {
        String ov5 = mtch.group(1);
        int mch = toInt(mtch.group(2));
        int refoff = position + toInt(mtch.group(2));
        int rdoff = toInt(mtch.group(2));
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
        if (rn >= 3) {
            cigarStr = DIG_M_END.matcher(cigarStr).replaceFirst(mch + "M" + rn + "S");
        }
    }

    private void captureMisSoftlyMS(Matcher mtch) {
        //prefix of CIGAR string before last matched sequence
        String ov5 = mtch.group(1);
        //length of matched sequence
        int mch = toInt(mtch.group(2));
        //length of soft-clipping
        int soft = toInt(mtch.group(3));
        //offset of soft-clipped sequence in the reference string (position + length of matched)
        int refoff = position + toInt(mtch.group(2));
        //offset of soft-clipped sequence in the read
        int rdoff = toInt(mtch.group(2));
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

    private boolean combineToCloseToOne(Matcher cm, boolean flag) {
        //If CIGAR string contains any letter except D, matched sequence, deletion and possibly insertion
        String op = cm.group(5);
        int g2 = toInt(cm.group(2));
        int g3 = toInt(cm.group(3));
        int g4 = toInt(cm.group(4));

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
                if (cm.groupCount() > 5) {
                    String istr = cm.group(6);
                    if (istr != null) { //If digit and after it Insertion is present after deletion,
                        // add its length to ilen
                        ilen += toInt(istr.substring(0, istr.length() - 1));
                    }
                }
            }
            cm = DIG_I_DIG_M_DIG_DI_DIGI.matcher(cigarStr);
            //Replace I-M-D-I? complex with deletion and insertion
            cigarStr = cm.replaceFirst(dlen + "D" + ilen + "I");
            flag = true;
        }
        return flag;
    }

    private boolean combineToCloseToCorrect(Matcher cm, boolean flag) {
        int g1 = toInt(cm.group(1));
        int g2 = toInt(cm.group(2));
        int g3 = toInt(cm.group(3));
        if (g2 <= 15) {
            String op = cm.group(4);
            //length of both deletions and matched sequence
            int dlen = g1 + g2;
            //matched sequence length
            int ilen = g2;
            if (op.equals("I")) {
                ilen += g3;
            } else if (op.equals("D")) {
                dlen += g3;
                if (cm.groupCount() > 4) {
                    //insertion string
                    String istr = cm.group(5);
                    if (istr != null) { //If insertion is present after 2nd deletion, add its length to $ilen
                        ilen += toInt(istr.substring(0, istr.length() - 1));
                    }
                }
            }
            cigarStr = cm.replaceFirst(dlen + "D" + ilen + "I");
            flag = true;
        }
        return flag;
    }

    private boolean threeIndels(jregex.Matcher threeIndelsMatcher, boolean flag) {
        int tslen = toInt(threeIndelsMatcher.group(5)) + toInt(threeIndelsMatcher.group(8));
        if ("I".equals(threeIndelsMatcher.group(4))) {
            tslen += toInt(threeIndelsMatcher.group(3));
        }
        if ("I".equals(threeIndelsMatcher.group(7))) {
            tslen += toInt(threeIndelsMatcher.group(6));
        }
        if ("I".equals(threeIndelsMatcher.group(10))) {
            tslen += toInt(threeIndelsMatcher.group(9));
        }

        int dlen = toInt(threeIndelsMatcher.group(5)) + toInt(threeIndelsMatcher.group(8));
        if ("D".equals(threeIndelsMatcher.group(4))) {
            dlen += toInt(threeIndelsMatcher.group(3));
        }
        if ("D".equals(threeIndelsMatcher.group(7))) {
            dlen += toInt(threeIndelsMatcher.group(6));
        }
        if ("D".equals(threeIndelsMatcher.group(10))) {
            dlen += toInt(threeIndelsMatcher.group(9));
        }

        int mid = toInt(threeIndelsMatcher.group(5)) + toInt(threeIndelsMatcher.group(8));
        String ov5 = threeIndelsMatcher.group(1);
        int refoff = position + toInt(threeIndelsMatcher.group(2));
        int rdoff = toInt(threeIndelsMatcher.group(2));
        int RDOFF = toInt(threeIndelsMatcher.group(2));
        int rm = toInt(threeIndelsMatcher.group(11));

        if (!ov5.isEmpty()) { //If the complex is not at start of CIGAR string
            refoff += sum(globalFind(ALIGNED_LENGTH_MND, ov5)); // reference position // - 975 -
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
            newCigarStr += dlen + "D" + rm + "M";
            if (dlen == 0) {
                RDOFF = RDOFF + rm;
                newCigarStr = RDOFF + "M";
            }
        } else {
            newCigarStr += dlen + "D" + tslen + "I" + rm + "M";
        }
        if (mid <= 15) {
            cigarStr = DIGM_D_DI_DIGM_D_DI_DIGM_DI_DIGM.matcher(cigarStr).replaceFirst(newCigarStr);
            flag = true;
        }
        return flag;
    }

    private boolean threeDeletions(jregex.Matcher threeDeletionsMatcher, boolean flag) {
        //length of both matched sequences and insertion
        int tslen = toInt(threeDeletionsMatcher.group(4)) + toInt(threeDeletionsMatcher.group(6));
        //length of deletions and internal matched sequences
        int dlen = toInt(threeDeletionsMatcher.group(3)) + toInt(threeDeletionsMatcher.group(4)) +
                toInt(threeDeletionsMatcher.group(5)) + toInt(threeDeletionsMatcher.group(6)) +
                toInt(threeDeletionsMatcher.group(7));

        //length of internal matched sequences
        int mid = toInt(threeDeletionsMatcher.group(4)) + toInt(threeDeletionsMatcher.group(6));

        //prefix of CIGAR string before the M-D-M-I-M-D complex
        String ov5 = threeDeletionsMatcher.group(1);

        //offset of first deletion in the reference sequence
        int refoff = position + toInt(threeDeletionsMatcher.group(2));

        //offset of first deletion in the read
        int rdoff = toInt(threeDeletionsMatcher.group(2));

        //offset of first deletion in the read corrected by possibly matching bases
        int RDOFF = toInt(threeDeletionsMatcher.group(2));

        int rm = toInt(threeDeletionsMatcher.group(8));

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

    private boolean twoDeletionsInsertionToComplex(jregex.Matcher mm, boolean flag) {
        //length of both matched sequences and insertion
        int tslen = toInt(mm.group(4)) + toInt(mm.group(5)) + toInt(mm.group(6));

        //length of deletions and internal matched sequences
        int dlen = toInt(mm.group(3)) + toInt(mm.group(4)) + toInt(mm.group(6)) + toInt(mm.group(7));

        //length of internal matched sequences
        int mid = toInt(mm.group(4)) + toInt(mm.group(6));

        //prefix of CIGAR string before the M-D-M-I-M-D complex
        String ov5 = mm.group(1);

        //offset of first deletion in the reference sequence
        int refoff = position + toInt(mm.group(2));

        //offset of first deletion in the read
        int rdoff = toInt(mm.group(2));

        //offset of first deletion in the read corrected by possibly matching bases
        int RDOFF = toInt(mm.group(2));

        int rm = toInt(mm.group(8));

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
            //If length of internal matched sequences is no more than 10, replace M-D-M-I-M-D complex with M-D-I
            cigarStr = D_M_D_DD_M_D_I_D_M_D_DD_prim.matcher(cigarStr).replaceFirst(newCigarStr);
            flag = true;
        }
        return flag;
    }

    private boolean beginDigitMNumberIorDNumberM(jregex.Matcher mm) {
        int tmid = toInt(mm.group(1));
        int mlen = toInt(mm.group(4));
        int tn = 0;
                    /*
                    Sequence and insertion/deletion are replaced with soft-clipping up to
                    1st matching base in last matched sequence
                    For deletion position is adjusted by deletion length
                     */
        int tslen = tmid + (mm.group(3).equals("I") ? toInt(mm.group(2)) : 0);
        position += tmid + (mm.group(3).equals("D") ? toInt(mm.group(2)) : 0);
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
