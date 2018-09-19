package com.astrazeneca.vardict.modules;

import com.astrazeneca.vardict.collection.Tuple;
import com.astrazeneca.vardict.variations.Variation;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import jregex.Replacer;

import java.util.*;
import java.util.regex.Matcher;

import static com.astrazeneca.vardict.Utils.*;
import static com.astrazeneca.vardict.variations.VariationUtils.*;
import static com.astrazeneca.vardict.collection.Tuple.tuple;
import static com.astrazeneca.vardict.data.Patterns.*;

/**
 * Utility methods for parsing CIGAR.
 */
public class CigarUtils {

    /**
     * Modify the CIGAR for potential mis-alignment for indels at the end of reads to softclipping
     * and let VarDict's algorithm to figure out indels
     * @param indel length of insert/delete
     * @param ref Map of reference sequence (key - position, value - base)
     * @param oPosition Position of first matched base in sequence
     * @param oCigar Original CIGAR string
     * @param querySeq Base sequence
     * @param queryQual read base quality
     * @param lowqual low quality in soft-clipped sequence from configuration
     * @return Tuple of (adjusted position of first matched base, modified CIGAR string)
     */
    public static Tuple.Tuple2<Integer, String> modifyCigar(int indel,
                                                            Map<Integer, Character> ref,
                                                            final int oPosition,
                                                            final String oCigar,
                                                            final String querySeq,
                                                            final String queryQual,
                                                            final int lowqual) {

        int position = oPosition;
        String cigarStr = oCigar;

        //flag is set to true if CIGAR string is modified and should be looked at again
        boolean flag = true;
        // if CIGAR starts with deletion cut it off
        jregex.Matcher mm = BEGIN_NUMBER_D.matcher(cigarStr);
        if ( mm.find()) {
            position += toInt(mm.group(1));
            Replacer r = BEGIN_NUMBER_D.replacer("");
            cigarStr = r.replace(cigarStr);
        }
        mm = END_NUMBER_D.matcher(cigarStr);
        if ( mm.find()) {
            Replacer r = END_NUMBER_D.replacer("");
            cigarStr = r.replace(cigarStr);
        }
        // replace insertion at the beggining and end with soft clipping
        mm  =  BEGIN_NUMBER_I.matcher(cigarStr);
        if (mm.find()) {
            Replacer r = BEGIN_NUMBER_I.replacer(toInt(mm.group(1))+ "S");
            cigarStr = r.replace(cigarStr);
        }
        mm  =  END_NUMBER_I.matcher(cigarStr);
        if (mm.find()) {
            Replacer r = END_NUMBER_I.replacer(toInt(mm.group(1))+ "S");
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
                position += mm.group(3).equals("D") ? toInt(mm.group(2))  : 0;
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
                while ((tn < mlen) && isHasAndNotEquals(querySeq.charAt(tslen + tn), ref, position + tn)) {
                    tn++;
                }
                tslen += tn;
                mlen -= tn;
                position += tn;
                cigarStr = BEGIN_DIGIT_M_NUMBER_IorD_NUMBER_M_.matcher(cigarStr).replaceFirst(tslen + "S" + mlen + "M");
                flag = true;
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
            jregex.Matcher threeIndelsMatcher  = threeIndelsPattern.matcher(cigarStr);
            jregex.Matcher threeDeletionsMatcher = threeDeletionsPattern.matcher(cigarStr);
            int tslen;
            int dlen = 0;
            int mid = 0;
            String ov5;
            int rdoff;
            int refoff;
            int RDOFF;
            int rm;
            if (mm.find()) { //If CIGAR string contains matched sequence, deletion, short (<=9 bases) matched
                // sequence, insertion, short (<=9 bases) matched sequence and deletion

                //length of both matched sequences and insertion
                tslen = toInt(mm.group(4)) + toInt(mm.group(5)) + toInt(mm.group(6));

                //length of deletions and internal matched sequences
                dlen = toInt(mm.group(3)) + toInt(mm.group(4)) + toInt(mm.group(6)) + toInt(mm.group(7));

                //length of internal matched sequences
                mid = toInt(mm.group(4)) + toInt(mm.group(6));

                //prefix of CIGAR string before the M-D-M-I-M-D complex
                ov5 = mm.group(1);

                //offset of first deletion in the reference sequence
                refoff = position + toInt(mm.group(2));

                //offset of first deletion in the read
                rdoff = toInt(mm.group(2));

                //offset of first deletion in the read corrected by possibly matching bases
                RDOFF = toInt(mm.group(2));

                rm = toInt(mm.group(8));

                if (!ov5.isEmpty()) { //If the complex is not at start of CIGAR string
                    refoff += sum(globalFind(ALIGNED_LENGTH_MND, ov5)); // reference position
                    rdoff += sum(globalFind(SOFT_CLIPPED, ov5)); // read position
                }

                //number of bases after refoff/rdoff that match in reference and read
                int rn = 0;
                while (rdoff + rn < querySeq.length() && isHasAndEquals(querySeq.charAt(rdoff + rn), ref, refoff + rn)) {
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
            } else if (threeDeletionsMatcher.find()) {
                //length of both matched sequences and insertion
                tslen = toInt(threeDeletionsMatcher.group(4)) + toInt(threeDeletionsMatcher.group(6));
                //length of deletions and internal matched sequences
                dlen =  toInt(threeDeletionsMatcher.group(3)) + toInt(threeDeletionsMatcher.group(4)) +
                        toInt(threeDeletionsMatcher.group(5)) + toInt(threeDeletionsMatcher.group(6)) +
                        toInt(threeDeletionsMatcher.group(7));

                //length of internal matched sequences
                mid = toInt(threeDeletionsMatcher.group(4)) + toInt(threeDeletionsMatcher.group(6));

                //prefix of CIGAR string before the M-D-M-I-M-D complex
                ov5 = threeDeletionsMatcher.group(1);

                //offset of first deletion in the reference sequence
                refoff = position + toInt(threeDeletionsMatcher.group(2));

                //offset of first deletion in the read
                rdoff = toInt(threeDeletionsMatcher.group(2));

                //offset of first deletion in the read corrected by possibly matching bases
                RDOFF = toInt(threeDeletionsMatcher.group(2));

                rm = toInt(threeDeletionsMatcher.group(8));

                if (!ov5.isEmpty()) { //If the complex is not at start of CIGAR string
                    refoff += sum(globalFind(ALIGNED_LENGTH_MND, ov5)); // reference position // - 937 -
                    rdoff += sum(globalFind(SOFT_CLIPPED, ov5)); // read position
                    // MND changes from 27.04.2018
                }

                //number of bases after refoff/rdoff that match in reference and read
                int rn = 0;
                while (rdoff + rn < querySeq.length() && isHasAndEquals(querySeq.charAt(rdoff + rn), ref, refoff + rn)) {
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
            } else if (threeIndelsMatcher.find()) {
                tslen = toInt(threeIndelsMatcher.group(5)) + toInt(threeIndelsMatcher.group(8));
                if ("I".equals(threeIndelsMatcher.group(4))) {
                    tslen += toInt(threeIndelsMatcher.group(3));
                }
                if ("I".equals(threeIndelsMatcher.group(7))) {
                    tslen += toInt(threeIndelsMatcher.group(6));
                }
                if ("I".equals(threeIndelsMatcher.group(10))) {
                    tslen += toInt(threeIndelsMatcher.group(9));
                }

                dlen = toInt(threeIndelsMatcher.group(5)) + toInt(threeIndelsMatcher.group(8));
                if ("D".equals(threeIndelsMatcher.group(4))) {
                    dlen += toInt(threeIndelsMatcher.group(3));
                }
                if ("D".equals(threeIndelsMatcher.group(7))) {
                    dlen += toInt(threeIndelsMatcher.group(6));
                }
                if ("D".equals(threeIndelsMatcher.group(10))) {
                    dlen += toInt(threeIndelsMatcher.group(9));
                }

                mid = toInt(threeIndelsMatcher.group(5)) + toInt(threeIndelsMatcher.group(8));
                ov5 = threeIndelsMatcher.group(1);
                refoff = position + toInt(threeIndelsMatcher.group(2));
                rdoff = toInt(threeIndelsMatcher.group(2));
                RDOFF = toInt(threeIndelsMatcher.group(2));
                rm = toInt(threeIndelsMatcher.group(11));

                if (!ov5.isEmpty()) { //If the complex is not at start of CIGAR string
                    refoff += sum(globalFind(ALIGNED_LENGTH_MND, ov5)); // reference position // - 975 -
                    rdoff += sum(globalFind(SOFT_CLIPPED, ov5)); // read position
                }

                int rn = 0;
                while (rdoff + rn < querySeq.length() && isHasAndEquals(querySeq.charAt(rdoff + rn), ref, refoff + rn)) {
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
                    cigarStr = DIGM_D_DI_DIGM_D_DI_DIGM_DI_DIGM.matcher(cigarStr).replaceFirst(newCigarStr);
                    flag = true; // - 995 -
                }
            }
            // Combine two close deletions (<15bp) into one correct from <10 to <15
            Matcher cm = DIG_D_DIG_M_DIG_DI_DIGI.matcher(cigarStr);
            if (cm.find()) { //If CIGAR string contains deletion, short (<= 9 bases) matched sequence, deletion and
                // possibly insertion
                int g1 = toInt(cm.group(1));
                int g2 = toInt(cm.group(2));
                int g3 = toInt(cm.group(3));
                if (g2 <= 15) {
                    String op = cm.group(4);
                    //length of both deletions and matched sequence
                    dlen = g1 + g2;
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
            }

            // Combine two close indels (<15bp) into one
            cm = NOTDIG_DIG_I_DIG_M_DIG_DI_DIGI.matcher(cigarStr);
            if (cm.find() && !"D".equals(cm.group(1))) {
                //If CIGAR string contains any letter except D, matched sequence, deletion and possibly insertion
                String op = cm.group(5);
                int g2 = toInt(cm.group(2));
                int g3 = toInt(cm.group(3));
                int g4 = toInt(cm.group(4));

                if (g3 <= 15) {
                    //length of matched sequence and deletion
                    dlen = g3;
                    //length of first insertion and matched sequence
                    int ilen = g2 + g3;
                    if (op.equals("I")) {
                        ilen += g4;
                    } else if (op.equals("D")) {
                        dlen += g4;
                        //last insertion string
                        if (cm.groupCount()>5) {
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
            }
            cm = DIG_D_DIG_D.matcher(cigarStr);
            if (cm.find()) {
                dlen = toInt(cm.group(1)) + toInt(cm.group(2));
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
            while (rn < soft && isHasAndEquals(ref, refoff + rn, querySeq, rdoff + rn)
                    && queryQual.charAt(rdoff + rn) - 33 > lowqual) {
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
                while (rn + 1 < soft && isHasAndEquals(ref, refoff + rn + 1, querySeq, rdoff + rn + 1) &&
                        queryQual.charAt(rdoff + rn + 1) - 33 > lowqual) {
                    rn++;
                    RN.add(ref.get(refoff + rn + 1)); // to my mind in perl value increasing here ? not adding el
                }
                int rn_nt = RN.size(); // don't adjust if homopolymer
                if ((rn > 4 && rn_nt > 1) || (isHasAndEquals(ref, refoff, querySeq, rdoff))) { //If more than 3 bases
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
                        if (!ref.containsKey(refoff - rrn - 1)) {
                            break;
                        }
                        if (rrn < rdoff && isHasAndNotEquals(ref, refoff - rrn - 1, querySeq, rdoff - rrn - 1)) {
                            rn = rrn + 1;
                            rmch = 0;
                        } else if (rrn < rdoff && isHasAndEquals(ref, refoff - rrn - 1, querySeq, rdoff - rrn - 1)) {
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
        } else if ((mtch = BEGIN_ANY_DIG_M_END.matcher(cigarStr)).find()) { // Make >=3 mismatches in the end as soft
            // clipping
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
                if (!ref.containsKey(refoff - rrn - 1)) {
                    break;
                }
                if (rrn < rdoff && isHasAndNotEquals(ref, refoff - rrn - 1, querySeq, rdoff - rrn - 1)) {
                    rn = rrn + 1;
                    rmch = 0;
                } else if (rrn < rdoff && isHasAndEquals(ref, refoff - rrn - 1, querySeq, rdoff - rrn - 1)) {
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

        mtch = DIG_S_DIG_M.matcher(cigarStr);
        if (mtch.find()) {
            //length of matched sequence
            int mch = toInt(mtch.group(2));
            //length of soft-clipping
            int soft = toInt(mtch.group(1));
            //number of bases before matched sequence that match in reference and read sequences
            int rn = 0;
            Set<Character> RN = new HashSet<>();
            while (rn < soft && isHasAndEquals(ref, position - rn - 1, querySeq, soft - rn - 1)
                    && queryQual.charAt(soft - rn - 1) - 33 > lowqual) {
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
                while (rn + 1 < soft && isHasAndEquals(ref, position - rn - 2, querySeq, soft - rn - 2) &&
                        queryQual.charAt(soft - rn - 2) - 33 > lowqual) {
                    rn++;
                    RN.add(ref.get(position - rn - 2));
                }
                int rn_nt = RN.size();  // Don't adjust for homopolymers
                if ((rn > 4 && rn_nt > 1) || isHasAndEquals(ref, position - 1, querySeq, soft - 1)) {
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
                        if (!ref.containsKey(position + rrn)) {
                            break;
                        }
                        if (isHasAndNotEquals(ref, position + rrn, querySeq, soft + rrn)) {
                            rn = rrn + 1;
                            rmch = 0;
                        } else if (isHasAndEquals(ref, position + rrn, querySeq, soft + rrn)) {
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
        } else if ((mtch = BEGIN_DIG_M.matcher(cigarStr)).find()) {// Make >=3 mismatches in the end as soft clipping
            int mch = toInt(mtch.group(1));
            int rn = 0;
            int rrn = 0;
            int rmch = 0;
            while (rrn < mch && rn < mch) {
                if (!ref.containsKey(position + rrn)) {
                    break;
                }
                if (isHasAndNotEquals(ref, position + rrn, querySeq, rrn)) {
                    rn = rrn + 1;
                    rmch = 0;
                } else if (isHasAndEquals(ref, position + rrn, querySeq, rrn)) {
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
        return tuple(position, cigarStr);
    }

    /**
     * Decrement variant counters.
     * @param variation reference variant to decrement    $vref
     * @param direction true if read has negative strand flag    $dir
     * @param readPosition read position    $rp
     * @param baseQuality base quality   $q
     * @param mappingBaseQuality bases's mapping quality    $Q
     * @param numberOfMismatches number of mismatches   $nm
     * @param goodq configuration based good quality
     */
    static void subCnt(Variation variation,
                              boolean direction,
                              int readPosition,
                              double baseQuality,
                              int mappingBaseQuality,
                              int numberOfMismatches,
                              double goodq) {
        // rp = read position, q = quality
        variation.cnt--;
        variation.decDir(direction);
        variation.pmean -= readPosition;
        variation.qmean -= baseQuality;
        variation.Qmean -= mappingBaseQuality;
        variation.nm -= numberOfMismatches;
        if (baseQuality >= goodq) {
            variation.hicnt--;
        } else {
            variation.locnt--;
        }
    }

    /**
     * Increment variant counters.
     * @param variation variant to increment    $vref
     * @param direction if read has negative strand flag    $dir
     * @param readPosition read position    $rp
     * @param baseQuality base quality   $q
     * @param mappingBaseQuality bases's mapping quality    $Q
     * @param numberOfMismatches number of mismatches   $nm
     * @param goodq configuration based good quality
     */
    static void addCnt(Variation variation,
                              boolean direction,
                              int readPosition,
                              double baseQuality,
                              int mappingBaseQuality,
                              int numberOfMismatches,
                              double goodq) {
        variation.cnt++;
        variation.incDir(direction);
        variation.pmean += readPosition;
        variation.qmean += baseQuality;
        variation.Qmean += mappingBaseQuality;
        variation.nm += numberOfMismatches;
        if (baseQuality >= goodq) {
            variation.hicnt++;
        } else {
            variation.locnt++;
        }
    }

    /**
     * Increase count for variation
     * @param counters map to increase count for variation
     * @param idx index in <code>counters</code> map
     * @param s string with variations
     */
    static void increment(Map<Integer, Map<String, Integer>> counters,
                                  int idx,
                                  String s) {
        Map<String, Integer> map = counters.get(idx);
        if (map == null) {
            map = new HashMap<>();
            counters.put(idx, map);
        }
        incCnt(map, s, 1);
    }

    static boolean isBEGIN_ATGC_AMP_ATGCs_END(String s) {
        if (s.length() > 2) {
            char ch1 = s.charAt(0);
            char ch2 = s.charAt(1);
            if (ch2 == '&' && isATGC(ch1)) {
                for (int i = 2; i < s.length(); i++) {
                    if (!isATGC(s.charAt(i))) {
                        return false;
                    }
                }
                return true;
            }
        }
        return false;
    }

    static boolean isATGC(char ch) {
        switch (ch) {
            case 'A':
            case 'T':
            case 'G':
            case 'C':
                return true;

            default:
                return false;
        }
    }

    /**
     * Find closest mismatches to combine with indels.
     * @param referencePosition position for reference  $refp
     * @param readPosition position for base sequence   $readp
     * @param cigarLength length of CIGAR   $mlen
     * @param querySequence base sequence   $rdseq
     * @param queryQuality base quality sequence    $qstr
     * @param reference reference bases by position     $ref
     * @param refCoverage reference coverage    $cov
     * @param vext configuration based Extension of bp
     * @param goodq configuration based good quality
     * @return Tuple of (offset, querySequence's substring, queryQuality substring, number of mismatches)
     */
    static Tuple.Tuple4<Integer, String, String, Integer> findOffset(int referencePosition,
                                                                            int readPosition,
                                                                            int cigarLength,
                                                                            String querySequence,
                                                                            String queryQuality,
                                                                            Map<Integer, Character> reference,
                                                                            Map<Integer, Integer> refCoverage,
                                                                            int vext,
                                                                            double goodq) {
        int offset = 0;
        String ss = "";
        String q = "";
        int tnm = 0;
        int vsn = 0;
        for (int vi = 0; vsn <= vext && vi < cigarLength; vi++) {
            if (querySequence.charAt(readPosition + vi) == 'N') {
                break;
            }
            if (queryQuality.charAt(readPosition + vi) - 33 < goodq) {
                break;
            }
            Character refCh = reference.get(referencePosition + vi);
            if (refCh != null) {
                char ch = querySequence.charAt(readPosition + vi);
                if (isNotEquals(ch, refCh)) {
                    offset = vi + 1;
                    tnm++;
                    vsn = 0;
                } else {
                    vsn++;
                }
            }
        }
        if (offset > 0) {
            ss = substr(querySequence, readPosition, offset);
            q = substr(queryQuality, readPosition, offset);
            for (int osi = 0; osi < offset; osi++) {
                incCnt(refCoverage, referencePosition + osi, 1);
            }
        }
        return tuple(offset, ss, q, tnm);
    }

    /**
     * Attempts to clean up a CIGAR string so that edge cases are avoided in the rest of the code.
     * Specifically this method will remove leading or trailing hard-clips, and then convert
     * leading or trailing insertions into soft-clips.
     */
    static void cleanupCigar(final SAMRecord rec) {
        if (rec.getCigar() != null) {
            final List<CigarElement> elems = new ArrayList<>(rec.getCigar().getCigarElements());

            // First leading elements
            {
                final ListIterator<CigarElement> iterator = elems.listIterator();
                boolean noMatchesYet = true;
                while (iterator.hasNext() && noMatchesYet) {
                    final CigarElement elem = iterator.next();
                    if (elem.getOperator() == CigarOperator.INSERTION) {
                        final CigarElement replacement = new CigarElement(elem.getLength(), CigarOperator.SOFT_CLIP);
                        iterator.set(replacement);
                    }
                    else if (elem.getOperator() == CigarOperator.HARD_CLIP) {
                        iterator.remove();
                    }
                    else if (elem.getOperator().consumesReadBases() && elem.getOperator().consumesReferenceBases()) {
                        noMatchesYet = false;
                    }
                }
            }
            // Then trailing elements
            {
                final ListIterator<CigarElement> iterator = elems.listIterator(elems.size());
                boolean noMatchesYet = true;
                while (iterator.hasPrevious() && noMatchesYet) {
                    final CigarElement elem = iterator.previous();
                    if (elem.getOperator() == CigarOperator.INSERTION) {
                        final CigarElement replacement = new CigarElement(elem.getLength(), CigarOperator.SOFT_CLIP);
                        iterator.set(replacement);
                    }
                    else if (elem.getOperator() == CigarOperator.HARD_CLIP) {
                        iterator.remove();
                    }
                    else if (elem.getOperator().consumesReadBases() && elem.getOperator().consumesReferenceBases()) {
                        noMatchesYet = false;
                    }
                }
            }
            // And lastly replace the cigar
            rec.setCigar(new Cigar(elems));
        }
    }

    static CigarOperator getCigarOperator(Cigar cigar,
                                                 int ci) {
        CigarOperator operator =  cigar.getCigarElement(ci).getOperator();
        // Treat insertions at the edge as soft-clipping
        if ((ci == 0 || ci == cigar.numCigarElements() - 1) && operator == CigarOperator.I) {
            operator = CigarOperator.S;
        }
        return operator;
    }

    static int getAlignedLength(Cigar readCigar) {
        int length = 0;
        for (CigarElement element : readCigar.getCigarElements()) {
            if (element.getOperator() == CigarOperator.M || element.getOperator() == CigarOperator.D) {
                length += element.getLength();
            }
        }
        return length;
    }

    static int getSoftClippedLength(Cigar readCigar) {
        int length = 0;
        for (CigarElement element : readCigar.getCigarElements()) {
            CigarOperator operator = element.getOperator();
            if (operator == CigarOperator.M || operator == CigarOperator.I || operator == CigarOperator.S) {
                length += element.getLength();
            }
        }
        return length;
    }

    static int getMatchInsertionLength(Cigar readCigar) {
        int length = 0;
        for (CigarElement element : readCigar.getCigarElements()) {
            if (element.getOperator() == CigarOperator.M || element.getOperator() == CigarOperator.I) {
                length += element.getLength();
            }
        }
        return length;
    }

    static int getInsertionDeletionLength(Cigar readCigar) {
        int length = 0;
        for (CigarElement element : readCigar.getCigarElements()) {
            if (element.getOperator() == CigarOperator.I || element.getOperator() == CigarOperator.D) {
                length += element.getLength();
            }
        }
        return length;
    }
}
