package com.astrazeneca.vardict.pipeline.modules;

import com.astrazeneca.utils.Tuple;
import com.astrazeneca.vardict.variations.SoftClip;
import com.astrazeneca.vardict.variations.Variation;
import com.astrazeneca.vardict.variations.VariationUtils;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import jregex.Replacer;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static com.astrazeneca.GlobalReadOnlyScope.instance;
import static com.astrazeneca.utils.Tuple.tuple;
import static com.astrazeneca.utils.Utils.*;
import static com.astrazeneca.utils.Utils.substr;
import static com.astrazeneca.vardict.VarDictLauncher.*;
import static com.astrazeneca.vardict.variations.VariationUtils.isHasAndEquals;
import static com.astrazeneca.vardict.variations.VariationUtils.isHasAndNotEquals;
import static com.astrazeneca.vardict.variations.VariationUtils.isNotEquals;

/**
 * Utility methods for parsing CIGAR.
 */
public class CigarUtils {
    /**
     * Regexp finds number followed by S followed by number and I or D at the start of string
     */
    private static final jregex.Pattern BEGIN_NUMBER_S_NUMBER_IorD = new jregex.Pattern("^(\\d+)S(\\d+)([ID])");
    /**
     * Regexp finds number followed by I or D followed by number followed by S followed by end of string
     */
    private static final jregex.Pattern NUMBER_IorD_NUMBER_S_END = new jregex.Pattern("(\\d+)([ID])(\\d+)S$");
    /**
     * Regexp finds number followed by S followed by number followed by M followed by number followed by I or D at the start of string
     */
    private static final jregex.Pattern BEGIN_NUMBER_S_NUMBER_M_NUMBER_IorD = new jregex.Pattern("^(\\d+)S(\\d+)M(\\d+)([ID])");
    /**
     * Regexp finds number followed by I or D followed by number followed by M followed by number followed by S followed by end of string
     */
    private static final jregex.Pattern NUMBER_IorD_NUMBER_M_NUMBER_S_END = new jregex.Pattern("(\\d+)([ID])(\\d+)M(\\d+)S$");
    /**
     * Regexp finds digit-M-number-(I or D)-number-M
     */
    private static final jregex.Pattern BEGIN_DIGIT_M_NUMBER_IorD_NUMBER_M = new jregex.Pattern("^(\\d)M(\\d+)([ID])(\\d+)M");
    private static final Pattern BEGIN_DIGIT_M_NUMBER_IorD_NUMBER_M_ = Pattern.compile("^\\dM\\d+[ID]\\d+M");
    /**
     * Regexp replaces number-(I or D)-number-M
     */
    private static final jregex.Pattern NUMBER_IorD_DIGIT_M_END = new jregex.Pattern("(\\d+)([ID])(\\d)M$");

    /**
     * Regexp finds number-M-number-D-digit-M-number-I-digit-M-number-D.
     * ^(.*?) captures everything before the M-D-M-I-M-D complex
     */
    private static final jregex.Pattern D_M_D_DD_M_D_I_D_M_D_DD = new jregex.Pattern("^(.*?)(\\d+)M(\\d+)D(\\d)M(\\d+)I(\\d)M(\\d+)D(\\d+)M");
    private static final Pattern D_M_D_DD_M_D_I_D_M_D_DD_prim = Pattern.compile("(\\d+)M(\\d+)D(\\d)M(\\d+)I(\\d)M(\\d+)D(\\d+)M");
    /**
     * Regexp finds number-D-digit-M-number-D and optionally number-I
     */
    private static final Pattern DIG_D_DIG_M_DIG_DI_DIGI = Pattern.compile("(\\d+)D(\\d)M(\\d+)([DI])(\\d+I)?");
    /**
     * Regexp findx number-I-digit-M-number and optionally number-I
     */
    private static final Pattern DIG_I_dig_M_DIG_DI_DIGI = Pattern.compile("(\\d+)I(\\d)M(\\d+)([DI])(\\d+I)?");

    /**
     * Regexp finds number-S-number-M at start of CIGAR string
     */
    private static final Pattern D_S_D_M = Pattern.compile("^(\\d+)S(\\d+)M");
    private static final Pattern D_M_D_S_END = Pattern.compile("\\d+M\\d+S$");
    /**
     * Regexp finds number-M-number-S at end of CIGAR string
     * ^(.*?) captures everything before M-S complex
     */
    private static final Pattern ANY_NUMBER_M_NUMBER_S_END = Pattern.compile("^(.*?)(\\d+)M(\\d+)S$");
    /**
     * Regexp finds number followed by D at the start of string
     */
    private static final jregex.Pattern BEGIN_NUMBER_D = new jregex.Pattern("^(\\d+)D");
    /**
     * Regexp finds number followed by I at the start of string
     */
    private static final jregex.Pattern BEGIN_NUMBER_I = new jregex.Pattern("^(\\d+)I");
    /**
     * Regexp finds number followed by I at the end of string
     */
    private static final jregex.Pattern END_NUMBER_I = new jregex.Pattern("(\\d+)I$");

    private static final jregex.Pattern SOFT_CLIPPED = new jregex.Pattern("(\\d+)[MIS]");
    /**
     * regexp finds numbers followed by M (matched) or D (deleted) in CIGAR string
     */
    private static final jregex.Pattern ALIGNED_LENGTH = new jregex.Pattern("(\\d+)[MD]");

    /**
     * Modify the CIGAR for potential mis-alignment for indels at the end of reads to softclipping and let VarDict's algorithm to figure out indels
     * @param indel length of insert/delete
     * @param reference map of reference sequence (key - position, value - base)
     * @param oPosition position of first matched base in sequence
     * @param oCigar original CIGAR string
     * @param querySeq base sequence
     * @param queryQual base quality sequence
     * @return Tuple of (adjusted position of first matched base, modified CIGAR string)
     */
    //perl version: 685
    // TODO: create new class CigarModified and devide this method for small private methods
    public static Tuple.Tuple2<Integer, String> modifyCigar(int indel,
                                                     Map<Integer, Character> reference,
                                                     int oPosition,
                                                     String oCigar,
                                                     String querySeq,
                                                     String queryQual) {
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
        // replace insertion at the beginning and end with soft clipping
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
                position += mm.group(3).equals("D") ? 2 : 0;
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
            mm = BEGIN_DIGIT_M_NUMBER_IorD_NUMBER_M.matcher(cigarStr); //If CIGAR starts with 1-9 bases long matched sequence, insertion or deletion and matched sequence
            if (mm.find()) {
                int tmid = toInt(mm.group(1));
                int mlen = toInt(mm.group(4));
                if (tmid <= 8) {
                    /*
                    If first matched sequence length is no more than 8, this sequence and insertion/deletion are
                    replaced with soft-clipping up to 1st matching base in last matched sequence
                    For deletion position is adjusted by deletion length
                     */
                    int tslen = tmid + (mm.group(3).equals("I") ? toInt(mm.group(2)) : 0);
                    position += tmid + (mm.group(3).equals("D") ? toInt(mm.group(2)) : 0);
                    int n0 = 0;
                    while (n0 < mlen
                            && isHasAndNotEquals(querySeq.charAt(tslen + n0), reference, position + n0)) {
                        n0++;
                    }
                    tslen += n0;
                    mlen -= n0;
                    position += n0;
                    //Regexp replaces digit-M-number-(I or D)-number-M with tslen + S + mlen + M
                    cigarStr = BEGIN_DIGIT_M_NUMBER_IorD_NUMBER_M_.matcher(cigarStr).replaceFirst(tslen + "S" + mlen + "M");
                    flag = true;
                }
            }
            mm = NUMBER_IorD_DIGIT_M_END.matcher(cigarStr);
            if (mm.find()) { //If CIGAR ends with insertion or deletion and 1-9 bases long matched sequence
                int tmid = toInt(mm.group(3));
                if (tmid <= 8) { //If matched sequence length is no more than 8, insertion/deletion and matched sequence are replaced with soft-clipping
                    String tslen = tmid + (mm.group(2).equals("I") ? toInt(mm.group(1)) : 0) + "S";
                    Replacer r = NUMBER_IorD_DIGIT_M_END.replacer(tslen);
                    cigarStr = r.replace(cigarStr);
                    flag = true;
                }
            }

            //perl version: 750
            // Combine two deletions and insertion into one complex if they are close
            mm = D_M_D_DD_M_D_I_D_M_D_DD.matcher(cigarStr);
            if (mm.find()) { //If CIGAR string contains matched sequence, deletion, short (<=9 bases) matched sequence, insertion, short (<=9 bases) matched sequence and deletion
                //length of internal matched sequences
                int mid = toInt(mm.group(4)) + toInt(mm.group(6));
                if (mid <= 10) {
                    //length of both matched sequences and insertion
                    int tslen = mid + toInt(mm.group(5));
                    //length of deletions and internal matched sequences
                    int dlen = toInt(mm.group(3)) + mid + toInt(mm.group(7));
                    //prefix of CIGAR string before the M-D-M-I-M-D complex
                    String ov5 = mm.group(1);
                    //offset of first deletion in the read
                    int rdoff = toInt(mm.group(2));
                    //offset of first deletion in the reference sequence
                    int refoff = position + rdoff;
                    //offset of first deletion in the read corrected by possibly matching bases
                    int RDOFF = rdoff;
                    int rm = toInt(mm.group(8));
                    if (!ov5.isEmpty()) { //If the complex is not at start of CIGAR string
                        rdoff += sum(globalFind(SOFT_CLIPPED, ov5)); // read position
                        refoff += sum(globalFind(ALIGNED_LENGTH, ov5)); // reference position
                    }
                    //number of bases after refoff/rdoff that match in reference and read
                    int rn = 0;
                    while (rdoff + rn < querySeq.length()
                            && isHasAndEquals(querySeq.charAt(rdoff + rn), reference, refoff + rn)) {
                        rn++;
                    }
                    RDOFF += rn;
                    dlen -= rn;
                    tslen -= rn;
                    String newCigarStr = RDOFF + "M";
                    if ( tslen <= 0 ) {
                        dlen -= tslen;
                        rm += tslen;
                        newCigarStr += dlen + "D" + rm + "M";
                    } else {
                        newCigarStr += dlen + "D" + tslen + "I" + rm + "M";
                    }
                    //If length of internal matched sequences is no more than 10, replace M-D-M-I-M-D complex with M-D-I
                    cigarStr = D_M_D_DD_M_D_I_D_M_D_DD_prim.matcher(cigarStr).replaceFirst(newCigarStr);
                    flag = true;
                }
            }
            //perl version: 786
            // Combine two close deletions (<10bp) into one
            Matcher cm = DIG_D_DIG_M_DIG_DI_DIGI.matcher(cigarStr);
            if (cm.find()) { //If CIGAR string contains deletion, short (<= 9 bases) matched sequence, deletion and possibly insertion
                int g2 = toInt(cm.group(2));
                int g3 = toInt(cm.group(3));
                String op = cm.group(4);

                //length of both deletions and matched sequence
                int dlen = toInt(cm.group(1)) + g2;
                //matched sequence length
                int ilen = g2;
                if (op.equals("I")) {
                    ilen += g3;
                } else { // op == "D"
                    dlen += g3;
                    //insertion string
                    String istr = cm.group(5);
                    if (istr != null) { //If insertion is present after 2nd deletion, add its length to $ilen
                        ilen += toInt(istr.substring(0, istr.length() - 1));
                    }
                }
                //Replace D-M-D-I? complex with deletion and insertion
                cigarStr = cm.replaceFirst(dlen + "D" + ilen + "I");
                flag = true;
            }

            //perl version: 798
            // Combine two close indels (<10bp) into one
            cm = DIG_I_dig_M_DIG_DI_DIGI.matcher(cigarStr);
            if (cm.find()) { //If CIGAR string contains insertion, short (<=9 bases) matched sequence, deletion and possibly insertion
                String op = cm.group(4);
                int g2 = toInt(cm.group(2));
                int g3 = toInt(cm.group(3));

                //length of matched sequence and deletion
                int dlen = g2;
                //length of first insertion and matched sequence
                int ilen = toInt(cm.group(1)) + g2;
                if (op.equals("I")) {
                    ilen += g3;
                } else { // op == "D"
                    dlen += g3;
                    //last insertion string
                    String istr = cm.group(5);
                    if (istr != null) { //If insertion is present after deletion, add its length to ilen
                        ilen += toInt(istr.substring(0, istr.length() - 1));
                    }
                }
                //Replace I-M-D-I? complex with deletion and insertion
                cigarStr = cm.replaceFirst(dlen + "D" + ilen + "I");
                flag = true;
            }
        }

        //perl version: 812
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
            int refoff = position + mch;
            //offset of soft-clipped sequence in the read
            int rdoff = mch;
            if (!ov5.isEmpty()) { //If prefix is present
                //Add all matched, insertion and soft-clipped lengths to read position
                rdoff += sum(globalFind(SOFT_CLIPPED, ov5)); // read position
                //Add all matched and deletion lengths to reference position
                refoff += sum(globalFind(ALIGNED_LENGTH, ov5)); // reference position
            }
            //number of bases after refoff/rdoff that match in reference and read sequences
            int rn = 0;
            Set<Character> RN = new HashSet<>();
            while( rn + 1 < soft
                    && isHasAndEquals(reference, refoff + rn + 1, querySeq, rdoff + rn + 1)
                    && queryQual.charAt(rdoff + rn + 1) - 33 > instance().conf.lowqual) {
                rn++;
                RN.add(reference.get(refoff + rn + 1));
            }
            int rn_nt = RN.size(); // don't adjust if homopolymer
            if ( (rn > 3 && rn_nt > 1) || (isHasAndEquals(reference, refoff, querySeq, rdoff))) { //If more than 3 bases match after refoff/rdoff or base at refoff/rdoff match
                mch += rn + 1;
                soft -= rn + 1;
                if (soft > 0) {
                    cigarStr = D_M_D_S_END.matcher(cigarStr).replaceFirst(mch + "M" + soft + "S");
                } else {
                    cigarStr = D_M_D_S_END.matcher(cigarStr).replaceFirst(mch + "M");
                }
            }
            else if (rn == 0) {
                while (rn < mch && isHasAndNotEquals(reference, refoff - rn - 1, querySeq, rdoff - rn - 1)) {
                    rn++;
                }
                if (rn > 0 && rn < mch) {
                    soft += rn;
                    mch -= rn;
                    cigarStr = D_M_D_S_END.matcher(cigarStr).replaceFirst(mch + "M" + soft + "S");
                }
            }
        }

        mtch = D_S_D_M.matcher(cigarStr);
        if (mtch.find()) {
            //length of matched sequence
            int mch = toInt(mtch.group(2));
            //length of soft-clipping
            int soft = toInt(mtch.group(1));
            //number of bases before matched sequence that match in reference and read sequences
            int rn = 0;
            Set<Character> RN = new HashSet<>();
            while (rn + 1 < soft && isHasAndEquals(reference, position - rn - 2, querySeq, soft - rn - 2)
                    && queryQual.charAt(soft - rn - 2) - 33 > instance().conf.lowqual) {
                rn++;
                RN.add(reference.get(position - rn - 2));
            }
            int rn_nt = RN.size();
            if ((rn > 3 && rn_nt > 1) || (isHasAndEquals(reference, position - 1, querySeq, soft -1 ))) {//If more than 3 bases match before matched sequence or last base of clipped sequence matches
                //Replace the S-M complex with either match or match-soft clip
                mch += rn + 1;
                soft -= rn + 1;
                if (soft > 0) {
                    cigarStr = mtch.replaceFirst(soft + "S" + mch + "M");
                } else {
                    cigarStr = mtch.replaceFirst(mch + "M");
                }
                position -= rn + 1;
                rn++;
            }
            if (rn == 0) {
                while (rn < mch && isHasAndNotEquals(reference, position + rn, querySeq, soft + rn)) {
                    rn++;
                }
                if (rn > 0 && rn < mch) {
                    soft += rn;
                    mch -= rn;
                    cigarStr = mtch.replaceFirst(soft + "S" + mch + "M");
                    position += rn;
                }
            }
        }
        return tuple(position, cigarStr);
    }

    /**
     * Decrement variant counters.
     * @param variation reference variant to decrement
     * @param direction if read has negative strand flag
     * @param readPosition read position
     * @param baseQuality base quality
     * @param mappingBaseQuality bases's mapping quality
     * @param numberOfMismatches number of mismatches
     */
    //perl version: 2141
    public static void subCnt(Variation variation,
                       boolean direction,
                       int readPosition,
                       double baseQuality,
                       int mappingBaseQuality,
                       int numberOfMismatches) {
        // referenceVariant dir read_position quality
        variation.varsCount--;
        variation.decCountForDir(direction);
        variation.meanPosition -= readPosition;
        variation.meanQuality -= baseQuality;
        variation.meanMappingQuality -= mappingBaseQuality;
        variation.numberOfMismatches -= numberOfMismatches;
        if (baseQuality >= instance().conf.goodq) {
            variation.highQualityReadsCount--;
        } else {
            variation.lowQualityReadsCount--;
        }
    }

    /**
     * Increment variant counters.
     * @param variation variant to increment
     * @param direction if read has negative strand flag
     * @param readPosition read position
     * @param baseQuality base quality
     * @param mappingBaseQuality bases's mapping quality
     * @param numberOfMismatches number of mismatches
     */
    //perl version: 550
    public static void addCnt(Variation variation,
                       boolean direction,
                       int readPosition,
                       double baseQuality,
                       int mappingBaseQuality,
                       int numberOfMismatches) {
        variation.varsCount++;
        variation.incCountForDir(direction);
        variation.meanPosition += readPosition;
        variation.meanQuality += baseQuality;
        variation.meanMappingQuality += mappingBaseQuality;
        variation.numberOfMismatches += numberOfMismatches;
        if (baseQuality >= instance().conf.goodq) {
            variation.highQualityReadsCount++;
        } else {
            variation.lowQualityReadsCount++;
        }
    }

    public static CigarOperator getCigarOperator(Cigar cigar, int ci) {
        CigarOperator operator =  cigar.getCigarElement(ci).getOperator();
        // Treat insertions at the edge as soft-clipping
        if ((ci == 0 || ci == cigar.numCigarElements() - 1) && operator == CigarOperator.I) {
            operator = CigarOperator.S;
        }
        return operator;
    }

    public static int getSoftClippedLength(Cigar readCigar) {
        int length = 0;
        for (CigarElement element : readCigar.getCigarElements()) {
            CigarOperator operator = element.getOperator();
            if (operator == CigarOperator.M || operator == CigarOperator.I || operator == CigarOperator.S) {
                length += element.getLength();
            }
        }
        return length;
    }

    public static int getMatchInsertionLength(Cigar readCigar) {
        int length = 0;
        for (CigarElement element : readCigar.getCigarElements()) {
            if (element.getOperator() == CigarOperator.M || element.getOperator() == CigarOperator.I) {
                length += element.getLength();
            }
        }
        return length;
    }

    public static int getInsertionDeletionLength(Cigar readCigar) {
        int length = 0;
        for (CigarElement element : readCigar.getCigarElements()) {
            if (element.getOperator() == CigarOperator.I || element.getOperator() == CigarOperator.D) {
                length += element.getLength();
            }
        }
        return length;
    }

    // TODO: regexp?
    public static boolean isBEGIN_ATGC_AMP_ATGCs_END(String s) {
        if (s.length() > 2) {
            char ch1 = s.charAt(0);
            char ch2 = s.charAt(1);
            if (ch2 == '&' && isATGC(ch1)) {
                for (int i = 2; i < s.length(); i++) {
                    if(!isATGC(s.charAt(i))) {
                        return false;
                    }
                }
                return true;
            }
        }
        return false;
    }

    private static boolean isATGC(char ch) {
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
     * Get {@link Variation} from {@link SoftClip#seq} field
     * @param softClip variation with soft clipped reads
     * @param idx index to get variation
     * @param ch key to get variation from {@link SoftClip#seq}
     * @return variation
     */
    public static Variation getVariationFromSeq(SoftClip softClip, int idx, Character ch) {
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
     * Find closest mismatches to combine with indels.
     * @param referencePosition position for reference
     * @param readPosition position for base sequence
     * @param cigarLength length of CIGAR
     * @param querySequence base sequence
     * @param queryQuality base quality sequence
     * @param reference reference bases by position
     * @param refCoverage reference coverage
     * @return Tuple of (offset, querySequence's substring, queryQuality substring, number of mismatches)
     */
    //perl version: 1633
    public static Tuple.Tuple4<Integer, String, String, Integer> findOffset(int referencePosition,
                                                                     int readPosition,
                                                                     int cigarLength,
                                                                     String querySequence,
                                                                     String queryQuality,
                                                                     Map<Integer, Character> reference,
                                                                     Map<Integer, Integer> refCoverage) {
        int offset = 0;
        String ss = "";
        String q = "";
        int tnm = 0;
        int vsn = 0;
        for (int vi = 0; vsn <=  instance().conf.vext && vi < cigarLength; vi++) {
            if (querySequence.charAt(readPosition + vi) == 'N') {
                break;
            }
            if (queryQuality.charAt(readPosition + vi) - 33 < instance().conf.goodq) {
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
                VariationUtils.incCnt(refCoverage, referencePosition + osi, 1);
            }
        }
        return tuple(offset, ss, q, tnm);
    }

    /**
     * Increase count for variation
     * @param counters map to increase count for variation
     * @param idx index in <code>counters</code> map
     * @param s string with variations
     */
    public static void increment(Map<Integer, Map<String, Integer>> counters, int idx, String s) {
        Map<String, Integer> map = counters.get(idx);
        if (map == null) {
            map = new HashMap<>();
            counters.put(idx, map);
        }
        VariationUtils.incCnt(map, s, 1);
    }

    public static int getAlignedLenght(Cigar readCigar) {
        int lenght = 0;
        for (CigarElement element : readCigar.getCigarElements()) {
            if (element.getOperator() == CigarOperator.M || element.getOperator() == CigarOperator.D) {
                lenght += element.getLength();
            }
        }
        return lenght;
    }
}
