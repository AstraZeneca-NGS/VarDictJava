package com.astrazeneca.vardict.data;

import java.util.regex.Pattern;

/**
 * Regex Patterns from all classes of VarDict stored in one place.
 */
public class Patterns {
    // SAMRecord patterns
    public static final jregex.Pattern MC_Z_NUM_S_ANY_NUM_S = new jregex.Pattern("\\d+S\\S*\\d+S");

    //Variation patterns
    public static final Pattern BEGIN_DIGITS = Pattern.compile("^(\\d+)");
    public static final Pattern UP_NUMBER_END = Pattern.compile("\\^(\\d+)$");
    public static final Pattern BEGIN_MINUS_NUMBER_ANY = Pattern.compile("^-\\d+(.*)");
    public static final Pattern BEGIN_MINUS_NUMBER_CARET = Pattern.compile("^-\\d+\\^");
    public static final Pattern BEGIN_MINUS_NUMBER = Pattern.compile("^-(\\d+)");
    public static final jregex.Pattern MINUS_NUM_NUM = new jregex.Pattern("-\\d\\d");
    public static final Pattern HASH_GROUP_CARET_GROUP = Pattern.compile("#(.+)\\^(.+)");

    //Sclip patterns
    public static final Pattern B_A7 = Pattern.compile("^.AAAAAAA");
    public static final Pattern B_T7 = Pattern.compile("^.TTTTTTT");

    //ATGC patterns
    public static final Pattern CARET_ATGNC = Pattern.compile("\\^([ATGNC]+)");
    public static final Pattern CARET_ATGC_END = Pattern.compile("\\^([ATGC]+)$");
    public static final Pattern AMP_ATGC = Pattern.compile("&([ATGC]+)");
    public static final Pattern BEGIN_PLUS_ATGC = Pattern.compile("^\\+([ATGC]+)");
    public static final Pattern HASH_ATGC = Pattern.compile("#([ATGC]+)");
    public static final Pattern ATGSs_AMP_ATGSs_END = Pattern.compile("(\\+[ATGC]+)&[ATGC]+$");
    public static final Pattern MINUS_NUMBER_AMP_ATGCs_END = Pattern.compile("(-\\d+)&[ATGC]+$");
    public static final Pattern MINUS_NUMBER_ATGNC_SV_ATGNC_END = Pattern.compile("^-\\d+\\^([ATGNC]+)<...\\d+>([ATGNC]+)$");
    public static final Pattern BEGIN_ATGC_END = Pattern.compile("^[ATGC]+$");

    //SV patterns
    public static final Pattern DUP_NUM = Pattern.compile("<dup(\\d+)");
    public static final Pattern DUP_NUM_ATGC = Pattern.compile("<dup(\\d+)>([ATGC]+)$");
    public static final Pattern INV_NUM = Pattern.compile("<inv(\\d+)");
    public static final Pattern SOME_SV_NUMBERS = Pattern.compile("<(...)\\d+>");
    public static final Pattern ANY_SV = Pattern.compile("<(...)>");

    //File and columns patterns
    public static final Pattern SAMPLE_PATTERN = Pattern.compile("([^\\/\\._]+).sorted[^\\/]*.bam");
    public static final Pattern SAMPLE_PATTERN2 = Pattern.compile("([^\\/]+)[_\\.][^\\/]*bam");
    public static final Pattern INTEGER_ONLY = Pattern.compile("^\\d+$");

    //CIGAR patterns
    /**
     * Regexp finds number followed by S followed by number and I or D at the start of string
     */
    public static final jregex.Pattern BEGIN_NUMBER_S_NUMBER_IorD = new jregex.Pattern("^(\\d+)S(\\d+)([ID])");
    /**
     * Regexp finds number followed by I or D followed by number followed by S followed by end of string
     */
    public static final jregex.Pattern NUMBER_IorD_NUMBER_S_END = new jregex.Pattern("(\\d+)([ID])(\\d+)S$");
    /**
     * Regexp finds number
     * followed by S followed by number followed by M followed by number followed by I or D at the start of string
     */
    public static final jregex.Pattern BEGIN_NUMBER_S_NUMBER_M_NUMBER_IorD =
            new jregex.Pattern("^(\\d+)S(\\d+)M(\\d+)([ID])");
    /**
     * Regexp finds number
     * followed by I or D followed by number followed by M followed by number followed by S followed by end of string
     */
    public static final jregex.Pattern NUMBER_IorD_NUMBER_M_NUMBER_S_END =
            new jregex.Pattern("(\\d+)([ID])(\\d+)M(\\d+)S$");
    /**
     * Regexp finds digit-M-number-(I or D)-number-M
     */
    public static final jregex.Pattern BEGIN_DIGIT_M_NUMBER_IorD_NUMBER_M =
            new jregex.Pattern("^(\\d)M(\\d+)([ID])(\\d+)M");

    public static final Pattern BEGIN_DIGIT_M_NUMBER_IorD_NUMBER_M_ = Pattern.compile("^\\dM\\d+[ID]\\d+M");
    /**
     * Regexp replaces number-(I or D)-number-M
     */
    public static final jregex.Pattern NUMBER_IorD_DIGIT_M_END = new jregex.Pattern("(\\d+)([ID])(\\d)M$");
    public static final jregex.Pattern NUMBER_IorD_NUMBER_M_END = new jregex.Pattern("(\\d+)([ID])(\\d+)M$");
    /**
     * Regexp finds number-M-number-D-digit-M-number-I-digit-M-number-D. ^(.*?) captures everything before the
     * M-D-M-I-M-D complex
     */
    public static final jregex.Pattern D_M_D_DD_M_D_I_D_M_D_DD =
            new jregex.Pattern("^(.*?)(\\d+)M(\\d+)D(\\d+)M(\\d+)I(\\d+)M(\\d+)D(\\d+)M");
    public static final Pattern D_M_D_DD_M_D_I_D_M_D_DD_prim =
            Pattern.compile("(\\d+)M(\\d+)D(\\d+)M(\\d+)I(\\d+)M(\\d+)D(\\d+)M");
    public static final jregex.Pattern threeDeletionsPattern =
            new jregex.Pattern("^(.*?)(\\d+)M(\\d+)D(\\d+)M(\\d+)D(\\d+)M(\\d+)D(\\d+)" + "M");
    public static final jregex.Pattern threeIndelsPattern =
            new jregex.Pattern("^(.*?)(\\d+)M(\\d+)([DI])(\\d+)M(\\d+)([DI])(\\d+)M(\\d+)([DI])(\\d+)M");
    public static final Pattern DIGM_D_DI_DIGM_D_DI_DIGM_DI_DIGM =
            Pattern.compile("\\d+M\\d+[DI]\\d+M\\d+[DI]\\d+M\\d+[DI]\\d+M");
    public static final Pattern DM_DD_DM_DD_DM_DD_DM = Pattern.compile("\\d+M\\d+D\\d+M\\d+D\\d+M\\d+D\\d+M");
    /**
     * Regexp finds number-D-digit-M-number-D and optionally number-I
     */
    public static final Pattern DIG_D_DIG_M_DIG_DI_DIGI = Pattern.compile("(\\d+)D(\\d+)M(\\d+)([DI])(\\d+I)?");
    /**
     * Regexp finds number-I-digit-M-number and optionally number-I
     */
    public static final Pattern DIG_I_DIG_M_DIG_DI_DIGI = Pattern.compile("(\\d+)I(\\d+)M(\\d+)([DI])(\\d+I)?");
    /**
     * Regexp finds not_digit-number-I-number-M-number-D and optionally number-I
     */
    public static final Pattern NOTDIG_DIG_I_DIG_M_DIG_DI_DIGI =
            Pattern.compile("(\\D)(\\d+)I(\\d+)M(\\d+)([DI])(\\d+I)?");
    /**
     * Regexp finds number-D-number-D
     */
    public static final Pattern DIG_D_DIG_D = Pattern.compile("(\\d+)D(\\d+)D");
    /**
     * Regexp finds number-I-number-I
     */
    public static final Pattern DIG_I_DIG_I = Pattern.compile("(\\d+)I(\\d+)I");
    public static final Pattern BEGIN_ANY_DIG_M_END = Pattern.compile("^(.*?)(\\d+)M$");

    public static final Pattern DIG_M_END = Pattern.compile("\\d+M$");
    public static final Pattern BEGIN_DIG_M = Pattern.compile("^(\\d+)M");
    /**
     * Regexp finds number-S-number-M at start of CIGAR string
     */
    public static final Pattern DIG_S_DIG_M = Pattern.compile("^(\\d+)S(\\d+)M");
    public static final Pattern DIG_M_DIG_S_END = Pattern.compile("\\d+M\\d+S$");
    /**
     * Regexp finds number-M-number-S at end of CIGAR string ^(.*?) captures everything before M-S complex
     */
    public static final Pattern ANY_NUMBER_M_NUMBER_S_END = Pattern.compile("^(.*?)(\\d+)M(\\d+)S$");
    /**
     * Regexp finds number followed by D at the start of string
     */
    public static final jregex.Pattern BEGIN_NUMBER_D = new jregex.Pattern("^(\\d+)D");
    public static final jregex.Pattern END_NUMBER_D = new jregex.Pattern("(\\d+)D$");
    /**
     * Regexp finds number followed by I at the start of string
     */
    public static final jregex.Pattern BEGIN_NUMBER_I = new jregex.Pattern("^(\\d+)I");
    /**
     * Regexp finds number followed by I at the end of string
     */
    public static final jregex.Pattern END_NUMBER_I = new jregex.Pattern("(\\d+)I$");

    /**
     * Regexp finds numbers followed by M (matched) or D (deleted) in CIGAR string
     */
    public static final jregex.Pattern ALIGNED_LENGTH_MND = new jregex.Pattern("(\\d+)[MND]");

    /**
     * The total aligned length, excluding soft-clipped bases and insertions
     */
    public static final jregex.Pattern ALIGNED_LENGTH_MD = new jregex.Pattern("(\\d+)[MD=X]");

    public static final jregex.Pattern SOFT_CLIPPED = new jregex.Pattern("(\\d+)[MIS]");
    public static final Pattern SA_CIGAR_D_S_5clip = Pattern.compile("^\\d\\d+S");
    public static final Pattern SA_CIGAR_D_S_5clip_GROUP = Pattern.compile("^(\\d\\d+)S");
    public static final jregex.Pattern SA_CIGAR_D_S_5clip_GROUP_Repl = new jregex.Pattern("^\\d+S");
    public static final Pattern SA_CIGAR_D_S_3clip = Pattern.compile("\\d\\dS$");
    public static final Pattern SA_CIGAR_D_S_3clip_GROUP = Pattern.compile("(\\d\\d+)S$");
    public static final jregex.Pattern SA_CIGAR_D_S_3clip_GROUP_Repl = new jregex.Pattern("\\d\\d+S$");
    public static final Pattern BEGIN_dig_dig_S_ANY_dig_dig_S_END = Pattern.compile("^\\d\\dS.*\\d\\dS$");
    public static final jregex.Pattern BEGIN_NUM_S_OR_BEGIN_NUM_H = new jregex.Pattern("^(\\d+)S|^\\d+H");
    public static final jregex.Pattern END_NUM_S_OR_NUM_H = new jregex.Pattern("(\\d+)S$|H$");

    //Exception patterns
    public static final Pattern UNABLE_FIND_CONTIG = Pattern.compile("Unable to find entry for contig");
    public static final Pattern WRONG_START_OR_END = Pattern.compile("Malformed query");
}
