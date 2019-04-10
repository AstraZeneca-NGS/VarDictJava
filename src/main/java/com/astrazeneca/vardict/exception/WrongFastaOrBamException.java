package com.astrazeneca.vardict.exception;


import java.util.Locale;

public class WrongFastaOrBamException extends RuntimeException {
    public final static String WrongFastaOrBamExceptionMeassage = "The name of this chromosome \"%s\" is missing in your" +
            " fasta file. Please be sure that chromosome names in BAM, fasta and BED are in correspondence " +
            "with each other and you use correct fasta for your BAM (can be checked in BAM header).";

    public WrongFastaOrBamException(String chr, Throwable e) {
            super(String.format(Locale.US, WrongFastaOrBamExceptionMeassage, chr), e);
    }
}