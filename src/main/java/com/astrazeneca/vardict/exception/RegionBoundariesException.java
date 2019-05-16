package com.astrazeneca.vardict.exception;


import java.util.Locale;

public class RegionBoundariesException extends RuntimeException {
    public final static String RegionBoundariesExceptionMessage = "The region %s:%d-%d is wrong. " +
            "We have problem while reading it, possible the start is after the end of the region or " +
            "the fasta doesn't contain this region.";

    public RegionBoundariesException(String chr, int start, int end, Throwable e) {
            super(String.format(Locale.US, RegionBoundariesExceptionMessage, chr, start, end) , e);
    }
}