package com.astrazeneca.vardict.exception;


public class RegionMissedSourceException extends RuntimeException {
    public final static String RegionSourceMissedMessage = "The required BED file or region missed, please, set it " +
            "with path to BED or with -R option.";

    public RegionMissedSourceException() {
            super(RegionSourceMissedMessage);
    }
}