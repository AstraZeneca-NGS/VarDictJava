package com.astrazeneca.vardict.printers;

public abstract class OutputVariant {
    protected String delimiter;

    // Common fields for all type of variants
    protected String sample = "";
    protected String gene;
    protected String chr = "";
    protected int startPosition;
    protected int endPosition;
    protected String refAllele = "";
    protected String varAllele = "";

    protected String leftSequence = "";
    protected String rightSequence = "";
    protected String region = "";
    protected String varType = "";
    protected String DEBUG = "";

    protected abstract String outputString(String delimiter);
}
