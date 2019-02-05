package com.astrazeneca.vardict.printers;

/**
 * Abstract class that can be inherited for different types of variants that will contain different number of fields.
 * Common fields for all type of variants are collected in this class. Inherited classes contain only specific fields.
 */
public abstract class OutputVariant {
    protected String delimiter = "\t";

    protected String sample = "";
    protected String gene;
    protected String chr = "";
    protected int startPosition;
    protected int endPosition;
    protected String refAllele = "";
    protected String varAllele = "";

    protected int shift3;
    protected double msi;
    protected int msint;

    protected String leftSequence = "";
    protected String rightSequence = "";
    protected String region = "";
    protected String varType = "";
    protected String DEBUG = "";

    /**
     * Set delimiter to print variants between fields. Default is <code>\t</code> (tab delimiter).
     * @param delimiter string contains delimiter
     */
    public void setDelimiter(String delimiter) {
        this.delimiter = delimiter;
    }
}
