package com.astrazeneca.vardict.data.scopedata;

/**
 * The data created after combine analysis in somatic mode. Used for determining a type of Variant.
 */
public class CombineAnalysisData {
    public int maxReadLength;
    public String type;

    public CombineAnalysisData(int maxReadLength, String type) {
        this.maxReadLength = maxReadLength;
        this.type = type;
    }
}
