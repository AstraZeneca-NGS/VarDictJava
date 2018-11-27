package com.astrazeneca.vardict.data.scopedata;

import com.astrazeneca.vardict.variations.Vars;

import java.util.Map;

/**
 * The data after ToVarsBuilder step in pipeline. Used for creating output variants in all modes of VarDict.
 * Contains max read lengths and map of positions-Vars.
 */

public class AlignedVarsData {
    public int maxReadLength;
    public Map<Integer, Vars> alignedVariants;

    public AlignedVarsData(int maxReadLength, Map<Integer, Vars> alignedVariants) {
        this.maxReadLength = maxReadLength;
        this.alignedVariants = alignedVariants;
    }
}
