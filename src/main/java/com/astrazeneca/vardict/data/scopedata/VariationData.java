package com.astrazeneca.vardict.data.scopedata;

import com.astrazeneca.vardict.collection.VariationMap;
import com.astrazeneca.vardict.data.SVStructures;
import com.astrazeneca.vardict.variations.Sclip;
import com.astrazeneca.vardict.variations.Variation;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;

/**
 * The data created after CigarParser step in pipeline. Used for process Variation in
 * Variation Realigner and Structural Variants analysis (realign it and searching for structural variants).
 */
public class VariationData {
    public Map<Integer, VariationMap<String, Variation>> nonInsertionVariants;
    public Map<Integer, VariationMap<String, Variation>> insertionVariants;
    public Map<Integer, Map<String, Integer>> positionToInsertionCount;
    public Map<Integer, Map<String, Integer>> positionToDeletionsCount;
    public SVStructures svStructures;
    public Map<Integer, Integer> refCoverage;
    public Map<Integer, Sclip> softClips5End;
    public Map<Integer, Sclip> softClips3End;
    public Integer maxReadLength;
    public Set<String> splice;
    public Map<Integer, Map<String, Integer>> mnp;
    public Map<String, int[]> spliceCount;
    public double duprate;

}

