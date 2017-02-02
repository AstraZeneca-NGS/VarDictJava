package com.astrazeneca.vardict.pipeline.data;

import com.astrazeneca.vardict.variations.SoftClip;
import com.astrazeneca.vardict.variations.Variation;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;

public class VariationData {

    public Map<Integer, Map<String, Variation>> nonInsertionVariants;
    public Map<Integer, Map<String, Variation>> insertionVariants;
    public Map<Integer, Map<String, Integer>> positionToInsertionCount;
    public Map<Integer, Map<String, Integer>> positionToDeletionsCount;
    public Map<Integer, Integer> refCoverage;
    public Map<Integer, SoftClip> softClips5End;
    public Map<Integer, SoftClip> softClips3End;
    public Integer maxReadLength;
    public Set<String> splice;
    public Map<Integer, Map<String, Integer>> mnp;
    public Map<String, int[]> spliceCount;

    public VariationData() {}

    public VariationData(Integer maxReadLength, Set<String> splice) {
        this.maxReadLength = maxReadLength;
        this.splice = splice;

        this.nonInsertionVariants = new HashMap<>();
        this.insertionVariants = new HashMap<>();
        this.positionToDeletionsCount = new HashMap<>();
        this.positionToInsertionCount = new HashMap<>();
        this.refCoverage = new HashMap<>();
        this.softClips3End = new HashMap<>();
        this.softClips5End = new HashMap<>();
        this.mnp = new HashMap<>();
        this.spliceCount = new HashMap<>();
    }
}

