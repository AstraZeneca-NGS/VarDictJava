package com.astrazeneca.vardict.data.scopedata;

import com.astrazeneca.vardict.collection.VariationMap;
import com.astrazeneca.vardict.variations.Sclip;
import com.astrazeneca.vardict.variations.Variation;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;

/**
 * Initial data to start VarDict pipelines. Create maps for variations, reference coverage and softclips.
 */
public class InitialData {
    public Map<Integer, VariationMap<String, Variation>> nonInsertionVariants;
    public Map<Integer, VariationMap<String, Variation>> insertionVariants;
    public Map<Integer, Integer> refCoverage;
    public Map<Integer, Sclip> softClips5End;
    public Map<Integer, Sclip> softClips3End;

    public InitialData() {
        this.nonInsertionVariants = new HashMap<>();
        this.insertionVariants = new HashMap<>();
        this.refCoverage = new HashMap<>();
        this.softClips3End = new HashMap<>();
        this.softClips5End = new HashMap<>();
    }

    public InitialData(Map<Integer, VariationMap<String, Variation>> nonInsertionVariants,
                       Map<Integer, VariationMap<String, Variation>> insertionVariants,
                       Map<Integer, Integer> refCoverage,
                       Map<Integer, Sclip> softClips3End,
                       Map<Integer, Sclip> softClips5End) {
        this.nonInsertionVariants = nonInsertionVariants;
        this.insertionVariants = insertionVariants;
        this.refCoverage = refCoverage;
        this.softClips3End = softClips3End;
        this.softClips5End = softClips5End;
    }
}

