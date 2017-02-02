package com.astrazeneca.vardict.pipeline.data;

import com.astrazeneca.vardict.variations.Variation;

import java.util.Map;

public class RealignedVariationData {

    public final Map<Integer, Map<String, Variation>> nonInsertionVariants;
    public final Map<Integer, Map<String, Variation>> insertionVariants;
    public final Map<Integer, Integer> refCoverage;
    public final Integer maxReadLength;

    public RealignedVariationData(Map<Integer, Map<String, Variation>> nonInsertionVariants,
                                  Map<Integer, Map<String, Variation>> insertionVariants,
                                  Map<Integer, Integer> refCoverage,
                                  Integer maxReadLength) {
        this.nonInsertionVariants = nonInsertionVariants;
        this.insertionVariants = insertionVariants;
        this.refCoverage = refCoverage;
        this.maxReadLength = maxReadLength;
    }
}
