package com.astrazeneca.vardict.data.scopedata;

import com.astrazeneca.vardict.collection.VariationMap;
import com.astrazeneca.vardict.data.CurrentSegment;
import com.astrazeneca.vardict.data.SVStructures;
import com.astrazeneca.vardict.variations.Sclip;
import com.astrazeneca.vardict.variations.Variation;

import java.util.List;
import java.util.Map;

/**
 * The data created after variation realigner and structural variants steps in pipeline.
 * Used for creating Variants from Variations in toVarsBuilder step.
 */
public class RealignedVariationData {
    public final Map<Integer, VariationMap<String, Variation>> nonInsertionVariants;
    public final Map<Integer, VariationMap<String, Variation>> insertionVariants;
    public final Map<Integer, Sclip> softClips5End;
    public final Map<Integer, Sclip> softClips3End;
    public final Map<Integer, Integer> refCoverage;
    public final Integer maxReadLength;
    public final SVStructures svStructures;
    public final double duprate;
    public final CurrentSegment CURSEG;
    public final Map<Integer, List<Sclip>> SOFTP2SV;
    public final Scope<VariationData> previousScope;

    public RealignedVariationData(Map<Integer, VariationMap<String, Variation>> nonInsertionVariants,
                                  Map<Integer, VariationMap<String, Variation>> insertionVariants,
                                  Map<Integer, Sclip> softClips3End,
                                  Map<Integer, Sclip> softClips5End,
                                  Map<Integer, Integer> refCoverage,
                                  Integer maxReadLength,
                                  double duprate,
                                  SVStructures svStructures,
                                  CurrentSegment CURSEG,
                                  Map<Integer, List<Sclip>> SOFTP2SV,
                                  Scope<VariationData> previousScope) {
        this.softClips3End = softClips3End;
        this.softClips5End = softClips5End;
        this.nonInsertionVariants = nonInsertionVariants;
        this.insertionVariants = insertionVariants;
        this.refCoverage = refCoverage;
        this.maxReadLength = maxReadLength;
        this.svStructures = svStructures;
        this.duprate = duprate;
        this.CURSEG = CURSEG;
        this.SOFTP2SV = SOFTP2SV;
        this.previousScope = previousScope;
    }
}
