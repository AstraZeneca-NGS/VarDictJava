package com.astrazeneca.vardict.printers;

import com.astrazeneca.vardict.data.Region;
import com.astrazeneca.vardict.variations.Variant;
import com.astrazeneca.vardict.variations.Vars;

import java.text.DecimalFormat;

import static com.astrazeneca.vardict.Utils.join;
import static com.astrazeneca.vardict.data.scopedata.GlobalReadOnlyScope.instance;

public class SomaticOutputVariant extends OutputVariant {
    private int var1totalCoverage;
    private int var1variantCoverage;
    private int var1refForwardCoverage;
    private int var1refReverseCoverage;
    private int var1variantForwardCount;
    private int var1variantReverseCount;
    private String var1genotype = "0";
    private double var1frequency;
    private String var1strandBiasFlag = "0";
    private double var1meanPosition;
    private int var1isAtLeastAt2Position;
    private double var1meanQuality;
    private int var1hasAtLeast2DiffQualities;
    private double var1meanMappingQuality;
    private double var1highQualityToLowQualityRatio;
    private double var1highQualityReadsFrequency;
    private double var1extraFrequency;
    private double var1nm;
    private double var1duprate;
    private String var1sv = "0";

    private int var2totalCoverage;
    private int var2variantCoverage;
    private int var2refForwardCoverage;
    private int var2refReverseCoverage;
    private int var2variantForwardCount;
    private int var2variantReverseCount;
    private String var2genotype = "0";
    private double var2frequency;
    private String var2strandBiasFlag = "0";
    private double var2meanPosition;
    private int var2isAtLeastAt2Position;
    private double var2meanQuality;
    private int var2hasAtLeast2DiffQualities;
    private double var2meanMappingQuality;
    private double var2highQualityToLowQualityRatio;
    private double var2highQualityReadsFrequency;
    private double var2extraFrequency;
    private double var2nm;
    private double var2duprate;
    private String var2sv = "0";

    private int shift3;
    private double msi;
    private int msint;
    private String varLabel = "";

    public SomaticOutputVariant(Variant common, Variant variant1, Variant variant2, Region region, String sv1,
                                String sv2, String varLabel) {
        this.sample = instance().sample;
        this.gene = region.gene;
        this.chr = region.chr;

        if (common!= null) {
            this.startPosition = common.startPosition;
            this.endPosition = common.endPosition;
            this.refAllele = common.refallele;
            this.varAllele = common.varallele;

            this.shift3 = common.shift3;
            this.msi = common.msi;
            this.msint = common.msint;
            this.leftSequence = common.leftseq.isEmpty() ? "0" : common.leftseq;
            this.rightSequence = common.rightseq.isEmpty() ? "0" : common.rightseq;
            this.varType = common.vartype;
            this.varLabel = varLabel;
            this.region = region.chr + ":" + region.start + "-" + region.end;
            this.DEBUG = common.DEBUG;
        }

        if (variant1 != null) {
            this.var1totalCoverage = variant1.totalPosCoverage;
            this.var1variantCoverage = variant1.positionCoverage;
            this.var1refForwardCoverage = variant1.refForwardCoverage;
            this.var1refReverseCoverage = variant1.refReverseCoverage;
            this.var1variantForwardCount = variant1.varsCountOnForward;
            this.var1variantReverseCount = variant1.varsCountOnReverse;
            this.var1genotype = variant1.genotype == null ? "0" : variant1.genotype;
            this.var1frequency = variant1.frequency;
            this.var1strandBiasFlag = variant1.strandBiasFlag == null ? "0" : variant1.strandBiasFlag;
            this.var1meanPosition = variant1.meanPosition;
            this.var1isAtLeastAt2Position = variant1.isAtLeastAt2Positions ? 1 : 0;
            this.var1meanQuality = variant1.meanQuality;
            this.var1hasAtLeast2DiffQualities = variant1.hasAtLeast2DiffQualities ? 1 : 0;
            this.var1meanMappingQuality = variant1.meanMappingQuality;
            this.var1highQualityToLowQualityRatio = variant1.highQualityToLowQualityRatio;
            this.var1highQualityReadsFrequency = variant1.highQualityReadsFrequency;
            this.var1extraFrequency = variant1.extraFrequency;
            this.var1nm = variant1.numberOfMismatches;
            this.var1duprate = variant1.duprate;
        }

        if (variant2 != null) {
            this.var2totalCoverage = variant2.totalPosCoverage;
            this.var2variantCoverage = variant2.positionCoverage;
            this.var2refForwardCoverage = variant2.refForwardCoverage;
            this.var2refReverseCoverage = variant2.refReverseCoverage;
            this.var2variantForwardCount = variant2.varsCountOnForward;
            this.var2variantReverseCount = variant2.varsCountOnReverse;
            this.var2genotype = variant2.genotype == null ? "0" : variant2.genotype;
            this.var2frequency = variant2.frequency;
            this.var2strandBiasFlag = variant2.strandBiasFlag == null ? "0" : variant2.strandBiasFlag;
            this.var2meanPosition = variant2.meanPosition;
            this.var2isAtLeastAt2Position = variant2.isAtLeastAt2Positions ? 1 : 0;
            this.var2meanQuality = variant2.meanQuality;
            this.var2hasAtLeast2DiffQualities = variant2.hasAtLeast2DiffQualities ? 1 : 0;
            this.var2meanMappingQuality = variant2.meanMappingQuality;
            this.var2highQualityToLowQualityRatio = variant2.highQualityToLowQualityRatio;
            this.var2highQualityReadsFrequency = variant2.highQualityReadsFrequency;
            this.var2extraFrequency = variant2.extraFrequency;
            this.var2nm = variant2.numberOfMismatches;
            this.var2duprate = variant2.duprate;
        }
        this.var1sv = sv1.equals("") ? "0" : sv1;
        this.var2sv = sv2.equals("") ? "0" : sv2;
    }

    @Override
    public String outputString(String delimiter) {
        this.delimiter = delimiter;
        // 55 columns
        String outputVariant = join(delimiter,
                sample,
                gene,
                chr,
                startPosition,
                endPosition,
                refAllele,
                varAllele,

                var1totalCoverage,
                var1variantCoverage,
                var1refForwardCoverage,
                var1refReverseCoverage,
                var1variantForwardCount,
                var1variantReverseCount,
                var1genotype,
                var1frequency == 0 ? 0 : new DecimalFormat("0.0000").format(var1frequency),
                var1strandBiasFlag,
                var1meanPosition == 0 ? 0 : new DecimalFormat("0.0").format(var1meanPosition),
                var1isAtLeastAt2Position,
                var1meanQuality == 0 ? 0 : new DecimalFormat("0.0").format(var1meanQuality),
                var1hasAtLeast2DiffQualities,
                var1meanMappingQuality == 0 ? 0 : new DecimalFormat("0.0").format(var1meanMappingQuality),
                var1highQualityToLowQualityRatio == 0 ? 0 : new DecimalFormat("0.000").format(var1highQualityToLowQualityRatio),
                var1highQualityReadsFrequency == 0 ? 0 : new DecimalFormat("0.0000").format(var1highQualityReadsFrequency),
                var1extraFrequency == 0 ? 0 : new DecimalFormat("0.0000").format(var1extraFrequency),
                var1nm > 0 ? new DecimalFormat("0.0").format(var1nm) : 0,

                var2totalCoverage,
                var2variantCoverage,
                var2refForwardCoverage,
                var2refReverseCoverage,
                var2variantForwardCount,
                var2variantReverseCount,
                var2genotype,
                var2frequency == 0 ? 0 : new DecimalFormat("0.0000").format(var2frequency),
                var2strandBiasFlag,
                var2meanPosition == 0 ? 0 : new DecimalFormat("0.0").format(var2meanPosition),
                var2isAtLeastAt2Position,
                var2meanQuality == 0 ? 0 : new DecimalFormat("0.0").format(var2meanQuality),
                var2hasAtLeast2DiffQualities,
                var2meanMappingQuality == 0 ? 0 : new DecimalFormat("0.0").format(var2meanMappingQuality),
                var2highQualityToLowQualityRatio == 0 ? 0 : new DecimalFormat("0.000").format(var2highQualityToLowQualityRatio),
                var2highQualityReadsFrequency == 0 ? 0 : new DecimalFormat("0.0000").format(var2highQualityReadsFrequency),
                var2extraFrequency == 0 ? 0 : new DecimalFormat("0.0000").format(var2extraFrequency),
                var2nm > 0 ? new DecimalFormat("0.0").format(var2nm) : 0,

                shift3,
                msi == 0 ? 0 : new DecimalFormat("0.000").format(msi),
                msint,
                leftSequence, rightSequence,
                region,
                varLabel,
                varType,
                var1duprate == 0 ? 0 : new DecimalFormat("0.0").format(var1duprate),
                var1sv,
                var2duprate == 0 ? 0 : new DecimalFormat("0.0").format(var2duprate),
                var2sv
        );
        if (instance().conf.debug) {
            outputVariant = join(delimiter, outputVariant, DEBUG);
        }
        return outputVariant;
    }
}
