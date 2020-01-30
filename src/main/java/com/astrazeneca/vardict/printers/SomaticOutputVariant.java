package com.astrazeneca.vardict.printers;

import com.astrazeneca.vardict.data.Region;
import com.astrazeneca.vardict.data.fishertest.FisherExact;
import com.astrazeneca.vardict.variations.Variant;

import java.text.DecimalFormat;

import static com.astrazeneca.vardict.Utils.getRoundedValueToPrint;
import static com.astrazeneca.vardict.Utils.join;
import static com.astrazeneca.vardict.data.scopedata.GlobalReadOnlyScope.instance;

/**
 * Variant created in Somatic mode. Must contains 55 total fields if no --fisher.
 */
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
    private double pvalue1;
    private String oddratio1 = "0";
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
    private double pvalue2;
    private String oddratio2 = "0";
    private double var2meanQuality;
    private int var2hasAtLeast2DiffQualities;
    private double var2meanMappingQuality;
    private double var2highQualityToLowQualityRatio;
    private double var2highQualityReadsFrequency;
    private double var2extraFrequency;
    private double var2nm;
    private double var2duprate;
    private String var2sv = "0";

    private String varLabel = "";

    private double pvalue;
    private String oddratio = "0";

    public SomaticOutputVariant(Variant beginVariant, Variant endVariant, Variant tumorVariant, Variant normalVariant,
                                Region region, String sv1, String sv2, String varLabel) {
        this.sample = instance().sample;
        this.gene = region.gene;
        this.chr = region.chr;

        if (beginVariant != null) {
            this.startPosition = beginVariant.startPosition;
            this.endPosition = beginVariant.endPosition;
            this.refAllele = beginVariant.refallele;
            this.varAllele = beginVariant.varallele;
            this.varType = beginVariant.vartype;
            this.DEBUG = beginVariant.DEBUG;
        }
        if (endVariant != null) {
            this.shift3 = endVariant.shift3;
            this.msi = endVariant.msi;
            this.msint = endVariant.msint;
            this.leftSequence = endVariant.leftseq.isEmpty() ? "0" : endVariant.leftseq;
            this.rightSequence = endVariant.rightseq.isEmpty() ? "0" : endVariant.rightseq;
        }

        if (tumorVariant != null) {
            this.var1totalCoverage = tumorVariant.totalPosCoverage;
            this.var1variantCoverage = tumorVariant.positionCoverage;
            this.var1refForwardCoverage = tumorVariant.refForwardCoverage;
            this.var1refReverseCoverage = tumorVariant.refReverseCoverage;
            this.var1variantForwardCount = tumorVariant.varsCountOnForward;
            this.var1variantReverseCount = tumorVariant.varsCountOnReverse;
            this.var1genotype = tumorVariant.genotype == null ? "0" : tumorVariant.genotype;
            this.var1frequency = tumorVariant.frequency;
            this.var1strandBiasFlag = tumorVariant.strandBiasFlag == null ? "0" : tumorVariant.strandBiasFlag;
            this.var1meanPosition = tumorVariant.meanPosition;
            this.var1isAtLeastAt2Position = tumorVariant.isAtLeastAt2Positions ? 1 : 0;
            this.var1meanQuality = tumorVariant.meanQuality;
            this.var1hasAtLeast2DiffQualities = tumorVariant.hasAtLeast2DiffQualities ? 1 : 0;
            this.var1meanMappingQuality = tumorVariant.meanMappingQuality;
            this.var1highQualityToLowQualityRatio = tumorVariant.highQualityToLowQualityRatio;
            this.var1highQualityReadsFrequency = tumorVariant.highQualityReadsFrequency;
            this.var1extraFrequency = tumorVariant.extraFrequency;
            this.var1nm = tumorVariant.numberOfMismatches;
            this.var1duprate = tumorVariant.duprate;
        }

        if (normalVariant != null) {
            this.var2totalCoverage = normalVariant.totalPosCoverage;
            this.var2variantCoverage = normalVariant.positionCoverage;
            this.var2refForwardCoverage = normalVariant.refForwardCoverage;
            this.var2refReverseCoverage = normalVariant.refReverseCoverage;
            this.var2variantForwardCount = normalVariant.varsCountOnForward;
            this.var2variantReverseCount = normalVariant.varsCountOnReverse;
            this.var2genotype = normalVariant.genotype == null ? "0" : normalVariant.genotype;
            this.var2frequency = normalVariant.frequency;
            this.var2strandBiasFlag = normalVariant.strandBiasFlag == null ? "0" : normalVariant.strandBiasFlag;
            this.var2meanPosition = normalVariant.meanPosition;
            this.var2isAtLeastAt2Position = normalVariant.isAtLeastAt2Positions ? 1 : 0;
            this.var2meanQuality = normalVariant.meanQuality;
            this.var2hasAtLeast2DiffQualities = normalVariant.hasAtLeast2DiffQualities ? 1 : 0;
            this.var2meanMappingQuality = normalVariant.meanMappingQuality;
            this.var2highQualityToLowQualityRatio = normalVariant.highQualityToLowQualityRatio;
            this.var2highQualityReadsFrequency = normalVariant.highQualityReadsFrequency;
            this.var2extraFrequency = normalVariant.extraFrequency;
            this.var2nm = normalVariant.numberOfMismatches;
            this.var2duprate = normalVariant.duprate;
        }

        this.varLabel = varLabel;
        this.region = region.chr + ":" + region.start + "-" + region.end;
        this.var1sv = sv1.equals("") ? "0" : sv1;
        this.var2sv = sv2.equals("") ? "0" : sv2;

        if (instance().conf.fisher) {
            calculateFisherSomatic(tumorVariant, normalVariant);
        }
    }

    private void calculateFisherSomatic(Variant tumorVariant, Variant normalVariant) {
        FisherExact fisher;

        if (tumorVariant != null) {
            fisher = new FisherExact(tumorVariant.refForwardCoverage, tumorVariant.refReverseCoverage,
                    tumorVariant.varsCountOnForward, tumorVariant.varsCountOnReverse);
        } else {
            fisher = new FisherExact(0, 0, 0, 0);
        }
        this.pvalue1 = fisher.getPValue();
        this.oddratio1 = fisher.getOddRatio();

        if (normalVariant != null) {
            fisher = new FisherExact(normalVariant.refForwardCoverage, normalVariant.refReverseCoverage,
                    normalVariant.varsCountOnForward, normalVariant.varsCountOnReverse);
        } else {
            fisher = new FisherExact(0, 0, 0, 0);
        }
        this.pvalue2 = fisher.getPValue();
        this.oddratio2 = fisher.getOddRatio();

        int tref = this.var1totalCoverage - this.var1variantCoverage;
        int rref = this.var2totalCoverage - this.var2variantCoverage;
        if (tref < 0) {
            tref = 0;
        }
        if (rref < 0) {
            rref = 0;
        }
        fisher = new FisherExact(this.var1variantCoverage, tref, this.var2variantCoverage, rref);
        double pvalue_greater = fisher.getPValueGreater();
        double pvalue_less = fisher.getPValueLess();

        if (pvalue_less < pvalue_greater) {
            this.pvalue = pvalue_less;
        } else {
            this.pvalue = pvalue_greater;
        }
        this.oddratio = fisher.getOddRatio();
    }

    @Override
    public String toString() {
        String outputVariant;
        if (!instance().conf.fisher) {
            outputVariant = create_somatic_variant_55columns();
        } else {
            outputVariant = create_somatic_variant_61columns();
        }
        if (instance().conf.debug) {
            outputVariant = join(delimiter, outputVariant, DEBUG);
        }
        return outputVariant;
    }

    /**
     * 61 columns: oddratio and pvalue added for each sample and total and format as it will be after using R script
     */
    private String create_somatic_variant_61columns() {
        String outputVariant;
        var1nm = var1nm > 0 ? var1nm : 0;
        var2nm = var2nm > 0 ? var2nm : 0;
        String msi_f = msi == 0 ? "0" : new DecimalFormat("0.000").format(msi);

        outputVariant = join(delimiter,
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
                getRoundedValueToPrint("0.0000", var1frequency),
                var1strandBiasFlag,
                getRoundedValueToPrint("0.0", var1meanPosition),
                var1isAtLeastAt2Position,
                getRoundedValueToPrint("0.0", var1meanQuality),
                var1hasAtLeast2DiffQualities,
                getRoundedValueToPrint("0.0", var1meanMappingQuality),
                getRoundedValueToPrint("0.000", var1highQualityToLowQualityRatio),
                getRoundedValueToPrint("0.0000", var1highQualityReadsFrequency),
                getRoundedValueToPrint("0.0000", var1extraFrequency),
                getRoundedValueToPrint("0.0", var1nm),
                getRoundedValueToPrint("0.00000", pvalue1),
                oddratio1,

                var2totalCoverage,
                var2variantCoverage,
                var2refForwardCoverage,
                var2refReverseCoverage,
                var2variantForwardCount,
                var2variantReverseCount,
                var2genotype,
                getRoundedValueToPrint("0.0000", var2frequency),
                var2strandBiasFlag,
                getRoundedValueToPrint("0.0", var2meanPosition),
                var2isAtLeastAt2Position,
                getRoundedValueToPrint("0.0", var2meanQuality),
                var2hasAtLeast2DiffQualities,
                getRoundedValueToPrint("0.0", var2meanMappingQuality),
                getRoundedValueToPrint("0.000", var2highQualityToLowQualityRatio),
                getRoundedValueToPrint("0.0000", var2highQualityReadsFrequency),
                getRoundedValueToPrint("0.0000", var2extraFrequency),
                getRoundedValueToPrint("0.0", var2nm),
                getRoundedValueToPrint("0.00000", pvalue2),
                oddratio2,

                shift3,
                msi_f,
                msint,
                leftSequence, rightSequence,
                region,
                varLabel,
                varType,
                getRoundedValueToPrint("0.0", var1duprate),
                var1sv,
                getRoundedValueToPrint("0.0", var2duprate),
                var2sv,
                getRoundedValueToPrint("0.00000", pvalue),
                oddratio
        );
        return outputVariant;
    }

    /**
     * 55 columns: not oddratio and pvalue columns
     */
    private String create_somatic_variant_55columns() {
        String outputVariant;
        outputVariant = join(delimiter,
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
        return outputVariant;
    }
}
