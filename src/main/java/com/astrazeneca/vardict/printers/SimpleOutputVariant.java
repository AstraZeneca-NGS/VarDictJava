package com.astrazeneca.vardict.printers;

import com.astrazeneca.vardict.data.Region;
import com.astrazeneca.vardict.data.fishertest.FisherExact;
import com.astrazeneca.vardict.variations.Variant;

import java.text.DecimalFormat;

import static com.astrazeneca.vardict.data.scopedata.GlobalReadOnlyScope.instance;
import static com.astrazeneca.vardict.Utils.join;
/**
 * Variant created in Simple mode. Must contain 36 total fields if no --fisher.
 */
public class SimpleOutputVariant extends OutputVariant {
    private int totalCoverage;
    private int variantCoverage;
    private int referenceForwardCount;
    private int referenceReverseCount;
    private int variantForwardCount;
    private int variantReverseCount;
    private String genotype = "";
    private double frequency;
    private String bias = "0;0";
    private double pmean;
    private int pstd;
    private double qual;
    private int qstd;
    private double mapq;
    private double qratio;
    private double hifreq;
    private double extrafreq;

    private double nm;
    private int hicnt;
    private int hicov;

    private double duprate;
    private int crispr;
    private String sv;

    private double pvalue;
    private String oddratio = "0";

    public SimpleOutputVariant(Variant variant, Region region, String sv, int position) {
        this.sample = instance().sample;
        this.gene = region.gene;
        this.chr = region.chr;

        if (variant == null) {
            this.startPosition = position;
            this.endPosition = position;
        } else {
            this.startPosition = variant.startPosition;
            this.endPosition = variant.endPosition;
            this.refAllele = variant.refallele;
            this.varAllele = variant.varallele;
            this.totalCoverage = variant.totalPosCoverage;
            this.variantCoverage = variant.positionCoverage;
            this.referenceForwardCount = variant.refForwardCoverage;
            this.referenceReverseCount = variant.refReverseCoverage;
            this.variantForwardCount = variant.varsCountOnForward;
            this.variantReverseCount = variant.varsCountOnReverse;
            this.genotype = variant.genotype == null ? "0" : variant.genotype;
            this.frequency = variant.frequency;
            this.bias = variant.strandBiasFlag;
            this.pmean = variant.meanPosition;
            this.pstd = variant.isAtLeastAt2Positions ? 1 : 0;
            this.qual = variant.meanQuality;
            this.qstd = variant.hasAtLeast2DiffQualities ? 1 : 0;
            this.mapq = variant.meanMappingQuality;
            this.qratio = variant.highQualityToLowQualityRatio;
            this.hifreq = variant.highQualityReadsFrequency;
            this.extrafreq = variant.extraFrequency;
            this.shift3 = variant.shift3;
            this.msi = variant.msi;
            this.msint = variant.msint;
            this.nm = variant.numberOfMismatches;
            this.hicnt = variant.hicnt;
            this.hicov = variant.hicov;
            this.leftSequence = variant.leftseq.isEmpty() ? "0" : variant.leftseq;
            this.rightSequence = variant.rightseq.isEmpty() ? "0" : variant.rightseq;
            this.varType = variant.vartype;
            this.duprate = variant.duprate;
            this.crispr = variant.crispr;
            this.DEBUG = variant.DEBUG;
        }

        if (instance().conf.fisher) {
            FisherExact fisher;
            if (variant != null) {
                fisher = new FisherExact(variant.refForwardCoverage, variant.refReverseCoverage,
                        variant.varsCountOnForward, variant.varsCountOnReverse);
            } else {
                fisher = new FisherExact(0, 0, 0, 0);
            }
            this.pvalue = fisher.getPValue();
            this.oddratio = fisher.getOddRatio();
        }

        this.region = region.chr + ":" + region.start + "-" + region.end;
        this.sv = sv.equals("") ? "0" : sv;
    }

    @Override
    public String toString() {
        String outputVariant;
        if (!instance().conf.fisher) {
            outputVariant = create_simple_variant_36columns();
        } else {
            outputVariant = create_simple_variant_38columns();
        }
        if (instance().conf.crisprCuttingSite != 0) {
            outputVariant = join(delimiter, outputVariant, crispr);
        }
        if (instance().conf.debug) {
            outputVariant = join(delimiter, outputVariant, DEBUG);
        }
        return outputVariant;
    }

    /**
     * 38 columns: oddratio and pvalue added + format as it will be after using R script
     */
    private String create_simple_variant_38columns() {
        String outputVariant;
        String frequency_f = frequency == Math.round(frequency)
                ? new DecimalFormat("0").format(frequency)
                : new DecimalFormat("0.0000").format(frequency).replaceAll("0+$", "");
        String pmean_f = pmean == Math.round(pmean)
                ? new DecimalFormat("0").format(pmean)
                : new DecimalFormat("0.0").format(pmean).replaceAll("0+$", "");
        String qual_f = qual == Math.round(qual)
                ? new DecimalFormat("0").format(qual)
                : new DecimalFormat("0.0").format(qual).replaceAll("0+$", "");
        String pvalue_f = pvalue == Math.round(pvalue)
                ? new DecimalFormat("0").format(pvalue)
                : new DecimalFormat("0.00000").format(pvalue).replaceAll("0+$", "");
        String mapq_f = mapq == Math.round(mapq)
                ? new DecimalFormat("0").format(mapq)
                : new DecimalFormat("0.0").format(mapq).replaceAll("0+$", "");
        String qratio_f = qratio == Math.round(qratio)
                ? new DecimalFormat("0").format(qratio)
                : new DecimalFormat("0.000").format(qratio).replaceAll("0+$", "");
        String hifreq_f = hifreq == 0
                ? "0"
                : new DecimalFormat("0.0000").format(hifreq);
        String extrafreq_f = extrafreq == Math.round(extrafreq)
                ? new DecimalFormat("0").format(extrafreq)
                :  new DecimalFormat("0.0000").format(extrafreq).replaceAll("0+$", "");
        String msi_f = msi == Math.round(msi)
                ? new DecimalFormat("0").format(msi)
                : new DecimalFormat("0.000").format(msi).replaceAll("0+$", "");
        nm = nm > 0 ? nm : 0;
        String nm_f = nm == 0
                ? "0"
                : new DecimalFormat("0.0").format(nm);
        String duprate_f = duprate == Math.round(duprate)
                ? new DecimalFormat("0").format(duprate)
                : new DecimalFormat("0.0").format(duprate).replaceAll("0+$", "");
        outputVariant = join(delimiter,
                sample,
                gene,
                chr,
                startPosition,
                endPosition,
                refAllele,
                varAllele,

                totalCoverage,
                variantCoverage,
                referenceForwardCount,
                referenceReverseCount,
                variantForwardCount,
                variantReverseCount,
                genotype,
                frequency_f,
                bias,
                pmean_f,
                pstd,
                qual_f,
                qstd,

                String.valueOf(pvalue_f),
                oddratio,

                mapq_f,
                qratio_f,
                hifreq_f,
                extrafreq_f,

                shift3,
                msi_f,
                msint,
                nm_f,
                hicnt,
                hicov,
                leftSequence, rightSequence,
                region,
                varType,
                duprate_f,
                sv
        );
        return outputVariant;
    }

    /**
     * 36 columns for simple mode variant: no columns for oddratio and pvalue
     */
    private String create_simple_variant_36columns() {
        String outputVariant;
        outputVariant = join(delimiter,
                sample,
                gene,
                chr,
                startPosition,
                endPosition,
                refAllele,
                varAllele,

                totalCoverage,
                variantCoverage,
                referenceForwardCount,
                referenceReverseCount,
                variantForwardCount,
                variantReverseCount,
                genotype,
                frequency == 0 ? 0 : new DecimalFormat("0.0000").format(frequency),
                bias,
                pmean == 0 ? 0 : new DecimalFormat("0.0").format(pmean),
                pstd,
                qual == 0 ? 0 : new DecimalFormat("0.0").format(qual),
                qstd,

                mapq == 0 ? 0 : new DecimalFormat("0.0").format(mapq),
                qratio == 0 ? 0 : new DecimalFormat("0.000").format(qratio),
                hifreq == 0 ? 0 : new DecimalFormat("0.0000").format(hifreq),
                extrafreq == 0 ? 0 : new DecimalFormat("0.0000").format(extrafreq),

                shift3,
                msi == 0 ? 0 : new DecimalFormat("0.000").format(msi),
                msint,
                nm > 0 ? new DecimalFormat("0.0").format(nm) : 0,
                hicnt,
                hicov,
                leftSequence, rightSequence,
                region,
                varType,
                duprate == 0 ? 0 : new DecimalFormat("0.0").format(duprate),
                sv
        );
        return outputVariant;
    }
}
