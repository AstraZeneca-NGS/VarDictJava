package com.astrazeneca.vardict.printers;

import com.astrazeneca.vardict.collection.Tuple;
import com.astrazeneca.vardict.data.Region;
import com.astrazeneca.vardict.data.fishertest.FisherExact;
import com.astrazeneca.vardict.variations.Variant;

import java.text.DecimalFormat;
import java.util.List;

import static com.astrazeneca.vardict.Utils.getRoundedValueToPrint;
import static com.astrazeneca.vardict.data.scopedata.GlobalReadOnlyScope.instance;
import static com.astrazeneca.vardict.Utils.join;

/**
 * Variant created in Amplicon mode. Must contain 38 total fields if no --fisher.
 */
public class AmpliconOutputVariant extends OutputVariant {
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
    private int goodVariantsCount;
    private int totalVariantsCount;
    private int noCoverage;
    private int ampliconFlag;

    private double pvalue;
    private String oddratio = "0";

    public AmpliconOutputVariant(Variant variant, Region region, List<Tuple.Tuple2<Variant, String>> goodVariants,
                                 List<Tuple.Tuple2<Variant, String>> badVariants, int position, int gvscnt,
                                 int noCov, boolean flag) {
        this.sample = instance().sample;
        this.gene = region.gene;
        this.chr = region.chr;

        if (variant == null) {
            this.startPosition = position;
            this.endPosition = position;
            this.region = region.chr + ":" + position + "-" + position;
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
            this.goodVariantsCount = gvscnt;
            this.totalVariantsCount = gvscnt + badVariants.size();
            this.noCoverage = noCov;
            this.ampliconFlag = flag ? 1 : 0;
            if (!goodVariants.isEmpty()) {
                this.region = goodVariants.get(0)._2;
            } else {
                this.region = region.chr + ":" + position + "-" + position;
            }

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

        if (instance().conf.debug && variant != null) {
            StringBuilder goodOutput = new StringBuilder();
            StringBuilder badOutput = new StringBuilder();
            for (int gvi = 0; gvi < goodVariants.size(); gvi++) {
                Tuple.Tuple2<Variant, String> tp = goodVariants.get(gvi);
                goodOutput.append("\tGood" + gvi + " " + join(" ", debugAmpVariant(" ", tp._1), tp._2));
            }
            for (int bvi = 0; bvi < badVariants.size(); bvi++) {
                Tuple.Tuple2<Variant, String> tp = badVariants.get(bvi);
                badOutput.append("\tBad" + bvi + " " + join(" ", debugAmpVariant(" ", tp._1), tp._2));
            }
            this.DEBUG = variant.DEBUG + goodOutput.append(badOutput).toString();
        }
    }

    @Override
    public String toString() {
        String outputVariant;
        if (!instance().conf.fisher) {
            outputVariant = create_amplicon_variant_38columns();
        } else {
            outputVariant = create_amplicon_variant_40columns();
        }
        if (instance().conf.debug) {
            outputVariant = join(delimiter, outputVariant, DEBUG);
        }
        return outputVariant;
    }

    /**
     * 40 columns: oddratio and pvalue added and format as it will be after using R script
     */
    private String create_amplicon_variant_40columns() {
        String outputVariant;
        String hifreq_f = hifreq == 0
                ? "0"
                : new DecimalFormat("0.0000").format(hifreq);
        nm = nm > 0 ? nm : 0;
        String nm_f = nm == 0
                ? "0"
                : new DecimalFormat("0.0").format(nm);

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
                getRoundedValueToPrint("0.0000", frequency),
                bias,
                getRoundedValueToPrint("0.0", pmean),
                pstd,
                getRoundedValueToPrint("0.0", qual),
                qstd,

                getRoundedValueToPrint("0.00000", pvalue),
                oddratio,

                getRoundedValueToPrint("0.0", mapq),
                getRoundedValueToPrint("0.000", qratio),
                hifreq_f,
                getRoundedValueToPrint("0.0000", extrafreq),

                shift3,
                getRoundedValueToPrint("0.000", msi),
                msint,
                nm_f,
                hicnt,
                hicov,
                leftSequence, rightSequence,
                region,
                varType,
                goodVariantsCount,
                totalVariantsCount,
                noCoverage,
                ampliconFlag
        );
        return outputVariant;
    }

    /**
     * 38 columns: no oddratio and pvalue columns
     */
    private String create_amplicon_variant_38columns() {
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
                goodVariantsCount,
                totalVariantsCount,
                noCoverage,
                ampliconFlag
        );
        return outputVariant;
    }

    /**
     * The line will appear in output if -D (debug mode) is enabled. It contains information about variant on each amplicon
     * @param delm delimiter to split the data, default is space.
     * @param variant initial variant to print
     * @return variant information
     */
    public String debugAmpVariant(String delm, Variant variant) {
        return join(delm,
                variant.totalPosCoverage,
                variant.positionCoverage,
                variant.refForwardCoverage,
                variant.refReverseCoverage,
                variant.varsCountOnForward,
                variant.varsCountOnReverse,
                variant.genotype == null ? "0" : variant.genotype,
                variant.frequency == 0 ? 0 : new DecimalFormat("0.0000").format(variant.frequency),
                variant.strandBiasFlag,
                new DecimalFormat("0.0").format(variant.meanPosition),
                variant.isAtLeastAt2Positions ? 1 : 0,
                new DecimalFormat("0.0").format(variant.meanQuality),
                variant.hasAtLeast2DiffQualities ? 1 : 0,
                new DecimalFormat("0.0").format(variant.meanMappingQuality),
                new DecimalFormat("0.000").format(variant.highQualityToLowQualityRatio),
                variant.highQualityReadsFrequency == 0 ? 0 : new DecimalFormat("0.0000").format(variant.highQualityReadsFrequency),
                variant.extraFrequency == 0 ? 0 : new DecimalFormat("0.0000").format(variant.extraFrequency)
        );
    }
}
