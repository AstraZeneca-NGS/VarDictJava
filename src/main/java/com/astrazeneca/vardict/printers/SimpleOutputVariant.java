package com.astrazeneca.vardict.printers;

import com.astrazeneca.vardict.data.Region;
import com.astrazeneca.vardict.variations.Variant;

import java.text.DecimalFormat;

import static com.astrazeneca.vardict.data.scopedata.GlobalReadOnlyScope.instance;
import static com.astrazeneca.vardict.Utils.join;

public class SimpleOutputVariant extends OutputVariant {
    public int totalCoverage;
    public int variantCoverage;
    public int referenceForwardCount;
    public int referenceReverseCount;
    public int variantForwardCount;
    public int variantReverseCount;
    public String genotype = "";
    public double frequency;
    public String bias = "0;0";
    public double pmean;
    public int pstd;
    public double qual;
    public int qstd;
    public double mapq;
    public double qratio;
    public double hifreq;
    public double extrafreq;

    public int shift3;
    public double msi;
    public int msint;
    public double nm;
    public int hicnt;
    public int hicov;

    public double duprate;
    public String sv;

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
            this.DEBUG = variant.DEBUG;
        }
        this.region = region.chr + ":" + region.start + "-" + region.end;
        this.sv = sv.equals("") ? "0" : sv;
    }

    @Override
    public String outputString(String delimiter) {
        this.delimiter = delimiter;
        // 36 columns
        String outputVariant = join(delimiter,
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
        if (instance().conf.debug) {
            outputVariant = join(delimiter, outputVariant, DEBUG);
        }
        return outputVariant;
    }
}
