package com.astrazeneca.vardict.postprocessmodules;

import com.astrazeneca.vardict.Configuration;
import com.astrazeneca.vardict.data.ReferenceResource;
import com.astrazeneca.vardict.collection.DirectThreadExecutor;
import com.astrazeneca.vardict.data.Reference;
import com.astrazeneca.vardict.data.Region;
import com.astrazeneca.vardict.data.scopedata.*;
import com.astrazeneca.vardict.printers.SomaticOutputVariant;
import com.astrazeneca.vardict.printers.VariantPrinter;
import com.astrazeneca.vardict.variations.Variant;
import com.astrazeneca.vardict.variations.Vars;

import java.util.*;
import java.util.function.BiConsumer;

import static com.astrazeneca.vardict.Utils.printExceptionAndContinue;
import static com.astrazeneca.vardict.data.scopedata.GlobalReadOnlyScope.instance;
import static com.astrazeneca.vardict.data.scopedata.GlobalReadOnlyScope.getMode;
import static com.astrazeneca.vardict.data.Patterns.MINUS_NUM_NUM;
import static com.astrazeneca.vardict.variations.VariationUtils.VarsType.var;
import static com.astrazeneca.vardict.variations.VariationUtils.VarsType.varn;
import static com.astrazeneca.vardict.variations.VariationUtils.*;

/**
 * Class for preparation of variants found in somatic (paired) analysis to the output
 */
public class SomaticPostProcessModule implements BiConsumer<Scope<AlignedVarsData>, Scope<AlignedVarsData>> {
    private static final String STRONG_SOMATIC = "StrongSomatic";
    private static final String SAMPLE_SPECIFIC = "SampleSpecific";
    private static final String DELETION = "Deletion";
    private static final String LIKELY_LOH = "LikelyLOH";
    private static final String GERMLINE = "Germline";
    private static final String STRONG_LOH = "StrongLOH";
    private static final String LIKELY_SOMATIC = "LikelySomatic";
    private static final String AF_DIFF = "AFDiff";
    private static final String FALSE = "FALSE";
    private static final String COMPLEX = "Complex";
    private static final String SNV = "SNV";
    private static final String EMPTY_STRING = "";

    private Integer maxReadLength = 0;

    private VariantPrinter variantPrinter;
    private ReferenceResource referenceResource;

    public SomaticPostProcessModule(ReferenceResource referenceResource, VariantPrinter variantPrinter) {
        this.referenceResource = referenceResource;
        this.variantPrinter = variantPrinter;
    }
    /**
     * Performs analysis of variants from both samples (BAM1 and BAM2).
     * Creates variant output for variants from aligned map and prints them.
     */
    @Override
    public void accept(Scope<AlignedVarsData> scopeFromBam2,
                       Scope<AlignedVarsData> scopeFromBam1) {

        Region region = scopeFromBam1.region;
        Set<String> splice = scopeFromBam1.splice;
        Map<Integer, Vars> variationsFromBam1 = scopeFromBam1.data.alignedVariants;
        Map<Integer, Vars> variationsFromBam2 = scopeFromBam2.data.alignedVariants;

        maxReadLength = Math.max(scopeFromBam1.maxReadLength, scopeFromBam2.maxReadLength);

        Set<Integer> allPositions = new HashSet<>(variationsFromBam1.keySet());
        allPositions.addAll(variationsFromBam2.keySet());
        List<Integer> variantPositions = new ArrayList<>(allPositions);
        Collections.sort(variantPositions);
        int lastPosition = 0;
        for (Integer position : variantPositions) {
            try {
                lastPosition = position;
                if (position < region.start || position > region.end) {
                    continue;
                }
                Vars v1 = variationsFromBam1.get(position);
                Vars v2 = variationsFromBam2.get(position);
                if (v1 == null && v2 == null) { // both samples have no coverage
                    continue;
                }
                if (v1 == null) { // no coverage for sample 1
                    callingForOneSample(v2, true, DELETION, region, splice);
                } else if (v2 == null) { // no coverage for sample 2
                    callingForOneSample(v1, false, SAMPLE_SPECIFIC, region, splice);
                } else { // both samples have coverage
                    callingForBothSamples(position, v1, v2, region, splice);
                }
            } catch (Exception exception) {
                printExceptionAndContinue(exception, "position", String.valueOf(lastPosition), region);
            }
        }
    }

    /**
     * Implements variations analysis from one sample, print out the result.
     * @param variants variations from one BAM
     * @param isFirstCover if the first calling
     * @param varLabel type of variation (LikelyLOH, LikelySomatic, Germline, AFDiff, StrongSomatic)
     * @param region region from BED file
     * @param splice set of strings representing introns in splice
     * */
    void callingForOneSample(Vars variants, boolean isFirstCover, String varLabel, Region region, Set<String> splice) {
        if (variants.variants.isEmpty()) {
            return;
        }
        for (Variant variant : variants.variants) {
            if (variant.refallele.equals(variant.varallele)) {
                continue;
            }
            SomaticOutputVariant outputVariant;
            variant.vartype = variant.varType();
            if (!variant.isGoodVar(variants.referenceVariant, variant.vartype, splice)) {
                continue;
            }
            if (variant.vartype.equals(COMPLEX)) {
                variant.adjComplex();
            }

            if (isFirstCover) {
                outputVariant = new SomaticOutputVariant(variant, variant, null, variant, region, "", variants.sv, varLabel);
                variantPrinter.print(outputVariant);
            } else {
                outputVariant = new SomaticOutputVariant(variant, variant, variant, null, region, variants.sv, "", varLabel);
                variantPrinter.print(outputVariant);
            }
        }
    }

    /**
     * Run analysis of variations from BAM1 and BAM2.
     * @param position position on which analysis is processed
     * @param v1 variations from BAM1 on target position
     * @param v2 variations from BAM2 on target position
     * @param region region from BED file
     * @param splice set of strings representing introns in splice
     */
    void callingForBothSamples(Integer position, Vars v1, Vars v2, Region region, Set<String> splice)  {
        if (v1.variants.isEmpty() && v2.variants.isEmpty()) {
            return;
        }
        if (v1.variants.size() > 0) {
            printVariationsFromFirstSample(position, v1, v2, region, splice);
        } else if (v2.variants.size() > 0) { // sample 1 has only reference
            printVariationsFromSecondSample(position, v1, v2, region, splice);
        }
    }
    /**
     * Analyse and print variations from BAM1 based on variations from BAM2.
     * @param position position on which analysis is processed
     * @param v1 variations from BAM1 on target position
     * @param v2 variations from BAM2 on target position
     * @param region region from BED file
     * @param splice set of strings representing introns in splice
     */
    private void printVariationsFromFirstSample(Integer position, Vars v1, Vars v2, Region region, Set<String> splice){
        int numberOfProcessedVariation = 0;
        while (numberOfProcessedVariation < v1.variants.size()
                && v1.variants.get(numberOfProcessedVariation).isGoodVar(v1.referenceVariant,
                v1.variants.get(numberOfProcessedVariation).varType(), splice)) {
            final Variant vref = v1.variants.get(numberOfProcessedVariation);
            if (vref.refallele.equals(vref.varallele)) {
                numberOfProcessedVariation++;
                continue;
            }
            final String nt = vref.descriptionString;
            vref.vartype = vref.varType();
            SomaticOutputVariant outputVariant;
            if (vref.vartype.equals(COMPLEX)) {
                vref.adjComplex();
            }
            Variant v2nt = getVarMaybe(v2, varn, nt);
            if (v2nt != null) {
                String type = determinateType(v2, vref, v2nt, splice);
                outputVariant = new SomaticOutputVariant(vref, v2nt, vref, v2nt, region, v1.sv, v2.sv, type);
                variantPrinter.print(outputVariant);
            } else { // sample 1 only, should be strong somatic
                Variant varForPrint = new Variant();
                if (!v2.variants.isEmpty()) {
                    Variant v2r = getVarMaybe(v2, var, 0);
                    int tcov = v2r != null && v2r.totalPosCoverage != 0 ? v2r.totalPosCoverage : 0;
                    int rfc = v2r != null && v2r.refForwardCoverage != 0 ? v2r.refForwardCoverage : 0;
                    int rrc = v2r != null && v2r.refReverseCoverage != 0 ? v2r.refReverseCoverage : 0;
                    varForPrint.totalPosCoverage = tcov;
                    varForPrint.refForwardCoverage = rfc;
                    varForPrint.refReverseCoverage = rrc;
                } else if (v2.referenceVariant != null) {
                    varForPrint = v2.referenceVariant;
                } else {
                    varForPrint = null;
                }
                String type = STRONG_SOMATIC;
                jregex.Matcher mm = MINUS_NUM_NUM.matcher(nt);
                if (!vref.vartype.equals(SNV) && (nt.length() > 10 || mm.find())) {
                    v2nt = new Variant();
                    v2.varDescriptionStringToVariants.put(nt, v2nt); // Ensure it's initialized before passing to combineAnalysis
                    if (vref.positionCoverage < instance().conf.minr + 3 && !nt.contains("<")) {
                        CombineAnalysisData tpl = combineAnalysis(vref, v2nt, region.chr, position, nt,
                                splice, maxReadLength);
                        maxReadLength = tpl.maxReadLength;
                        String newtype = tpl.type;
                        if ("FALSE".equals(newtype)) {
                            numberOfProcessedVariation++;
                            continue;
                        }
                        if (newtype.length() > 0) {
                            type = newtype;
                        }
                    }
                }
                if (type.equals(STRONG_SOMATIC)) {
                    outputVariant = new SomaticOutputVariant(vref, vref, vref, varForPrint, region, v1.sv, v2.sv, STRONG_SOMATIC);
                    variantPrinter.print(outputVariant);
                } else {
                    outputVariant = new SomaticOutputVariant(vref, vref, vref, v2nt, region, v1.sv, v2.sv, type);
                    variantPrinter.print(outputVariant);
                }
            }
            numberOfProcessedVariation++;
        }
        if (numberOfProcessedVariation == 0) {
            if (v2.variants.isEmpty()) {
                return;
            }
            for (Variant v2var : v2.variants) {
                SomaticOutputVariant outputVariant;
                v2var.vartype = v2var.varType();
                if (!v2var.isGoodVar(v2.referenceVariant, v2var.vartype, splice)) {
                    continue;
                }
                // potential LOH
                String nt = v2var.descriptionString;
                Variant v1nt = getVarMaybe(v1, varn, nt);
                if (v1nt != null) {
                    if (v1nt.refallele.equals(v1nt.varallele)) {
                        continue;
                    }
                    String type = v1nt.frequency < instance().conf.lofreq ? LIKELY_LOH : GERMLINE;
                    if (COMPLEX.equals(v2var.vartype)) {
                        v1nt.adjComplex();
                    }

                    v1nt.vartype = v1nt.varType();
                    outputVariant = new SomaticOutputVariant(v1nt, v2var, v1nt, v2var, region, v1.sv, v2.sv, type);
                    variantPrinter.print(outputVariant);
                } else {
                    if (v2var.refallele.equals(v2var.varallele)) {
                        continue;
                    }
                    Variant v1var = getVarMaybe(v1, var, 0);
                    int tcov = v1var != null && v1var.totalPosCoverage != 0 ? v1var.totalPosCoverage : 0;

                    Variant v1ref = v1.referenceVariant;
                    int fwd = v1ref != null ? v1ref.varsCountOnForward : 0;
                    int rev = v1ref != null ? v1ref.varsCountOnReverse : 0;

                    String genotype = v1var != null ? v1var.genotype :
                            (v1ref != null ? v1ref.descriptionString + "/" + v1ref.descriptionString : "N/N");

                    if (COMPLEX.equals(v2var.vartype)) {
                        v2var.adjComplex();
                    }

                    Variant varForPrint = new Variant();
                    varForPrint.totalPosCoverage = tcov;
                    varForPrint.refForwardCoverage = fwd;
                    varForPrint.refReverseCoverage = rev;
                    varForPrint.genotype = genotype;

                    outputVariant = new SomaticOutputVariant(v2var, v2var, varForPrint, v2var, region, "", v2.sv, STRONG_LOH);
                    variantPrinter.print(outputVariant);
                }
            }
        }
    }

    /**
     * Analyse and print variations from BAM2 based on variations from BAM1.
     * @param position position on which analysis is processed
     * @param v1 variations from BAM1 on target position
     * @param v2 variations from BAM2 on target position
     * @param region region from BED file
     * @param splice set of strings representing introns in splice
     */
    private void printVariationsFromSecondSample(Integer position, Vars v1, Vars v2, Region region, Set<String> splice){
        for (Variant v2var : v2.variants) {
            if (v2var.refallele.equals(v2var.varallele)) {
                continue;
            }
            v2var.vartype = v2var.varType();
            if (!v2var.isGoodVar(v2.referenceVariant, v2var.vartype, splice)) {
                continue;
            }
            // potential LOH
            String descriptionString = v2var.descriptionString;
            String type = STRONG_LOH;
            Variant v1nt = v1.varDescriptionStringToVariants.computeIfAbsent(descriptionString, k -> new Variant());
            v1nt.positionCoverage = 0;
            String newType = EMPTY_STRING;
            jregex.Matcher mm = MINUS_NUM_NUM.matcher(descriptionString);
            if (v2.varDescriptionStringToVariants.get(descriptionString).positionCoverage < instance().conf.minr + 3
                    && !descriptionString.contains("<") && (descriptionString.length() > 10 || mm.find())) {
                CombineAnalysisData tpl = combineAnalysis(
                        v2.varDescriptionStringToVariants.get(descriptionString),
                        v1nt,
                        region.chr,
                        position,
                        descriptionString,
                        splice,
                        maxReadLength);
                maxReadLength = tpl.maxReadLength;
                newType = tpl.type;
                if (FALSE.equals(newType)) {
                    continue;
                }
            }
            Variant varForPrint;
            if (newType.length() > 0) {
                type = newType;
                varForPrint = v1nt;
            } else {
                Variant v1ref = v1.referenceVariant;
                if (v1ref != null) {
                    varForPrint = v1ref;
                } else {
                    varForPrint = null;
                }
            }
            if (COMPLEX.equals(v2var.vartype)) {
                v2var.adjComplex();
            }

            SomaticOutputVariant outputVariant = new SomaticOutputVariant(v2var, v2var, varForPrint, v2var, region, "", v2.sv, type);
            variantPrinter.print(outputVariant);
        }
    }

    /**
     * Analyse two variations and return their type.
     * @param variants variations from BAM2
     * @param standardVariant a variation to compare with
     * @param variantToCompare a variation to be compared
     * @param splice set of strings representing introns in splice
     * @return type of variation (LikelyLOH, LikelySomatic, Germline, AFDiff, StrongSomatic)
     * */
     String determinateType(Vars variants, Variant standardVariant, Variant variantToCompare, Set<String> splice) {
        String type;
        if (variantToCompare.isGoodVar(variants.referenceVariant, standardVariant.vartype, splice)) {
            if (standardVariant.frequency > (1 - instance().conf.lofreq) && variantToCompare.frequency < 0.8d && variantToCompare.frequency > 0.2d) {
                type = LIKELY_LOH;
            } else {
                if (variantToCompare.frequency < instance().conf.lofreq || variantToCompare.positionCoverage <= 1) {
                    type = LIKELY_SOMATIC;
                } else {
                    type = GERMLINE;
                }
            }
        } else {
            if (variantToCompare.frequency < instance().conf.lofreq || variantToCompare.positionCoverage <= 1) {
                type = LIKELY_SOMATIC;
            } else {
                type = AF_DIFF;
            }
        }
        if (variantToCompare.isNoise() && standardVariant.vartype.equals(SNV)) {
            type = STRONG_SOMATIC;
        }
        return type;
    }


    /**
     * Taken a likely somatic indels and see whether combine two bam files still support somatic status. This is mainly for Indels
     * that softclipping overhang is too short to positively being called in one bam file, but can be called in the other bam file,
     * thus creating false positives
     * @param variant1 variant 1
     * @param variant2 variant 2
     * @param chrName chromosome
     * @param position position
     * @param descriptionString description string of variant
     * @param splice set of strings representing introns in splice
     * @param maxReadLength max read length
     * @return (new <code>maxReadLength</code>, "FALSE" | "")
     */
      CombineAnalysisData combineAnalysis(Variant variant1, Variant variant2,
                                          String chrName, int position,
                                          String descriptionString, Set<String> splice,
                                          int maxReadLength)  {
        Configuration config = instance().conf;
        if (config.y) {
            System.err.printf("Start Combine %s %s\n", position, descriptionString);
        }
        // Don't do it for structural variants
        if (variant1.endPosition - variant1.startPosition > instance().conf.SVMINLEN) {
            return new CombineAnalysisData(maxReadLength, EMPTY_STRING);
        }

        Region region = new Region(chrName, variant1.startPosition - maxReadLength, variant1.endPosition + maxReadLength, "");
        Reference ref = referenceResource.getReference(region);

        Scope<InitialData> currentScope = new Scope<>(config.bam.getBam1() + ":" + config.bam.getBam2(),
                  region, ref, referenceResource, maxReadLength, splice,
                  variantPrinter, new InitialData());
        AlignedVarsData tpl = getMode().pipeline(currentScope, new DirectThreadExecutor()).join().data;

        maxReadLength = tpl.maxReadLength;
        Map<Integer, Vars> vars = tpl.alignedVariants;
        Variant vref = getVarMaybe(vars, position, varn, descriptionString);
        if (vref != null) {
            if (config.y) {
                System.err.printf("Combine: 1: %s comb: %s\n", variant1.positionCoverage, vref.positionCoverage);
            }
            if (vref.positionCoverage - variant1.positionCoverage >= config.minr) {
                variant2.totalPosCoverage = vref.totalPosCoverage - variant1.totalPosCoverage;
                if (variant2.totalPosCoverage < 0)
                    variant2.totalPosCoverage = 0;

                variant2.positionCoverage = vref.positionCoverage - variant1.positionCoverage;
                if (variant2.positionCoverage < 0)
                    variant2.positionCoverage = 0;

                variant2.refForwardCoverage = vref.refForwardCoverage - variant1.refForwardCoverage;
                if (variant2.refForwardCoverage < 0)
                    variant2.refForwardCoverage = 0;

                variant2.refReverseCoverage = vref.refReverseCoverage - variant1.refReverseCoverage;
                if (variant2.refReverseCoverage < 0)
                    variant2.refReverseCoverage = 0;

                variant2.varsCountOnForward = vref.varsCountOnForward - variant1.varsCountOnForward;
                if (variant2.varsCountOnForward < 0)
                    variant2.varsCountOnForward = 0;

                variant2.varsCountOnReverse = vref.varsCountOnReverse - variant1.varsCountOnReverse;
                if (variant2.varsCountOnReverse < 0)
                    variant2.varsCountOnReverse = 0;

                if (variant2.positionCoverage != 0) {
                    variant2.meanPosition = (vref.meanPosition * vref.positionCoverage - variant1.meanPosition * variant1.positionCoverage) / variant2.positionCoverage;
                    variant2.meanQuality = (vref.meanQuality * vref.positionCoverage - variant1.meanQuality * variant1.positionCoverage) / variant2.positionCoverage;
                    variant2.meanMappingQuality = (vref.meanMappingQuality * vref.positionCoverage - variant1.meanMappingQuality * variant1.positionCoverage) / variant2.positionCoverage;
                    variant2.highQualityReadsFrequency = (vref.highQualityReadsFrequency * vref.positionCoverage - variant1.highQualityReadsFrequency * variant1.positionCoverage) / variant2.positionCoverage;
                    variant2.extraFrequency = (vref.extraFrequency * vref.positionCoverage - variant1.extraFrequency * variant1.positionCoverage) / variant2.positionCoverage;
                    variant2.numberOfMismatches = (vref.numberOfMismatches * vref.positionCoverage - variant1.numberOfMismatches * variant1.positionCoverage) / variant2.positionCoverage;
                } else {
                    variant2.meanPosition = 0;
                    variant2.meanQuality = 0;
                    variant2.meanMappingQuality = 0;
                    variant2.highQualityReadsFrequency = 0;
                    variant2.extraFrequency = 0;
                    variant2.numberOfMismatches = 0;
                }
                variant2.isAtLeastAt2Positions = true;
                variant2.hasAtLeast2DiffQualities = true;

                if (variant2.totalPosCoverage <= 0) {
                    return new CombineAnalysisData(maxReadLength, FALSE);
                }

                variant2.frequency = variant2.positionCoverage / (double)variant2.totalPosCoverage;
                variant2.highQualityToLowQualityRatio = variant1.highQualityToLowQualityRatio; // Can't back calculate and should be inaccurate
                variant2.genotype = vref.genotype;
                variant2.strandBiasFlag = strandBias(variant2.refForwardCoverage, variant2.refReverseCoverage) + ";" +
                        strandBias(variant2.varsCountOnForward, variant2.varsCountOnReverse);
                return new CombineAnalysisData(maxReadLength, GERMLINE);
            } else if (vref.positionCoverage < variant1.positionCoverage - 2) {
                if (config.y) {
                    System.err.printf("Combine produce less: %s %s %s %s %s\n", chrName, position, descriptionString, vref.positionCoverage, variant1.positionCoverage);
                }
                return new CombineAnalysisData(maxReadLength, FALSE);
            } else {
                return new CombineAnalysisData(maxReadLength, EMPTY_STRING);
            }
        }
        return new CombineAnalysisData(maxReadLength, FALSE);
    }
}
