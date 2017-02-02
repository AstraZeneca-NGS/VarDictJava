package com.astrazeneca.vardict.pipeline.postprocessmodules;

import com.astrazeneca.utils.DirectThreadExecutor;
import com.astrazeneca.utils.Tuple;
import com.astrazeneca.vardict.*;
import com.astrazeneca.vardict.modes.AbstractMode;
import com.astrazeneca.vardict.pipeline.data.Scope;
import com.astrazeneca.vardict.variations.Variant;
import com.astrazeneca.vardict.variations.Vars;

import java.io.IOException;
import java.io.PrintStream;
import java.util.*;
import java.util.function.BiConsumer;

import static com.astrazeneca.GlobalReadOnlyScope.instance;
import static com.astrazeneca.utils.Tuple.tuple;
import static com.astrazeneca.utils.Utils.join;
import static com.astrazeneca.utils.Utils.round;
import static com.astrazeneca.vardict.variations.VariationUtils.VarsType.var;
import static com.astrazeneca.vardict.variations.VariationUtils.VarsType.varn;
import static com.astrazeneca.vardict.variations.VariationUtils.getVarMaybe;
import static com.astrazeneca.vardict.variations.VariationUtils.isGoodVar;
import static com.astrazeneca.vardict.variations.VariationUtils.strandBias;
import static java.lang.String.format;

/**
 * Class for paired sample variant calling.
 */
public class SomaticPostProcessModule implements BiConsumer<Scope<Tuple.Tuple2<Integer, Map<Integer, Vars>>>,
        Scope<Tuple.Tuple2<Integer, Map<Integer, Vars>>>> {

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
    public static final String EMPTY_STRING = "";

    private Integer maxReadLength = 0;

    /**
     * Performs analysis of variations from both samples (BAM1 and BAM2).
     * @return maximum read length
     * @throws IOException
     */
    //perl version: 368
    @Override
    public void accept(Scope<Tuple.Tuple2<Integer, Map<Integer, Vars>>> scopeFromBam2,
                       Scope<Tuple.Tuple2<Integer, Map<Integer, Vars>>> scopeFromBam1) {

        Region region = scopeFromBam1.region;
        PrintStream out = scopeFromBam1.out;
        Set<String> splice = scopeFromBam1.splice;
        Map<Integer, Vars> variationsFromBam1 = scopeFromBam1.data._2;
        Map<Integer, Vars> variationsFromBam2 = scopeFromBam2.data._2;

        maxReadLength = Math.max(scopeFromBam1.maxReadLength, scopeFromBam2.maxReadLength);

        Set<Integer> allPositions = new HashSet<>(variationsFromBam1.keySet());
        allPositions.addAll(variationsFromBam2.keySet());
        List<Integer> variantPositions = new ArrayList<>(allPositions);
        Collections.sort(variantPositions);

        //perl version: 375
        for (Integer position : variantPositions) {
            Vars v1 = variationsFromBam1.get(position);
            Vars v2 = variationsFromBam2.get(position);
            if (v1 == null && v2 == null) { // both samples have no coverage
                continue;
            }
            if (v1 == null) { // no coverage for sample 1
                callingForOneSample(v2, true, DELETION, region, splice, out);
            } else if (v2 == null) { // no coverage for sample 2
                callingForOneSample(v1, false, SAMPLE_SPECIFIC, region, splice, out);
            } else { // both samples have coverage
                try {
                    callingForBothSamples(position, v1, v2, region, splice, out, variationsFromBam2);
                } catch (IOException e) {
                    throw new RuntimeException(e);
                }
            }
        }
    }

    // perl version: 389
    private void callingForBothSamples(Integer position, Vars v1, Vars v2, Region region, Set<String> splice, PrintStream out,
                                       Map<Integer, Vars> variationsFromBam2) throws IOException {
        if (v1.variants.isEmpty() && v2.variants.isEmpty()) {
            return;
        }
        if (v1.variants.size() > 0) {
            printVariationsFromFirstSample(position, v1, v2, region, splice, out, variationsFromBam2);
        } else if (v2.variants.size() > 0) { // sample 1 has only reference
            printVariationsFromSecondSample(position, v1, v2, region, splice, out);
        }
    }

    /**
     * Analyse variations from BAM2 based on variations from BAM1.
     * @param position position on which analysis is processed
     * @param v1 variations from BAM1 on target position
     * @param v2 variations from BAM2 on target position
     * */
    //perl version: 456
    private void printVariationsFromSecondSample(Integer position, Vars v1, Vars v2, Region region, Set<String> splice, PrintStream out) throws IOException {
        final Variant.Type vartype;
        Variant v2var = v2.variants.get(0);
        vartype = v2var.getType();
        if (!isGoodVar(v2var, v2.referenceVariant, vartype, splice)) {
            return;
        }
        // potential LOH
        String descriptionString = v2var.descriptionString;
        String type = STRONG_LOH;
        Variant v1nt = v1.varDescriptionStringToVariants.computeIfAbsent(descriptionString, k -> new Variant());
        v1nt.positionCoverage = 0;
        String newType = EMPTY_STRING;
        if (v2.varDescriptionStringToVariants.get(descriptionString).positionCoverage < instance().conf.minr + 3) {
            Tuple.Tuple2<Integer, String> tpl = combineAnalysis(
                    v2.varDescriptionStringToVariants.get(descriptionString),
                    v1nt,
                    region.chr,
                    position,
                    descriptionString,
                    splice,
                    maxReadLength);
            maxReadLength = tpl._1;
            newType = tpl._2;
            if (FALSE.equals(newType)) {
                return;
            }
        }
        String th1;
        if (newType.length() > 0) {
            type = newType;
            th1 = join("\t", v1nt.addDelimiter("\t"), format("%.1f", v1nt.numberOfMismatches));
        } else {
            Variant v1ref = v1.referenceVariant;
            th1 = join("\t", v1ref.addDelimiter("\t"), format("%.1f", v1ref.numberOfMismatches));
        }
        if (vartype == Variant.Type.complex) {
            v2var.adjComplex();
        }
        String info = join("\t", th1, v2var.addDelimiter("\t"), format("%.1f", v2var.numberOfMismatches));
        out.println(constructVariationReportString(type, vartype, v2var, info, region));
    }

    /**
     * Analyse variations from BAM1 based on variations from BAM2.
     * @param position position on which analysis is processed
     * @param v1 variations from BAM1 on target position
     * @param v2 variations from BAM2 on target position
     */
    //perl version: 392
    private void printVariationsFromFirstSample(Integer position, Vars v1, Vars v2, Region region, Set<String> splice, PrintStream out,
                                                Map<Integer, Vars> variationsFromBam2) throws IOException {
        Variant.Type vartype;
        int numberOfProcessedVariation = 0;
        for (Variant vref : v1.variants) {
            //perl version: 409
            if (!isGoodVar(vref, v1.referenceVariant, vref.getType(), splice)) {
                break;
            }
            final String nt = vref.descriptionString;
            vartype = vref.getType();
            if (vartype == Variant.Type.complex) {
                vref.adjComplex();
            }
            Variant v2nt = getVarMaybe(v2, varn, nt);
            if (v2nt != null) {
                String type = determinateType(v2, vartype, vref, v2nt, splice);
                String info = join("\t",
                        vref.addDelimiter("\t"), format("%.1f", vref.numberOfMismatches),
                        v2nt.addDelimiter("\t"), format("%.1f", v2nt.numberOfMismatches)
                );
                out.println(constructVariationReportString(type, vartype, vref, v2nt, info, region));
            } else { // sample 1 only, should be strong somatic
                String type = STRONG_SOMATIC;
                v2nt = new Variant();
                v2.varDescriptionStringToVariants.put(nt, v2nt); // Ensure it's initialized before passing to combineAnalysis
                if (vref.positionCoverage < instance().conf.minr + 3) {
                    Tuple.Tuple2<Integer, String> tpl = combineAnalysis(
                            vref,
                            v2nt,
                            region.chr,
                            position,
                            nt,
                            splice,
                            maxReadLength
                    );
                    maxReadLength = tpl._1;
                    String newtype = tpl._2;
                    if (FALSE.equals(newtype)) {
                        numberOfProcessedVariation++;
                        continue;
                    }
                    if (newtype.length() > 0) {
                        type = newtype;
                    }
                }
                if (type.equals(STRONG_SOMATIC)) {
                    String tvf;
                    if (v2.referenceVariant == null) {
                        Variant v2m = getVarMaybe(v2, var, 0);
                        int tcov = v2m != null && v2m.totalPosCoverage != 0 ? v2m.totalPosCoverage : 0;
                        tvf = join("\t", tcov, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
                    } else {
                        Variant v2ref = variationsFromBam2.get(position).referenceVariant;
                        tvf = join("\t", v2ref.addDelimiter("\t"), format("%.1f", v2ref.numberOfMismatches));
                    }

                    String info = join("\t", vref.addDelimiter("\t"), format("%.1f", vref.numberOfMismatches), tvf);
                    out.println(constructVariationReportString(type, vartype, vref, vref, info, region));
                } else {
                    String info = join("\t",
                            vref.addDelimiter("\t"), format("%.1f", vref.numberOfMismatches),
                            v2nt.addDelimiter("\t"), format("%.1f", v2nt.numberOfMismatches)
                    );
                    out.println(constructVariationReportString(type, vartype, vref, info, region));
                }
            }
            numberOfProcessedVariation++;
        }

        if (numberOfProcessedVariation == 0) {
            //perl version: 440
            if (v2.variants.isEmpty()) {
                return;
            }
            Variant v2var = v2.variants.get(0);
            vartype = v2var.getType();
            if (!isGoodVar(v2var, v2.referenceVariant, vartype, splice)) {
                return;
            }
            // potentail LOH
            Variant v1nt = getVarMaybe(v1, varn, v2var.descriptionString);
            if (v1nt != null) {
                String type = v1nt.frequency < instance().conf.lofreq ? LIKELY_LOH : GERMLINE;
                if (vartype == Variant.Type.complex) {
                    v1nt.adjComplex();
                }
                String info = join("\t",
                        v1nt.addDelimiter("\t"), format("%.1f", v1nt.numberOfMismatches),
                        v2var.addDelimiter("\t"), format("%.1f", v2var.numberOfMismatches)
                );
                out.println(constructVariationReportString(type, vartype, v1nt, v2var, info, region));
            } else {
                String th1;
                Variant v1ref = v1.referenceVariant;
                if (v1ref != null) {
                    th1 = join("\t", v1ref.addDelimiter("\t"), format("%.1f", v1ref.numberOfMismatches));
                } else {
                    th1 = join("\t", v1.variants.get(0).totalPosCoverage, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
                }
                if (vartype == Variant.Type.complex) {
                    v2var.adjComplex();
                }
                String info = join("\t", th1, v2var.addDelimiter("\t"), format("%.1f", v2var.numberOfMismatches));
                out.println(constructVariationReportString(STRONG_LOH, vartype, v2var, v2var, info, region));
            }
        }
    }

    /**
     * Analyse two variations and return their type.
     * @param variants variations fro BAM2
     * @param vartype variation type (SNV, Insertion, Deletion, Complex)
     * @param etalonVariant a variation to compare with
     * @param variantToCompare a variation to be compared
     * @return type of variation (LikelyLOH, LikelySomatic, Germline, AFDiff, StrongSomatic)
     * */
    //perl version: 413
    private String determinateType(Vars variants, Variant.Type vartype, Variant etalonVariant, Variant variantToCompare, Set<String> splice) {
        String type;
        final Configuration config = instance().conf;
        if (isGoodVar(variantToCompare, variants.referenceVariant, vartype, splice)) {
            if (etalonVariant.frequency > (1 - config.lofreq) && variantToCompare.frequency < 0.8d && variantToCompare.frequency > 0.2d) {
                type = LIKELY_LOH;
            } else {
                if (variantToCompare.frequency < config.lofreq || variantToCompare.positionCoverage <= 1) {
                    type = LIKELY_SOMATIC;
                } else {
                    type = GERMLINE;
                }
            }
        } else {
            if (variantToCompare.frequency < config.lofreq || variantToCompare.positionCoverage <= 1) {
                type = LIKELY_SOMATIC;
            } else {
                type = AF_DIFF;
            }
        }
        if (variantToCompare.isNoise(config.goodq, config.lofreq) && vartype == Variant.Type.SNV) {
            type = STRONG_SOMATIC;
        }
        return type;
    }

    /**
     * Implements variations analysis from one sample, print out the result.
     * @param variants variations from one BAM
     * @param isFirstCover if the first calling
     * @param varLabel type of variation (LikelyLOH, LikelySomatic, Germline, AFDiff, StrongSomatic)
     * */
    //perl version: 379, 384
    private void callingForOneSample(Vars variants, boolean isFirstCover, String varLabel, Region region, Set<String> splice, PrintStream out) {
        Variant.Type variationType;
        if (variants.variants.isEmpty()) {
            return;
        }
        Variant var = variants.variants.get(0);
        variationType = var.getType();
        if (!isGoodVar(var, variants.referenceVariant, variationType, splice)) {
            return;
        }
        if (variationType == Variant.Type.complex) {
            var.adjComplex();
        }

        Object[] infoToReport = isFirstCover
                ? new Object[] {
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                        var.addDelimiter("\t"), format("%.1f", var.numberOfMismatches)
                }
                : new Object[] {
                        var.addDelimiter("\t"), format("%.1f", var.numberOfMismatches),
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
                };
            String info = join("\t", infoToReport);
        out.println(constructVariationReportString(varLabel, variationType, var, info, region));
    }

    /**
     * Prepare to print the result of paired sample variant calling.
     * @param varLabel type of variation (LikelyLOH, LikelySomatic, Germline, AFDiff, StrongSomatic)
     * @param vartype variation type (SNV, Insertion, Deletion, Complex)
     * @param variant variation which info will be printed out
     * @param info some additional information regarding variation
     * @return String to print out.
     */
    private String constructVariationReportString(String varLabel, Variant.Type vartype, Variant variant, String info,
                                                  Region region) {
        return constructVariationReportString(varLabel, vartype, variant, variant, info, region);
    }

    /**
     * Prepare to print the result of paired sample variant calling.
     * @param varLabel type of variation (LikelyLOH, LikelySomatic, Germline, AFDiff, StrongSomatic)
     * @param vartype variation type (SNV, Insertion, Deletion, Complex)
     * @param variant1 variation which info will be printed out
     * @param variant2 variation which info will be printed out
     * @param info some additional information regarding variation
     * @return String to print out.
     */
    private String constructVariationReportString(String varLabel, Variant.Type vartype,
                                                  Variant variant1, Variant variant2, String info, Region region) {
        return join("\t",
                instance().sample, region.gene, region.chr,
                variant1.startPosition, variant1.endPosition, variant1.refAllele, variant1.varAllele,
                info,
                variant2.shift3, variant2.msi == 0 ? 0 : format("%.3f", variant2.msi), variant2.msint, variant2.precedingRefSequence, variant2.followingRefSequence,
                region.chr + ":" + region.start + "-" + region.end, varLabel, vartype
        );
    }

    /**
     * Taken a likely somatic indels and see whether combine two bam files still support somatic status. This is mainly for Indels
     * that softclipping overhang is too short to positively being called in one bam file, but can be called in the other bam file,
     * thus creating false positives
     *
     * @param variant1 variant 1
     * @param variant2 variant 2
     * @param chrName chromosome
     * @param position position
     * @param nucleotide nucleotide
     *
     * @param splice set of strings representing introns in splice
     * @param maxReadLength max read length
     * @return (new <code>maxReadLength</code>, "FALSE" | "")
     * @throws IOException
     */
    //perl version: 476
    static Tuple.Tuple2<Integer, String> combineAnalysis(Variant variant1,
                                                         Variant variant2,
                                                         String chrName,
                                                         int position,
                                                         String nucleotide,
                                                         Set<String> splice,
                                                         int maxReadLength) throws IOException {
        Configuration config = instance().conf;
        Map<String, Integer> chrsLengths = instance().chrLengths;
        if (config.y) {
            System.err.printf("Start Combine %s %s\n", position, nucleotide);
        }
        Region region = new Region(chrName, variant1.startPosition - maxReadLength, variant1.endPosition + maxReadLength, "");
        Map<Integer, Character> ref = ReferenceResource.getReference(region);

        Tuple.Tuple2<Integer, Map<Integer, Vars>> tpl = AbstractMode.pipeline(
                config.bam.getBam1() + ":" + config.bam.getBam2(),
                region,
                ref,
                maxReadLength,
                splice,
                instance().out,
                new DirectThreadExecutor()
        ).join().data;

        maxReadLength = tpl._1;
        Map<Integer, Vars> vars = tpl._2;
        Variant vref = getVarMaybe(vars, position, varn, nucleotide);
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

                variant2.meanPosition = (vref.meanPosition * vref.positionCoverage - variant1.meanPosition * variant1.positionCoverage) / variant2.positionCoverage;
                variant2.meanQuality = (vref.meanQuality * vref.positionCoverage - variant1.meanQuality * variant1.positionCoverage) / variant2.positionCoverage;
                variant2.meanMappingQuality = (vref.meanMappingQuality * vref.positionCoverage - variant1.meanMappingQuality * variant1.positionCoverage) / variant2.positionCoverage;
                variant2.highQualityReadsFrequency = (vref.highQualityReadsFrequency * vref.positionCoverage - variant1.highQualityReadsFrequency * variant1.positionCoverage) / variant2.positionCoverage;
                variant2.extraFrequency = (vref.extraFrequency * vref.positionCoverage - variant1.extraFrequency * variant1.positionCoverage) / variant2.positionCoverage;
                variant2.numberOfMismatches = (vref.numberOfMismatches * vref.positionCoverage - variant1.numberOfMismatches * variant1.positionCoverage) / variant2.positionCoverage;

                variant2.isAtLeastAt2Positions = true;
                variant2.hasAtLeast2DiffQualities = true;

                if (variant2.totalPosCoverage <= 0) {
                    return tuple(maxReadLength, FALSE);
                }

                variant2.frequency = round(variant2.positionCoverage / (double)variant2.totalPosCoverage, 4);
                variant2.highQualityToLowQualityRatio = variant1.highQualityToLowQualityRatio; // Can't back calculate and should be inaccurate
                variant2.genotype = vref.genotype;
                variant2.strandBiasFlag = strandBias(variant2.refForwardCoverage, variant2.refReverseCoverage, config.bias, config.minb) + ";" + strandBias(variant2.varsCountOnForward, variant2.varsCountOnReverse, config.bias, config.minb);
                return tuple(maxReadLength, GERMLINE);
            } else if (vref.positionCoverage < variant1.positionCoverage - 2) {
                if (config.y) {
                    System.err.printf("Combine produce less: %s %s %s %s %s\n", chrName, position, nucleotide, vref.positionCoverage, variant1.positionCoverage);
                }
                return tuple(maxReadLength, FALSE);
            } else {
                return tuple(maxReadLength, EMPTY_STRING);
            }

        }
        return tuple(maxReadLength, FALSE);
    }
}
