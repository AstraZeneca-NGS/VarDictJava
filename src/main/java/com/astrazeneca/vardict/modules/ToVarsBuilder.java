package com.astrazeneca.vardict.modules;

import com.astrazeneca.vardict.Configuration;
import com.astrazeneca.vardict.Utils;
import com.astrazeneca.vardict.collection.VariationMap;
import com.astrazeneca.vardict.data.scopedata.AlignedVarsData;
import com.astrazeneca.vardict.data.scopedata.RealignedVariationData;
import com.astrazeneca.vardict.data.Region;
import com.astrazeneca.vardict.data.scopedata.Scope;
import com.astrazeneca.vardict.variations.Variant;
import com.astrazeneca.vardict.variations.Variation;
import com.astrazeneca.vardict.variations.Vars;

import java.time.LocalDateTime;
import java.util.regex.Pattern;
import java.util.regex.Matcher;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static com.astrazeneca.vardict.data.scopedata.GlobalReadOnlyScope.instance;
import static com.astrazeneca.vardict.Utils.*;
import static com.astrazeneca.vardict.data.Patterns.*;
import static com.astrazeneca.vardict.variations.VariationUtils.*;
import static java.lang.Math.abs;

/**
 * This step of pipeline takes already realigned variations (as an intermediate structures) and create Aligned
 * variants data contains information about reference variants, maps of vars on positions and list of description string
 * and variants in region.
 */
public class ToVarsBuilder implements Module<RealignedVariationData, AlignedVarsData> {
    // Left 20 and right 20 bases in reference to output in VCF file
    public static final int REF_20_BASES = 20;
    // Number of bases where to check for MSIs
    public static final int REF_30_BASES = 30;
    public static final int REF_50_BASES = 50;
    public static final int REF_70_BASES = 70;

    private Region region;
    private Map<Integer, Integer> refCoverage;
    private Map<Integer, VariationMap<String, Variation>> insertionVariants;
    private Map<Integer, VariationMap<String, Variation>> nonInsertionVariants;
    private Map<Integer, Character> ref;
    private Double duprate;

    // Map of IUPAC ambiguity codes that we can observe in reference.
    // By VCF 4.3 specification they aren't allowed and must be reduced to a first alphabetically base.
    private Map<String, String> IUPAC_AMBIGUITY_CODES = Stream.of(new String[][] {
            {"M","A"},
            {"R","A"},
            {"W","A"},
            {"S","C"},
            {"Y","C"},
            {"K","G"},
            {"V","A"},
            {"H","A"},
            {"D","A"},
            {"B","C"},
    }).collect(Collectors.toMap(key -> key[0], key -> key[1]));

    public Map<Integer, VariationMap<String, Variation>> getInsertionVariants() {
        return insertionVariants;
    }

    public Map<Integer, VariationMap<String, Variation>> getNonInsertionVariants() {
        return nonInsertionVariants;
    }

    private void initFromScope(Scope<RealignedVariationData> scope) {
        this.ref = scope.regionRef.referenceSequences;
        this.region = scope.region;
        this.refCoverage = scope.data.refCoverage;
        this.insertionVariants = scope.data.insertionVariants;
        this.nonInsertionVariants = scope.data.nonInsertionVariants;
        this.duprate = scope.data.duprate;
    }

    /**
     * Collect variants from variations and fill map of aligned variants on position
     * @param scope realigned variation data to process (contains filled insertionVariations and non-insertion variations maps)
     * @return object contains max read lengths and map of positions-Vars
     */
    @Override
    public Scope<AlignedVarsData> process(Scope<RealignedVariationData> scope) {
        initFromScope(scope);
        Configuration config = instance().conf;

        if (config.y) {
            System.err.printf("Current segment: %s:%d-%d \n", region.chr, region.start, region.end);
        }
        //the variant structure
        Map<Integer, Vars> alignedVariants = new HashMap<>();
        int lastPosition = 0;
        //Loop over positions
        for (Map.Entry<Integer, VariationMap<String, Variation>> entH : getNonInsertionVariants().entrySet()) {
            try {
                final int position = entH.getKey();
                lastPosition = position;
                VariationMap<String, Variation> varsAtCurPosition = entH.getValue();

                // skip position if there are no variants on position (both insertion and non-insertion)
                if (varsAtCurPosition.isEmpty() && !getInsertionVariants().containsKey(position)) {
                    continue;
                }

                //Skip if there are no structural variants on position or if the delete duplication option is on
                if (varsAtCurPosition.sv == null || instance().conf.deleteDuplicateVariants) {
                    //Skip if start position is outside region of interest
                    if (position < region.start || position > region.end) {
                        continue;
                    }
                }

                //Skip position if it has no coverage (except SVs)
                if (varsAtCurPosition.sv == null && !refCoverage.containsKey(position)) {
                    continue;
                }

                if (isTheSameVariationOnRef(position, varsAtCurPosition)) {
                    continue;
                }

                if (!refCoverage.containsKey(position) || refCoverage.get(position) == 0) { // ignore when there's no coverage
                    System.err.printf("Error tcov: %s %d %d %d %s\n",
                            region.chr, position, region.start, region.end, varsAtCurPosition.sv.type);
                    continue;
                }

                //total position coverage
                int totalPosCoverage = refCoverage.get(position);

                //position coverage by high-quality reads
                final int hicov = calcHicov(getInsertionVariants().get(position), varsAtCurPosition);

                //array of all variants for the position
                List<Variant> var = new ArrayList<>();
                List<String> keys = new ArrayList<>(varsAtCurPosition.keySet());
                Collections.sort(keys);

                //temporary array used for debugging
                List<String> debugLines = new ArrayList<>();

                createVariant(duprate, alignedVariants, position, varsAtCurPosition, totalPosCoverage, var, debugLines, keys, hicov);
                totalPosCoverage = createInsertion(duprate, position, totalPosCoverage, var, debugLines, hicov);
                sortVariants(var);

                double maxfreq = collectVarsAtPosition(alignedVariants, position, var);

                if (!config.doPileup && maxfreq <= config.freq && instance().ampliconBasedCalling == null) {
                    if (!config.bam.hasBam2()) {
                        alignedVariants.remove(position);
                        continue;
                    }
                }
                //if reference variant has frequency more than $FREQ, set genotype1 to reference variant
                //else if non-reference variant exists, set genotype1 to first (largest quality*coverage) variant
                //otherwise set genotype1 to reference variant

                Vars variationsAtPos = getOrPutVars(alignedVariants, position);

                collectReferenceVariants(position, totalPosCoverage, variationsAtPos, debugLines);
            } catch (Exception exception) {
                printExceptionAndContinue(exception, "position", String.valueOf(lastPosition), region);
            }
        }

        if (config.y) {
            System.err.println("TIME: Finish preparing vars:" + LocalDateTime.now());
        }

        return new Scope<>(
                scope,
                new AlignedVarsData(scope.maxReadLength, alignedVariants));
    }

    /**
     * Checks if there is only reference variant on position and pileup/amplicon mode/somatic mode
     * @param position start position of the variant
     * @param varsAtCurPosition map of description strings on variations
     * @return true if no pileup/amplicon mode/somatic mode enabled and only reference variant is on position
     */
    private boolean isTheSameVariationOnRef(int position, VariationMap<String, Variation> varsAtCurPosition) {
        Set<String> vk = new HashSet<>(varsAtCurPosition.keySet());
        if (getInsertionVariants().containsKey(position)) {
            vk.add("I");
        }
        if (vk.size() == 1 && ref.containsKey(position) && vk.contains(ref.get(position).toString())) {
            // ignore if only reference were seen and no pileup to avoid computation
            if (!instance().conf.doPileup && !instance().conf.bam.hasBam2() && instance().ampliconBasedCalling == null) {
                return true;
            }
        }
        return false;
    }

    /**
     * For deletion variant calculates shift3 (number of bases to be shifted to 3' for deletions due to alternative alignment),
     * msi (which indicates Microsatellite instability) and msint (MicroSattelite unit length in base pairs) fields.
     * @param position position to seek in reference to get left sequence of variant
     * @param dellen length of deletion part to get them from reference
     * @return MSI object.
     */
    private MSI proceedVrefIsDeletion(int position, int dellen) {
        //left 70 bases in reference sequence
        String leftseq = joinRef(ref, (Math.max(position - REF_70_BASES, 1)), position - 1); // left 10 nt
        int chr0 = getOrElse(instance().chrLengths, region.chr, 0);
        //right 70 + dellen bases in reference sequence
        String tseq = joinRef(ref, position, Math.min(position + dellen + REF_70_BASES, chr0));

        //Try to adjust for microsatellite instability
        MSI msiData = findMSI(substr(tseq, 0, dellen), substr(tseq, dellen), leftseq);
        double msi = msiData.msi;
        int shift3 = msiData.shift3;
        String msint = msiData.msint;

        MSI msiDataWIthoutLeft = findMSI(leftseq, substr(tseq, dellen), null);
        double tmsi = msiDataWIthoutLeft.msi;
        String tmsint = msiDataWIthoutLeft.msint;
        if (msi < tmsi) {
            msi = tmsi;
            // Don't change shift3
            msint = tmsint;
        }
        if (msi <= shift3 / (double) dellen) {
            msi = shift3 / (double) dellen;
        }
        return new MSI(msi, shift3, msint);
    }

    /**
     * For insertion variant calculates shift3 (number of bases to be shifted to 3' for deletions due to alternative alignment),
     * msi (which indicates Microsatellite instability) and msint (MicroSatellite unit length in base pairs) fields.
     * @param position position to seek in reference to get left sequence of variant
     * @param vn variant description string
     * @return tuple of msi, shift3 and msint.
     */
    private MSI proceedVrefIsInsertion(int position, String vn) {
        //variant description string without first symbol '+'
        String tseq1 = vn.substring(1);
        //left 50 bases in reference sequence
        String leftseq = joinRef(ref, Math.max(position - REF_50_BASES, 1), position); // left 10 nt
        int x = getOrElse(instance().chrLengths, region.chr, 0);
        //right 70 bases in reference sequence
        String tseq2 = joinRef(ref, position + 1, (Math.min(position + REF_70_BASES, x)));

        MSI msiData = findMSI(tseq1, tseq2, leftseq);
        double msi = msiData.msi;
        int shift3 = msiData.shift3;
        String msint = msiData.msint;

        //Try to adjust for microsatellite instability
        MSI msiDataWithoutLeft = findMSI(leftseq, tseq2, null);
        double tmsi = msiDataWithoutLeft.msi;
        String tmsint = msiDataWithoutLeft.msint;
        if (msi < tmsi) {
            msi = tmsi;
            // Don't change shift3
            msint = tmsint;
        }
        if (msi <= shift3 / (double)tseq1.length()) {
            msi = shift3 / (double)tseq1.length();
        }
        return new MSI(msi, shift3, msint);
    }

    /**
     * Collect variants to map of position to Vars and fill the data about reference variant and description
     * strings on position
     * @param alignedVariants map to be filled
     * @param position position to be put in map
     * @param var list of variants to be sorted in place
     * @return maxfreq on position
     */
    private double collectVarsAtPosition(Map<Integer, Vars> alignedVariants, int position, List<Variant> var) {
        double maxfreq = 0;
        for (Variant tvar : var) {
            //If variant description string is 1-char base and it matches reference base at this position
            if (tvar.descriptionString.equals(String.valueOf(ref.get(position)))) {
                //this is a reference variant
                getOrPutVars(alignedVariants, position).referenceVariant = tvar;
            } else {
                //append variant to VAR and put it to VARN with key tvar.n (variant description string)
                getOrPutVars(alignedVariants, position).variants.add(tvar);
                getOrPutVars(alignedVariants, position).varDescriptionStringToVariants.put(tvar.descriptionString, tvar);
                if (tvar.frequency > maxfreq) {
                    maxfreq = tvar.frequency;
                }
            }
        }
        return maxfreq;
    }

    /**
     * Sort the list of variants by meanQuality product positionCoverage in descending order
     * @param var list of variants to be sorted in place
     */
    private void sortVariants(List<Variant> var) {
        //sort variants by product of quality and coverage
        Collections.sort(var, new Comparator<Variant>() {
            @Override
            public int compare(Variant o1, Variant o2) {
                int compare = Double.compare(o2.meanQuality * o2.positionCoverage, o1.meanQuality * o1.positionCoverage);
                if (compare != 0) {
                    return compare;
                }
                return o1.descriptionString.compareTo(o2.descriptionString);
            }
        });
    }

    /**
     * Creates insertion variant from each variation on this position and adds it to the list of aligned variants.
     * @param duprate duplication rate of insertion
     * @param position position of the insertion start
     * @param totalPosCoverage current total position coverage (total depth) that can be updated by extra counts
     * @param var list of variants at the position to be updated
     * @param debugLines list of debug lines to be updated
     * @param hicov position coverage by high quality reads
     * @return updated total position coverage
     */
    int createInsertion(double duprate, int position, int totalPosCoverage,
                        List<Variant> var, List<String> debugLines, int hicov) {
        //Handle insertions separately
        Map<String, Variation> insertionVariations = getInsertionVariants().get(position);
        if (insertionVariations != null) {
            List<String> insertionDescriptionStrings = new ArrayList<>(insertionVariations.keySet());
            Collections.sort(insertionDescriptionStrings);
            //Loop over insertion variants
            for (String descriptionString : insertionDescriptionStrings) {
                if (descriptionString.contains("&") && refCoverage.containsKey(position + 1)) {
                    totalPosCoverage = refCoverage.get(position + 1);
                }
                // String n = entV.getKey();
                Variation cnt = insertionVariations.get(descriptionString);
                //count of variants in forward strand
                int fwd = cnt.getDir(false);
                //count of variants in reverse strand
                int rev = cnt.getDir(true);
                //strand bias flag (0, 1 or 2)
                int bias = strandBias(fwd, rev);
                //mean base quality for variant
                double vqual = roundHalfEven("0.0", cnt.meanQuality / cnt.varsCount); // base quality
                //mean mapping quality for variant
                double mq = roundHalfEven("0.0", cnt.meanMappingQuality / (double)cnt.varsCount); // mapping quality
                //number of high-quality reads for variant
                int hicnt = cnt.highQualityReadsCount;
                //number of low-quality reads for variant
                int locnt = cnt.lowQualityReadsCount;

                // Also commented in Perl
                // hicov += hicnt;

                //adjust position coverage if variant count is more than position coverage and no more than
                // position coverage + extracnt
                int ttcov = totalPosCoverage;
                if (cnt.varsCount > totalPosCoverage && cnt.extracnt != 0 &&cnt.varsCount - totalPosCoverage < cnt.extracnt) {
                    ttcov = cnt.varsCount;
                }

                if (ttcov < cnt.varsCount) {
                    ttcov = cnt.varsCount;
                    if (refCoverage.containsKey(position + 1) && ttcov < refCoverage.get(position + 1) - cnt.varsCount) {
                        ttcov = refCoverage.get(position + 1);
                        // Adjust the reference
                        Variation variantNextPosition = getVariationMaybe(getNonInsertionVariants(), position + 1, ref.get(position + 1));
                        if (variantNextPosition != null) {
                            variantNextPosition.varsCountOnForward -= fwd;
                            variantNextPosition.varsCountOnReverse -= rev;
                        }
                    }
                    totalPosCoverage = ttcov;
                }

                Variant tvref = new Variant();
                if (hicov < hicnt) {
                    hicov = hicnt;
                }
                tvref.descriptionString = descriptionString;
                tvref.positionCoverage = cnt.varsCount;
                tvref.varsCountOnForward = fwd;
                tvref.varsCountOnReverse = rev;
                tvref.strandBiasFlag = String.valueOf(bias);
                tvref.frequency = Utils.roundHalfEven("0.0000", cnt.varsCount / (double) ttcov);
                tvref.meanPosition = Utils.roundHalfEven("0.0", cnt.meanPosition / (double) cnt.varsCount);
                tvref.isAtLeastAt2Positions = cnt.pstd;
                tvref.meanQuality = vqual;
                tvref.hasAtLeast2DiffQualities = cnt.qstd;
                tvref.meanMappingQuality = mq;
                tvref.highQualityToLowQualityRatio = hicnt / (locnt != 0 ? locnt : 0.5d);
                tvref.highQualityReadsFrequency = hicov > 0 ? hicnt / (double)hicov : 0;
                tvref.extraFrequency = cnt.extracnt != 0 ? cnt.extracnt / (double)ttcov : 0;
                tvref.shift3 = 0;
                tvref.msi = 0;
                tvref.numberOfMismatches = Utils.roundHalfEven("0.0", cnt.numberOfMismatches / (double)cnt.varsCount);
                tvref.hicnt = hicnt;
                tvref.hicov = hicov;
                tvref.duprate = duprate;

                var.add(tvref);
                if (instance().conf.debug) {
                    tvref.debugVariantsContentInsertion(debugLines, descriptionString);
                }
            }
        }
        return totalPosCoverage;
    }

    /**
     * Creates non-insertion variant from each variation on this position and adds it to the list of aligned variants.
     * @param duprate duplication rate of non-insertion variant
     * @param alignedVars map of variants for position to get SV info
     * @param position position of the non-insertion variant start
     * @param nonInsertionVariations map of description strings and intermediate variations
     * @param totalPosCoverage current total position coverage (total depth) that can be updated by extra counts
     * @param var list of variants at the position to be updated
     * @param debugLines list of debug lines to be updated
     * @param keys sorted list of variant description strings
     * @param hicov position coverage by high quality reads
     */
    void createVariant(double duprate, Map<Integer, Vars> alignedVars, int position,
                       VariationMap<String, Variation> nonInsertionVariations, int totalPosCoverage, List<Variant> var,
                       List<String> debugLines, List<String> keys, int hicov) {
        //Loop over all variants found for the position except insertions
        for (String descriptionString : keys) {
            if (descriptionString.equals("SV")) {
                VariationMap.SV sv = nonInsertionVariations.sv;
                getOrPutVars(alignedVars, position).sv = sv.splits + "-" + sv.pairs + "-" + sv.clusters;
                continue;
            }
            Variation cnt = nonInsertionVariations.get(descriptionString);
            if (cnt.varsCount == 0) { //Skip variant if it does not have count
                continue;
            }
            //count of variants in forward strand
            int fwd = cnt.getDir(false);
            //count of variants in reverse strand
            int rev = cnt.getDir(true);
            //strand bias flag (0, 1 or 2)
            int bias = strandBias(fwd, rev);
            //$vqual mean base quality for variant
            double baseQuality = roundHalfEven("0.0", cnt.meanQuality / cnt.varsCount); // base quality
            //$mq mean mapping quality for variant
            double mappinqQuality = roundHalfEven("0.0", cnt.meanMappingQuality / (double) cnt.varsCount);
            //number of high-quality reads for variant
            int hicnt = cnt.highQualityReadsCount;
            //number of low-quality reads for variant
            int locnt = cnt.lowQualityReadsCount;
            /**
             * Condition:
             # 1). cnt.cnt > tcov                         - variant count is more than position coverage
             # 2). cnt.cnt - tcov < cnt.extracnt          - variant count is no more than position coverage + extracnt
             */
            int ttcov = totalPosCoverage;
            if (cnt.varsCount > totalPosCoverage && cnt.extracnt > 0 && cnt.varsCount - totalPosCoverage < cnt.extracnt) { //adjust position coverage if condition holds
                ttcov = cnt.varsCount;
            }

            //create variant record
            Variant tvref = new Variant();
            tvref.descriptionString = descriptionString;
            tvref.positionCoverage = cnt.varsCount;
            tvref.varsCountOnForward = fwd;
            tvref.varsCountOnReverse = rev;
            tvref.strandBiasFlag = String.valueOf(bias);
            tvref.frequency = Utils.roundHalfEven("0.0000", cnt.varsCount / (double) ttcov);
            tvref.meanPosition = Utils.roundHalfEven("0.0", cnt.meanPosition / (double) cnt.varsCount);
            tvref.isAtLeastAt2Positions = cnt.pstd;
            tvref.meanQuality = baseQuality;
            tvref.hasAtLeast2DiffQualities = cnt.qstd;
            tvref.meanMappingQuality = mappinqQuality;
            tvref.highQualityToLowQualityRatio = hicnt / (locnt != 0 ? locnt : 0.5d);
            tvref.highQualityReadsFrequency = hicov > 0 ? hicnt / (double) hicov : 0;
            tvref.extraFrequency = cnt.extracnt != 0 ? cnt.extracnt / (double) ttcov : 0;
            tvref.shift3 = 0;
            tvref.msi = 0;
            tvref.numberOfMismatches = Utils.roundHalfEven("0.0", cnt.numberOfMismatches / (double) cnt.varsCount);
            tvref.hicnt = hicnt;
            tvref.hicov = hicov;
            tvref.duprate = duprate;

            //append variant record
            var.add(tvref);
            if (instance().conf.debug) {
                tvref.debugVariantsContentSimple(debugLines, descriptionString);
            }
        }
    }

    /**
     * Adjust variant negative counts of fields FWD, REV, RFC, RRC to zeros and print the information message to console
     * @param p start position of variant
     * @param vref variant to adjust
     */
    private void adjustVariantCounts(int p, Variant vref) {
        String message = "column in variant on position: " + p + " " + vref.refallele + "->" +
                vref.varallele + " was negative, adjusted to zero.";

        if (vref.refForwardCoverage < 0 ) {
            vref.refForwardCoverage = 0;
            System.err.println("Reference forward count " + message);
        }
        if (vref.refReverseCoverage < 0) {
            vref.refReverseCoverage = 0;
            System.err.println("Reference reverse count " + message);
        }
        if (vref.varsCountOnForward < 0) {
            vref.varsCountOnForward = 0;
            System.err.println("Variant forward count " + message);
        }
        if (vref.varsCountOnReverse < 0 ) {
            vref.varsCountOnReverse = 0;
            System.err.println("Variant reverse count " + message);
        }
    }

    /**
     * Calculates coverage by high quality reads on position
     * @param insertionVariations Map of description string on insertion cariation for this position
     * @param nonInsertionVariations Map of description string on non-insertion variation for this position
     * @return coverage for high quality reads
     */
    private int calcHicov(VariationMap<String, Variation> insertionVariations,
                          VariationMap<String, Variation> nonInsertionVariations) {
        int hicov = 0;
        for (Map.Entry<String, Variation> descVariantEntry : nonInsertionVariations.entrySet()) {
            if (descVariantEntry.getKey().equals("SV") || descVariantEntry.getKey().startsWith("+")) {
                continue;
            }
            hicov += descVariantEntry.getValue().highQualityReadsCount;
        }
        if (insertionVariations != null) {
            for (Variation variation : insertionVariations.values()) {
                //hicov += variation.highQualityReadsCount;
            }
        }
        return hicov;
    }

    /**
     * Find microsatellite instability
     * Tandemly repeated short sequence motifs ranging from 1– 6(8 in our case) base pairs are called microsatellites.
     * Other frequently used terms for these DNA regions are simple sequences or short tandem repeats (STRs)
     * @param tseq1 variant description string
     * @param tseq2 right 70 bases in reference sequence
     * @param left left 50 bases in reference sequence
     * @return MSI
     */
    private MSI findMSI(String tseq1, String tseq2, String left) {

        //Number of nucleotides in microsattelite
        int nmsi = 1;
        //Number of bases to be shifted to 3 prime
        int shift3 = 0;
        String maxmsi = "";
        double msicnt = 0;
        while (nmsi <= tseq1.length() && nmsi <= 6) {
            //Microsattelite nucleotide sequence; trim nucleotide(s) from the end
            String msint = substr(tseq1, -nmsi, nmsi);
            Pattern pattern = Pattern.compile("((" + msint + ")+)$");
            Matcher mtch = pattern.matcher(tseq1);
            String msimatch = "";
            if (mtch.find()) {
                msimatch = mtch.group(1);
            }
            if (left != null && !left.isEmpty()) {
                mtch = pattern.matcher(left + tseq1);
                if (mtch.find()) {
                    msimatch = mtch.group(1);
                }
            }
            double curmsi = msimatch.length() / (double)nmsi;
            mtch = Pattern.compile("^((" + msint + ")+)").matcher(tseq2);
            if (mtch.find()) {
                curmsi += mtch.group(1).length() / (double)nmsi;
            }
            if (curmsi > msicnt) {
                maxmsi = msint;
                msicnt = curmsi;
            }
            nmsi++;
        }

        String tseq = tseq1 + tseq2;
        while (shift3 < tseq2.length() && tseq.charAt(shift3) == tseq2.charAt(shift3)) {
            shift3++;
        }
        return new MSI(msicnt, shift3, maxmsi);
    }

    /**
     * Fill the information about reference, genotype, refallele and varallele to the variant.
     * @param position position to get data from reference
     * @param totalPosCoverage total coverage on position
     * @param variationsAtPos information about variant on position to fill reference variant data
     * @param debugLines list of debug lines to fill DEBUG field in variant
     */
    private void collectReferenceVariants(int position, int totalPosCoverage, Vars variationsAtPos, List<String> debugLines) {
        int referenceForwardCoverage = 0; // $rfc
        int referenceReverseCoverage = 0; // $rrc
        //description string for reference or best non-reference variant
        String genotype1;

        if (variationsAtPos.referenceVariant != null && variationsAtPos.referenceVariant.frequency >= instance().conf.freq) {
            genotype1 = variationsAtPos.referenceVariant.descriptionString;
        } else if (variationsAtPos.variants.size() > 0) {
            genotype1 = variationsAtPos.variants.get(0).descriptionString;
        } else {
            genotype1 = variationsAtPos.referenceVariant.descriptionString;
        }
        if (variationsAtPos.referenceVariant != null) {
            referenceForwardCoverage = variationsAtPos.referenceVariant.varsCountOnForward;
            referenceReverseCoverage = variationsAtPos.referenceVariant.varsCountOnReverse;
        }

        if (genotype1.startsWith("+")) {
            Matcher mm = DUP_NUM.matcher(genotype1);
            if (mm.find()) {
                genotype1 = "+" + (Configuration.SVFLANK + toInt(mm.group(1)));
            }
            else {
                genotype1 = "+" + (genotype1.length() - 1);
            }
        }
        //description string for any other variant
        String genotype2;

        if (totalPosCoverage > refCoverage.get(position) && getNonInsertionVariants().containsKey(position + 1)
                && ref.containsKey(position + 1)
                && getNonInsertionVariants().get(position + 1).containsKey(ref.get(position + 1).toString())) {
            Variation tpref = getVariationMaybe(getNonInsertionVariants(), position + 1, ref.get(position + 1));
            referenceForwardCoverage = tpref.varsCountOnForward;
            referenceReverseCoverage = tpref.varsCountOnReverse;
        }
        List<Integer> positionsForChangedRefVariant = new ArrayList<>();
        // only reference reads are observed.
        if (variationsAtPos.variants.size() > 0) { //Condition: non-reference variants are found
            //Loop over non-reference variants
            for (Variant vref : variationsAtPos.variants) {
                //vref - variant reference
                String genotype1current = genotype1;
                genotype2 = vref.descriptionString;
                if (genotype2.startsWith("+")) {
                    genotype2 = "+" + (genotype2.length() - 1);
                }
                //variant description string
                final String descriptionString = vref.descriptionString; //$vn
                //length of deletion in variant (0 if no deletion)
                int deletionLength = 0; //$dellen
                Matcher matcher = BEGIN_MINUS_NUMBER.matcher(descriptionString);
                if (matcher.find()) {
                    deletionLength = toInt(matcher.group(1));
                }
                //effective position (??): p + dellen - 1 for deletion, p otherwise
                int endPosition = position;
                if (descriptionString.startsWith("-")) {
                    endPosition = position + deletionLength - 1;
                }
                //reference sequence for variation (to be written to .vcf file)
                String refallele = "";
                //variant sequence (to be written to .vcf file)
                String varallele;

                // how many bp can a deletion be shifted to 3 prime
                //3' shift (integer) for MSI adjustment
                int shift3 = 0;
                double msi = 0;
                String msint = "";

                int startPosition = position;

                MSI refVariantMsi;
                //if variant is an insertion
                if (descriptionString.startsWith("+")) {
                    //If no '&' and '#' symbols are found in variant string
                    //These symbols are in variant if a matched sequence follows insertion
                    if (!descriptionString.contains("&") && !descriptionString.contains("#") && !descriptionString.contains("<dup")) {
                        refVariantMsi = proceedVrefIsInsertion(position, descriptionString);
                        msi = refVariantMsi.msi;
                        shift3 = refVariantMsi.shift3;
                        msint = refVariantMsi.msint;
                    }
                    //Shift position to 3' if -3 option is set
                    if (instance().conf.moveIndelsTo3) {
                        startPosition += shift3;
                        endPosition += shift3;
                    }
                    //reference allele is 1 base
                    refallele = ref.containsKey(position) ? ref.get(position).toString() : "";
                    //variant allele is reference base concatenated with insertion
                    varallele = refallele + descriptionString.substring(1);
                    if (varallele.length() > instance().conf.SVMINLEN) {
                        endPosition += varallele.length();
                        varallele = "<DUP>";
                    }
                    Matcher mm = DUP_NUM.matcher(varallele);
                    if (mm.find()) {
                        int dupCount = toInt(mm.group(1));
                        endPosition = startPosition + (2 * Configuration.SVFLANK + dupCount) - 1;
                        genotype2 = "+" + (2 * Configuration.SVFLANK + dupCount);
                        varallele = "<DUP>";
                    }
                } else if (descriptionString.startsWith("-")) { //deletion variant
                    Matcher matcherINV = INV_NUM.matcher(descriptionString);
                    Matcher matcherStartMinusNum = BEGIN_MINUS_NUMBER_CARET.matcher(descriptionString);

                    if (deletionLength < instance().conf.SVMINLEN) {
                        //variant allele is in the record
                        //remove '-' and number from beginning of variant string
                        varallele = descriptionString.replaceFirst("^-\\d+", "");

                        refVariantMsi = proceedVrefIsDeletion(position, deletionLength);
                        msi = refVariantMsi.msi;
                        shift3 = refVariantMsi.shift3;
                        msint = refVariantMsi.msint;

                        if (matcherINV.find()) {
                            varallele = "<INV>";
                            genotype2 = "<INV" + deletionLength + ">";
                        }
                    } else if (matcherStartMinusNum.find()) {
                        varallele = "<INV>";
                        genotype2 = "<INV" + deletionLength + ">";
                    } else {
                        varallele = "<DEL>";
                    }
                    //If no matched sequence or indel follows the variant
                    if (!descriptionString.contains("&") && !descriptionString.contains("#") && !descriptionString.contains("^")) {
                        //Shift position to 3' if -3 option is set
                        if (instance().conf.moveIndelsTo3) {
                            startPosition += shift3;
                        }
                        //variant allele is 1 base from reference string preceding p
                        if (!varallele.equals("<DEL>")) {
                            varallele = ref.containsKey(position - 1) ? ref.get(position - 1).toString() : "";
                        }
                        //prepend same base to reference allele
                        refallele = ref.containsKey(position - 1) ? ref.get(position - 1).toString() : "";
                        startPosition--;
                    }
                    Matcher mm = SOME_SV_NUMBERS.matcher(descriptionString);
                    if (mm.find()) {
                        refallele = ref.containsKey(position) ? ref.get(position).toString() : "";
                    }
                    else if (deletionLength < instance().conf.SVMINLEN) {
                        //append dellen bases from reference string to reference allele
                        refallele += joinRef(ref, position, position + deletionLength - 1);
                    }
                } else { //Not insertion/deletion variant. SNP or MNP
                    //Find MSI adjustment
                    String tseq1 = joinRef(ref, Math.max(position - REF_30_BASES, 1), position + 1);
                    int chr0 = getOrElse(instance().chrLengths, region.chr, 0);
                    String tseq2 = joinRef(ref, position + 2, Math.min(position + REF_70_BASES, chr0));

                    MSI msiData = findMSI(tseq1, tseq2, null);
                    msi = msiData.msi;
                    shift3 = msiData.shift3;
                    msint = msiData.msint;
                    //reference allele is 1 base from reference sequence
                    refallele = ref.containsKey(position) ? ref.get(position).toString() : "";
                    //variant allele is same as description string
                    varallele = descriptionString;
                }

                Matcher mtch = AMP_ATGC.matcher(descriptionString);
                if (mtch.find()) { //If variant is followed by matched sequence
                    //following matching sequence
                    String extra = mtch.group(1);
                    //remove '&' symbol from variant allele
                    varallele = varallele.replaceFirst("&", "");
                    //append length(extra) bases from reference sequence to reference allele and genotype1
                    String tch = joinRef(ref, endPosition + 1, endPosition + extra.length());
                    refallele += tch;
                    genotype1current += tch;

                    //Adjust position
                    endPosition += extra.length();

                    mtch = AMP_ATGC.matcher(varallele);
                    if (mtch.find()) {
                        String vextra = mtch.group(1);
                        varallele = varallele.replaceFirst("&", "");
                        tch = joinRef(ref, endPosition + 1, endPosition + vextra.length());
                        refallele += tch;
                        genotype1current += tch;
                        endPosition += vextra.length();
                    }

                    //If description string starts with '+' sign, remove it from reference and variant alleles
                    if (descriptionString.startsWith("+")) {
                        refallele = refallele.substring(1);
                        varallele = varallele.substring(1);
                        startPosition++;
                    }

                    if (varallele.equals("<DEL>") && refallele.length() >= 1) {
                        refallele = ref.containsKey(startPosition) ? ref.get(startPosition).toString() : "";
                        if (refCoverage.containsKey(startPosition - 1)) {
                            totalPosCoverage = refCoverage.get(startPosition - 1);
                        }
                        if (vref.positionCoverage > totalPosCoverage ){
                            totalPosCoverage = vref.positionCoverage;
                        }
                        vref.frequency = vref.positionCoverage / (double) totalPosCoverage;
                    }
                }

                //If variant is followed by short matched sequence and insertion/deletion
                mtch = HASH_GROUP_CARET_GROUP.matcher(descriptionString);
                if (mtch.find()) {
                    String matchedSequence = mtch.group(1); //$mseq
                    //insertion/deletion tail
                    String tail = mtch.group(2);

                    //adjust position by length of matched sequence
                    endPosition += matchedSequence.length();

                    //append bases from reference sequence to reference allele
                    refallele += joinRef(ref, endPosition - matchedSequence.length() + 1, endPosition);

                    //If tail is a deletion
                    mtch = BEGIN_DIGITS.matcher(tail);
                    if (mtch.find()) {
                        //append (deletion length) bases from reference sequence to reference allele
                        int deletion = toInt(mtch.group(1)); //$d
                        refallele += joinRef(ref, endPosition + 1, endPosition + deletion);

                        //shift position by deletion length
                        endPosition += deletion;
                    }
                    //clean special symbols from alleles
                    varallele = varallele.replaceFirst("#", "").replaceFirst("\\^(\\d+)?", "");

                    //replace '#' with 'm' and '^' with 'i' in genotypes
                    genotype1current = genotype1current.replaceFirst("#", "m").replaceFirst("\\^", "i");
                    genotype2 = genotype2.replaceFirst("#", "m").replaceFirst("\\^", "i");
                }
                mtch = CARET_ATGNC.matcher(descriptionString); // for deletion followed directly by insertion in novolign
                if (mtch.find()) {
                    //remove '^' sign from varallele
                    varallele = varallele.replaceFirst("\\^", "");

                    //replace '^' sign with 'i' in genotypes
                    genotype1current = genotype1current.replaceFirst("\\^", "i");
                    genotype2 = genotype2.replaceFirst("\\^", "i");
                }

                // Perform adjustment to get as close to CRISPR site as possible
                int cutSite = instance().conf.crisprCuttingSite;
                if (cutSite != 0 && refallele.length() > 1 && varallele.length() > 1 ) { // fix 5' for complex in CRISPR mode
                    int n = 0;
                    while (refallele.length() > n + 1 && varallele.length() > n + 1 && substr(refallele, n, 1).equals(substr(varallele, n, 1))) {
                        n++;
                    }
                    if (n != 0) {
                        startPosition += n;
                        refallele = substr(refallele, n);
                        varallele = substr(varallele, n);
                    }
		        // Let adjComplex to take care the 3'
                }

                if (cutSite != 0 && (refallele.length() != varallele.length() && substr(refallele, 0, 1).equals(substr(varallele, 0, 1)))) {
                    if (!(startPosition == cutSite || endPosition == cutSite)) {
                        int n = 0;
                        int dis = Math.min(abs(cutSite - startPosition), abs(cutSite - endPosition));
                        if (startPosition < cutSite) {
                            while(startPosition + n < cutSite && n < shift3 && endPosition + n != cutSite) {
                                n++;
                            }
                            if (abs(startPosition + n - cutSite) > dis && abs(endPosition + n - cutSite) > dis ) {
                                n = 0;
                                // Don't move if it makes it worse
                            }
                        }
                        if (endPosition < cutSite && n == 0) {
                            if (abs(endPosition - cutSite) <= abs(startPosition - cutSite) ) {
                                while(endPosition + n < cutSite && n < shift3) {
                                    n++;
                                }
                            }
                        }
                        if (n > 0) {
                            startPosition += n;
                            endPosition += n;
                            refallele = "";
                            for (int i = startPosition; i <= endPosition; i++) {
                                refallele += ref.get(i);
                            }
                            String tva = "";
                            if (refallele.length() < varallele.length()) { // Insertion
                                tva = substr(varallele, 1);
                                if (tva.length() > 1) {
                                    int ttn = n % tva.length();
                                    if (ttn != 0) {
                                        tva = substr(tva, ttn) + substr(tva, 0, ttn);
                                    }
                                }
                            }
                            varallele = ref.get(startPosition) + tva;
                        }
                        vref.crispr = n;
                    }
                }

                //preceding reference sequence
                vref.leftseq = joinRef(ref, Math.max(startPosition - REF_20_BASES, 1), startPosition - 1); // left 20 nt
                int chr0 = getOrElse(instance().chrLengths, region.chr, 0);
                //following reference sequence
                vref.rightseq = joinRef(ref, endPosition + 1, Math.min(endPosition + REF_20_BASES, chr0)); // right 20 nt
                //genotype description string
                String genotype = genotype1current + "/" + genotype2;
                //remove '&' and '#' symbols from genotype string
                //replace '^' symbol with 'i' in genotype string
                genotype = genotype
                        .replace("&", "")
                        .replace("#", "")
                        .replace("^", "i");
                //convert extrafreq, freq, hifreq, msi fields to strings
                vref.extraFrequency = roundHalfEven("0.0000", vref.extraFrequency);
                vref.frequency = roundHalfEven("0.0000", vref.frequency);
                vref.highQualityReadsFrequency = roundHalfEven("0.0000", vref.highQualityReadsFrequency);
                vref.msi = roundHalfEven("0.000", msi);
                vref.msint = msint.length();
                vref.shift3 = shift3;
                vref.startPosition = startPosition;
                vref.endPosition = endPosition;
                vref.refallele = validateRefallele(refallele);
                vref.varallele = varallele;
                vref.genotype = genotype;
                vref.totalPosCoverage = totalPosCoverage;
                vref.refForwardCoverage = referenceForwardCoverage;
                vref.refReverseCoverage = referenceReverseCoverage;

                //bias is [0-2];[0-2] where first flag is for reference, second for variant
                //if reference variant is not found, first flag is 0
                if (variationsAtPos.referenceVariant != null) {
                    vref.strandBiasFlag = variationsAtPos.referenceVariant.strandBiasFlag + ";" + vref.strandBiasFlag;
                } else {
                    vref.strandBiasFlag = "0;" + vref.strandBiasFlag;
                }

                adjustVariantCounts(position, vref);
                if (startPosition != position && instance().conf.doPileup) {
                    positionsForChangedRefVariant.add(position);
                }
                constructDebugLines(debugLines, vref);
            }
            //TODO: It is a "lazy" solution because current logic in realignment methods can't be changed simply for --nosv option
            if (instance().conf.disableSV) {
                variationsAtPos.variants.removeIf(vref -> ANY_SV.matcher(vref.varallele).find());
            }
        } else if (variationsAtPos.referenceVariant != null ) {  //no variant reads are detected
            Variant vref = variationsAtPos.referenceVariant;
            updateRefVariant(position, totalPosCoverage, vref, debugLines, referenceForwardCoverage, referenceReverseCoverage);
        } else {
            variationsAtPos.referenceVariant = new Variant();
        }

        // Update reference variants if there were indels and start position were changed, or we work with amplicons
        // so after update ref variants can be output in pileup
        if (variationsAtPos.referenceVariant != null && instance().conf.doPileup &&
                (positionsForChangedRefVariant.contains(position)|| instance().ampliconBasedCalling != null)) {
            Variant vref = variationsAtPos.referenceVariant;
            updateRefVariant(position, totalPosCoverage, vref, debugLines, referenceForwardCoverage, referenceReverseCoverage);
        }
    }

    private void updateRefVariant(int position, int totalPosCoverage, Variant vref, List<String> debugLines,
                                  int referenceForwardCoverage, int referenceReverseCoverage) {
        vref.totalPosCoverage = totalPosCoverage;
        vref.positionCoverage = 0;
        vref.frequency = 0;
        vref.refForwardCoverage = referenceForwardCoverage;
        vref.refReverseCoverage = referenceReverseCoverage;
        vref.varsCountOnForward = 0;
        vref.varsCountOnReverse = 0;
        vref.msi = 0;
        vref.msint = 0;
        if (vref.strandBiasFlag.indexOf(';') == -1) {
            vref.strandBiasFlag += ";0";
        }
        vref.shift3 = 0;
        vref.startPosition = position;
        vref.endPosition = position;
        vref.highQualityReadsFrequency = roundHalfEven("0.0000", vref.highQualityReadsFrequency);
        String referenceBase = ref.containsKey(position) ? ref.get(position).toString() : ""; // $r
        //both refallele and varallele are 1 base from reference string
        vref.refallele = validateRefallele(referenceBase);
        vref.varallele = validateRefallele(referenceBase);
        vref.genotype = referenceBase + "/" + referenceBase;
        vref.leftseq = "";
        vref.rightseq = "";
        vref.duprate = duprate;

        constructDebugLines(debugLines, vref);
    }

    /**
     * Construct DEBUG lines for the variant
     */
    private void constructDebugLines(List<String> debugLines, Variant vref) {

        if (instance().conf.debug) {
            StringBuilder sb = new StringBuilder();
            for (String str : debugLines) {
                if (sb.length() > 0) {
                    sb.append(" & ");
                }
                sb.append(str);
            }
            vref.DEBUG = sb.toString();
        }
    }

    /**
     * Validate reference allele according to VCF 4.3 specification in case if IUPAC ambiguity codes are present
     * in reference.
     * @param refallele sequence of reference bases that covers variant
     * @return reference allele sequence where IUPAC ambuguity bases are changed to the one that is
     * first alphabetically.
     */
    String validateRefallele(String refallele) {
        for (int i = 0; i < refallele.length(); i++) {
            String refBase = substr(refallele, i, 1);
            if (IUPAC_AMBIGUITY_CODES.containsKey(refBase)) {
                refallele = refallele.replaceFirst(refBase, IUPAC_AMBIGUITY_CODES.get(refBase));
            }
        }
        return refallele;
    }

    /**
     * Microsatellite instability
     * Tandemly repeated short sequence motifs ranging from 1– 6(8 in our case) base pairs are called microsatellites.
     */
    class MSI {
        /**
         * MSI count
         */
        double msi;

        /**
         * No. of bases to be shifted to 3 prime for deletions due to alternative alignment
         */
        int shift3;
        /**
         * MicroSattelite unit length in base pairs
         */
        String msint;

        public MSI(double msi, int shift3, String msint) {
            this.msi = msi;
            this.shift3 = shift3;
            this.msint = msint;
        }
    }
}
