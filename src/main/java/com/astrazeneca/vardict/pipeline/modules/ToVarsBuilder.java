package com.astrazeneca.vardict.pipeline.modules;

import com.astrazeneca.utils.Tuple;
import com.astrazeneca.vardict.Configuration;
import com.astrazeneca.vardict.Region;
import com.astrazeneca.vardict.pipeline.data.RealignedVariationData;
import com.astrazeneca.vardict.pipeline.data.Scope;
import com.astrazeneca.vardict.variations.Variant;
import com.astrazeneca.vardict.variations.Variation;
import com.astrazeneca.vardict.variations.VariationUtils;
import com.astrazeneca.vardict.variations.Vars;

import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static com.astrazeneca.GlobalReadOnlyScope.instance;
import static com.astrazeneca.utils.Tuple.tuple;
import static com.astrazeneca.utils.Utils.*;
import static java.lang.String.format;

/**
 * Class for creation of variant structure
 */
public class ToVarsBuilder implements Module<RealignedVariationData, Tuple.Tuple2<Integer, Map<Integer, Vars>>> {

    private Region region;
    private Map<Integer, Character> reference;
    private List<String> debugLines;

    /**
     * 14	     * @param region region of interest
     * 15	     * @param reference part of reference sequence (key - position, value - base)
     * 16	     * @param chromosomesLengths map of chromosome lengths
     * 17	     * @param config VarDict configuration
     * 18
     */
    public ToVarsBuilder() {
        this.debugLines = new ArrayList<>();
    }

    public void setReference(Map<Integer, Character> reference) {
        this.reference = reference;
    }

    @Override
    public Scope<Tuple.Tuple2<Integer, Map<Integer, Vars>>> process(Scope<RealignedVariationData> scope) {

        initFromScope(scope);

        //the variant structure
        Map<Integer, Vars> alignedVariants = new HashMap<>();
        //perl version: 1344
        //Loop over positions
        for (Map.Entry<Integer, Map<String, Variation>> entH : scope.data.nonInsertionVariants.entrySet()) {
            final int position = entH.getKey();
            Map<String, Variation> varsAtCurPosition = entH.getValue();

            if (varsAtCurPosition.isEmpty()) {
                continue;
            }
            //perl version: 1346
            //Skip if position is outside region of interest
            if (position < region.start || position > region.end) {
                continue;
            }

            //perl version: 1361
            //skip position if it has no coverage
            if (!scope.data.refCoverage.containsKey(position) || scope.data.refCoverage.get(position) == 0) {
                continue;
            }
            //perl version: 1351
            if (isTheSameVariationOnRef(scope.data.insertionVariants, position, varsAtCurPosition))
                continue;
            //perl version: 1357
            //total position coverage
            final int totalPosCoverage = scope.data.refCoverage.get(position);
            if (totalPosCoverage == 0) { // ignore when there's no coverage
                continue;
            }

            //position coverage by high-quality reads
            final int hicov = calcHicov(varsAtCurPosition, scope.data.insertionVariants.get(position));
            //perl version: 1370
            //array of all variants for the position
            List<Variant> var = new ArrayList<>();
            //temporary array used for debugging
            List<String> tmp = new ArrayList<>();

            Configuration config = instance().conf;

            //perl version: 1387
            createVariants(varsAtCurPosition, totalPosCoverage, hicov, var, tmp);
            createInsertionVariants(scope.data.insertionVariants, position, totalPosCoverage, hicov, var, tmp);
            //perl version: 1404
            sortVariants(var);

            //perl version: 1405
            double maxfreq = collectVarsAtPosition(alignedVariants, position, var);
            //perl version: 1415
            if (!config.doPileup && maxfreq <= config.freq && config.ampliconBasedCalling == null) {
                if (!config.bam.hasBam2()) {
                    alignedVariants.remove(position);
                    continue;
                }
            }
            Vars variationsAtPos = getOrPutVars(alignedVariants, position);
            //perl version: 1423
            collectReferenceVariants(position, totalPosCoverage, variationsAtPos);
        }

        return new Scope<>(
                scope,
                tuple(scope.data.maxReadLength, alignedVariants));
    }

    private void initFromScope(Scope<RealignedVariationData> scope) {
        this.reference = scope.regionRef;
        this.region = scope.region;
    }

    /**
     * Collect reference variants method
     *
     * @param variationsAtPos  variants at position
     * @param position         position of variants
     * @param totalPosCoverage total position coverage
     */
    private void collectReferenceVariants(int position, int totalPosCoverage, Vars variationsAtPos) {
        // Make sure the first strandBiasFlag is always for the reference nucleotide
        //perl version: 1423
        //reference forward strand coverage
        int refForwardCoverage = 0;
        //reference reverse strand coverage
        int refReverseCoverage = 0;
        //description string for reference or best non-reference variant
        String genotype1 = "";

        //if reference variant has frequency more than $FREQ, set genotype1 to reference variant
        //else if non-reference variant exists, set genotype1 to first (largest quality*coverage) variant
        //otherwise set genotype1 to reference variant
        if (variationsAtPos.referenceVariant != null) {
            if (variationsAtPos.referenceVariant.frequency >= instance().conf.freq) {
                genotype1 = variationsAtPos.referenceVariant.descriptionString;
            } else if (variationsAtPos.variants.size() > 0) {
                genotype1 = variationsAtPos.variants.get(0).descriptionString;
            }
            refForwardCoverage = variationsAtPos.referenceVariant.varsCountOnForward;
            refReverseCoverage = variationsAtPos.referenceVariant.varsCountOnReverse;
        } else if (variationsAtPos.variants.size() > 0) {
            genotype1 = variationsAtPos.variants.get(0).descriptionString;
        }
        //perl version: 1425
        if (genotype1.startsWith("+")) {
            genotype1 = "+" + (genotype1.length() - 1);
        }

        //perl version: 1432
        // only reference reads are observed.
        if (variationsAtPos.variants.size() > 0) { //Condition: non-reference variants are found
            //Loop over non-reference variants
            for (Variant nonRefVariant : variationsAtPos.variants) {
                genotype1 = proceedNonRefVariant(position, totalPosCoverage, variationsAtPos, refForwardCoverage, refReverseCoverage, genotype1, nonRefVariant);
            }
        } else {
            if (variationsAtPos.referenceVariant == null)
                variationsAtPos.referenceVariant = new Variant();
        }
        //perl version: 1605
        if (variationsAtPos.referenceVariant != null) {
            Variant refVariation = variationsAtPos.referenceVariant; //no variant reads are detected.
            refVariation.totalPosCoverage = totalPosCoverage;
            refVariation.positionCoverage = 0;
            refVariation.frequency = 0;
            refVariation.refForwardCoverage = refForwardCoverage;
            refVariation.refReverseCoverage = refReverseCoverage;
            refVariation.varsCountOnForward = 0;
            refVariation.varsCountOnReverse = 0;
            refVariation.msi = 0;
            refVariation.msint = 0;
            refVariation.strandBiasFlag += ";0";
            refVariation.shift3 = 0;
            refVariation.startPosition = position;
            refVariation.endPosition = position;
            refVariation.highQualityReadsFrequency = refVariation.highQualityReadsFrequency;
            String referenceValueAtPos = reference.containsKey(position) ? reference.get(position).toString() : "";
            //both refAllele and varAllele are 1 base from reference string
            refVariation.refAllele = referenceValueAtPos;
            refVariation.varAllele = referenceValueAtPos;
            refVariation.genotype = referenceValueAtPos + "/" + referenceValueAtPos;
            refVariation.precedingRefSequence = "";
            refVariation.followingRefSequence = "";
            setDebugField(refVariation);
        }
    }

    /**
     * Processing non-reference variants method:
     * fill field of non reference variation such as bias, genotype, totalCoverage and other
     *
     * @param position           position of variants
     * @param totalPosCoverage   total position coverage
     * @param refForwardCoverage reference forward strand coverage
     * @param refReverseCoverage reference reverse strand coverage
     * @param genotype1          description string for reference or best non-reference variant
     * @param variationsAtPos    variants at current position
     * @param nonRefVariant      non-reference variants
     */
    private String proceedNonRefVariant(int position, int totalPosCoverage, Vars variationsAtPos, int refForwardCoverage, int refReverseCoverage, String genotype1, Variant nonRefVariant) {
        //perl version: 1436
        String genotype2 = nonRefVariant.descriptionString; // description string for any other variant
        //perl version: 1437
        if (genotype2.startsWith("+")) {
            genotype2 = "+" + (genotype2.length() - 1);
        }
        //perl version: 1438
        //variant description string
        final String varDescriptionStr = nonRefVariant.descriptionString;
        //length of deletion in variant (0 if no deletion)
        int deletionLength = 0;
        //perl version: 1439
        Matcher matcher = PatternsUtils.START_DIG.matcher(varDescriptionStr);
        if (matcher.find()) {
            deletionLength = toInt(matcher.group(1));
        }
        //perl version: 1440
        //effective position (??): position + deletionLength - 1 for deletion, position otherwise
        int effectivePosition = position;
        if (varDescriptionStr.startsWith("-")) {
            effectivePosition = position + deletionLength - 1;
        }
        //perl version: 1441
        //reference sequence for variation (to be written to .vcf file)
        String refVariationSequence = "";
        //variant sequence (to be written to .vcf file)
        String variantSequence = "";
        //perl version: 1442
        // how many bp can a deletion be shifted to 3 prime
        //3' shift (integer) for MSI adjustment
        int shift3 = 0;
        double msi = 0;
        String msint = "";
        //perl version: 1443
        int startPosition = position;
        Tuple.Tuple3<Double, Integer, String> tupleRefVariant = Tuple.tuple(msi, shift3, msint);
        //perl version: 1444
        //if variant is an insertion
        if (varDescriptionStr.startsWith("+")) {
            //perl version: 1445
            //If no '&' and '#' symbols are found in variant string
            //These symbols are in variant if a matched sequence follows insertion
            if (!varDescriptionStr.contains("&") && !varDescriptionStr.contains("#")) {
                tupleRefVariant = proceedVrefIsInsertion(position, varDescriptionStr);
            }
            //perl version: 1463
            //Shift position to 3' if -3 option is set
            if (instance().conf.moveIndelsTo3) {
                startPosition += tupleRefVariant._2;
                effectivePosition += tupleRefVariant._2;
            }
            //perl version: 1467
            //reference allele is 1 base
            refVariationSequence = reference.containsKey(position) ? reference.get(position).toString() : "";
            //variant allele is reference base concatenated with insertion
            variantSequence = refVariationSequence + varDescriptionStr.substring(1);
        }
        //perl version: 1470
        else if (varDescriptionStr.startsWith("-")) { //deletion variant
            //variant allele is in the record
            //remove '-' and number from beginning of variant string
            variantSequence = varDescriptionStr.replaceFirst("^-\\d+", "");
            //perl version: 1473
            tupleRefVariant = proceedVariantDeletion(position, deletionLength);
            //perl version: 1479
            //If no matched sequence or indel follows the variant
            if (!varDescriptionStr.contains("&") && !varDescriptionStr.contains("#") && !varDescriptionStr.contains("^")) {
                //perl version: 1486
                //Shift position to 3' if -3 option is set
                if (instance().conf.moveIndelsTo3) {
                    startPosition += tupleRefVariant._2;
                }
                //perl version: 1489
                //variant allel is 1 base from reference string preceding position
                variantSequence = reference.containsKey(position - 1) ? reference.get(position - 1).toString() : "";
                //prepend same base to reference allele
                refVariationSequence = variantSequence;
                startPosition--;
            }
            //perl version: 1493
            //append deletionLength bases from reference string to reference allele
            refVariationSequence += VariationUtils.joinRef(reference, position, position + deletionLength - 1);
        }
        //perl version: 1497
        else { //Not insertion/deletion variant. SNP or MNP
            //perl version: 1498
            tupleRefVariant = findMSIAdjustment(position);
            //reference allele is 1 base from reference sequence
            refVariationSequence = reference.containsKey(position) ? reference.get(position).toString() : "";
            //variant allele is same as description string
            variantSequence = varDescriptionStr;
        }

        genotype1 = matchingNonRefVariant(totalPosCoverage, refForwardCoverage, refReverseCoverage, genotype1, nonRefVariant, genotype2, varDescriptionStr, effectivePosition, refVariationSequence, variantSequence, startPosition, tupleRefVariant);
        //perl version: 1601
        //strandBiasFlag is [0-2];[0-2] where first flag is for reference, second for variant
        //if reference variant is not found, first flag is 0
        if (variationsAtPos.referenceVariant != null) {
            nonRefVariant.strandBiasFlag = variationsAtPos.referenceVariant.strandBiasFlag + ";" + nonRefVariant.strandBiasFlag;
        } else {
            nonRefVariant.strandBiasFlag = "0;" + nonRefVariant.strandBiasFlag;
        }
        setDebugField(nonRefVariant);
        return genotype1;
    }

    /**
     * Method for finding of variant is followed by matched sequence
     *
     * @param totalPosCoverage     - total position coverage
     * @param refForwardCoverage   - reference variant forward strand coverage
     * @param refReverseCoverage   - reference variant reverse strand coverage
     * @param genotype1            - description string for reference or best non-reference variant
     * @param nonRefVariant        - non-reference variants
     * @param genotype2            - description string for any other variant
     * @param varDescriptionStr    - variant description string
     * @param effectivePosition    - effective position (??): position + deletionLength - 1 for deletion, position otherwise
     * @param refVariationSequence - reference sequence for variation
     * @param variantSequence      - variant sequence (to be written to .vcf file)
     * @param startPosition        - sum of shift3
     * @param tupleRefVariant      - tuple(msi, shift3, msint),
     *                             where msi: Microsatellite instability,
     *                             shift3: number of bases to be shifted to 3 prime,
     *                             msint: microsattelite nucleotide sequence; trim nucleotide(s) from the end
     * @return genotype1 - description string for reference or best non-reference variant
     */
    private String matchingNonRefVariant(int totalPosCoverage, int refForwardCoverage, int refReverseCoverage, String genotype1, Variant nonRefVariant, String genotype2, String varDescriptionStr, int effectivePosition, String refVariationSequence, String variantSequence, int startPosition, Tuple.Tuple3<Double, Integer, String> tupleRefVariant) {
        //perl version: 1503
        Matcher mtch = PatternsUtils.AMP_ATGC.matcher(varDescriptionStr);
        if (mtch.find()) { //If variant is followed by matched sequence
            //following matching sequence
            String extra = mtch.group(1);
            //remove '&' symbol from variant allele
            variantSequence = variantSequence.replaceFirst("&", "");
            //append length(extra) bases from reference sequence to reference allele and genotype1
            String tch = VariationUtils.joinRef(reference, effectivePosition + 1, effectivePosition + extra.length());
            refVariationSequence += tch;
            genotype1 += tch;
            //perl version: 1510
            //Adjust position
            effectivePosition += extra.length();
            //perl version: 1511
            mtch = PatternsUtils.AMP_ATGC.matcher(variantSequence);
            if (mtch.find()) {
                String vextra = mtch.group(1);
                variantSequence = variantSequence.replaceFirst("&", "");
                tch = VariationUtils.joinRef(reference, effectivePosition + 1, effectivePosition + vextra.length());
                refVariationSequence += tch;
                genotype1 += tch;
                effectivePosition += vextra.length();
            }
            //perl version: 1520
            //If description string starts with '+' sign, remove it from reference and variant alleles
            if (varDescriptionStr.startsWith("+")) {
                refVariationSequence = refVariationSequence.substring(1);
                variantSequence = variantSequence.substring(1);
                startPosition++;
            }
        }
        //perl version: 1526
        //If variant is followed by short matched sequence and insertion/deletion
        mtch = PatternsUtils.HASH_GROUP_CARET_GROUP.matcher(varDescriptionStr);
        if (mtch.find()) {
            //matched sequence
            String matchedSequence = mtch.group(1);
            //insertion/deletion tail
            String tail = mtch.group(2);

            //adjust position by length of matched sequence
            effectivePosition += matchedSequence.length();

            //append bases from reference sequence to reference allele
            refVariationSequence += VariationUtils.joinRef(reference, effectivePosition - matchedSequence.length() + 1, effectivePosition);

            //If tail is a deletion
            mtch = PatternsUtils.BEGIN_DIGITS.matcher(tail);
            if (mtch.find()) {
                //append (deletion length) bases from reference sequence to reference allele
                int d = toInt(mtch.group(1));
                refVariationSequence += VariationUtils.joinRef(reference, effectivePosition + 1, effectivePosition + d);

                //shift position by deletion length
                effectivePosition += d;
            }

            //clean special symbols from alleles
            variantSequence = variantSequence.replaceFirst("#", "").replaceFirst("\\^(\\d+)?", "");
            //perl version: 1539
            //replace '#' with 'm' and '^' with 'i' in genotypes
            genotype1 = genotype1.replaceFirst("#", "m").replaceFirst("\\^", "i");
            genotype2 = genotype2.replaceFirst("#", "m").replaceFirst("\\^", "i");
        }
        //perl version: 1544
        mtch = PatternsUtils.CARET_ATGNC.matcher(varDescriptionStr); // for deletion followed directly by insertion in novolign
        if (mtch.find()) {
            //remove '^' sign from varAllele
            variantSequence = variantSequence.replaceFirst("\\^", "");
            //perl version: 1510
            //replace '^' sign with 'i' in genotypes
            genotype1 = genotype1.replaceFirst("\\^", "i");
            genotype2 = genotype2.replaceFirst("\\^", "i");
        }

        //preceding reference sequence
        nonRefVariant.precedingRefSequence = VariationUtils.joinRef(reference, startPosition - 20 < 1 ? 1 : startPosition - 20, startPosition - 1); // left 20 nt
        int chr0 = getOrElse(instance().chrLengths, region.chr, 0);
        //following reference sequence
        nonRefVariant.followingRefSequence = VariationUtils.joinRef(reference, effectivePosition + 1, effectivePosition + 20 > chr0 ? chr0 : effectivePosition + 20); // right 20 nt
        //perl version: 1549
        //genotype description string
        String genotype = genotype1 + "/" + genotype2;
        //remove '&' and '#' symbols from genotype string
        //replace '^' symbol with 'i' in genotype string
        genotype = genotype.replace("&", "").replace("#", "").replace("^", "i");
        //convert extraFrequency, frequency, highQualityReadsFrequency, msi fields to strings
        nonRefVariant.extraFrequency = nonRefVariant.extraFrequency;
        //perl version: 1588
        nonRefVariant.frequency = round(nonRefVariant.frequency, 4);
        nonRefVariant.highQualityReadsFrequency = round(nonRefVariant.highQualityReadsFrequency, 4);
        nonRefVariant.msi = tupleRefVariant._1;
        nonRefVariant.msint = tupleRefVariant._3.length();
        nonRefVariant.shift3 = tupleRefVariant._2;
        nonRefVariant.startPosition = startPosition;
        final Matcher pnMtch = PatternsUtils.PLUS_NUMBER.matcher(genotype);
        if (pnMtch.find()) {
            nonRefVariant.endPosition = startPosition + toInt(pnMtch.group(1));
        } else {
            nonRefVariant.endPosition = effectivePosition;
        }
        nonRefVariant.refAllele = refVariationSequence;
        nonRefVariant.varAllele = variantSequence;
        nonRefVariant.genotype = genotype;
        nonRefVariant.totalPosCoverage = totalPosCoverage;
        nonRefVariant.refForwardCoverage = refForwardCoverage;
        nonRefVariant.refReverseCoverage = refReverseCoverage;
        return genotype1;
    }

    /**
     * Find MSI adjustment (Not insertion/deletion variant)
     *
     * @param position position of variants
     * @return tuple(msi, shift3, msint),
     * where msi: Microsatellite instability,
     * shift3: number of bases to be shifted to 3 prime,
     * msint: microsattelite nucleotide sequence; trim nucleotide(s) from the end
     */
    private Tuple.Tuple3<Double, Integer, String> findMSIAdjustment(int position) {
        //Find MSI adjustment
        String tseq1 = VariationUtils.joinRef(reference, position - 30 > 1 ? position - 30 : 1, position + 1);
        int chr0 = getOrElse(instance().chrLengths, region.chr, 0);
        String tseq2 = VariationUtils.joinRef(reference, position + 2, position + 70 > chr0 ? chr0 : position + 70);

        Tuple.Tuple3<Double, Integer, String> tpl = findMSI(tseq1, tseq2, null);
        return tpl;
    }

    /**
     * Processing deletion variant method
     *
     * @param position       position of variants
     * @param deletionLength length of deletion
     * @return Tuple.tuple(msi, shift3, msint),
     * where msi: Microsatellite instability,
     * shift3: number of bases to be shifted to 3 prime,
     * msint: microsattelite nucleotide sequence; trim nucleotide(s) from the end
     */
    private Tuple.Tuple3<Double, Integer, String> proceedVariantDeletion(int position, int deletionLength) {
        //perl version: 1473
        //left 70 bases in reference sequence
        String leftseq = VariationUtils.joinRef(reference, (position - 70 > 1 ? position - 70 : 1), position - 1); // left 10 nt
        int chr0 = getOrElse(instance().chrLengths, region.chr, 0);
        //perl version: 1474
        //right 70 + deletionLength bases in reference sequence
        String tseq = VariationUtils.joinRef(reference, position, position + deletionLength + 70 > chr0 ? chr0 : position + deletionLength + 70);

        //Try to adjust for microsatellite instability
        Tuple.Tuple3<Double, Integer, String> tpl = findMSI(substr(tseq, 0, deletionLength), substr(tseq, deletionLength), leftseq);
        double msi = tpl._1;
        int shift3 = tpl._2;
        String msint = tpl._3;

        tpl = findMSI(leftseq, substr(tseq, deletionLength), leftseq);
        double tmsi = tpl._1;
        int tshift3 = tpl._2;
        String tmsint = tpl._3;
        //perl version: 1477
        if (msi < tmsi) {
            msi = tmsi;
            shift3 = tshift3;
            msint = tmsint;
        }
        if (msi <= shift3 / (double) deletionLength) {
            msi = shift3 / (double) deletionLength;
        }
        return Tuple.tuple(msi, shift3, msint);
    }

    /**
     * Proceesing non-reference variant if variant is an insertion
     *
     * @param position          position of variants
     * @param varDescriptionStr variant description string
     * @return tuple(msi, shift3, msint),
     * where msi: Microsatellite instability,
     * shift3: number of bases to be shifted to 3 prime,
     * msint: microsattelite nucleotide sequence; trim nucleotide(s) from the end
     */
    private Tuple.Tuple3<Double, Integer, String> proceedVrefIsInsertion(int position, String varDescriptionStr) {
        //variant description string without first symbol '+'
        String tseq1 = varDescriptionStr.substring(1);
        //perl version: 1448
        //left 50 bases in reference sequence
        String leftseq = VariationUtils.joinRef(reference, position - 50 > 1 ? position - 50 : 1, position); // left 10 nt
        int x = getOrElse(instance().chrLengths, region.chr, 0);
        //perl version: 1449
        //right 70 bases in reference sequence
        String tseq2 = VariationUtils.joinRef(reference, position + 1, (position + 70 > x ? x : position + 70));
        //perl version: 1450
        Tuple.Tuple3<Double, Integer, String> tpl = findMSI(tseq1, tseq2, leftseq);
        double msi = tpl._1;
        int shift3 = tpl._2;
        String msint = tpl._3;

        //Try to adjust for microsatellite instability
        tpl = findMSI(leftseq, tseq2, null);
        double tmsi = tpl._1;
        int tshift3 = tpl._2;
        String tmsint = tpl._3;
        //perl version: 1452
        if (msi < tmsi) {
            msi = tmsi;
            shift3 = tshift3;
            msint = tmsint;
        }
        if (msi <= shift3 / (double) tseq1.length()) {
            msi = shift3 / (double) tseq1.length();
        }
        return Tuple.tuple(msi, shift3, msint);
    }

    /**
     * The method of determining the same variation on reference
     *
     * @param insertionVariants variants with insertion
     * @param v                 variants at current position
     * @param position          position of variants
     * @return false - if not the same variation on reference, true - if the same
     */
    private boolean isTheSameVariationOnRef(Map<Integer, Map<String, Variation>> insertionVariants, int position, Map<String, Variation> v) {
        Set<String> vk = new HashSet<String>(v.keySet());
        if (insertionVariants.containsKey(position)) {
            vk.add("I");
        }
        //perl version: 1352
        Configuration config = instance().conf;
        if (vk.size() == 1 && reference.containsKey(position) && vk.contains(reference.get(position).toString())) {
            if (!config.doPileup && !config.bam.hasBam2() && config.ampliconBasedCalling == null) { // ignore if only reference were seen and no pileup to avoid computation
                return true;
            }
        }
        return false;
    }

    /**
     * Collect variants at position
     *
     * @param vars     - map where to collect Vars
     * @param position position of variants
     * @param variants array of all variants for the position
     * @return maximum variant frequency
     */
    //perl version: 1415
    private double collectVarsAtPosition(Map<Integer, Vars> vars, int position, List<Variant> variants) {
        double maxfreq = 0;
        for (Variant tvar : variants) {
            //If variant description string is 1-char base and it matches reference base at this position
            if (tvar.descriptionString.equals(String.valueOf(reference.get(position)))) {
                //this is a reference variant
                getOrPutVars(vars, position).referenceVariant = tvar;
            } else {
                //append variant to VAR and put it to VARN with key tvar.descriptionString (variant description string)
                getOrPutVars(vars, position).variants.add(tvar);
                getOrPutVars(vars, position).varDescriptionStringToVariants.put(tvar.descriptionString, tvar);
                if (tvar.frequency > maxfreq) {
                    maxfreq = tvar.frequency;
                }
            }
        }
        return maxfreq;
    }

    /**
     * Form Variant#DEBUG field from debug information
     *
     * @param vref reference variant
     */
    private void setDebugField(Variant vref) {
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
     * Sort variants by product of quality and coverage
     *
     * @param variantsForPosition array of all variants for the position
     */
    private void sortVariants(List<Variant> variantsForPosition) {
        //perl version: 1404
        //sort variants by product of quality and coverage
        Collections.sort(variantsForPosition, new Comparator<Variant>() {
            @Override
            public int compare(Variant o1, Variant o2) {
                return Double.compare(o2.meanQuality * o2.positionCoverage, o1.meanQuality * o1.positionCoverage);
            }
        });
    }

    /**
     * Create variants
     *
     * @param insertionVariants - insertion variants at position
     * @param position          - position of variants
     * @param totalCovarage     - total position coverage
     * @param hicov             - position coverage by high-quality reads
     * @param var               - array of all variants for the position
     * @param tmp               - temporary array used for debugging
     */
    private static void createInsertionVariants(Map<Integer, Map<String, Variation>> insertionVariants, int position, int totalCovarage, int hicov, List<Variant> var, List<String> tmp) {
        //perl version: 1393
        //Handle insertions separately
        Map<String, Variation> iv = insertionVariants.get(position);
        if (iv != null) {
            List<String> ikeys = new ArrayList<>(iv.keySet());
            Collections.sort(ikeys);
            // for (Entry<String, Variation> entV : iv.entrySet()) {
            //Loop over insertion variants
            for (String n : ikeys) {
                // String descriptionString = entV.getKey();
                Variation cnt = iv.get(n);
                //count of variants in forward strand
                int fwd = cnt.getCountForDir(false);
                //count of variants in reverse strand
                int rev = cnt.getCountForDir(true);
                //strand strandBiasFlag flag (0, 1 or 2)
                int bias = strandBias(fwd, rev, instance().conf.bias, instance().conf.minb);
                //mean base quality for variant
                double vqual = round(cnt.meanQuality / cnt.varsCount, 1); // base quality
                //mean mapping quality for variant
                double mq = cnt.meanMappingQuality / (double) cnt.varsCount; // mapping quality
                //number of high-quality reads for variant
                int hicnt = cnt.highQualityReadsCount;
                //number of low-quality reads for variant
                int locnt = cnt.lowQualityReadsCount;

//                    highQPosCoverage += highQualityReadsCount;

                //adjust position coverage if variant count is more than position coverage and no more than poistion coverage + extraCount
                int ttcov = totalCovarage;
                if (cnt.varsCount > totalCovarage && cnt.varsCount - totalCovarage < cnt.extraCount) {
                    ttcov = cnt.varsCount;
                }

                Variant tvref = createVariantRecord(hicov, n, cnt, fwd, rev, bias, vqual, mq, hicnt, locnt, ttcov);

                var.add(tvref);
                if (instance().conf.debug) {
                    tmp.add("I" + n
                            + ":" + (fwd + rev)
                            + ":F-" + fwd
                            + ":R-" + rev
                            + ":" + format("%.3f", tvref.frequency)
                            + ":" + tvref.strandBiasFlag
                            + ":" + tvref.meanPosition
                            + ":" + tvref.isAtLeastAt2Positions
                            + ":" + vqual
                            + ":" + tvref.hasAtLeast2DiffQualities
                            + ":" + format("%.3f", tvref.highQualityReadsFrequency)
                            + ":" + tvref.meanMappingQuality
                            + ":" + tvref.highQualityToLowQualityRatio);
                }

            }

        }
    }

    /**
     * Create variants
     *
     * @param variants      - variants at current position
     * @param totalCovarage - total position coverage
     * @param hicov         - position coverage by high-quality reads
     * @param var           - array of all variants for the position
     * @param tmp           - temporary array used for debugging
     */
    //perl version: 1370
    private static void createVariants(Map<String, Variation> variants, int totalCovarage, int hicov, List<Variant> var, List<String> tmp) {
        List<String> keys = new ArrayList<>(variants.keySet());
        //perl version: 1404
        //sort variants by product of quality and coverage
        Collections.sort(keys);

        //perl version: 1374
        //Loop over all varints found for the position except insertions
        for (String n : keys) {
            // for (Entry<String, Variation> entV : v.entrySet()) {
            // String descriptionString = entV.getKey();
            // Variation varsCount = entV.getValue();
            Variation cnt = variants.get(n);
            if (cnt.varsCount == 0) { //Skip variant if it does not have count
                continue;
            }
            //count of variants in forward strand
            int fwd = cnt.getCountForDir(false);
            //count of variants in reverse strand
            int rev = cnt.getCountForDir(true);

            //perl version: 1375
            //strand strandBiasFlag flag (0, 1 or 2)
            int bias = strandBias(fwd, rev, instance().conf.bias, instance().conf.minb);
            //mean base quality for variant
            double vqual = round(cnt.meanQuality / cnt.varsCount, 1); // base quality
            //mean mapping quality for variant
            double mq = cnt.meanMappingQuality / (double) cnt.varsCount; // mapping quality
            //number of high-quality reads for variant
            int hicnt = cnt.highQualityReadsCount;
            //number of low-quality reads for variant
            int locnt = cnt.lowQualityReadsCount;
            /**
             * Condition:
             # 1). varsCount.varsCount > totalPosCoverage                         - variant count is more than position coverage
             # 2). varsCount.varsCount - totalPosCoverage < varsCount.extraCount          - variant count is no more than poistion coverage + extraCount
             */
            //perl version: 1379
            int ttcov = totalCovarage;
            if (cnt.varsCount > totalCovarage && cnt.varsCount - totalCovarage < cnt.extraCount) { //adjust position coverage if condition holds
                ttcov = cnt.varsCount;
            }
            //perl version: 1380
            //create variant record
            Variant tvref = createVariantRecord(hicov, n, cnt, fwd, rev, bias, vqual, mq, hicnt, locnt, ttcov);


            //perl version: 1388
            //append variant record
            var.add(tvref);
            if (instance().conf.debug) { //debugging output
                tmp.add(n
                        + ":" + (fwd + rev)
                        + ":F-" + fwd
                        + ":R-" + rev
                        + ":" + format("%.3f", tvref.frequency)
                        + ":" + tvref.strandBiasFlag
                        + ":" + tvref.meanPosition
                        + ":" + tvref.isAtLeastAt2Positions
                        + ":" + vqual
                        + ":" + tvref.hasAtLeast2DiffQualities
                        + ":" + format("%.3f", tvref.highQualityReadsFrequency)
                        + ":" + tvref.meanMappingQuality
                        + ":" + tvref.highQualityToLowQualityRatio);
            }

        }
    }

    //perl version: 1380
    //create variant record
    private static Variant createVariantRecord(int hicov, String n, Variation cnt, int fwd, int rev, int bias, double vqual, double mq, int hicnt, int locnt, double ttcov) {
        Variant tvref = new Variant();
        tvref.descriptionString = n;
        tvref.positionCoverage = cnt.varsCount;
        tvref.varsCountOnForward = fwd;
        tvref.varsCountOnReverse = rev;
        tvref.strandBiasFlag = String.valueOf(bias);
        tvref.frequency = cnt.varsCount / ttcov;
        tvref.meanPosition = round(cnt.meanPosition / (double) cnt.varsCount, 1);
        tvref.isAtLeastAt2Positions = cnt.isAtLeastAt2Positions;
        tvref.meanQuality = vqual;
        tvref.hasAtLeast2DiffQualities = cnt.hasAtLeast2DiffQualities;
        tvref.meanMappingQuality = mq;
        tvref.highQualityToLowQualityRatio = hicnt / (locnt != 0 ? locnt : 0.5d);
        tvref.highQualityReadsFrequency = hicov > 0 ? hicnt / (double) hicov : 0;
        tvref.extraFrequency = cnt.extraCount != 0 ? cnt.extraCount / ttcov : 0;
        tvref.shift3 = 0;
        tvref.msi = 0;
        tvref.numberOfMismatches = cnt.numberOfMismatches / (double) cnt.varsCount;
        tvref.highQualityReadsCount = hicnt;
        tvref.highQPosCoverage = hicov;
        return tvref;
    }

    //perl version: 1364
    private static int calcHicov(Map<String, Variation> v, Map<String, Variation> iv) {
        int hicov = 0;
        for (Variation vr : v.values()) {
            hicov += vr.highQualityReadsCount;
        }
        if (iv != null) {
            for (Variation vr : iv.values()) {
                hicov += vr.highQualityReadsCount;
            }
        }
        return hicov;
    }

    /**
     * Find microsatellite instability
     * Tandemly repeated short sequence motifs ranging from 1– 6(8 in our case) base pairs are called microsatellites.
     * Other frequently used terms for these DNA regions are simple sequences or short tandem repeats (STRs)
     *
     * @param tseq1
     * @param tseq2
     * @param left
     * @return Tuple of (MSI count, No. of bases to be shifted to 3 prime for deletions due to alternative alignment, MicroSattelite unit length in base pairs)
     */
    //perl version: 1721
    static Tuple.Tuple3<Double, Integer, String> findMSI(String tseq1, String tseq2, String left) {

        //number of nucleotides in microsattelite
        int nmsi = 1;
        //number of bases to be shifted to 3 prime
        int shift3 = 0;
        String maxmsi = "";
        double msicnt = 0;
        while (nmsi <= tseq1.length() && nmsi <= 8) {
            //microsattelite nucleotide sequence; trim nucleotide(s) from the end
            String msint = substr(tseq1, -nmsi);
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
            //perl version: 1731
            double curmsi = msimatch.length() / (double) nmsi;
            mtch = Pattern.compile("^((" + msint + ")+)").matcher(tseq2);
            if (mtch.find()) {
                curmsi += mtch.group(1).length() / (double) nmsi;
            }
            if (curmsi > msicnt) {
                maxmsi = msint;
                msicnt = curmsi;
            }
            nmsi++;
        }
        //perl version: 1739
        String tseq = tseq1 + tseq2;
        //perl version: 1740
        while (shift3 < tseq2.length() && tseq.charAt(shift3) == tseq2.charAt(shift3)) {
            shift3++;
        }

        return tuple(msicnt, shift3, maxmsi);
    }

    /**
     * Calculate strand strandBiasFlag flag
     *
     * @param fwd  Variant count for forward strand
     * @param rev  Variant count for reverse strand
     * @param bias
     * @param minb minimum reads for strandBiasFlag
     * @return 0 - small total count, only one of strands, 1 - strand strandBiasFlag, 2 - no strand strandBiasFlag
     */
    //perl version: 2775
    static int strandBias(int fwd, int rev, double bias, int minb) {

        if (fwd + rev <= 12) { // using p=0.01, because prop.test(1,12) = 0.01
            return fwd * rev > 0 ? 2 : 0;
        }

        return (fwd / (double) (fwd + rev) >= bias && rev / (double) (fwd + rev) >= bias && fwd >= minb && rev >= minb) ? 2 : 1;
    }

    //Should be replaces with ap.putIfAbsent(position, new Vars());
    static Vars getOrPutVars(Map<Integer, Vars> map, int position) {
        Vars vars = map.get(position);
        if (vars == null) {
            vars = new Vars();
            map.put(position, vars);
        }
        return vars;
    }
}
