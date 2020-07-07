package com.astrazeneca.vardict.postprocessmodules;

import com.astrazeneca.vardict.collection.Tuple;
import com.astrazeneca.vardict.data.Region;
import com.astrazeneca.vardict.printers.AmpliconOutputVariant;
import com.astrazeneca.vardict.printers.VariantPrinter;
import com.astrazeneca.vardict.variations.Variant;
import com.astrazeneca.vardict.variations.Vars;

import java.util.*;

import static com.astrazeneca.vardict.Utils.printExceptionAndContinue;
import static com.astrazeneca.vardict.data.scopedata.GlobalReadOnlyScope.instance;
import static com.astrazeneca.vardict.collection.Tuple.tuple;
import static java.lang.String.format;

/**
 * Class for preparation of variants found in amplicon analysis to the output
 */
public class AmpliconPostProcessModule {
    private final static Comparator<Variant> VAR_TCOV_COMPARATOR = (o1, o2) -> Integer.compare(o2.totalPosCoverage, o1.totalPosCoverage);
    private final static Comparator<Tuple.Tuple2<Variant, String>> GVS_COMPARATOR = (a, b) -> Double.compare(b._1.frequency, a._1.frequency);

    /**
     * Method processed Vars for each amplicon on position.
     * @param rg region
     * @param vars result of calling toVarsBuilder
     * @param ampliconsOnPositions map of position =&gt; (list of (region number, region))
     * @param splice set of strings representing spliced regions
     * @param variantPrinter specific variant printer where output must be print
     */
    public void process(Region rg, List<Map<Integer, Vars>> vars,
                        Map<Integer, List<Tuple.Tuple2<Integer, Region>>> ampliconsOnPositions,
                        final Set<String> splice, VariantPrinter variantPrinter) {

        List<Integer> positions = new ArrayList<>(ampliconsOnPositions.keySet());
        Collections.sort(positions);
        int lastPosition = 0;
        for (Integer position : positions) {
            try {
                lastPosition = position;
                final List<Tuple.Tuple2<Integer, Region>> ampliconRegions = ampliconsOnPositions.get(position);
                // good variants
                List<Tuple.Tuple2<Variant, String>> gvs = new ArrayList<>();
                //reference variants
                List<Variant> ref = new ArrayList<>();

                List<Variant> vrefList = new ArrayList<>();
                //good amplicon
                Set<String> goodmap = new HashSet<>();
                List<Integer> vcovs = new ArrayList<>();
                //Contains list of good variants on each amplicon in position
                Map<Integer, List<Variant>> goodVariantsOnAmp = new LinkedHashMap<>();

                //DNA sequencing coverage
                int nocov = 0;
                //max DNA sequencing coverage (max depth)
                int maxcov = 0;
                double maxaf = 0;

                for (Tuple.Tuple2<Integer, Region> amps : ampliconRegions) {
                    final int ampliconNumber = amps._1;
                    //chromosome name
                    final String chr = amps._2.chr;
                    //start index
                    final int start = amps._2.start;
                    //end index
                    final int end = amps._2.end;

                    Vars vtmp = vars.get(ampliconNumber).get(position);
                    List<Variant> variantsOnAmplicon = vtmp == null ? null : vtmp.variants;
                    Variant refAmplicon = vtmp == null ? null : vtmp.referenceVariant;

                    if (variantsOnAmplicon != null && !variantsOnAmplicon.isEmpty()) {
                        List<Variant> goodVars = new ArrayList<>();
                        for (Variant tv : variantsOnAmplicon) {
                            vcovs.add(tv.totalPosCoverage);
                            if (tv.totalPosCoverage > maxcov) {
                                maxcov = tv.totalPosCoverage;
                            }
                            // vartype may take values SNV (Single Nucleotide Variant),
                            // Complex (or MNV (Multiple Nucleotide Variant)), Insertion, Deletion
                            tv.vartype = tv.varType();
                            if (tv.isGoodVar(refAmplicon, tv.vartype, splice)) {
                                gvs.add(tuple(tv, chr + ":" + start + "-" + end));
                                goodVars.add(tv);
                                goodVariantsOnAmp.put(ampliconNumber, goodVars);
                                if (tv.frequency > maxaf) {
                                    maxaf = tv.frequency;
                                }
                                goodmap.add(format("%s-%s-%s", ampliconNumber, tv.refallele, tv.varallele));
                            }
                        }
                    } else if (refAmplicon != null) {
                        vcovs.add(refAmplicon.totalPosCoverage);
                    } else {
                        vcovs.add(0);
                    }
                    if (refAmplicon != null) {
                        ref.add(refAmplicon);
                    }
                }

                //Depth (coverage) in DNA sequencing refers to the number of times a nucleotide is read during the sequencing process.
                // Coverage is the average number of reads representing a given nucleotide in the reconstructed sequence.
                for (int t : vcovs) {
                    if (t < maxcov / (double) 50) { //The amplicon that has depth less than 1/50 of the max depth will be considered not working and thus not used.
                        nocov++;
                    }
                }

                if (gvs.size() > 1) {
                    Collections.sort(gvs, GVS_COMPARATOR);
                }
                if (ref.size() > 1) {
                    Collections.sort(ref, VAR_TCOV_COMPARATOR);
                }

                if (gvs.isEmpty()) { // Only reference
                    if (instance().conf.doPileup) {
                        if (!ref.isEmpty()) {
                            //vref = ref.get(0);
                            vrefList.add(ref.get(0));
                        } else {
                            AmpliconOutputVariant outputVariant = new AmpliconOutputVariant(null, rg, null, null, position, 0, nocov, false);
                            variantPrinter.print(outputVariant);
                            continue;
                        }
                    } else {
                        continue;
                    }
                } else {
                    fillVrefList(gvs, vrefList);
                }
                boolean flag = isAmpBiasFlag(goodVariantsOnAmp);

                List<Tuple.Tuple2<Variant, String>> goodVariants = gvs;
                for (int i = 0; i < vrefList.size(); i++) {
                    Variant vref = vrefList.get(i);
                    if (flag) { // different good variants detected in different amplicons
                        String gdnt = gvs.get(0)._1.descriptionString;
                        List<Tuple.Tuple2<Variant, String>> gcnt = new ArrayList<>();
                        for (Tuple.Tuple2<Integer, Region> amps : ampliconRegions) {
                            final Vars vtmp = vars.get(amps._1).get(position);
                            final Variant variant = vtmp == null ? null : vtmp.varDescriptionStringToVariants.get(gdnt);
                            if (variant != null && variant.isGoodVar(vtmp.referenceVariant, null, splice)) {
                                gcnt.add(tuple(variant, amps._2.chr + ":" + amps._2.start + "-" + amps._2.end));
                            }
                        }
                        if (gcnt.size() == gvs.size()) {
                            flag = false;
                        }
                        Collections.sort(gcnt, GVS_COMPARATOR);
                        goodVariants = gcnt;
                    }
                    int initialGvscnt = countVariantOnAmplicons(vref, goodVariantsOnAmp);
                    int currentGvscnt = initialGvscnt;
                    //bad variants
                    List<Tuple.Tuple2<Variant, String>> badVariants = new ArrayList<>();
                    if (initialGvscnt != ampliconRegions.size() || flag) {
                        for (Tuple.Tuple2<Integer, Region> amps : ampliconRegions) {
                            int amp = amps._1;
                            Region reg = amps._2;
                            if (goodmap.contains(format("%s-%s-%s", amp, vref.refallele, vref.varallele))) {
                                continue;
                            }
                            if (instance().conf.doPileup && vref.refallele.equals(vref.varallele)) {
                                continue;
                            }
                            if (vref.startPosition >= reg.insertStart && vref.endPosition <= reg.insertEnd) {
                                String regStr = reg.chr + ":" + reg.start + "-" + reg.end;

                                if (vars.get(amp).containsKey(position) && vars.get(amp).get(position).variants.size() > 0) {
                                    badVariants.add(tuple(vars.get(amp).get(position).variants.get(0), regStr));
                                } else if (vars.get(amp).containsKey(position) && vars.get(amp).get(position).referenceVariant != null) {
                                    badVariants.add(tuple(vars.get(amp).get(position).referenceVariant, regStr));
                                } else {
                                    badVariants.add(tuple(null, regStr));
                                }
                            } else if ((vref.startPosition < reg.insertEnd && reg.insertEnd < vref.endPosition)
                                    || (vref.startPosition < reg.insertStart && reg.insertStart < vref.endPosition)) { // the variant overlap with amplicon's primer
                                if (currentGvscnt > 1)
                                    currentGvscnt--;
                            }
                        }
                    }
                    if (flag && currentGvscnt < initialGvscnt) {
                        flag = false;
                    }
                    vref.vartype = vref.varType();
                    if (vref.vartype.equals("Complex")) {
                        vref.adjComplex();
                    }
                    AmpliconOutputVariant outputVariant = new AmpliconOutputVariant(vref, rg, goodVariants, badVariants, position, currentGvscnt, nocov, flag);
                    variantPrinter.print(outputVariant);
                }
            } catch (Exception exception) {
                printExceptionAndContinue(exception, "position", String.valueOf(lastPosition), rg);
            }
        }
    }
    /**
     * Count amplicons where good variant appears
     * @param vref variant to check
     * @param goodVariantsOnAmp map of amplicons and lists of variants
     * @return number of amplicons
     */
    private int countVariantOnAmplicons(Variant vref, Map<Integer, List<Variant>> goodVariantsOnAmp) {
        int gvscnt = 0;
        for (Map.Entry<Integer, List<Variant>> entry: goodVariantsOnAmp.entrySet()) {
            List<Variant> variants = entry.getValue();
            for(Variant variant : variants) {
                if (variant.refallele.equals(vref.refallele) && variant.varallele.equals(vref.varallele)) {
                    gvscnt++;
                }
            }
        }
        return gvscnt;
    }
    /**
     * If variant with the same varallele and refallele is already added to output list, skip it.
     * The variant with the biggest frequency will be added.
     * Variant must be skipped to avoid duplicates because identical variants can be on different amplicons.
     * @param gvs good variants per start-end
     * @param vrefList list of variants on all amplicons in position
     */
    private void fillVrefList(List<Tuple.Tuple2<Variant, String>> gvs, List<Variant> vrefList) {
        for (Tuple.Tuple2<Variant, String> goodVariant : gvs) {
            boolean variantWasAdded = false;
            Variant variantToAdd = goodVariant._1;
            for (Variant var : vrefList) {
                if (var.varallele.equals(goodVariant._1.varallele) && var.refallele.equals(goodVariant._1.refallele)) {
                    variantWasAdded = true;
                }
            }
            if (!variantWasAdded) vrefList.add(variantToAdd);
        }
    }
    /**
     * Determine if different amplicon contains different variants and set AMPBIAS flag.
     * Check if each amplicon on position contains identical lists of Variants. If some variants are absent between
     * amplicons, or differ by variant description string, returns true.
     * @param goodVariantsOnAmp map amplicons on list of its variants.
     * @return true if variants are differ, false if not
     */
    private boolean isAmpBiasFlag(Map<Integer, List<Variant>> goodVariantsOnAmp) {
        if (goodVariantsOnAmp.isEmpty()) return false;
        List<Integer> ampliconList = new ArrayList<>(goodVariantsOnAmp.keySet());
        Collections.sort(ampliconList);
        int ampliconLength = ampliconList.size() - 1;

        for (int i = 0; i < ampliconLength; i++) {
            int currentAmplicon = ampliconList.get(i);
            int nextAmplicon = ampliconList.get(i + 1);
            List<Variant> goodVariantListCurrentAmp = goodVariantsOnAmp.get(currentAmplicon);
            List<Variant> goodVariantListNextAmp = goodVariantsOnAmp.get(nextAmplicon);

            if (goodVariantListNextAmp == null || goodVariantListCurrentAmp.size() != goodVariantListNextAmp.size()) {
                return true;
            }
            Collections.sort(goodVariantListCurrentAmp, VAR_TCOV_COMPARATOR);
            Collections.sort(goodVariantListNextAmp, VAR_TCOV_COMPARATOR);
            for (int j = 0; j < goodVariantListCurrentAmp.size(); j++) {
                Variant var1 = goodVariantListCurrentAmp.get(j);
                Variant var2 = goodVariantListNextAmp.get(j);
                if (!var1.descriptionString.equals(var2.descriptionString)) {
                    return true;
                }
            }
        }
        return false;
    }
}
