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
     * @param positions map of positon =&gt; (list of (region number, region))
     * @param splice set of strings representing spliced regions
     * @param variantPrinter specific variant printer where output must be print
     */
    public void process(Region rg, List<Map<Integer, Vars>> vars, Map<Integer, List<Tuple.Tuple2<Integer, Region>>> positions,
                               final Set<String> splice, VariantPrinter variantPrinter) {

        List<Integer> pp = new ArrayList<>(positions.keySet());
        Collections.sort(pp);
        int lastPosition = 0;
        for (Integer p : pp) {
            try {
                lastPosition = p;
                final List<Tuple.Tuple2<Integer, Region>> v = positions.get(p);
                // good variants
                List<Tuple.Tuple2<Variant, String>> gvs = new ArrayList<>();
                //reference variants
                List<Variant> ref = new ArrayList<>();
                double maxaf = 0;
                // vartype may take values SNV (Single Nucleotide Variant), Complex (or MNV (Multiple Nucleotide Variant)), Insertion, Deletion
                Variant vref;
                List<Variant> vrefList = new ArrayList<>();
                //DNA sequencing coverage
                int nocov = 0;
                //max DNA sequencing coverage (max depth)
                int maxcov = 0;
                //good amplicon
                Set<String> goodmap = new HashSet<>();
                List<Integer> vcovs = new ArrayList<>();
                //Contains list of good variants on each amplicon in position
                Map<Integer, List<Variant>> goodVariantsOnAmp = new LinkedHashMap<>();
                for (Tuple.Tuple2<Integer, Region> amps : v) {
                    final int amp = amps._1;
                    //chromosome name
                    final String chr = amps._2.chr;
                    //start index
                    final int S = amps._2.start;
                    //end index
                    final int E = amps._2.end;

                    Vars vtmp = vars.get(amp).get(p);
                    List<Variant> l = vtmp == null ? null : vtmp.variants;
                    Variant refAmpP = vtmp == null ? null : vtmp.referenceVariant;
                    if (l != null && !l.isEmpty()) {
                        List<Variant> goodVars = new ArrayList<>();
                        for (Variant tv : l) {
                            vcovs.add(tv.totalPosCoverage);
                            if (tv.totalPosCoverage > maxcov) {
                                maxcov = tv.totalPosCoverage;
                            }
                            tv.vartype = tv.varType();
                            if (tv.isGoodVar(refAmpP, tv.vartype, splice)) {
                                gvs.add(tuple(tv, chr + ":" + S + "-" + E));
                                goodVars.add(tv);
                                goodVariantsOnAmp.put(amp, goodVars);
                                if (tv.frequency > maxaf) {
                                    maxaf = tv.frequency;
                                }
                                goodmap.add(format("%s-%s-%s", amp, tv.refallele, tv.varallele));
                            }
                        }
                    } else if (refAmpP != null) {
                        vcovs.add(refAmpP.totalPosCoverage);
                    } else {
                        vcovs.add(0);
                    }
                    if (refAmpP != null) {
                        ref.add(refAmpP);
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

                if (gvs.isEmpty()) { // Only referenece
                    if (instance().conf.doPileup) {
                        if (!ref.isEmpty()) {
                            //vref = ref.get(0);
                            vrefList.add(ref.get(0));
                        } else {
                            AmpliconOutputVariant outputVariant = new AmpliconOutputVariant(null, rg, null, null, p, 0, nocov, false);
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
                    vref = vrefList.get(i);
                    if (flag) { // different good variants detected in different amplicons
                        String gdnt = gvs.get(0)._1.descriptionString;
                        List<Tuple.Tuple2<Variant, String>> gcnt = new ArrayList<>();
                        for (Tuple.Tuple2<Integer, Region> amps : v) {
                            final Vars vtmp = vars.get(amps._1).get(p);
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
                    if (initialGvscnt != v.size() || flag) {
                        for (Tuple.Tuple2<Integer, Region> amps : v) {
                            int amp = amps._1;
                            Region reg = amps._2;
                            if (goodmap.contains(format("%s-%s-%s", amp, vref.refallele, vref.varallele))) {
                                continue;
                            }
                            // my $tref = $vars[$amp]->{ $p }->{ VAR }->[0]; ???
                            if (vref.startPosition >= reg.insertStart && vref.endPosition <= reg.insertEnd) {

                                String regStr = reg.chr + ":" + reg.start + "-" + reg.end;

                                if (vars.get(amp).containsKey(p) && vars.get(amp).get(p).variants.size() > 0) {
                                    badVariants.add(tuple(vars.get(amp).get(p).variants.get(0), regStr));
                                } else if (vars.get(amp).containsKey(p) && vars.get(amp).get(p).referenceVariant != null) {
                                    badVariants.add(tuple(vars.get(amp).get(p).referenceVariant, regStr));
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
                    AmpliconOutputVariant outputVariant = new AmpliconOutputVariant(vref, rg, goodVariants, badVariants, p, currentGvscnt, nocov, flag);
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
    private static int countVariantOnAmplicons(Variant vref, Map<Integer, List<Variant>> goodVariantsOnAmp) {
        int gvscnt = 0;
        for (Map.Entry<Integer, List<Variant>> entry: goodVariantsOnAmp.entrySet()) {
            List<Variant> variants = entry.getValue();
            for(Variant variant : variants) {
                if (variant.equals(vref)) {
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
    private static void fillVrefList(List<Tuple.Tuple2<Variant, String>> gvs, List<Variant> vrefList) {
        for (Tuple.Tuple2<Variant, String> goodVariant : gvs) {
            boolean variantWasAdded = false;
            for (Variant var : vrefList) {
                if (var.varallele.equals(goodVariant._1.varallele)
                        && var.refallele.equals(goodVariant._1.refallele)) {
                    variantWasAdded = true;
                }
            }
            if (!variantWasAdded) vrefList.add(goodVariant._1);
        }
    }
    /**
     * Determine if different amplicon contains different variants and set AMPBIAS flag.
     * Check if each amplicon on position contains identical lists of Variants. If some variants are absent between
     * amplicons, or differ by variant description string, returns true.
     * @param goodVariantsOnAmp map amplicons on list of its variants.
     * @return true if variants are differ, false if not
     */
    private static boolean isAmpBiasFlag(Map<Integer, List<Variant>> goodVariantsOnAmp) {
        if (goodVariantsOnAmp.isEmpty()) return false;
        for (int i = goodVariantsOnAmp.keySet().iterator().next(); i < goodVariantsOnAmp.keySet().size() - 1; i++) {
            List<Variant> goodVariantListFirst = goodVariantsOnAmp.get(i);
            List<Variant> goodVariantListSecond = goodVariantsOnAmp.get(i + 1);

            if (goodVariantListSecond == null || goodVariantListFirst.size() != goodVariantListSecond.size()) {
                return true;
            }
            Collections.sort(goodVariantListFirst, VAR_TCOV_COMPARATOR);
            Collections.sort(goodVariantListSecond, VAR_TCOV_COMPARATOR);
            for (int j = 0; j < goodVariantListFirst.size(); j++) {
                Variant var1 = goodVariantListFirst.get(j);
                Variant var2 = goodVariantListSecond.get(j);
                if (!var1.descriptionString.equals(var2.descriptionString)) {
                    return true;
                }
            }
        }
        return false;
    }

}
