package com.astrazeneca.vardict.pipeline.postprocessmodules;

import com.astrazeneca.utils.Tuple;
import com.astrazeneca.vardict.Region;
import com.astrazeneca.vardict.pipeline.data.Scope;
import com.astrazeneca.vardict.pipeline.modules.ToVarsBuilder;
import com.astrazeneca.vardict.variations.Variant;
import com.astrazeneca.vardict.variations.VariationUtils;
import com.astrazeneca.vardict.variations.Vars;

import java.io.PrintStream;
import java.util.*;

import static com.astrazeneca.GlobalReadOnlyScope.instance;
import static com.astrazeneca.utils.Tuple.tuple;
import static com.astrazeneca.utils.Utils.join;
import static java.lang.String.format;

public class AmpliconPostProcessModule {

    private final static Comparator<Variant> VAR_TCOV_COMPARATOR = (o1, o2) -> Integer.compare(o2.totalPosCoverage, o1.totalPosCoverage);


    private final static Comparator<Tuple.Tuple2<Variant, String>> GVS_COMPARATOR = (a, b) -> Double.compare(b._1.frequency, a._1.frequency);

    /**
     * amplicon variant calling
     *
     * @param rg
     *            region
     * @param vars result of {@link ToVarsBuilder#process(Scope)} calling
     * @param positions map of positon =&gt; (list of (region number, region))
     * @param splice set of strings representing spliced regions
     */
    public static void process(Region rg, List<Map<Integer, Vars>> vars, Map<Integer, List<Tuple.Tuple2<Integer, Region>>> positions,
                               final Set<String> splice, PrintStream out) {

        List<Integer> pp = new ArrayList<>(positions.keySet());
        Collections.sort(pp);
        for (Integer p : pp) {

            final List<Tuple.Tuple2<Integer, Region>> v = positions.get(p);

            // good variants
            List<Tuple.Tuple2<Variant, String>> gvs = new ArrayList<>();
            //reference variants
            List<Variant> ref = new ArrayList<>();
            String nt = null;
            double maxaf = 0;
            //vartype may take values SNV (Single Nucleotide Variant), Complex (or MNV (Multiple Nucleotide Variant)), Insertion, Deletion
            Variant.Type vartype;
            boolean flag = false;
            Variant vref;
            //DNA sequencing coverage
            int nocov = 0;
            //max DNA sequencing coverage (max depth)
            int maxcov = 0;
            //good amplicon
            Set<String> goodmap = new HashSet<>();
            List<Integer> vcovs = new ArrayList<>();
            //amps map of amplicons.
            //An amplicon is a piece of DNA or RNA that is the source and/or product of
            //natural or artificial amplification or replication events.
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
                    Variant tv = l.get(0);
                    vcovs.add(tv.totalPosCoverage);
                    if (tv.totalPosCoverage > maxcov) {
                        maxcov = tv.totalPosCoverage;
                    }
                    vartype = tv.getType();
                    if (VariationUtils.isGoodVar(tv, refAmpP, vartype, splice)) {
                        gvs.add(tuple(tv, chr + ":" + S + "-" + E));
                        if (nt != null && !tv.descriptionString.equals(nt)) {
                            flag = true;
                        }
                        if (tv.frequency > maxaf) {
                            maxaf = tv.frequency;
                            nt = tv.descriptionString;
                        }
                        goodmap.add(format("%s-%s-%s", amp, tv.refAllele, tv.varAllele));
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
                if (t < maxcov / 50) { //The amplicon that has depth less than 1/50 of the max depth will be considered not working and thus not used.
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
                        vref = ref.get(0);
                    } else {
                        out.println(
                                join("\t", instance().sample, rg.gene, rg.chr, p, p, "", "", 0, 0, 0, 0, 0, 0, "", 0,
                                "0;0", 0, 0, 0, 0, 0, "", 0, 0, 0, 0, 0, 0, "", "", 0, 0,
                                rg.chr + ":", +p + "-" + p, "", 0, 0, 0, 0));
                        continue;
                    }
                } else {
                    continue;
                }
            } else {
                vref = gvs.get(0)._1;
            }
            if (flag) { // different good variants detected in different amplicons
                String gdnt = gvs.get(0)._1.descriptionString;
                List<Tuple.Tuple2<Variant, String>> gcnt = new ArrayList<>();
                for (Tuple.Tuple2<Integer, Region> amps : v) {
                    final Vars vtmp = vars.get(amps._1).get(p);
                    final Variant variant = vtmp == null ? null : vtmp.varDescriptionStringToVariants.get(gdnt);
                    if (variant != null && VariationUtils.isGoodVar(variant, vtmp.referenceVariant, null, splice)) {
                        gcnt.add(tuple(variant, amps._2.chr + ":" + amps._2.start + "-" + amps._2.end));
                    }
                }
                if (gcnt.size() == gvs.size()) {
                    flag = false;
                }
                Collections.sort(gcnt, GVS_COMPARATOR);
                gvs = gcnt;
            }

            //bad variants
            List<Tuple.Tuple2<Variant, String>> badv = new ArrayList<>();
            int gvscnt = gvs.size();
            if (gvscnt != v.size() || flag) {
                for (Tuple.Tuple2<Integer, Region> amps : v) {
                    int amp = amps._1;
                    Region reg = amps._2;
                    if (goodmap.contains(format("%s-%s-%s", amp, vref.refAllele, vref.varAllele))) {
                        continue;
                    }
                    // my $tref = $vars[$amp]->{ $p }->{ VAR }->[0]; ???
                    if (vref.startPosition >= reg.variantsStart && vref.endPosition <= reg.variantsEnd) {

                        String regStr = reg.chr + ":" + reg.start + "-" + reg.end;

                        if (vars.get(amp).containsKey(p) && vars.get(amp).get(p).variants.size() > 0) {
                            badv.add(tuple(vars.get(amp).get(p).variants.get(0), regStr));
                        } else if (vars.get(amp).containsKey(p) && vars.get(amp).get(p).referenceVariant != null) {
                            badv.add(tuple(vars.get(amp).get(p).referenceVariant, regStr));
                        } else {
                            badv.add(tuple(null, regStr));
                        }
                    } else if ((vref.startPosition < reg.variantsEnd && reg.variantsEnd < vref.endPosition)
                            || (vref.startPosition < reg.variantsStart && reg.variantsStart < vref.endPosition)) { // the variant overlap with amplicon's primer
                        if (gvscnt > 1)
                            gvscnt--;
                    }
                }
            }
            if (flag && gvscnt < gvs.size()) {
                flag = false;
            }
            vartype = vref.getType();
            if (vartype == Variant.Type.complex) {
                vref.adjComplex();
            }
            out.print(join("\t", instance().sample, rg.gene, rg.chr,
                    vref.addDelimiterExtended("\t"), gvs.get(0)._2, vartype, gvscnt, gvscnt + badv.size(), nocov, flag ? 1 : 0));
            if (instance().conf.debug) {
                out.print("\t" + vref.DEBUG);
            }
            if (instance().conf.debug) {
                for (int gvi = 0; gvi < gvs.size(); gvi++) {
                    Tuple.Tuple2<Variant, String> tp = gvs.get(gvi);
                    out.print("\tGood" + gvi + " " + join(" ", tp._1.addDelimiter(" "), tp._2));
                }
                for (int bvi = 0; bvi < badv.size(); bvi++) {
                    Tuple.Tuple2<Variant, String> tp = badv.get(bvi);
                    out.print("\tBad" + bvi + " " + join(" ", tp._1.addDelimiter(" "), tp._2));
                }
            }
            out.println();
        }
    }
}
