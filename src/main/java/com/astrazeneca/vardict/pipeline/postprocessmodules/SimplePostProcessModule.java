package com.astrazeneca.vardict.pipeline.postprocessmodules;

import com.astrazeneca.utils.Tuple;
import com.astrazeneca.vardict.Configuration;
import com.astrazeneca.vardict.pipeline.data.Scope;
import com.astrazeneca.vardict.variations.Variant;
import com.astrazeneca.vardict.variations.VariationUtils;
import com.astrazeneca.vardict.variations.Vars;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.function.Consumer;

import static com.astrazeneca.GlobalReadOnlyScope.instance;
import static com.astrazeneca.utils.Utils.join;
import static com.astrazeneca.vardict.variations.VariationUtils.VarsType.ref;
import static java.lang.String.format;

public class SimplePostProcessModule implements Consumer<Scope<Tuple.Tuple2<Integer, Map<Integer, Vars>>>> {
    /**
     * Single sample mode variant calling
     */
    //perl version: 305
    @Override
    public void accept(Scope<Tuple.Tuple2<Integer, Map<Integer, Vars>>> mapScope) {
        for (int p = mapScope.region.start; p <= mapScope.region.end; p++) {
            Configuration conf = instance().conf;
            String sample = instance().sample;

            List<Variant.Type> vts = new ArrayList<>();
            List<Variant> vrefs = new ArrayList<>();
            if (!mapScope.data._2.containsKey(p) || mapScope.data._2.get(p).variants.isEmpty()) {
                if (!conf.doPileup) {
                    continue;
                }
                Variant vref = VariationUtils.getVarMaybe(mapScope.data._2, p, ref);
                if (vref == null) {
                    mapScope.out.println(join("\t", sample, mapScope.region.gene, mapScope.region.chr, p, p,
                            "", "", 0, 0, 0, 0, 0, 0, "", 0, "0;0", 0, 0, 0, 0, 0, "", 0, 0, 0, 0, 0, 0, "", "", 0, 0,
                            mapScope.region.chr + ":" + mapScope.region.start + "-" + mapScope.region.end, ""
                    ));
                    continue;
                }
                vts.add(Variant.Type.noInfo);
                vrefs.add(vref);
            } else {
                List<Variant> vvar = mapScope.data._2.get(p).variants;
                Variant rref = VariationUtils.getVarMaybe(mapScope.data._2, p, ref);
                for (Variant vref : vvar) {
                    if (vref.refAllele.contains("N")) {
                        continue;
                    }
                    Variant.Type vartype = vref.getType();
                    if (!VariationUtils.isGoodVar(vref, rref, vartype, mapScope.splice)) {
                        if (!conf.doPileup) {
                            continue;
                        }
                    }

                    vts.add(vartype);
                    vrefs.add(vref);
                }
            }
            for (int vi = 0; vi < vts.size(); vi++) {
                Variant.Type vartype = vts.get(vi);
                Variant vref = vrefs.get(vi);
                if (vartype == Variant.Type.complex) {
                    vref.adjComplex();
                }
                mapScope.out.println(join("\t", sample, mapScope.region.gene, mapScope.region.chr,
                        vref.startPosition, vref.endPosition, vref.refAllele, vref.varAllele,
                        vref.totalPosCoverage,
                        vref.positionCoverage,
                        vref.refForwardCoverage,
                        vref.refReverseCoverage,
                        vref.varsCountOnForward,
                        vref.varsCountOnReverse,
                        vref.genotype,
                        format("%.4f", vref.frequency),
                        vref.strandBiasFlag,
                        format("%.1f", vref.meanPosition),
                        vref.isAtLeastAt2Positions ? 1 : 0,
                        format("%.1f", vref.meanQuality),
                        vref.hasAtLeast2DiffQualities ? 1 : 0,
                        format("%.1f", vref.meanMappingQuality),
                        format("%.3f", vref.highQualityToLowQualityRatio),
                        vref.highQualityReadsFrequency == 0 ? 0 : format("%.4f", vref.highQualityReadsFrequency),
                        vref.extraFrequency == 0 ? 0 : format("%.4f", vref.extraFrequency),
                        vref.shift3,
                        vref.msi == 0 ? 0 : format("%.3f", vref.msi),
                        vref.msint,
                        format("%.1f", vref.numberOfMismatches),
                        vref.highQualityReadsCount,
                        vref.highQPosCoverage,
                        vref.precedingRefSequence.isEmpty() ? "0" : vref.precedingRefSequence,
                        vref.followingRefSequence.isEmpty() ? "0" : vref.followingRefSequence,
                        mapScope.region.chr + ":" + mapScope.region.start + "-" + mapScope.region.end, vartype
                ));
                if (conf.debug) {
                    mapScope.out.println("\t" + vref.DEBUG);
                }
            }

        }
    }
}
