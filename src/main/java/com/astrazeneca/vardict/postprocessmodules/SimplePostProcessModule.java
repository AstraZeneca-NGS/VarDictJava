package com.astrazeneca.vardict.postprocessmodules;

import com.astrazeneca.vardict.Configuration;
import com.astrazeneca.vardict.data.scopedata.AlignedVarsData;
import com.astrazeneca.vardict.data.scopedata.Scope;
import com.astrazeneca.vardict.variations.Variant;
import com.astrazeneca.vardict.printers.SimpleOutputVariant;
import com.astrazeneca.vardict.printers.VariantPrinter;
import com.astrazeneca.vardict.variations.Vars;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.function.Consumer;

import static com.astrazeneca.vardict.Utils.printExceptionAndContinue;
import static com.astrazeneca.vardict.data.scopedata.GlobalReadOnlyScope.instance;

/**
 * Class for preparation of variants found in simple analysis to the output
 */
public class SimplePostProcessModule implements Consumer<Scope<AlignedVarsData>> {
    private VariantPrinter variantPrinter;

    public SimplePostProcessModule(VariantPrinter variantPrinter) {
        this.variantPrinter = variantPrinter;
    }

    /**
     * Single sample mode variant calling. Creates variant output for variants from aligned map and prints them.
     */
    @Override
    public void accept(Scope<AlignedVarsData> mapScope) {
        int lastPosition = 0;
        for (Map.Entry<Integer, Vars> ent : mapScope.data.alignedVariants.entrySet()) {
            try {
                int position = ent.getKey();
                lastPosition = position;
                Vars variantsOnPosition = ent.getValue();

                Configuration conf = instance().conf;

                List<Variant> vrefs = new ArrayList<>();
                if (variantsOnPosition.sv.isEmpty()) {
                    if (position < mapScope.region.start || position > mapScope.region.end) {
                        continue;
                    }
                }

                if (variantsOnPosition.variants.isEmpty()) {
                    if (!conf.doPileup) {
                        continue;
                    }
                    Variant vref = variantsOnPosition.referenceVariant;
                    if (vref == null) {
                        SimpleOutputVariant outputVariant = new SimpleOutputVariant(vref, mapScope.region, variantsOnPosition.sv, position);
                        variantPrinter.print(outputVariant);
                        continue;
                    }
                    vref.vartype = "";
                    vrefs.add(vref);
                } else {
                    List<Variant> vvar = variantsOnPosition.variants;
                    for (Variant vref : vvar) {
                        if (vref.refallele.contains("N")) {
                            continue;
                        }
                        if (vref.refallele.equals(vref.varallele)) {
                            if (!conf.doPileup) {
                                continue;
                            }
                        }
                        if (vref.startPosition != position && conf.doPileup && vvar.size() == 1) {
                            Variant refVar = variantsOnPosition.referenceVariant;
                            if (refVar == null) {
                                SimpleOutputVariant outputVariant = new SimpleOutputVariant(refVar, mapScope.region, variantsOnPosition.sv, position);
                                variantPrinter.print(outputVariant);
                                refVar = new Variant();
                            }
                            refVar.vartype = "";
                            vrefs.add(refVar);
                        }
                        vref.vartype = vref.varType();
                        if (!vref.isGoodVar(variantsOnPosition.referenceVariant, vref.vartype, mapScope.splice)) {
                            if (!conf.doPileup) {
                                continue;
                            }
                        }
                        vrefs.add(vref);
                    }
                }
                for (int vi = 0; vi < vrefs.size(); vi++) {
                    Variant vref = vrefs.get(vi);
                    if ("Complex".equals(vref.vartype)) {
                        vref.adjComplex();
                    }
                    if (instance().conf.crisprCuttingSite == 0) {
                        vref.crispr = 0;
                    }
                    SimpleOutputVariant outputVariant = new SimpleOutputVariant(vref, mapScope.region, variantsOnPosition.sv, position);
                    variantPrinter.print(outputVariant);
                }
            } catch (Exception exception) {
                printExceptionAndContinue(exception, "position", String.valueOf(lastPosition), mapScope.region);
            }
        }
    }
}
