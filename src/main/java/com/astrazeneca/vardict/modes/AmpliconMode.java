package com.astrazeneca.vardict.modes;

import com.astrazeneca.vardict.data.Reference;
import com.astrazeneca.vardict.data.ReferenceResource;
import com.astrazeneca.vardict.collection.ConcurrentHashSet;
import com.astrazeneca.vardict.collection.DirectThreadExecutor;
import com.astrazeneca.vardict.collection.Tuple;
import com.astrazeneca.vardict.data.Region;
import com.astrazeneca.vardict.data.scopedata.AlignedVarsData;
import com.astrazeneca.vardict.data.scopedata.InitialData;
import com.astrazeneca.vardict.data.scopedata.Scope;
import com.astrazeneca.vardict.postprocessmodules.AmpliconPostProcessModule;
import com.astrazeneca.vardict.printers.VariantPrinter;
import com.astrazeneca.vardict.variations.Vars;

import java.io.ByteArrayOutputStream;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.*;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutionException;

import static com.astrazeneca.vardict.Utils.join;
import static com.astrazeneca.vardict.data.scopedata.GlobalReadOnlyScope.instance;
import static com.astrazeneca.vardict.collection.Tuple.tuple;

/**
 * Amplicon variant calling
 * Mode starting Amplicon analysis: BAM file is processed by amplicons in sequence, then received Vars by each amplicon processed in
 * AmpliconPostProcessModule.
 *
 * An amplicon is a piece of DNA or RNA that is the source and/or product of
 * natural or artificial amplification or replication events.
 * */
public class AmpliconMode extends AbstractMode {

    public AmpliconMode(List<List<Region>> segments, ReferenceResource referenceResource) {
        super(segments, referenceResource);
        printHeader();
    }

    /**
     * In not parallel mode each region will be processed in sequence.
     */
    @Override
    public void notParallel() {
        VariantPrinter variantPrinter = VariantPrinter.createPrinter(instance().printerTypeOut);

        for (List<Region> regions : segments) {
            Map<Integer, List<Tuple.Tuple2<Integer, Region>>> pos = new HashMap<>();
            int ampliconNumber = 0;
            Region currentRegion = regions.get(0);
            final Set<String> splice = new HashSet<>();
            List<Map<Integer, Vars>> vars = new ArrayList<>();
            for (Region region : regions) {
                currentRegion = region;
                for (int p = region.insertStart; p <= region.insertEnd; p++) {
                    List<Tuple.Tuple2<Integer, Region>> list = pos.computeIfAbsent(p, k -> new ArrayList<>());
                    list.add(tuple(ampliconNumber, region));
                }
                Scope<InitialData> initialScope = new Scope<>(instance().conf.bam.getBam1(), region,
                        tryToGetReference(region), referenceResource, 0, splice,
                        variantPrinter, new InitialData());
                CompletableFuture<Scope<AlignedVarsData>> pipeline = pipeline(initialScope,
                        new DirectThreadExecutor());

                AlignedVarsData data = pipeline.join().data;
                vars.add(data.alignedVariants);
                ampliconNumber++;
            }
            new AmpliconPostProcessModule().process(currentRegion, vars, pos, splice, variantPrinter);
        }
    }

    /**
     * In parallel mode workers are created for each region and are processed in parallel.
     */
    @Override
    protected AbstractParallelMode createParallelMode() {
        return new AbstractParallelMode() {
            @Override
            void produceTasks() throws InterruptedException, ExecutionException {
                for (List<Region> regions : segments) {
                    Map<Integer, List<Tuple.Tuple2<Integer, Region>>> pos = new HashMap<>();
                    int j = 0;
                    Region currentRegion = regions.get(0);
                    final Set<String> splice = new ConcurrentHashSet<>();
                    List<CompletableFuture<Scope<AlignedVarsData>>> workers = new ArrayList<>(regions.size() - 1);

                    for (Region region : regions) {
                        currentRegion = region;
                        for (int p = region.insertStart; p <= region.insertEnd; p++) {
                            List<Tuple.Tuple2<Integer, Region>> list = pos.computeIfAbsent(p, k -> new ArrayList<>());
                            list.add(tuple(j, region));
                        }
                        VariantPrinter variantPrinter = VariantPrinter.createPrinter(instance().printerTypeOut);
                        Reference reference = tryToGetReference(region);
                        Scope<InitialData> initialScope = new Scope<>(instance().conf.bam.getBam1(), region,
                                reference, referenceResource, 0, splice,
                                variantPrinter, new InitialData());

                        CompletableFuture<Scope<AlignedVarsData>> pipeline = pipeline(initialScope, executor);
                        workers.add(pipeline);
                        j++;
                    }

                    List<Map<Integer, Vars>> vars = new ArrayList<>();
                    for (CompletableFuture<Scope<AlignedVarsData>> future : workers) {
                        vars.add(future.join().data.alignedVariants);
                    }

                    Region lastRegion = currentRegion;
                    CompletableFuture<OutputStream> processAmpliconOutput = CompletableFuture
                            .supplyAsync(() -> {
                                OutputStream baos = new ByteArrayOutputStream();
                                PrintStream out = new PrintStream(baos);
                                VariantPrinter variantPrinter = VariantPrinter.createPrinter(instance().printerTypeOut);
                                variantPrinter.setOut(out);
                                new AmpliconPostProcessModule().process(lastRegion, vars, pos, splice, variantPrinter);
                                out.close();
                                return baos;
                                }, executor)
                            .exceptionally(ex -> {
                                stopVardictWithException(lastRegion, ex);
                                throw new RuntimeException(ex);
                            });
                    toPrint.add(processAmpliconOutput);
                }
                toPrint.put(AbstractMode.LAST_SIGNAL_FUTURE);
            }
        };
    }

    @Override
    public void printHeader() {
        if (instance().conf.printHeader) {
            String header = join("\t","Sample", "Gene", "Chr", "Start", "End", "Ref", "Alt", "Depth", "AltDepth", "RefFwdReads",
                    "RefRevReads", "AltFwdReads", "AltRevReads", "Genotype", "AF", "Bias", "PMean", "PStd",
                    "QMean", "QStd", "MQ", "Sig_Noise", "HiAF", "ExtraAF", "shift3", "MSI", "MSI_NT", "NM",
                    "HiCnt", "HiCov", "5pFlankSeq", "3pFlankSeq", "Seg", "VarType", "GoodVarCount", "TotalVarCount",
                    "Nocov", "Ampflag");
            System.out.println(header);
        }
    }
}
