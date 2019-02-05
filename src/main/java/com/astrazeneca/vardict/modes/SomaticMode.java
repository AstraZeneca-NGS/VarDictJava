package com.astrazeneca.vardict.modes;

import com.astrazeneca.vardict.data.ReferenceResource;
import com.astrazeneca.vardict.collection.ConcurrentHashSet;
import com.astrazeneca.vardict.collection.DirectThreadExecutor;
import com.astrazeneca.vardict.collection.Tuple;
import com.astrazeneca.vardict.data.Reference;
import com.astrazeneca.vardict.data.Region;
import com.astrazeneca.vardict.data.scopedata.Scope;
import com.astrazeneca.vardict.postprocessmodules.SomaticPostProcessModule;
import com.astrazeneca.vardict.printers.VariantPrinter;
import com.astrazeneca.vardict.variations.Vars;

import java.io.ByteArrayOutputStream;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;

import static com.astrazeneca.vardict.Utils.join;
import static com.astrazeneca.vardict.data.scopedata.GlobalReadOnlyScope.instance;

/**
 * Mode starting Somatic analysis: tumor and normal BAM files are processed in sequence by regions, then received Vars processed in
 * SomaticPostProcessModule..
 */
public class SomaticMode extends AbstractMode {

    public SomaticMode(List<List<Region>> segments, ReferenceResource referenceResource) {
        super(segments, referenceResource);
        printHeader();
    }

    @Override
    public void notParallel() {
        VariantPrinter variantPrinter = VariantPrinter.createPrinter(instance().printerTypeOut);

        for (List<Region> list : segments) {
            for (Region region : list) {

                final Set<String> splice = new ConcurrentHashSet<>();
                Reference ref = referenceResource.getReference(region);
                CompletableFuture<Scope<Tuple.Tuple2<Integer, Map<Integer, Vars>>>> bam1VariationsFuture = pipeline(instance().conf.bam.getBam1(),
                        region, ref, referenceResource, 0, splice, variantPrinter, new DirectThreadExecutor());

                Scope<Tuple.Tuple2<Integer, Map<Integer, Vars>>> bam1Variations = bam1VariationsFuture.join();

                CompletableFuture<Scope<Tuple.Tuple2<Integer, Map<Integer, Vars>>>> bam2VariationFuture = pipeline(instance().conf.bam.getBam2(),
                        region, ref, referenceResource, bam1Variations.maxReadLength, splice, variantPrinter, new DirectThreadExecutor());

                CompletableFuture<Void> somaticProcessOutput = bam2VariationFuture.thenAcceptBoth(bam1VariationsFuture, new SomaticPostProcessModule(referenceResource, variantPrinter));

                somaticProcessOutput.join();
            }
        }
    }

    @Override
    protected AbstractParallelMode createParallelMode() {
        return new AbstractParallelMode() {
            @Override
            void produceTasks() throws InterruptedException, ExecutionException {
                for (List<Region> list : segments) {
                    for (Region region : list) {
                        final Set<String> splice = new ConcurrentHashSet<>();
                        Reference ref1 = referenceResource.getReference(region);
                        Future<OutputStream> f2 = executor.submit(new SomdictWorker(region, splice, ref1));
                        toPrint.put(f2);
                    }
                }
                toPrint.put(AbstractMode.LAST_SIGNAL_FUTURE);
            }
        };
    }

    private class SomdictWorker implements Callable<OutputStream> {
        private final Region region;
        private final Set<String> splice;
        final Reference ref;

        public SomdictWorker(Region region, Set<String> splice, Reference ref) {
            this.region = region;
            this.splice = splice;
            this.ref = ref;
        }

        @Override
        public OutputStream call() {
            ByteArrayOutputStream baos = new ByteArrayOutputStream();
            PrintStream out = new PrintStream(baos);
            VariantPrinter variantPrinter = VariantPrinter.createPrinter(instance().printerTypeOut);
            variantPrinter.setOut(out);

            CompletableFuture<Scope<Tuple.Tuple2<Integer, Map<Integer, Vars>>>> bam1VariationFuture =
                    pipeline(instance().conf.bam.getBam1(),
                            region, ref, referenceResource, 0, splice, variantPrinter, new DirectThreadExecutor());
            CompletableFuture<Scope<Tuple.Tuple2<Integer, Map<Integer, Vars>>>> bam2VariationFuture =
                    pipeline(instance().conf.bam.getBam2(),
                            region, ref, referenceResource, 0, splice, variantPrinter, new DirectThreadExecutor());

            CompletableFuture<Void> somaticProcessOutput = bam2VariationFuture.thenAcceptBoth(bam1VariationFuture,
                    new SomaticPostProcessModule(referenceResource, variantPrinter));
            somaticProcessOutput.join();
            out.close();

            return baos;
        }
    }

    @Override
    public void printHeader() {
        if (instance().conf.printHeader) {
            System.out.println(join("\t",
                    "Sample", "Gene", "Chr", "Start", "End", "Ref", "Alt",
                    "Depth", "AltDepth", "RefFwdReads", "RefRevReads", "AltFwdReads", "AltRevReads", "Genotype", "AF",
                    "Bias", "PMean", "PStd", "QMean", "QStd", "MQ", "Sig_Noise", "HiAF", "ExtraAF", "NM",
                    "Depth", "AltDepth", "RefFwdReads", "RefRevReads", "AltFwdReads", "AltRevReads", "Genotype", "AF",
                    "Bias", "PMean", "PStd", "QMean", "QStd", "MQ", "Sig_Noise", "HiAF", "ExtraAF", "NM",
                    "shift3", "MSI", "MSI_NT", "5pFlankSeq", "3pFlankSeq", "Seg", "VarLabel", "VarType",
                    "Duprate1", "SV_info1", "Duprate2", "SV_info2"));
        }
    }
}
