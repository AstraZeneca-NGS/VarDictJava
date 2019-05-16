package com.astrazeneca.vardict.modes;

import com.astrazeneca.vardict.data.ReferenceResource;
import com.astrazeneca.vardict.collection.ConcurrentHashSet;
import com.astrazeneca.vardict.collection.DirectThreadExecutor;
import com.astrazeneca.vardict.data.Reference;
import com.astrazeneca.vardict.data.Region;
import com.astrazeneca.vardict.data.scopedata.AlignedVarsData;
import com.astrazeneca.vardict.data.scopedata.InitialData;
import com.astrazeneca.vardict.data.scopedata.Scope;
import com.astrazeneca.vardict.postprocessmodules.SomaticPostProcessModule;
import com.astrazeneca.vardict.printers.VariantPrinter;

import java.io.ByteArrayOutputStream;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.List;
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

    /**
     * In not parallel mode each region will be processed in sequence.
     */
    @Override
    public void notParallel() {
        VariantPrinter variantPrinter = VariantPrinter.createPrinter(instance().printerTypeOut);

        for (List<Region> list : segments) {
            for (Region region : list) {
                final Set<String> splice = new ConcurrentHashSet<>();
                Reference ref = tryToGetReference(region);
                processBothBamsInPipeline(variantPrinter, region, splice, ref);
            }
        }
    }

    /**
     * In parallel mode somatic workers are created for each region and are processed in parallel.
     */
    @Override
    protected AbstractParallelMode createParallelMode() {
        return new AbstractParallelMode() {
            @Override
            void produceTasks() throws InterruptedException, ExecutionException {
                for (List<Region> list : segments) {
                    for (Region region : list) {
                        final Set<String> splice = new ConcurrentHashSet<>();
                        Reference ref1 = tryToGetReference(region);
                        Future<OutputStream> f2 = executor.submit(new SomdictWorker(region, splice, ref1));
                        toPrint.put(f2);
                    }
                }
                toPrint.put(AbstractMode.LAST_SIGNAL_FUTURE);
            }
        };
    }

    /**
     * Class needed for somatic parallel mode. Each worker will process pipeline for region on both BAM files
     * (tumor and normal).
     */
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
            processBothBamsInPipeline(variantPrinter, region, splice, ref);
            out.close();

            return baos;
        }
    }

    private void processBothBamsInPipeline(VariantPrinter variantPrinter, Region region, Set<String> splice, Reference ref) {
        Scope<InitialData> initialScope1 = new Scope<>(
                instance().conf.bam.getBam1(),
                region, ref, referenceResource,
                0,
                splice,
                variantPrinter,
                new InitialData());
        CompletableFuture<Scope<AlignedVarsData>> bam1VariationsFuture = pipeline(initialScope1, new DirectThreadExecutor());

        Scope<AlignedVarsData> bam1Variations = bam1VariationsFuture.join();

        Scope<InitialData> initialScope2 = new Scope<>(
                instance().conf.bam.getBam2(),
                region, ref, referenceResource,
                bam1Variations.maxReadLength,
                splice,
                variantPrinter,
                new InitialData());
        CompletableFuture<Scope<AlignedVarsData>> bam2VariationFuture = pipeline(initialScope2, new DirectThreadExecutor());

        CompletableFuture<Void> somaticProcessOutput = bam2VariationFuture
                .thenAcceptBoth(bam1VariationsFuture, new SomaticPostProcessModule(referenceResource, variantPrinter))
                .exceptionally(ex -> {
                    stopVardictWithException(region, ex);
                    throw new RuntimeException(ex);
                });
        somaticProcessOutput.join();
    }

    @Override
    public void printHeader() {
        if (instance().conf.printHeader) {
            String header = join("\t",
                    "Sample", "Gene", "Chr", "Start", "End", "Ref", "Alt",
                    "Depth", "AltDepth", "RefFwdReads", "RefRevReads", "AltFwdReads", "AltRevReads", "Genotype", "AF",
                    "Bias", "PMean", "PStd", "QMean", "QStd", "MQ", "Sig_Noise", "HiAF", "ExtraAF", "NM",
                    "Depth", "AltDepth", "RefFwdReads", "RefRevReads", "AltFwdReads", "AltRevReads", "Genotype", "AF",
                    "Bias", "PMean", "PStd", "QMean", "QStd", "MQ", "Sig_Noise", "HiAF", "ExtraAF", "NM",
                    "shift3", "MSI", "MSI_NT", "5pFlankSeq", "3pFlankSeq", "Seg", "VarLabel", "VarType",
                    "Duprate1", "SV_info1", "Duprate2", "SV_info2");
            System.out.println(header);
        }
    }
}
