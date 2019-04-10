package com.astrazeneca.vardict.modes;

import com.astrazeneca.vardict.data.ReferenceResource;
import com.astrazeneca.vardict.collection.DirectThreadExecutor;
import com.astrazeneca.vardict.data.Reference;
import com.astrazeneca.vardict.data.Region;
import com.astrazeneca.vardict.data.scopedata.AlignedVarsData;
import com.astrazeneca.vardict.data.scopedata.InitialData;
import com.astrazeneca.vardict.data.scopedata.Scope;
import com.astrazeneca.vardict.postprocessmodules.SimplePostProcessModule;
import com.astrazeneca.vardict.printers.VariantPrinter;

import java.io.ByteArrayOutputStream;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.HashSet;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.Future;

import static com.astrazeneca.vardict.Utils.join;
import static com.astrazeneca.vardict.data.scopedata.GlobalReadOnlyScope.instance;

/**
 * Mode starting Simple analysis: BAM file is processed in sequence by regions, then received Vars processed in
 * SimplePostProcessModule.
 */
public class SimpleMode extends AbstractMode {

    public SimpleMode(List<List<Region>> segments, ReferenceResource referenceResource) {
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
                processBamInPipeline(region, variantPrinter);
            }
        }
    }

    /**
     * In parallel mode workers are created for each region and are processed in parallel.
     */
    @Override
    protected AbstractParallelMode createParallelMode() {
        return new AbstractParallelMode() {
            @Override
            void produceTasks() throws InterruptedException {
                for (List<Region> list : segments) {
                    for (Region region : list) {
                        Future<OutputStream> submit = executor.submit(new VardictWorker(region));
                        toPrint.put(submit);
                    }
                }
                toPrint.put(LAST_SIGNAL_FUTURE);
            }
        };
    }

    /**
     * Class needed for simple parallel mode. Each worker will process pipeline for region.
     */
    private class VardictWorker implements Callable<OutputStream> {
        private Region region;

        public VardictWorker(Region region) {
            super();
            this.region = region;
        }

        @Override
        public OutputStream call() {
            ByteArrayOutputStream baos = new ByteArrayOutputStream();
            PrintStream out = new PrintStream(baos);
            VariantPrinter variantPrinter = VariantPrinter.createPrinter(instance().printerTypeOut);
            variantPrinter.setOut(out);
            processBamInPipeline(region, variantPrinter);
            out.close();
            return baos;
        }
    }

    /**
     * For each segment and region starts the pipeline.
     * @param region region from BED file/-R option to process.
     * @param out variant printer used for output
     */
    private void processBamInPipeline(Region region, VariantPrinter out) {
        Reference reference = tryToGetReference(region);
        Scope<InitialData> initialScope = new Scope<>(instance().conf.bam.getBam1(), region,
                reference, referenceResource, 0, new HashSet<>(),
                out, new InitialData());

        CompletableFuture<Scope<AlignedVarsData>> pipeline = pipeline(initialScope, new DirectThreadExecutor());
        CompletableFuture<Void> simpleProcessOutput = pipeline
                .thenAccept(new SimplePostProcessModule(out))
                .exceptionally(ex -> {
                    stopVardictWithException(region, ex);
                    throw new RuntimeException(ex);
                });
        simpleProcessOutput.join();
    }

    @Override
    public void printHeader() {
        if (instance().conf.printHeader) {
            String header = join("\t",
                    "Sample", "Gene", "Chr", "Start", "End", "Ref", "Alt", "Depth", "AltDepth", "RefFwdReads",
                    "RefRevReads", "AltFwdReads", "AltRevReads", "Genotype", "AF", "Bias", "PMean", "PStd",
                    "QMean", "QStd", "MQ", "Sig_Noise", "HiAF", "ExtraAF", "shift3", "MSI", "MSI_NT", "NM",
                    "HiCnt", "HiCov", "5pFlankSeq", "3pFlankSeq", "Seg", "VarType", "Duprate", "SV_info");
            if (instance().conf.crisprCuttingSite != 0) {
                header = join("\t", header, "CRISPR");
            }
            System.out.println(header);
        }
    }
}
