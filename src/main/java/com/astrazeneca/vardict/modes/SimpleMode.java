package com.astrazeneca.vardict.modes;

import com.astrazeneca.vardict.collection.Tuple;
import com.astrazeneca.vardict.data.ReferenceResource;
import com.astrazeneca.vardict.collection.DirectThreadExecutor;
import com.astrazeneca.vardict.data.Reference;
import com.astrazeneca.vardict.data.Region;
import com.astrazeneca.vardict.data.scopedata.Scope;
import com.astrazeneca.vardict.postprocessmodules.SimplePostProcessModule;
import com.astrazeneca.vardict.printers.VariantPrinter;
import com.astrazeneca.vardict.variations.Vars;

import java.io.ByteArrayOutputStream;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
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

    @Override
    public void notParallel() {
        VariantPrinter variantPrinter = VariantPrinter.createPrinter(instance().printerTypeOut);

        for (List<Region> list : segments) {
            for (Region region : list) {
                simple(region, variantPrinter);
            }
        }
    }

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

    private void simple(Region region, VariantPrinter out) {
        Reference ref = referenceResource.getReference(region);
        CompletableFuture<Scope<Tuple.Tuple2<Integer, Map<Integer, Vars>>>> pipeline = pipeline(
                instance().conf.bam.getBam1(),
                region,
                ref,
                referenceResource,
                0,
                new HashSet<>(),
                out,
                new DirectThreadExecutor()
        );
        CompletableFuture<Void> simpleProcessOutput = pipeline.thenAccept(new SimplePostProcessModule(out));
        simpleProcessOutput.join();

    }

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
            simple(region, variantPrinter);

            out.close();
            return baos;
        }

    }

    @Override
    public void printHeader() {
        if (instance().conf.printHeader) {
            System.out.println(join("\t",
                    "Sample", "Gene", "Chr", "Start", "End", "Ref", "Alt", "Depth", "AltDepth", "RefFwdReads",
                    "RefRevReads", "AltFwdReads", "AltRevReads", "Genotype", "AF", "Bias", "PMean", "PStd",
                    "QMean", "QStd", "MQ", "Sig_Noise", "HiAF", "ExtraAF", "shift3", "MSI", "MSI_NT", "NM",
                    "HiCnt", "HiCov", "5pFlankSeq", "3pFlankSeq", "Seg", "VarType", "Duprate", "SV_info"));
        }
    }

}
