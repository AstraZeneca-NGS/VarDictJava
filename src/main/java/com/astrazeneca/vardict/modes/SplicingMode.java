package com.astrazeneca.vardict.modes;

import com.astrazeneca.vardict.collection.DirectThreadExecutor;
import com.astrazeneca.vardict.data.Reference;
import com.astrazeneca.vardict.data.ReferenceResource;
import com.astrazeneca.vardict.data.Region;
import com.astrazeneca.vardict.data.scopedata.InitialData;
import com.astrazeneca.vardict.data.scopedata.Scope;
import com.astrazeneca.vardict.data.scopedata.VariationData;
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
 * Mode starting only splicing output without running all the Vardict pipeline: BAM file is processed in sequence by regions,
 * then printed only information about splice.
 */
public class SplicingMode extends AbstractMode {
    public SplicingMode(List<List<Region>> segments, ReferenceResource referenceResource) {
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
                processRegion(region, variantPrinter);
            }
        }
    }

    /**
     * For each segment and region starts the pipeline only for the part of outputting splice counts.
     * @param region region from BED file/-R option to process.
     * @param out variant printer used for output
     */
    private void processRegion(Region region, VariantPrinter out) {
        Reference reference = tryToGetReference(region);
        Scope<InitialData> initialScope = new Scope<>(instance().conf.bam.getBam1(), region,
                reference, referenceResource, 0, new HashSet<>(),
                out, new InitialData());

        CompletableFuture<Scope<VariationData>> pipeline = splicingPipeline(initialScope, new DirectThreadExecutor());
        pipeline.join();
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
            processRegion(region, variantPrinter);
            out.close();
            return baos;
        }
    }

    @Override
    public void printHeader() {
        if (instance().conf.printHeader) {
            String header = join("\t","Sample", "Chr", "Intron", "Intron count");
            System.out.println(header);
        }
    }
}
