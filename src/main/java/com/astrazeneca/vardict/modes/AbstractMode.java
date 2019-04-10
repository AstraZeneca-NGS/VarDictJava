package com.astrazeneca.vardict.modes;

import com.astrazeneca.vardict.data.ReferenceResource;
import com.astrazeneca.vardict.data.*;
import com.astrazeneca.vardict.data.scopedata.AlignedVarsData;
import com.astrazeneca.vardict.data.scopedata.InitialData;
import com.astrazeneca.vardict.data.scopedata.Scope;
import com.astrazeneca.vardict.data.scopedata.VariationData;
import com.astrazeneca.vardict.modules.*;
import com.astrazeneca.vardict.printers.VariantPrinter;

import java.io.OutputStream;
import java.util.List;
import java.util.concurrent.*;

import static com.astrazeneca.vardict.data.scopedata.GlobalReadOnlyScope.instance;

/**
 * Abstract Mode of VarDict. Provide interfaces for possible modes (not-parallel and parallel) and typical pipelines.
 */
public abstract class AbstractMode {
    /**
     * CompletableFuture which used to trigger the executor that it was the last Future in the chain. After adding it,
     * the variants will be printed.
     */
    final static CompletableFuture<OutputStream> LAST_SIGNAL_FUTURE = CompletableFuture.completedFuture(null);

    protected List<List<Region>> segments;
    ReferenceResource referenceResource;

    public AbstractMode(List<List<Region>> segments, ReferenceResource referenceResource) {
        this.segments = segments;
        this.referenceResource = referenceResource;
    }

    /**
     * Starts the typical pipeline of VarDict on each region. It parse SAM/BAM file, for each record parse CIGAR,
     * modify it if needed, create Variations, realign variations, search for Structural Variants and create map of
     * aligned variants ready for output preparation.
     * @param initialDataScope initial data for pipeline. contains data about BAM, region, reference
     * @param executor current Executor for parallel/single mode
     * @return object contains map of aligned variants
     */
    public CompletableFuture<Scope<AlignedVarsData>> pipeline(Scope<InitialData> initialDataScope,
                                                                     Executor executor) {
        return CompletableFuture.supplyAsync(
                () -> new SAMFileParser().process(initialDataScope), executor)
                .thenApply(new CigarParser(false)::process)
                .thenApply(new VariationRealigner()::process)
                .thenApply(new StructuralVariantsProcessor()::process)
                .thenApply(new ToVarsBuilder()::process)
                .exceptionally(ex -> {
                    stopVardictWithException(initialDataScope.region, ex);
                    throw new RuntimeException(ex);
                });
    }

    /**
     * Method for stopping VarDict with Exception. Number of exceptions set by Configuration.MAX_EXCEPTION_COUNT can be
     * skipped, but after this limit program must stop even in parallel mode (where Executor may hang on Exception).
     * @param region region where Exception occurs
     * @param ex initial exception
     */
    void stopVardictWithException(Region region, Throwable ex) {
        System.err.println("Critical exception occurs on region: "
                + region.chr +":" + region.start + "-" + region.end + ", program will be stopped.");
        ex.printStackTrace();
        System.exit(1);
    }

    /**
     * Starts the partial pipeline of VarDict on each region needed while searching structural variants on extended regions.
     * It parse SAM/BAM file, for each record parse CIGAR, modify it if needed and create Variations.
     * @param currentScope current data for pipeline. contains data about BAM, region, reference and already
     *                     filled maps of variations and softclips
     * @param executor current Executor for parallel/single mode
     * @return object contains variation data (updated maps of variations and softclips).
     */
    public CompletableFuture<Scope<VariationData>> partialPipeline(Scope<InitialData> currentScope,
                                                                          Executor executor) {
        return CompletableFuture.supplyAsync(
                () -> new SAMFileParser().process(currentScope), executor)
                .thenApply(new CigarParser(true)::process);
    }

    /**
     * Starts pipeline for splicing mode: output only information about splice counts (N in CIGAR reads).
     * @param currentScope current data for pipeline. contains data about BAM, region, reference and already
     *                     filled maps of variations and softclips
     * @param executor current Executor for parallel/single mode
     * @return empty object for variation data.
     */
    public CompletableFuture<Scope<VariationData>> splicingPipeline(Scope<InitialData> currentScope,
                                                                          Executor executor) {
        return CompletableFuture.supplyAsync(
                () -> new SAMFileParser().process(currentScope), executor)
                .thenApply(new CigarParser(false)::process);
    }

    public abstract void notParallel();

    public void parallel() {
        createParallelMode().process();
    }

    protected abstract AbstractParallelMode createParallelMode();

    /**
     * Abstract class for parallel modes of VarDict. Initializes executor, starts tasks and prints the variants.
     * The tasks producer must be overriding in the childs.
     */
    protected static abstract class AbstractParallelMode {
        static final int CAPACITY = 10;

        final ExecutorService executor = Executors.newFixedThreadPool(instance().conf.threads);
        final BlockingQueue<Future<OutputStream>> toPrint = new LinkedBlockingQueue<>(CAPACITY);

        void process() {
            executor.submit(() -> {
                try {
                    produceTasks();
                } catch (InterruptedException | ExecutionException e) {
                    throw new RuntimeException(e);
                }
            });
            try {
                while (true) {
                    Future<OutputStream> wrk = toPrint.take();
                    if (wrk == LAST_SIGNAL_FUTURE) {
                        break;
                    }
                    VariantPrinter variantPrinter = VariantPrinter.createPrinter(instance().printerTypeOut);
                    variantPrinter.print(wrk.get());
                }
            } catch (InterruptedException | ExecutionException e) {
                throw new RuntimeException(e);
            }

            executor.shutdown();
        }

        abstract void produceTasks() throws InterruptedException, ExecutionException;
    }

    /**
     * Print header to output with option -h. Each mode creates own string for header.
     */
    public abstract void printHeader();

    public Reference tryToGetReference(Region region) {
        Reference reference = new Reference();
        try {
            reference = referenceResource.getReference(region);
        } catch (Exception ex) {
            stopVardictWithException(region, ex);
        }
        return reference;
    }
}
