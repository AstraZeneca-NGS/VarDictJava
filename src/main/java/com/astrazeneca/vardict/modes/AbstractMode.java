package com.astrazeneca.vardict.modes;

import com.astrazeneca.vardict.data.ReferenceResource;
import com.astrazeneca.vardict.collection.Tuple;
import com.astrazeneca.vardict.data.*;
import com.astrazeneca.vardict.data.scopedata.InitialData;
import com.astrazeneca.vardict.data.scopedata.Scope;
import com.astrazeneca.vardict.data.scopedata.VariationData;
import com.astrazeneca.vardict.modules.*;
import com.astrazeneca.vardict.printers.VariantPrinter;
import com.astrazeneca.vardict.variations.Vars;

import java.io.OutputStream;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.*;

import static com.astrazeneca.vardict.data.scopedata.GlobalReadOnlyScope.instance;

/**
 * Abstract Mode of VarDict. Provide interfaces for possible modes (not-parallel and parallel) and typical pipelines.
 */
public abstract class AbstractMode {

    final static CompletableFuture<OutputStream> LAST_SIGNAL_FUTURE = CompletableFuture.completedFuture(null);

    protected List<List<Region>> segments;
    ReferenceResource referenceResource;

    public AbstractMode(List<List<Region>> segments, ReferenceResource referenceResource) {
        this.segments = segments;
        this.referenceResource = referenceResource;
    }

    public static CompletableFuture<Scope<Tuple.Tuple2<Integer, Map<Integer, Vars>>>> pipeline(String bam,
                                                                                               Region region,
                                                                                               Reference ref,
                                                                                               ReferenceResource referenceResource,
                                                                                               int maxReadLength,
                                                                                               Set<String> splice,
                                                                                               VariantPrinter out,
                                                                                               Executor executor) {

        return CompletableFuture.supplyAsync(
                () -> new SAMFileParser().process(
                        createInitialPipelineScope(bam, region, ref, referenceResource, maxReadLength, splice, out, new InitialData())
                ), executor)
                .thenApply(new CigarParser(false)::process)
                .thenApply(new VariationRealigner()::process)
                .thenApply(new StructuralVariantsProcessor()::process)
                .thenApply(new ToVarsBuilder()::process);
    }

    public static CompletableFuture<Scope<VariationData>> partialPipeline(InitialData initialData,
                                                                          String bam,
                                                                          Region region,
                                                                          Reference ref,
                                                                          ReferenceResource referenceResource,
                                                                          int maxReadLength,
                                                                          Set<String> splice,
                                                                          VariantPrinter out,
                                                                          Executor executor) {


        return CompletableFuture.supplyAsync(
                () -> new SAMFileParser().process(
                        createInitialPipelineScope(bam, region, ref, referenceResource, maxReadLength, splice, out, initialData)
                ), executor)
                .thenApply(new CigarParser(true)::process);
    }

    private static Scope<InitialData> createInitialPipelineScope(String bam,
                                                                 Region region,
                                                                 Reference ref,
                                                                 ReferenceResource referenceResource,
                                                                 int maxReadLength,
                                                                 Set<String> splice,
                                                                 VariantPrinter out,
                                                                 InitialData initialData) {
        Scope<InitialData> scope;

        scope = new Scope<>(bam, region, (ref == null ? referenceResource.getReference(region) : ref), referenceResource,
                maxReadLength, splice, out, initialData);
        return scope;
    }

    public abstract void notParallel();

    public void parallel() {
        createParallelMode().process();
    }

    protected abstract AbstractParallelMode createParallelMode();

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

    public abstract void printHeader();

}
