package com.astrazeneca.vardict.modes;

import com.astrazeneca.utils.Tuple;
import com.astrazeneca.vardict.ReferenceResource;
import com.astrazeneca.vardict.Region;
import com.astrazeneca.vardict.pipeline.data.Scope;
import com.astrazeneca.vardict.pipeline.modules.CigarParser;
import com.astrazeneca.vardict.pipeline.modules.SAMFileParser;
import com.astrazeneca.vardict.pipeline.modules.ToVarsBuilder;
import com.astrazeneca.vardict.pipeline.modules.VariationRealigner;
import com.astrazeneca.vardict.variations.Vars;

import java.io.OutputStream;
import java.io.PrintStream;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.*;

import static com.astrazeneca.GlobalReadOnlyScope.instance;

public abstract class AbstractMode {

    final static CompletableFuture<OutputStream> LAST_SIGNAL_FUTURE = CompletableFuture.completedFuture(null);

    protected List<List<Region>> segments;

    public AbstractMode(List<List<Region>> segments) {
        this.segments = segments;
    }

    public static CompletableFuture<Scope<Tuple.Tuple2<Integer, Map<Integer, Vars>>>> pipeline(String bam,
                                                                                               Region region,
                                                                                               Map<Integer, Character> ref,
                                                                                               int maxReadLength,
                                                                                               Set<String> splice,
                                                                                               PrintStream out,
                                                                                               Executor executor) {


        return CompletableFuture.supplyAsync(
                () -> new SAMFileParser().process(
                        createInitialPipelineScope(bam, region, ref, maxReadLength, splice, out)
                ), executor)
                .thenApply(new CigarParser()::process)
                .thenApply(new VariationRealigner()::process)
                .thenApply(new ToVarsBuilder()::process);
    }

    private static Scope<Object> createInitialPipelineScope(String bam, Region region, Map<Integer, Character> ref, int maxReadLength, Set<String> splice, PrintStream out) {
        Scope<Object> scope;
        scope = new Scope<>(bam, region, ref == null ? ReferenceResource.getReference(region) : ref, maxReadLength, splice, out, null);
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
                    System.out.print(wrk.get());
                }
            } catch (InterruptedException | ExecutionException e) {
                throw new RuntimeException(e);
            }

            executor.shutdown();
        }

        abstract void produceTasks() throws InterruptedException, ExecutionException;
    }
}
