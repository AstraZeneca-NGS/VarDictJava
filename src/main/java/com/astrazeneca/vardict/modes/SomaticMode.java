package com.astrazeneca.vardict.modes;

import com.astrazeneca.utils.ConcurrentHashSet;
import com.astrazeneca.utils.DirectThreadExecutor;
import com.astrazeneca.utils.Tuple;
import com.astrazeneca.vardict.ReferenceResource;
import com.astrazeneca.vardict.Region;
import com.astrazeneca.vardict.pipeline.data.Scope;
import com.astrazeneca.vardict.pipeline.postprocessmodules.SomaticPostProcessModule;
import com.astrazeneca.vardict.variations.Vars;

import java.io.ByteArrayOutputStream;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.*;

import static com.astrazeneca.GlobalReadOnlyScope.instance;

public class SomaticMode extends AbstractMode {

    public SomaticMode(List<List<Region>> segments) {
        super(segments);
    }

    @Override
    public void notParallel() {
        for (List<Region> list : segments) {
            for (Region region : list) {

                final Set<String> splice = new ConcurrentHashSet<>();
                Map<Integer, Character> ref = ReferenceResource.getReference(region);

                CompletableFuture<Scope<Tuple.Tuple2<Integer, Map<Integer, Vars>>>> bam1VariationsFuture = AbstractMode.pipeline(instance().conf.bam.getBam1(),
                        region, ref, 0, splice, instance().out, new DirectThreadExecutor());

                Scope<Tuple.Tuple2<Integer, Map<Integer, Vars>>> bam1Variations = bam1VariationsFuture.join();

                CompletableFuture<Scope<Tuple.Tuple2<Integer, Map<Integer, Vars>>>> variations = AbstractMode.pipeline(instance().conf.bam.getBam2(),
                        region, ref, bam1Variations.maxReadLength, splice, instance().out, new DirectThreadExecutor());

                variations.thenAcceptBoth(bam1VariationsFuture, new SomaticPostProcessModule()).join();
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
                        Map<Integer, Character> ref1 = ReferenceResource.getReference(region);
                        Future<OutputStream> f2 = executor.submit(new SomdictWorker(region, splice, ref1));
                        toPrint.put(f2);
                    }
                }
                toPrint.put(AbstractMode.LAST_SIGNAL_FUTURE);
            }
        };
    }

    private static class SomdictWorker implements Callable<OutputStream> {
        private final Region region;
        private final Set<String> splice;

        final Map<Integer, Character> ref;

        public SomdictWorker(Region region, Set<String> splice, Map<Integer, Character> ref) {
            this.region = region;
            this.splice = splice;
            this.ref = ref;
        }

        @Override
        public OutputStream call() throws Exception {
            ByteArrayOutputStream baos = new ByteArrayOutputStream();
            PrintStream out = new PrintStream(baos);

            CompletableFuture<Scope<Tuple.Tuple2<Integer, Map<Integer, Vars>>>> bam1VariationFuture = AbstractMode.pipeline(instance().conf.bam.getBam1(),
                    region, ref, 0, splice, out, new DirectThreadExecutor());

            AbstractMode.pipeline(instance().conf.bam.getBam2(), region, ref, 0, splice, out, new DirectThreadExecutor())
                    .thenAcceptBoth(bam1VariationFuture, new SomaticPostProcessModule())
                    .join();

            out.close();
            return baos;
        }

    }
}
