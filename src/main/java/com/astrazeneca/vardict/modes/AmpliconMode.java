package com.astrazeneca.vardict.modes;

import com.astrazeneca.utils.ConcurrentHashSet;
import com.astrazeneca.utils.DirectThreadExecutor;
import com.astrazeneca.utils.Tuple;
import com.astrazeneca.vardict.ReferenceResource;
import com.astrazeneca.vardict.Region;
import com.astrazeneca.vardict.pipeline.data.Scope;
import com.astrazeneca.vardict.pipeline.postprocessmodules.AmpliconPostProcessModule;
import com.astrazeneca.vardict.variations.Vars;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.util.*;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutionException;

import static com.astrazeneca.GlobalReadOnlyScope.instance;
import static com.astrazeneca.utils.Tuple.tuple;

public class AmpliconMode extends AbstractMode {

    public AmpliconMode(List<List<Region>> segments) {
        super(segments);
    }

    @Override
    public void notParallel() {
        for (List<Region> regions : segments) {
            Map<Integer, List<Tuple.Tuple2<Integer, Region>>> pos = new HashMap<>();
            int j = 0;
            Region rg = null;
            final Set<String> splice = new HashSet<>();
            List<Map<Integer, Vars>> vars = new ArrayList<>();
            for (Region region : regions) {
                rg = region; // ??
                for (int p = region.variantsStart; p <= region.variantsEnd; p++) {
                    List<Tuple.Tuple2<Integer, Region>> list = pos.computeIfAbsent(p, k -> new ArrayList<>());
                    list.add(tuple(j, region));
                }

                Tuple.Tuple2<Integer, Map<Integer, Vars>> data = pipeline(
                        instance().conf.bam.getBam1(),
                        region,
                        ReferenceResource.getReference(region),
                        0,
                        null,
                        instance().out,
                        new DirectThreadExecutor()
                ).join().data;

                vars.add(data._2);
                j++;
            }

            AmpliconPostProcessModule.process(rg, vars, pos, splice, instance().out);
        }
    }

    @Override
    protected AbstractParallelMode createParallelMode() {
        return new AbstractParallelMode() {
            @Override
            void produceTasks() throws InterruptedException, ExecutionException {
                for (List<Region> regions : segments) {
                    Map<Integer, List<Tuple.Tuple2<Integer, Region>>> pos = new HashMap<>();
                    int j = 0;
                    Region rg = null;
                    List<CompletableFuture<Scope<Tuple.Tuple2<Integer, Map<Integer, Vars>>>>> workers = new ArrayList<>(regions.size() - 1);
                    final Set<String> splice = new ConcurrentHashSet<>();
                    for (Region region : regions) {
                        rg = region; // ??
                        for (int p = region.variantsStart; p <= region.variantsEnd; p++) {
                            List<Tuple.Tuple2<Integer, Region>> list = pos.computeIfAbsent(p, k -> new ArrayList<>());
                            list.add(tuple(j, region));
                        }

                        workers.add(pipeline(instance().conf.bam.getBam1(),
                                region, null, 0, splice, null, executor));
                        j++;

                    }

                    List<Map<Integer, Vars>> vars = new ArrayList<>();
                    for (CompletableFuture<Scope<Tuple.Tuple2<Integer, Map<Integer, Vars>>>> future : workers) {
                        vars.add(future.get().data._2);
                    }

                    Region lastRegion = rg;
                    toPrint.add(CompletableFuture.supplyAsync(() -> {
                                ByteArrayOutputStream baos = new ByteArrayOutputStream();
                                PrintStream out = new PrintStream(baos);
                                AmpliconPostProcessModule.process(lastRegion, vars, pos, splice, out);
                                out.close();
                                return baos;
                            }, executor)
                    );
                }
                toPrint.put(AbstractMode.LAST_SIGNAL_FUTURE);
            }
        };
    }
}
