package com.astrazeneca.vardict.modes;

import com.astrazeneca.utils.DirectThreadExecutor;
import com.astrazeneca.vardict.ReferenceResource;
import com.astrazeneca.vardict.Region;
import com.astrazeneca.vardict.pipeline.postprocessmodules.SimplePostProcessModule;

import java.io.ByteArrayOutputStream;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.concurrent.*;

import static com.astrazeneca.GlobalReadOnlyScope.instance;

public class SimpleMode extends AbstractMode {

    public SimpleMode(List<List<Region>> segments) {
        super(segments);
    }

    @Override
    public void notParallel() {
        for (List<Region> list : segments) {
            for (Region region : list) {
                simple(region, instance().out);
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
                        toPrint.put(executor.submit(new SimpleMode.VardictWorker(region)));
                    }
                }
                toPrint.put(LAST_SIGNAL_FUTURE);
            }
        };
    }


    private static void simple(Region region, PrintStream out) {
        Map<Integer, Character> ref = ReferenceResource.getReference(region);
        AbstractMode.pipeline(
                instance().conf.bam.getBam1(),
                region,
                ref,
                0,
                new HashSet<>(),
                out,
                new DirectThreadExecutor()
        )
                .thenAccept(new SimplePostProcessModule())
                .join();
    }

    private static class VardictWorker implements Callable<OutputStream> {

        private Region region;

        public VardictWorker(Region region) {
            super();
            this.region = region;
        }

        @Override
        public OutputStream call() throws Exception {
            ByteArrayOutputStream baos = new ByteArrayOutputStream();
            PrintStream out = new PrintStream(baos);

            SimpleMode.simple(region, out);

            out.close();
            return baos;
        }

    }
}
