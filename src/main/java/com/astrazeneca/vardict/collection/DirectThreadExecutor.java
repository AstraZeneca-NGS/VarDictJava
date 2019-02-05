package com.astrazeneca.vardict.collection;

import java.util.concurrent.Executor;

public class DirectThreadExecutor implements Executor {
    @Override
    public void execute(Runnable command) {
        command.run();
    }
}
