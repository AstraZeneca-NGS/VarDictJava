package com.astrazeneca.utils;

import java.util.concurrent.Executor;

public class DirectThreadExecutor implements Executor {
    @Override
    public void execute(Runnable command) {
        command.run();
    }
}
