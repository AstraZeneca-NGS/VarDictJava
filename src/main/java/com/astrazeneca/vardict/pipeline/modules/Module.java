package com.astrazeneca.vardict.pipeline.modules;

import com.astrazeneca.vardict.pipeline.data.Scope;

@FunctionalInterface
public interface Module<T, R> {

    Scope<R> process(Scope<T> scope);
}
