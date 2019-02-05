package com.astrazeneca.vardict.modules;

import com.astrazeneca.vardict.data.scopedata.Scope;

/**
 * Functional interface for all Modules of Vardict (they can be the steps of pipeline in AbstractMode).
 * @param <T> means input data needed on step (module)
 * @param <R> means output data that step (module) produces
 */
@FunctionalInterface
public interface Module<T, R> {

    Scope<R> process(Scope<T> scope);
}
