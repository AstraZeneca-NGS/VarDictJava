package com.astrazeneca.vardict.modules;

import com.astrazeneca.vardict.data.scopedata.Scope;

@FunctionalInterface
public interface Module<T, R> {

    Scope<R> process(Scope<T> scope);
}
