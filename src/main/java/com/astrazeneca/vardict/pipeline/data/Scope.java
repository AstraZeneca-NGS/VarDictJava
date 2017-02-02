package com.astrazeneca.vardict.pipeline.data;

import com.astrazeneca.vardict.Region;

import java.io.PrintStream;
import java.util.Map;
import java.util.Set;

public class Scope<T> {

    public final String bam;
    public final Region region;
    public final Map regionRef;
    public final int maxReadLength;
    public final Set<String> splice;
    public final PrintStream out;

    public final T data;

    public Scope(String bam, Region region, Map regionRef, int maxReadLength, Set<String> splice, PrintStream out, T data) {
        this.bam = bam;
        this.region = region;
        this.regionRef = regionRef;
        this.maxReadLength = maxReadLength;
        this.splice = splice;
        this.out = out;
        this.data = data;
    }

    public Scope(Scope inheritableScope, T data) {
        this(inheritableScope.bam, inheritableScope.region, inheritableScope.regionRef, inheritableScope.maxReadLength, inheritableScope.splice, inheritableScope.out, data);
    }
}
