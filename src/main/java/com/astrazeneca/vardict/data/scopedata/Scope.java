package com.astrazeneca.vardict.data.scopedata;

import com.astrazeneca.vardict.data.ReferenceResource;
import com.astrazeneca.vardict.data.Reference;
import com.astrazeneca.vardict.data.Region;
import com.astrazeneca.vardict.printers.VariantPrinter;

import java.util.Set;

/**
 * Common scope of data must be storing between steps of VarDict pipeline.
 * @param <T> data of current step of pipeline
 */
public class Scope<T> {

    public final String bam;
    public final Region region;
    public final Reference regionRef;
    public final ReferenceResource referenceResource;
    public final int maxReadLength;
    public final Set<String> splice;
    public final VariantPrinter out;

    public final T data;

    public Scope(String bam, Region region, Reference regionRef, ReferenceResource referenceResource, int maxReadLength,
                 Set<String> splice, VariantPrinter out, T data) {
        this.bam = bam;
        this.region = region;
        this.regionRef = regionRef;
        this.referenceResource = referenceResource;
        this.maxReadLength = maxReadLength;
        this.splice = splice;
        this.data = data;
        this.out = out;
    }

    public Scope(Scope<?> inheritableScope, T data) {
        this(inheritableScope.bam, inheritableScope.region, inheritableScope.regionRef, inheritableScope.referenceResource,
                inheritableScope.maxReadLength, inheritableScope.splice, inheritableScope.out, data);
    }
}
