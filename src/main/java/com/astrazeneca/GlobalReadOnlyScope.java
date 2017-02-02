package com.astrazeneca;

import com.astrazeneca.vardict.Configuration;

import java.io.PrintStream;
import java.util.Map;

public class GlobalReadOnlyScope {

    private volatile static GlobalReadOnlyScope instance;

    public static GlobalReadOnlyScope instance() {
        return instance;
    }

    public static synchronized void init(Configuration conf, Map<String, Integer> chrLengths, String sample, String samplem, String ampliconBasedCalling) {
        if (instance != null) {
            throw new IllegalStateException("GlobalReadOnlyScope was already initialized. Must be initialized only once.");
        }
        instance = new GlobalReadOnlyScope(conf, chrLengths, sample, samplem, ampliconBasedCalling);
    }

    public final Configuration conf;
    public final Map<String, Integer> chrLengths;
    public final String sample;
    public final String samplem;
    public final String ampliconBasedCalling;
    public final PrintStream out;

    public GlobalReadOnlyScope(Configuration conf, Map<String, Integer> chrLengths, String sample, String samplem, String ampliconBasedCalling) {
        this.conf = conf;
        this.chrLengths = chrLengths;
        this.sample = sample;
        this.samplem = samplem;
        this.ampliconBasedCalling = ampliconBasedCalling;
        this.out = System.out;
    }
}
