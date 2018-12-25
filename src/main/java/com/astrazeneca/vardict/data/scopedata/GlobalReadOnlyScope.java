package com.astrazeneca.vardict.data.scopedata;

import com.astrazeneca.vardict.Configuration;
import com.astrazeneca.vardict.printers.PrinterType;

import java.util.Map;

/**
 * Global scope of the VarDict. Contains configuration that must be available from all the classes and methods.
 * Must be initialized only once. Clear method created only for testing purposes.
 */
public class GlobalReadOnlyScope {

    private volatile static GlobalReadOnlyScope instance;

    public static GlobalReadOnlyScope instance() {
        return instance;
    }

    public static synchronized void init(Configuration conf, Map<String, Integer> chrLengths, String sample, String samplem,
                                         String ampliconBasedCalling, Map<String, Integer> adaptorForward, Map<String, Integer> adaptorReverse) {
        if (instance != null) {
            throw new IllegalStateException("GlobalReadOnlyScope was already initialized. Must be initialized only once.");
        }
        instance = new GlobalReadOnlyScope(conf, chrLengths, sample, samplem, ampliconBasedCalling, adaptorForward, adaptorReverse);
    }

    /**
     * TEST usage only
     */
    public static synchronized void clear(){
        instance = null;
    }

    public final Configuration conf;
    public final Map<String, Integer> chrLengths;
    public final String sample;
    public final String samplem;
    public final String ampliconBasedCalling;
    public final PrinterType printerTypeOut;
    public final Map<String, Integer> adaptorForward;
    public final Map<String, Integer> adaptorReverse;

    public GlobalReadOnlyScope(Configuration conf, Map<String, Integer> chrLengths, String sample, String samplem,
                               String ampliconBasedCalling, Map<String, Integer> adaptorForward, Map<String, Integer> adaptorReverse) {
        this.conf = conf;
        this.chrLengths = chrLengths;
        this.sample = sample;
        this.samplem = samplem;
        this.ampliconBasedCalling = ampliconBasedCalling;
        this.printerTypeOut = conf.printerType;
        this.adaptorForward = adaptorForward;
        this.adaptorReverse = adaptorReverse;
    }
}
