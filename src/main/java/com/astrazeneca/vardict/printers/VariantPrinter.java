package com.astrazeneca.vardict.printers;

import java.io.OutputStream;
import java.io.PrintStream;

/**
 * Universal class for print Variants and OutputStreams. The "out" stream can be set by implementing
 * new type in enum PrinterType, adding it to method createPrinter() and creating new child VariantPrinter class
 * extending this class. The choice of variant printer can be made by parameter --defaultPrinter (PrinterType.OUT by default)
 */
public abstract class VariantPrinter {
    protected PrintStream out;
    private String delimiter = "\t";

    /**
     * Prints special variant output structure to the set output
     * @param variant output variant structure
     */
    public void print(OutputVariant variant) {
        out.println(variant.outputString(delimiter));
    }

    /**
     * Prints output stream to the set output. Usually used to print variants in multi-threading mode after they
     * were prepared from the queue.
     * @param outputStream output stream to print variants
     */
    public void print(OutputStream outputStream) {
        out.print(outputStream);
    }

    /**
     * Set out stream to the parameter. Usually used to update output in multi-threading mode to prepare variants for the queue.
     * @param printStream print stream where to save variants
     */
    public void setOut(PrintStream printStream) {
        out = printStream;
    }

    public PrintStream getOut() {
        return out;
    }

    /**
     * Set delimiter to print variants between fields. Default is <code>\t</code> (tab delimiter).
     * @param delimiter string contains delimiter
     */
    public void setDelimiter(String delimiter) {
        this.delimiter = delimiter;
    }

    /**
     * Factory method for creating needed printer classes for each printer type set in configuration.
     * @param type needed type (usually from instance)
     * @return created specific VariantPrinter
     */
    public static VariantPrinter createPrinter(PrinterType type) {
        switch(type) {
            case OUT: return new SystemOutVariantPrinter();
            default:  return new SystemOutVariantPrinter();
            }
    }
}
