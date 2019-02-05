package com.astrazeneca.vardict.printers;

/**
 * Printer types for further extending of possible types to rule them through command line. The needed variant printer
 * will be created for each type of PrinterType.
 */
public enum PrinterType {
    OUT("OUT"),
    ERR("ERR");

    private final String out;

    PrinterType(String out) {
        this.out = out;
    }
}
