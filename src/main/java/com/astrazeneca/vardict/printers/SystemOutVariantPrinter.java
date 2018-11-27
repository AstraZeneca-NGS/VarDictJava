package com.astrazeneca.vardict.printers;

/**
 * Standard output for variant printer (will print to STDOUT).
 */
public class SystemOutVariantPrinter extends VariantPrinter {
    public SystemOutVariantPrinter() {
        out = System.out;
    }
}
