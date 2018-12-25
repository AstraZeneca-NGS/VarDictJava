package com.astrazeneca.vardict.printers;

/**
 * Standard output for variant printer (will print to STDOUT).
 */
public class SystemErrVariantPrinter extends VariantPrinter {
    public SystemErrVariantPrinter() {
        out = System.err;
    }
}
