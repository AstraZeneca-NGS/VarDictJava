package com.astrazeneca.vardict;

public class Region {
    final String chr;
    final int start;
    final int end;
    final String gene;
    final int istart;
    final int iend;

    public Region(String chr, int start, int end, String gene) {
        this(chr, start, end, gene, 0, 0);
    }

    public Region(String chr, int start, int end, String gene, int istart, int iend) {
        this.chr = chr;
        this.start = start;
        this.end = end;
        this.gene = gene;
        this.istart = istart;
        this.iend = iend;
    }

    @Override
    public String toString() {
        return "Region [chr=" + chr + ", start=" + start + ", end=" + end + ", gene=" + gene + ", istart=" + istart + ", iend=" + iend + "]";
    }

}