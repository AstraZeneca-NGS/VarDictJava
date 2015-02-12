package com.astrazeneca.vardict;

/**
 * Class for holding region from BED file
 */
public class Region {
    /**
     * Chromosome name
     */
    final String chr;

    /**
     * Region start position
     */
    final int start;

    /**
     * Region end position
     */
    final int end;

    /**
     * Gene name
     */
    final String gene;

    /**
     * Position to start looking for variants (for amplicon based calling)
     */
    final int istart;

    /**
     * Position to end looking for variants (for amplicon based calling)
     */
    final int iend;

    /**
     * Constructor for 4-column BED file line
     * @param chr chromosome name
     * @param start start position
     * @param end end position
     * @param gene gene name
     */
    public Region(String chr, int start, int end, String gene) {
        this(chr, start, end, gene, 0, 0);
    }

    /**
     * Constructor for 8-column BED file line
     * @param chr chromosome name
     * @param start start position
     * @param end end position
     * @param gene gene name
     * @param istart Position to start looking for variants (for amplicon based calling)
     * @param iend Position to end looking for variants (for amplicon based calling)
     */
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