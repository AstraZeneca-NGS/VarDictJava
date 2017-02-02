package com.astrazeneca.vardict;

/**
 * Class for holding region from BED file
 */
public class Region {
    /**
     * Chromosome name
     */
    public final String chr;

    /**
     * Region start position
     */
    public final int start;

    /**
     * Region end position
     */
    public final int end;

    /**
     * Gene name
     */
    public final String gene;

    /**
     * Position to start looking for variants (for amplicon based calling)
     */
    public final int variantsStart;

    /**
     * Position to end looking for variants (for amplicon based calling)
     */
    public final int variantsEnd;

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
     * @param variantsStart Position to start looking for variants (for amplicon based calling)
     * @param variantsEnd Position to end looking for variants (for amplicon based calling)
     */
    public Region(String chr, int start, int end, String gene, int variantsStart, int variantsEnd) {
        this.chr = chr;
        this.start = start;
        this.end = end;
        this.gene = gene;
        this.variantsStart = variantsStart;
        this.variantsEnd = variantsEnd;
    }

    @Override
    public String toString() {
        return "Region [chr=" + chr + ", start=" + start + ", end=" + end + ", gene=" + gene + ", variantsStart=" + variantsStart + ", variantsEnd=" + variantsEnd + "]";
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        Region region = (Region) o;

        if (start != region.start) return false;
        if (end != region.end) return false;
        if (variantsStart != region.variantsStart) return false;
        if (variantsEnd != region.variantsEnd) return false;
        if (!chr.equals(region.chr)) return false;
        return gene.equals(region.gene);

    }
}