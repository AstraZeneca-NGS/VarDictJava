package com.astrazeneca.vardict.data;


import java.util.Objects;

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
    public final int insertStart;

    /**
     * Position to end looking for variants (for amplicon based calling)
     */
    public final int insertEnd;

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
     * @param insertStart Position to start looking for variants (for amplicon based calling)
     * @param insertEnd Position to end looking for variants (for amplicon based calling)
     */
    public Region(String chr, int start, int end, String gene, int insertStart, int insertEnd) {
        this.chr = chr;
        this.start = start;
        this.end = end;
        this.gene = gene;
        this.insertStart = insertStart;
        this.insertEnd = insertEnd;
    }

    /**
     * Method returns new Region created on the base of Region from parameter. Used for extended regions
     * in StructuralVariantsProcessor and VariationRealigner
     * @param region region to create new from
     * @param changedStart new start of region
     * @param changedEnd new end of region
     * @return created Region
     */
    public static Region newModifiedRegion(Region region, int changedStart, int changedEnd) {
        return new Region(region.chr, changedStart, changedEnd, region.gene, region.insertStart, region.insertEnd);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        Region region = (Region) o;
        return start == region.start &&
                end == region.end &&
                insertStart == region.insertStart &&
                insertEnd == region.insertEnd &&
                Objects.equals(chr, region.chr) &&
                Objects.equals(gene, region.gene);
    }

    @Override
    public int hashCode() {
        return Objects.hash(chr, start, end, gene, insertStart, insertEnd);
    }

    @Override
    public String toString() {
        return "Region [chr=" + chr + ", start=" + start + ", end=" + end + ", gene=" + gene + ", istart=" + insertStart + ", iend=" + insertEnd + "]";
    }

    public String printRegion(){
        return chr + ":" + start + "-" + end;
    }

}