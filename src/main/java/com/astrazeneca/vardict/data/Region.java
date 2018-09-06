package com.astrazeneca.vardict.data;

import java.util.Map;

import static com.astrazeneca.vardict.Utils.correctChr;
import static com.astrazeneca.vardict.Utils.toInt;

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
    public final int istart;

    /**
     * Position to end looking for variants (for amplicon based calling)
     */
    public final int iend;

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

    public static Region newModifiedRegion(Region region, int changedStart, int changedEnd) {
        return new Region(region.chr, changedStart, changedEnd, region.gene);
    }

    /**
     * Create region from command-line option
     * @param region Region string from command line
     * @param numberNucleotideToExtend number of nucleotides to extend region
     * @param chrs map of chromosome lengths
     * @param zeroBased Are positions zero-based or 1-based
     * @return Region
     */
    public static Region buildRegion(String region, final int numberNucleotideToExtend, Map<String, Integer> chrs,
                                      final boolean zeroBased) {
        String[] split = region.split(":");
        String chr = split[0];
        chr = correctChr(chrs, chr);
        String gene = split.length < 3 ? chr : split[2];
        String[] range = split[1].split("-");
        int start = toInt(range[0].replaceAll(",", ""));
        int end = range.length < 2 ? start : toInt(range[1].replaceAll(",", ""));
        start -= numberNucleotideToExtend;
        end += numberNucleotideToExtend;
        if (zeroBased && start < end) {
            start++;
        }
        if (start > end)
            start = end;

        return new Region(chr, start, end, gene);
    }


    @Override
    public String toString() {
        return "Region [chr=" + chr + ", start=" + start + ", end=" + end + ", gene=" + gene + ", istart=" + istart + ", iend=" + iend + "]";
    }

}