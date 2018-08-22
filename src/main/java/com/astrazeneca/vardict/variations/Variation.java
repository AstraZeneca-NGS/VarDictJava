package com.astrazeneca.vardict.variations;

/**
 * Intermediate variant structure
 */
public class Variation {
    /**
     * Variant count
     */
    public int cnt;

    /**
     * Variant count on forward strand
     */
    public int dirPlus;

    /**
     * Variant count on reverse strand
     */
    public int dirMinus;

    /**
     * Sum of variant positions in read
     */
    public double pmean;

    /**
     * Sum of base qualities for variant
     */
    public double qmean;

    /**
     * Sum of mapping qualities for variant
     */
    public double Qmean;

    /**
     * Sum of number of mismatches for variant
     */
    public double nm;

    /**
     * Number of low-quality reads with the variant
     */
    public int locnt;

    /**
     * Number of high-quality reads with the variant
     */
    public int hicnt;

    /**
     * Flags that is true when variant is found in at least 2 different positions
     */
    public boolean pstd;

    /**
     * Flags that is 1 when variant is read with at least 2 different qualities
     */
    public boolean qstd;

    /**
     * Position in read for previous instance of this variant (used for pstd)
     */
    public int pp;

    /**
     * Base quality for previous instance of this variant (used for qstd)
     */
    public double pq;

    /**
     * Adjusted count for indels due to local realignment
     */
    public int extracnt;

    /**
     * Increment count for direction
     * @param dir false for forward strand, true for reverse strand
     */
    public void incDir(boolean dir) {
        if (dir)
            this.dirMinus++;
        else
            this.dirPlus++;
    }

    /**
     * Decrement count for direction
     * @param dir false for forward strand, true for reverse strand
     */
    public void decDir(boolean dir) {
        if (dir)
            this.dirMinus--;
        else
            this.dirPlus--;
    }

    /**
     * Get variant count for direction
     * @param dir false for forward strand, true for reverse strand
     * @return variant count
     */
    public int getDir(boolean dir) {
        if (dir)
            return this.dirMinus;
        return this.dirPlus;
    }

    /**
     * $addDir
     * Add count for direction
     * @param dir false for forward strand, true for reverse strand
     */
    public void addDir(boolean dir, int add) {
        if (dir)
            this.dirMinus += add;
        else
            this.dirPlus += add;
    }

    /**
     * $subDir
     * Subtract count for direction
     * @param dir false for forward strand, true for reverse strand
     */
    public void subDir(boolean dir, int sub) {
        if (dir)
            this.dirMinus -= sub;
        else
            this.dirPlus -= sub;
    }

}
