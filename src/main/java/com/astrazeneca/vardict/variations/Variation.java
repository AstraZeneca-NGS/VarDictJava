package com.astrazeneca.vardict.variations;

/**
 * Intermediate variant structure
 */
public class Variation {
    /**
     * Variant count
     */
    public int varsCount;

    //TODO: incapsulate those vars counts
    /**
     * Variant count on forward strand
     */
    public int varsCountOnForward;

    /**
     * Variant count on reverse strand
     */
    public int varsCountOnReverse;

    /**
     * Sum of variant positions in read
     */
    public int meanPosition;

    /**
     * Sum of base qualities for variant
     */
    public double meanQuality;

    /**
     * Sum of mapping qualities for variant
     */
    public int meanMappingQuality;

    /**
     * Sum of number of mismatches for variant
     */
    //TODO: in Variation it is int
    public int numberOfMismatches;

    /**
     * Number of low-quality reads with the variant
     */
    public int lowQualityReadsCount;

    /**
     * Number of high-quality reads with the variant
     */
    public int highQualityReadsCount;

    /**
     * Flags that is true when variant is found in at least 2 different positions
     */
    public boolean isAtLeastAt2Positions;

    /**
     * Flags that is 1 when variant is read with at least 2 different qualities
     */
    public boolean hasAtLeast2DiffQualities;

    /**
     * Position in read for previous instance of this variant (used for isAtLeastAt2Positions)
     */
    public int previousPosition;

    /**
     * Base quality for previous instance of this variant (used for hasAtLeast2DiffQualities)
     */
    public double previousQuality;

    /**
     * Adjusted count for indels due to local realignment
     */
    public int extraCount;

    //perl version: 2152
    public void rmCnt(Variation tv) {
        adjVariant(tv, 1);
    }

    public void adjVariant(Variation variant, double adjustmentFactor) {
        varsCount -= adjustmentFactor * variant.varsCount;
        highQualityReadsCount -= adjustmentFactor * variant.highQualityReadsCount;
        lowQualityReadsCount -= adjustmentFactor * variant.lowQualityReadsCount;
        meanPosition -= adjustmentFactor * variant.meanPosition;
        meanQuality -= adjustmentFactor * variant.meanQuality;
        meanMappingQuality -= adjustmentFactor * variant.meanMappingQuality;
        subCountForDir(true, (int)(adjustmentFactor * variant.getCountForDir(true)));
        subCountForDir(false, (int)(adjustmentFactor * variant.getCountForDir(false)));
        correctCnt();
    }

    /**
     * correct counts for negative values
     */
    public void correctCnt() {
        if (this.varsCount < 0)
            this.varsCount = 0;
        if (this.highQualityReadsCount < 0)
            this.highQualityReadsCount = 0;
        if (this.lowQualityReadsCount < 0)
            this.lowQualityReadsCount = 0;
        if (this.meanPosition < 0)
            this.meanPosition = 0;
        if (this.meanQuality < 0)
            this.meanQuality = 0;
        if (this.meanMappingQuality < 0)
            this.meanMappingQuality = 0;
        if (this.getCountForDir(true) < 0)
            this.addCountForDir(true, -this.getCountForDir(true));
        if (this.getCountForDir(false) < 0)
            this.addCountForDir(false, -this.getCountForDir(false));
    }

    /**
     * Increment count for direction
     * @param dir true for forward strand, false for reverse strand
     */
    public void incCountForDir(boolean dir) {
        if (dir)
            this.varsCountOnReverse++;
        else
            this.varsCountOnForward++;
    }

    /**
     * Decrement count for direction
     * @param dir true for forward strand, false for reverse strand
     */
    public void decCountForDir(boolean dir) {
        if (dir)
            this.varsCountOnReverse--;
        else
            this.varsCountOnForward--;
    }

    /**
     * Get variant count for direction
     * @param dir true for forward strand, false for reverse strand
     * @return variant count
     */
    public int getCountForDir(boolean dir) {
        if (dir)
            return this.varsCountOnReverse;
        return this.varsCountOnForward;
    }

    /**
     * Add count for direction
     * @param dir true for forward strand, false for reverse strand
     */
    public void addCountForDir(boolean dir, int add) {
        if (dir)
            this.varsCountOnReverse += add;
        else
            this.varsCountOnForward += add;
    }

    /**
     * Subtract count for direction
     * @param dir true for forward strand, false for reverse strand
     */
    public void subCountForDir(boolean dir, int sub) {
        if (dir)
            this.varsCountOnReverse -= sub;
        else
            this.varsCountOnForward -= sub;
    }

}