package com.astrazeneca.vardict.variations;

/**
 * Intermediate variant structure
 */
public class Variation {
    /**
     * Variant count $cnt
     */
    public int varsCount;

    /**
     * Variant count on forward strand $dirPlus
     */
    public int varsCountOnForward;

    /**
     * Variant count on reverse strand $dirMinus
     */
    public int varsCountOnReverse;

    /**
     * Sum of variant positions in read $pmean
     */
    public double meanPosition;

    /**
     * Sum of base qualities for variant $qmean
     */
    public double meanQuality;

    /**
     * Sum of mapping qualities for variant $Qmean
     */
    public double meanMappingQuality;

    /**
     * Sum of number of mismatches for variant  $nm
     */
    public double numberOfMismatches;

    /**
     * Number of low-quality reads with the variant $locnt
     */
    public int lowQualityReadsCount;

    /**
     * Number of high-quality reads with the variant $hicnt
     */
    public int highQualityReadsCount;

    /**
     * Flags that is true when variant is found in at least 2 different positions $pstd
     */
    public boolean pstd;

    /**
     * Flags that is 1 when variant is read with at least 2 different qualities $qstd
     */
    public boolean qstd;

    /**
     * Position in read for previous instance of this variant (used for pstd) $pp
     */
    public int pp;

    /**
     * Base quality for previous instance of this variant (used for qstd)  $pq
     */
    public double pq;

    /**
     * Adjusted count for indels due to local realignment $extracnt
     */
    public int extracnt;

    /**
     * Increment count for direction
     * @param dir false for forward strand, true for reverse strand
     */
    public void incDir(boolean dir) {
        if (dir)
            this.varsCountOnReverse++;
        else
            this.varsCountOnForward++;
    }

    /**
     * Decrement count for direction
     * @param dir false for forward strand, true for reverse strand
     */
    public void decDir(boolean dir) {
        if (dir)
            this.varsCountOnReverse--;
        else
            this.varsCountOnForward--;
    }

    /**
     * Get variant count for direction
     * @param dir false for forward strand, true for reverse strand
     * @return variant count
     */
    public int getDir(boolean dir) {
        if (dir)
            return this.varsCountOnReverse;
        return this.varsCountOnForward;
    }

    /**
     * Add count for direction
     * @param dir false for forward strand, true for reverse strand
     * @param add amount of counts need to be added in the specific direction
     */
    public void addDir(boolean dir, int add) {
        if (dir)
            this.varsCountOnReverse += add;
        else
            this.varsCountOnForward += add;
    }

    /**
     * Subtract count for direction
     * @param dir false for forward strand, true for reverse strand
     * @param sub amount of counts need to be subtracted in the specific direction
     */
    public void subDir(boolean dir, int sub) {
        if (dir)
            this.varsCountOnReverse -= sub;
        else
            this.varsCountOnForward -= sub;
    }

    @Override
    public String toString() {
        return "Variation{" +
                "varsCount=" + varsCount +
                ", varsCountOnForward=" + varsCountOnForward +
                ", varsCountOnReverse=" + varsCountOnReverse +
                ", meanPosition=" + meanPosition +
                ", meanQuality=" + meanQuality +
                ", meanMappingQuality=" + meanMappingQuality +
                ", numberOfMismatches=" + numberOfMismatches +
                ", lowQualityReadsCount=" + lowQualityReadsCount +
                ", highQualityReadsCount=" + highQualityReadsCount +
                ", pstd=" + pstd +
                ", qstd=" + qstd +
                ", pp=" + pp +
                ", pq=" + pq +
                ", extracnt=" + extracnt +
                '}';
    }
}
