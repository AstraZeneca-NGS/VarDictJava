package com.astrazeneca.vardict.variations;

import static com.astrazeneca.utils.Utils.substr;
import static java.lang.String.format;

/**
 * Class for holding variant structure
 */
public class Variant extends Variation {

    /**
     * Variant description string
     * Description string format:
     * 1). single letter                   - for SNPs
     * 2). + sequence                      - for insertions
     * 3). - number                        - for deletions
     * 4). ... # sequence                  - for insertion/deletion variants followed by short matched sequence
     * 5). ... ^ sequence                  - followed by insertion
     * 6). ... ^ number                    - followed by deletion
     * 7). ... &amp; sequence                  - for insertion/deletion variants followed by matched sequence
     */
    public String descriptionString;

    /**
     * Position coverage
     */
    public int positionCoverage;

    /**
     * Strand bias flag (0, 1 or 2)
     */
    public String strandBiasFlag = "0";

    /**
     * Variant frequency
     */
    public double frequency;

    /**
     * Mean variant position in the read
     */
    //TODO: in Variation it is int
    public double meanPosition;

    /**
     * Mean mapping quality for variant
     */
    //TODO: in Variation it is int
    public double meanMappingQuality;

    /**
     * Ratio of high-quality reads to low-quality reads
     */
    public double highQualityToLowQualityRatio;

    /**
     * Variant frequency for high-quality reads
     */
    public double highQualityReadsFrequency;

    /**
     * Adjusted allele frequency for indels due to local realignment
     */
    public double extraFrequency;

    /**
     * No. of bases to be shifted to 3' for deletions due to alternative alignment
     */
    public int shift3;

    /**
     * msi. &gt; 1 indicates Microsatellite instability
     */
    public double msi;

    /**
     * MicroSattelite unit length in base pairs
     */
    public int msint;

    /**
     * Average number of mismatches for reads containing variant
     */
    public double numberOfMismatches;

    /**
     * Position coverage by high quality reads
     */
    public int highQPosCoverage;

    /**
     * Preceding reference sequence
     */
    public String precedingRefSequence;

    /**
     * Following reference sequence
     */
    public String followingRefSequence;

    /**
     * Start position
     */
    public int startPosition;

    /**
     * End position
     */
    public int endPosition;

    /**
     * Reference variant reverse strand coverage
     */
    public int refReverseCoverage;

    /**
     * Reference variant forward strand coverage
     */
    public int refForwardCoverage;

    /**
     * Total position coverage
     */
    public int totalPosCoverage;

    /**
     * Genotype description string
     */
    public String genotype;

    /**
     * Variant allele (to be written to .vcf file)
     */
    public String varAllele;

    /**
     * Reference allele for variation (to be written to .vcf file)
     */
    public String refAllele;

    /**
     * Debug information
     */
    public String DEBUG;

    public String addDelimiter(String stringToAdd) {
        return String.valueOf(totalPosCoverage) + stringToAdd +
                positionCoverage + stringToAdd +
                refForwardCoverage + stringToAdd +
                refReverseCoverage + stringToAdd +
                varsCountOnForward + stringToAdd +
                varsCountOnReverse + stringToAdd +
                genotype + stringToAdd +
                (frequency == 0 ? 0 : format("%.4f", frequency)) + stringToAdd +
                strandBiasFlag + stringToAdd +
                format("%.1f", meanPosition) + stringToAdd +
                (isAtLeastAt2Positions ? 1 : 0) + stringToAdd +
                format("%.1f", meanQuality) + stringToAdd +
                (hasAtLeast2DiffQualities ? 1 : 0) + stringToAdd +
                format("%.1f", meanMappingQuality) + stringToAdd +
                format("%.3f", highQualityToLowQualityRatio) + stringToAdd +
                (highQualityReadsFrequency == 0 ? 0 : format("%.4f", highQualityReadsFrequency)) + stringToAdd +
                (extraFrequency == 0 ? 0 : format("%.4f", extraFrequency));
    }

    public String addDelimiterExtended(String stringToAdd) {
        return String.valueOf(startPosition) + stringToAdd +
                endPosition + stringToAdd +
                refAllele + stringToAdd +
                varAllele + stringToAdd +
                addDelimiter(stringToAdd) + stringToAdd +
                shift3 + stringToAdd +
                (msi == 0 ? 0 : format("%.3f", msi)) + stringToAdd +
                msint + stringToAdd +
                format("%.1f", numberOfMismatches) + stringToAdd +
                highQualityReadsCount + stringToAdd +
                highQPosCoverage + stringToAdd +
                precedingRefSequence + stringToAdd +
                followingRefSequence;
    }

    /**
     * A variance is considered noise if the quality is below <code>threshold</code> and there're no more than 3 reads
     *
     * @param threshold quality threshold
     * @param lofreq The minimun alelle frequency allowed in normal for a somatic mutation
     * @return Returns <tt>true</tt> if variance is considered noise if the quality is below <code>threshold</code> and there're no more
     *         than 3 reads
     */
    //perl version: 509
    public boolean isNoise(double threshold, double lofreq) {
        final double qual = meanQuality;
        if (((qual < 4.5d || (qual < 12 && !hasAtLeast2DiffQualities)) && positionCoverage <= 3)
                || (qual < threshold && frequency < 2 * lofreq && positionCoverage <= 1)) {

            totalPosCoverage -= positionCoverage;
            positionCoverage = 0;
            varsCountOnForward = 0;
            varsCountOnReverse = 0;
            frequency = 0;
            highQualityReadsFrequency = 0;

            return true;

        }
        return false;
    }

    /**
     * Adjust the complex variant
     */
    //perl version: 345
    public void adjComplex() {
        String refAllele = this.refAllele;
        String varAllele = this.varAllele;
        int n = 0;
        while (refAllele.length() - n > 1 && varAllele.length() - n > 1 && refAllele.charAt(n) == varAllele.charAt(n)) {
            n++;
        }
        if (n > 0) {
            startPosition += n;
            this.refAllele = substr(refAllele, n);
            this.varAllele = substr(varAllele, n);
            precedingRefSequence += substr(refAllele, 0, n);
            precedingRefSequence = substr(precedingRefSequence, n);
        }
        refAllele = this.refAllele;
        varAllele = this.varAllele;
        n = 1;
        while (refAllele.length() - n > 0 && varAllele.length() - n > 0 && substr(refAllele, -n, 1).equals(substr(varAllele, -n, 1))) {
            n++;
        }
        if (n > 1) {
            endPosition -= n - 1;
            this.refAllele = substr(refAllele, 0, 1 - n);
            this.varAllele = substr(varAllele, 0, 1 - n);
            followingRefSequence = substr(refAllele, 1 - n, n - 1) + substr(followingRefSequence, 0, 1 - n);
        }

    }

    /**
     * Find variant type based on variant sequence field <code>refAllele</code> and
     * <code>varAllele</code>
     * @return variant type
     */
    //perl version: 291
    public Type getType() {
        if (refAllele.length() == 1 && varAllele.length() == 1) {
            return Type.SNV;
        } else if (refAllele.length() == 0 || varAllele.length() == 0) { //issue #19
            return Type.complex;
        } else if (refAllele.charAt(0) != varAllele.charAt(0)) {
            return Type.complex;
        } else if (refAllele.length() == 1 && varAllele.length() > 1 && varAllele.startsWith(refAllele)) {
            return Type.insertion;
        } else if (refAllele.length() > 1 && varAllele.length() == 1 && refAllele.startsWith(varAllele)) {
            return Type.deletion;
        }
        return Type.complex;
    }

    //vartype may take values SNV (Single Nucleotide Variant), Complex (or MNV (Multiple Nucleotide Variant)), Insertion, Deletion
    public enum Type {
        SNV("SNV"), complex("Complex"), insertion("Insertion"), deletion("Deletion"), noInfo("");
        private String stringName;

        Type(String stringName) {
            this.stringName = stringName;
        }

        @Override
        public String toString(){
            return stringName;
        }

    }
}
