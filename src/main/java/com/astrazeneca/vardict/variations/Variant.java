package com.astrazeneca.vardict.variations;

import java.text.DecimalFormat;
import java.util.List;
import java.util.Set;
import java.util.regex.Matcher;

import static com.astrazeneca.vardict.data.scopedata.GlobalReadOnlyScope.instance;
import static com.astrazeneca.vardict.data.Patterns.ANY_SV;
import static com.astrazeneca.vardict.Utils.substr;

/**
 * Class for holding variant structure
 */
public class Variant {

    /**
     * Variant description string $n
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
     * Position coverage $ cov
     */
    public int positionCoverage;

    /**
     * Forward strand reads count for variant $fwd
     */
    public int varsCountOnForward;

    /**
     * Reverse strand reads count for variant  $rev
     */
    public int varsCountOnReverse;

    /**
     * Strand bias flag (0, 1 or 2) $bias
     */
    public String strandBiasFlag = "0";

    /**
     * Variant frequency $freq
     */
    public double frequency;

    /**
     * Mean variant position in the read $pmean
     */
    public double meanPosition;

    /**
     * Flag that is true when variant is found in at least 2 different positions $pstd
     */
    public boolean isAtLeastAt2Positions;

    /**
     * Mean base quality for variant $pmean
     */
    public double meanQuality;

    /**
     * Flag that is true when variant is read with at least 2 different qualities $qstd
     */
    public boolean hasAtLeast2DiffQualities;

    /**
     * Mean mapping quality for variant $mapq
     */
    public double meanMappingQuality;

    /**
     * Ratio of high-quality reads to low-quality reads $qratio
     */
    public double highQualityToLowQualityRatio;

    /**
     * Variant frequency for high-quality reads $hifreq
     */
    public double highQualityReadsFrequency;

    /**
     * Adjusted allele frequency for indels due to local realignment $extrafreq
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
     * Average number of mismatches for reads containing variant $nm
     */
    public double numberOfMismatches;

    /**
     * Number of high-quality reads with the variant
     */
    public int hicnt;

    /**
     * Position coverage by high quality reads
     */
    public int hicov;

    /**
     * Preceding reference sequence
     */
    public String leftseq = "";

    /**
     * Following reference sequence
     */
    public String rightseq = "";

    /**
     * Start position $sp
     */
    public int startPosition;

    /**
     * End position $ep
     */
    public int endPosition;

    /**
     * Reference variant forward strand coverage $rfc
     */
    public int refReverseCoverage;

    /**
     * Reference variant reverse strand coverage $rrc
     */
    public int refForwardCoverage;

    /**
     * Total position coverage $tcov
     */
    public int totalPosCoverage;

    /**
     * Duplication rate
     */
    public double duprate;

    /**
     * Genotype description string
     */
    public String genotype;

    /**
     * Variant allele (to be written to .vcf file)
     */
    public String varallele = "";

    /**
     * Reference allele for variation (to be written to .vcf file)
     */
    public String refallele = "";

    public String vartype = "";

    /**
     * Debug information
     */
    public String DEBUG = "";

    /**
     * CRISPR information
     */
    public int crispr;

    /**
     * A variant is considered noise if the quality is below <code>goodq</code> and
     * there're no more than 3 reads
     * @return Returns true if variance is considered noise if the quality is below <code>goodq</code>
     * and there're no more than 3 reads in coverage
     */
    public boolean isNoise() {
        final double qual = this.meanQuality;
        if (((qual < 4.5d || (qual < 12 && !this.hasAtLeast2DiffQualities)) && this.positionCoverage <= 3)
                || (qual < instance().conf.goodq
                && this.frequency < 2 * instance().conf.lofreq
                && this.positionCoverage <= 1)) {

            this.totalPosCoverage -= this.positionCoverage;
            this.positionCoverage = 0;
            this.varsCountOnForward = 0;
            this.varsCountOnReverse = 0;
            this.frequency = 0;
            this.highQualityReadsFrequency = 0;

            return true;

        }
        return false;
    }

    /**
     * Adjust the complex variant
     */
    public void adjComplex() {
        String refAllele = this.refallele;
        String varAllele = this.varallele;
        if (varAllele.charAt(0) == '<') {
            return;
        }
        int n = 0;
        while (refAllele.length() - n > 1
                && varAllele.length() - n > 1
                && refAllele.charAt(n) == varAllele.charAt(n)) {
            n++;
        }
        if (n > 0) {
            this.startPosition += n;
            this.refallele = substr(refAllele, n);
            this.varallele = substr(varAllele, n);
            this.leftseq += substr(refAllele, 0, n);
            this.leftseq = substr(this.leftseq, n);
        }
        refAllele = this.refallele;
        varAllele = this.varallele;
        n = 1;
        while (refAllele.length() - n > 0
                && varAllele.length() - n > 0
                && substr(refAllele, -n, 1).equals(substr(varAllele, -n, 1))) {
            n++;
        }
        if (n > 1) {
            this.endPosition -= n - 1;
            this.refallele = substr(refAllele, 0, 1 - n);
            this.varallele = substr(varAllele, 0, 1 - n);
            this.rightseq = substr(refAllele, 1 - n, n - 1) + substr(this.rightseq, 0, 1 - n);
        }
    }

    public void debugVariantsContentSimple(List<String> tmp, String n) {
        StringBuilder sb = debugVariantsContent(n);
        tmp.add(sb.toString());
    }

    public void debugVariantsContentInsertion(List<String> tmp, String n) {
        StringBuilder sb = new StringBuilder();
        sb.append("I");
        sb.append(debugVariantsContent(n));
        tmp.add(sb.toString());
    }

    private StringBuilder debugVariantsContent(String n) {
        StringBuilder sb = new StringBuilder();
        sb.append(n
                + ":" + (varsCountOnForward + varsCountOnReverse)
                + ":F-" + varsCountOnForward
                + ":R-" + varsCountOnReverse
                + ":" + new DecimalFormat("0.0000").format(frequency)
                + ":" + strandBiasFlag
                + ":" + new DecimalFormat("0.0").format(meanPosition)
                + ":" + (isAtLeastAt2Positions ? "1" : "0")
                + ":" + new DecimalFormat("0.0").format(meanQuality)
                + ":" + (hasAtLeast2DiffQualities ? "1" : "0")
                + ":" + new DecimalFormat("0.0000").format(highQualityReadsFrequency)
                + ":" + meanMappingQuality
                + ":" + new DecimalFormat("0.000").format(highQualityToLowQualityRatio));
        return sb;
    }

    /**
     * $varType
     * Find variant type based on variant sequence field <code>refAllele</code> and
     * <code>varAllele</code>
     * @return variant type
     */
    public String varType() {
        Matcher mm = ANY_SV.matcher(varallele);
        if (refallele.equals(varallele) && refallele.length() == 1) {
            return "";
        } else if (refallele.length() == 1 && varallele.length() == 1) {
            return "SNV";
        } else if (mm.find()) {
            return mm.group(1);
        } else if (refallele.length() == 0 || varallele.length() == 0) { //issue #19
            return "Complex";
        } else if (refallele.charAt(0) != varallele.charAt(0)) {
            return "Complex";
        } else if (refallele.length() == 1 && varallele.length() > 1
                && varallele.startsWith(refallele)) {
            return "Insertion";
        } else if (refallele.length() > 1 && varallele.length() == 1
                && refallele.startsWith(varallele)) {
            return "Deletion";
        }
        return "Complex";
    }

    /**
     * Returns true whether a variant meet specified criteria
     * @param referenceVar reference variant    $rref
     * @param type Type of variant
     * @param splice set of strings representing introns in splice
     * @return true if variant meet specified criteria
     */
    public boolean isGoodVar(Variant referenceVar, String type,
                             Set<String> splice) {
        if (this == null || this.refallele == null || this.refallele.isEmpty()) {
            return false;
        }
        if (type == null || type.isEmpty()) {
            type = varType();
        }
        if (frequency < instance().conf.freq
                || hicnt < instance().conf.minr
                || meanPosition < instance().conf.readPosFilter
                || meanQuality < instance().conf.goodq) {
            return false;
        }

        if (referenceVar != null && referenceVar.hicnt > instance().conf.minr && frequency < 0.25d) {
            //The reference allele has much better mean mapq than var allele, thus likely false positives
            double d = meanMappingQuality + refallele.length() + varallele.length();
            double f = (1 + d) / (referenceVar.meanMappingQuality + 1);
            if ((d - 2 < 5 && referenceVar.meanMappingQuality > 20)
                    || f < 0.25d) {
                return false;
            }
        }

        if (type.equals("Deletion") && splice.contains(startPosition + "-" + endPosition)) {
            return false;
        }
        if (highQualityToLowQualityRatio < instance().conf.qratio) {
            return false;
        }
        if (frequency > 0.30d) {
            return true;
        }
        if (meanMappingQuality < instance().conf.mapq) {
            return false;
        }
        if (msi >= 15 && frequency <= instance().conf.monomerMsiFrequency && msint == 1) {
            return false;
        }
        if (msi >= 12 && frequency <= instance().conf.nonMonomerMsiFrequency && msint > 1) {
            return false;
        }
        if (strandBiasFlag.equals("2;1") && frequency < 0.20d) {
            if (type == null || type.equals("SNV") || (refallele.length() < 3 && varallele.length() < 3)) {
                return false;
            }
        }
        return true;
    }

    @Override
    public String toString() {
        return "Variant{" +
                "descriptionString='" + descriptionString + '\'' +
                ", positionCoverage=" + positionCoverage +
                ", varsCountOnForward=" + varsCountOnForward +
                ", varsCountOnReverse=" + varsCountOnReverse +
                ", strandBiasFlag='" + strandBiasFlag + '\'' +
                ", frequency=" + frequency +
                ", meanPosition=" + meanPosition +
                ", isAtLeastAt2Positions=" + isAtLeastAt2Positions +
                ", meanQuality=" + meanQuality +
                ", hasAtLeast2DiffQualities=" + hasAtLeast2DiffQualities +
                ", meanMappingQuality=" + meanMappingQuality +
                ", highQualityToLowQualityRatio=" + highQualityToLowQualityRatio +
                ", highQualityReadsFrequency=" + highQualityReadsFrequency +
                ", extraFrequency=" + extraFrequency +
                ", shift3=" + shift3 +
                ", msi=" + msi +
                ", msint=" + msint +
                ", numberOfMismatches=" + numberOfMismatches +
                ", hicnt=" + hicnt +
                ", hicov=" + hicov +
                ", leftseq='" + leftseq + '\'' +
                ", rightseq='" + rightseq + '\'' +
                ", startPosition=" + startPosition +
                ", endPosition=" + endPosition +
                ", refReverseCoverage=" + refReverseCoverage +
                ", refForwardCoverage=" + refForwardCoverage +
                ", totalPosCoverage=" + totalPosCoverage +
                ", duprate=" + duprate +
                ", genotype='" + genotype + '\'' +
                ", varallele='" + varallele + '\'' +
                ", refallele='" + refallele + '\'' +
                ", vartype='" + vartype + '\'' +
                ", crispr='" + crispr + '\'' +
                ", DEBUG='" + DEBUG + '\'' +
                '}';
    }
}
