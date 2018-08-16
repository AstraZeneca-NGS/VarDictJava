package com.astrazeneca.vardict.variations;

import java.util.regex.Matcher;

import static com.astrazeneca.vardict.data.Patterns.ANY_SV;
import static com.astrazeneca.vardict.Utils.substr;

/**
 * Class for holding variant structure
 */
public class Variant {

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
    public String n;

    /**
     * Position coverage
     */
    public int cov;

    /**
     * Forward strand reads count for variant
     */
    public int fwd;

    /**
     * Reverse strand reads count for variant
     */
    public int rev;

    /**
     * Strand bias flag (0, 1 or 2)
     */
    public String bias = "0";

    /**
     * Variant frequency
     */
    public double freq;

    /**
     * Mean variant position in the read
     */
    public double pmean;

    /**
     * Flag that is true when variant is found in at least 2 different positions
     */
    public boolean pstd;

    /**
     * Mean base quality for variant
     */
    public double qual;

    /**
     * Flag that is 1 when variant is read with at least 2 different qualities
     */
    public boolean qstd;

    /**
     * Mean mapping quality for variant
     */
    public double mapq;

    /**
     * Ratio of high-quality reads to low-quality reads
     */
    public double qratio;

    /**
     * Variant frequency for high-quality reads
     */
    public double hifreq;

    /**
     * Adjusted allele frequency for indels due to local realignment
     */
    public double extrafreq;

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
    public double nm;

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
    public String leftseq;

    /**
     * Following reference sequence
     */
    public String rightseq;

    /**
     * Start position
     */
    public int sp;

    /**
     * End position
     */
    public int ep;

    /**
     * Reference variant forward strand coverage
     */
    public int rrc;

    /**
     * Reference variant reverse strand coverage
     */
    public int rfc;

    /**
     * Total position coverage
     */
    public int tcov;

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
    public String varallele;

    /**
     * Reference allele for variation (to be written to .vcf file)
     */
    public String refallele;

    /**
     * Debug information
     */
    public String DEBUG;


    /**
     * A variant is considered noise if the quality is below <code>goodq</code> and
     * there're no more than 3 reads
     * @param goodq quality threshold
     * @param lofreq The minimum allele frequency allowed in normal for a somatic mutation
     * @return Returns <tt>true</tt> if variance is considered noise if the quality is below <code>goodq</code>
     * and there're no more than 3 reads in coverage
     */
    public boolean isNoise(double goodq,
                           double lofreq) {
        final double qual = this.qual;
        if (((qual < 4.5d || (qual < 12 && !this.qstd)) && this.cov <= 3)
                || (qual < goodq && this.freq < 2 * lofreq && this.cov <= 1)) {

            this.tcov -= this.cov;
            this.cov = 0;
            this.fwd = 0;
            this.rev = 0;
            this.freq = 0;
            this.hifreq = 0;

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
            this.sp += n;
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
            this.ep -= n - 1;
            this.refallele = substr(refAllele, 0, 1 - n);
            this.varallele = substr(varAllele, 0, 1 - n);
            this.rightseq = substr(refAllele, 1 - n, n - 1) + substr(this.rightseq, 0, 1 - n);
        }
    }

    /**
     * $varType
     * Find variant type based on variant sequence field <code>refAllele</code> and
     * <code>varAllele</code>
     * @return variant type
     */
    public String varType() {
        Matcher mm = ANY_SV.matcher(varallele);
        if (refallele.length() == 1 && varallele.length() == 1) {
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
}
