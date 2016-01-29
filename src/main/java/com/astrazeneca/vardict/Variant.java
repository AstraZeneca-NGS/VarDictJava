package com.astrazeneca.vardict;

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
    String n;

    /**
     * Position coverage
     */
    int cov;

    /**
     * Forward strand reads count for variant
     */
    int fwd;

    /**
     * Reverse strand reads count for variant
     */
    int rev;

    /**
     * Strand bias flag (0, 1 or 2)
     */
    String bias = "0";

    /**
     * Variant frequency
     */
    double freq;

    /**
     * Mean variant position in the read
     */
    double pmean;

    /**
     * Flag that is true when variant is found in at least 2 different positions
     */
    boolean pstd;

    /**
     * Mean base quality for variant
     */
    double qual;

    /**
     * Flag that is 1 when variant is read with at least 2 different qualities
     */
    boolean qstd;

    /**
     * Mean mapping quality for variant
     */
    double mapq;

    /**
     * Ratio of high-quality reads to low-quality reads
     */
    double qratio;

    /**
     * Variant frequency for high-quality reads
     */
    double hifreq;

    /**
     * Adjusted allele frequency for indels due to local realignment
     */
    double extrafreq;

    /**
     * No. of bases to be shifted to 3' for deletions due to alternative alignment
     */
    int shift3;

    /**
     * msi. &gt; 1 indicates Microsatellite instability
     */
    double msi;

    /**
     * MicroSattelite unit length in base pairs
     */
    int msint;

    /**
     * Average number of mismatches for reads containing variant
     */
    double nm;

    /**
     * Number of high-quality reads with the variant
     */
    int hicnt;

    /**
     * Position coverage by high quality reads
     */
    int hicov;

    /**
     * Preceding reference sequence
     */
    String leftseq;

    /**
     * Following reference sequence
     */
    String rightseq;

    /**
     * Start position
     */
    int sp;

    /**
     * End position
     */
    int ep;

    /**
     * Reference variant forward strand coverage
     */
    int rrc;

    /**
     * Reference variant reverse strand coverage
     */
    int rfc;

    /**
     * Total position coverage
     */
    int tcov;

    /**
     * Genotype description string
     */
    String genotype;

    /**
     * Variant allele (to be written to .vcf file)
     */
    String varallele;

    /**
     * Reference allele for variation (to be written to .vcf file)
     */
    String refallele;

    /**
     * Debug information
     */
    String DEBUG;
}
