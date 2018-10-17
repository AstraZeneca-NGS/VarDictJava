package com.astrazeneca.vardict.variations;

import com.astrazeneca.vardict.data.Region;

import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.List;
import java.util.regex.Matcher;

import static com.astrazeneca.vardict.Utils.join;
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

    public String vartype;

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

    public void callingSimpleVariant(String sample, Region region, PrintStream out, String sv) {
        out.println(join("\t",
                joinSimpleVariant(sample, region),
                region.chr + ":" + region.start + "-" + region.end,
                vartype,
                duprate == 0 ? 0 : new DecimalFormat("0.0").format(duprate),
                sv.equals("") ? 0 : sv
        ));
    }

    public String joinSimpleVariant(String sample, Region region) {
        return join("\t",
                outputStartOfVariantLine(sample, region),

                joinVar2("\t"),

                shift3,
                msi == 0 ? 0 : new DecimalFormat("0.000").format(msi),
                msint,
                nm > 0 ? new DecimalFormat("0.0").format(nm) : 0,
                hicnt,
                hicov,
                leftseq.isEmpty() ? "0" : leftseq,
                rightseq.isEmpty() ? "0" : rightseq
        );
    }

    public static void callingEmptySimpleVariant(Region region, String sample, PrintStream out,
                                                 int position) {
        out.println(join("\t",
                sample, region.gene, region.chr,
                position,
                position,
                "", "", 0, 0, 0, 0, 0, 0, "", 0, "0;0", 0, 0,
                0, 0, 0, "", 0, 0, 0, 0, 0, 0, "", "", 0, 0,
                region.chr + ":" + region.start + "-" + region.end,
                "", 0, 0
        ));
    }

    public void callingOneSample(Region region, String sample,
                                 PrintStream out, String tvf, String type, String sv, boolean isFirstCover) {
        String infoVariantWithNM = isFirstCover
                ? join("\t", tvf, joinVariantWithNM(this))
                : join("\t", joinVariantWithNM(this), tvf);

        String infoDuprateSV = isFirstCover
                ? join("\t", 0, 0,
                duprate == 0 ? 0 : new DecimalFormat("0.0").format(duprate), sv.equals("") ? 0 : sv)
                : join("\t",
                duprate == 0 ? 0 : new DecimalFormat("0.0").format(duprate), sv.equals("") ? 0 : sv, 0, 0);

        out.println(join("\t",
                outputStartOfVariantLine(sample, region),
                infoVariantWithNM,
                outputMiddleOfVariantLine(this),
                outputEndOfVariantLine(region, type),
                infoDuprateSV)
        );
    }

    public void constructBothSamplesWithoutSecondVariant(Region region, String sample,
                                                         PrintStream out, String info, String type,
                                                         Variant variant2,
                                                         String v1sv, String v2sv) {
        out.println(join("\t",
                outputStartOfVariantLine(sample, region),
                info,
                outputMiddleOfVariantLine(this),
                outputEndOfVariantLine(region, type),
                outputEndOfVariantDuprate(v1sv, v2sv, variant2))
        );
    }

    public void constructBothSamplesWithZeroDuprateSecondVariant(Region region, String sample,
                                                         PrintStream out, String info, String type,
                                                         String v2sv) {
        out.println(join("\t",
                outputStartOfVariantLine(sample, region),
                info,
                outputMiddleOfVariantLine(this),
                outputEndOfVariantLine(region, type),
                0, 0,
                duprate == 0 ? 0 : new DecimalFormat("0.000").format(duprate),
                v2sv.equals("") ? 0 : v2sv
        ));
    }

    public void constructBothSamplesWithSecondVariant(Region region, String sample, PrintStream out,
                                                      String info, Variant variant2, String type,
                                                      String v1sv, String v2sv) {
        out.println(join("\t",
                outputStartOfVariantLine(sample, region),
                info,
                outputMiddleOfVariantLine(variant2),
                outputEndOfVariantLine(region, type),
                outputEndOfVariantDuprate(v1sv, v2sv, variant2))
        );
    }

    private String outputStartOfVariantLine(String sample, Region region) {
        return join("\t",
                sample,
                region.gene,
                region.chr,
                sp,
                ep,
                refallele,
                varallele
        );
    }

    private String outputMiddleOfVariantLine(Variant variant) {
        return join("\t",
                variant.shift3,
                variant.msi == 0 ? 0 : new DecimalFormat("0.000").format(variant.msi),
                variant.msint,
                variant.leftseq.isEmpty() ? "0" : leftseq,
                variant.rightseq.isEmpty() ? "0" : rightseq
        );
    }

    private String outputEndOfVariantLine(Region region, String type) {
        return join("\t",
                region.chr + ":" + region.start + "-" + region.end,
                type,
                vartype
        );
    }

    private String outputEndOfVariantDuprate(String v1sv, String v2sv, Variant variant2) {
        return join("\t",
                duprate == 0 ? 0 : new DecimalFormat("0.0").format(duprate),
                v1sv.equals("") ? 0 : v1sv,
                variant2.duprate == 0 ? 0 : new DecimalFormat("0.0").format(variant2.duprate),
                v2sv.equals("") ? 0 : v2sv
        );
    }

    public String joinVar2(String delm) {
        return join(delm,
                tcov,
                cov,
                rfc,
                rrc,
                fwd,
                rev,
                genotype == null ? "0" : genotype,
                freq == 0 ? 0 : new DecimalFormat("0.0000").format(freq),
                bias,
                new DecimalFormat("0.0").format(pmean),
                pstd ? 1 : 0,
                new DecimalFormat("0.0").format(qual),
                qstd ? 1 : 0,
                new DecimalFormat("0.0").format(mapq),
                new DecimalFormat("0.000").format(qratio),
                hifreq == 0 ? 0 : new DecimalFormat("0.0000").format(hifreq),
                extrafreq == 0 ? 0 : new DecimalFormat("0.0000").format(extrafreq)
        );
    }

    public static String joinVariantWithNM(Variant variant) {
        return join("\t",
                variant.joinVar2("\t"),
                variant.nm > 0 ? new DecimalFormat("0.0").format(variant.nm) : 0
        );
    }

    public static String joinEmptyVariantWithTcov(int tcov) {
        return join("\t", tcov, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
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
                + ":" + (fwd + rev)
                + ":F-" + fwd
                + ":R-" + rev
                + ":" + new DecimalFormat("0.0000").format(freq)
                + ":" + bias
                + ":" + new DecimalFormat("0.0").format(pmean)
                + ":" + (pstd ? "1" : "0")
                + ":" + new DecimalFormat("0.0").format(qual)
                + ":" + (qstd ? "1" : "0")
                + ":" + new DecimalFormat("0.0000").format(hifreq)
                + ":" + mapq
                + ":" + new DecimalFormat("0.000").format(qratio));
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
