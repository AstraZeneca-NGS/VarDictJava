package com.astrazeneca.vardict.modules;

import com.astrazeneca.vardict.Configuration;
import com.astrazeneca.vardict.Utils;
import com.astrazeneca.vardict.collection.Tuple;
import com.astrazeneca.vardict.collection.VariationMap;
import com.astrazeneca.vardict.data.Reference;
import com.astrazeneca.vardict.data.Region;
import com.astrazeneca.vardict.variations.Sclip;
import com.astrazeneca.vardict.variations.Variant;
import com.astrazeneca.vardict.variations.Variation;
import com.astrazeneca.vardict.variations.Vars;

import java.io.IOException;
import java.time.LocalDateTime;
import java.util.regex.Pattern;
import java.util.regex.Matcher;
import java.util.*;

import static com.astrazeneca.vardict.Utils.*;
import static com.astrazeneca.vardict.collection.Tuple.tuple;
import static com.astrazeneca.vardict.data.Patterns.*;
import static com.astrazeneca.vardict.modules.SAMFileParser.parseSAM;
import static com.astrazeneca.vardict.variations.VariationUtils.*;

public class ToVarsBuilder {

    /**
     * The current processing segment (chr, start, end)
     */
    public static final ThreadLocal<Tuple.Tuple3<String, Integer, Integer>> CURSEG =
            ThreadLocal.withInitial(()-> new Tuple.Tuple3<>("", 0, 0));

    /**
     * Read BAM files and create variant structure
     * @param region region of interest
     * @param bam BAM file names (':' delimiter)
     * @param reference Reference object contains part of reference sequence (key - position, value - base)
     * @param chrs map of chromosome lengths
     * @param sample sample name
     * @param SPLICE set of strings representing spliced regions
     * @param ampliconBasedCalling string of maximum_distance:minimum_overlap for amplicon based calling
     * @param Rlen maximum read length
     * @param conf VarDict configuration
     * @return Tuple of (maximum read length, variant structure)
     * @throws IOException
     */
    public static Tuple.Tuple2<Integer, Map<Integer, Vars>> toVars(Region region,
                                                            String bam,
                                                            Reference reference,
                                                            Map<String, Integer> chrs,
                                                            String sample,
                                                            Set<String> SPLICE,
                                                            String ampliconBasedCalling,
                                                            int Rlen,
                                                            Configuration conf) throws IOException {
        Map<Integer, VariationMap<String, Variation>> hash = new HashMap<>(); // variations map
        Map<Integer, VariationMap<String, Variation>> iHash = new HashMap<>(); // insertions map
        Map<Integer, Integer> cov = new HashMap<>();
        Map<Integer, Sclip> sclip3 = new HashMap<>(); // soft clipped at 3'
        Map<Integer, Sclip> sclip5 = new HashMap<>(); // soft clipped at 5'
        boolean svflag = false;
        CURSEG.set(new Tuple.Tuple3<>(region.chr, region.start, region.end));

        Tuple.Tuple5<Map<Integer, VariationMap<String, Variation>>, Map<Integer, VariationMap<String, Variation>>,
                Map<Integer, Integer>, Integer, Double> parseTpl =
                parseSAM(region, bam, chrs, sample, SPLICE, ampliconBasedCalling, Rlen, reference,
                        conf, hash,
                        iHash, cov, sclip3, sclip5, svflag);

        hash = parseTpl._1;
        iHash = parseTpl._2;
        cov = parseTpl._3;
        Rlen = parseTpl._4;
        double duprate = parseTpl._5;
        Map<Integer, Character> ref = reference.referenceSequences;

        if (conf.y) {
            System.err.printf("Current segment: %s:%d-%d \n", region.chr, region.start, region.end);
        }

        //the variant structure
        Map<Integer, Vars> vars = new HashMap<>();
        //Loop over positions
        for (Map.Entry<Integer, VariationMap<String, Variation>> entH : hash.entrySet()) {
            final int p = entH.getKey();
            VariationMap<String, Variation> v = entH.getValue();

            if (v.isEmpty()) {
                continue;
            }

            //Skip if there are no structural variants on position or if the delete duplication option is on
            if (v.sv == null || conf.deleteDuplicateVariants) {
                //Skip if start position is outside region of interest
                if (p < region.start || p > region.end) {
                    continue;
                }
            }

            //Skip position if it has no coverage (except SVs)
            if (v.sv == null && !cov.containsKey(p)) {
                continue;
            }

            Set<String> vk = new HashSet<String>(v.keySet());
            if (iHash.containsKey(p)) {
                vk.add("I");
            }
            if (vk.size() == 1 && ref.containsKey(p) && vk.contains(ref.get(p).toString())) {
                // ignore if only reference were seen and no pileup to avoid computation
                if (!conf.doPileup && !conf.bam.hasBam2() && conf.ampliconBasedCalling == null) {
                    continue;
                }
            }

            if (!cov.containsKey(p) || cov.get(p) == 0) { // ignore when there's no coverage
                System.err.printf("Error tcov: %s %d %d %d %s\n",
                        region.chr, p, region.start, region.end, v.sv.type);
                continue;
            }

            //total position coverage
            int tcov = cov.get(p);

            //array of all variants for the position
            List<Variant> var = new ArrayList<>();
            //temporary array used for debugging
            List<String> tmp = new ArrayList<>();
            List<String> keys = new ArrayList<>(v.keySet());
            Collections.sort(keys);

            //position coverage by high-quality reads
            final int hicov = calcHicov(iHash.get(p), v);

            //Loop over all variants found for the position except insertions
            for (String n : keys) {
                if (n.equals("SV")) {
                    VariationMap.SV sv = v.sv;
                    getOrPutVars(vars, p).sv = sv.splits + "-" + sv.pairs + "-" + sv.clusters;
                    continue;
                }
                Variation cnt = v.get(n);
                if (cnt.cnt == 0) { //Skip variant if it does not have count
                    continue;
                }
                //count of variants in forward strand
                int fwd = cnt.getDir(false);
                //count of variants in reverse strand
                int rev = cnt.getDir(true);
                //strand bias flag (0, 1 or 2)
                int bias = strandBias(fwd, rev, conf.bias, conf.minb);
                //mean base quality for variant
                double vqual = roundHalfEven("0.0", cnt.qmean / cnt.cnt); // base quality
                //mean mapping quality for variant
                double mq = roundHalfEven("0.0", cnt.Qmean / (double) cnt.cnt);
                //number of high-quality reads for variant
                int hicnt = cnt.hicnt;
                //number of low-quality reads for variant
                int locnt = cnt.locnt;
                /**
                 * Condition:
                 # 1). cnt.cnt > tcov                         - variant count is more than position coverage
                 # 2). cnt.cnt - tcov < cnt.extracnt          - variant count is no more than position coverage + extracnt
                 */
                int ttcov = tcov;
                if (cnt.cnt > tcov && cnt.extracnt > 0 && cnt.cnt - tcov < cnt.extracnt) { //adjust position coverage if condition holds
                    ttcov = cnt.cnt;
                }

                //create variant record
                Variant tvref = new Variant();
                tvref.n = n;
                tvref.cov = cnt.cnt;
                tvref.fwd = fwd;
                tvref.rev = rev;
                tvref.bias = String.valueOf(bias);
                tvref.freq = cnt.cnt / (double) ttcov;
                tvref.pmean = cnt.pmean / (double) cnt.cnt;
                tvref.pstd = cnt.pstd;
                tvref.qual = vqual;
                tvref.qstd = cnt.qstd;
                tvref.mapq = mq;
                tvref.qratio = hicnt / (locnt != 0 ? locnt : 0.5d);
                tvref.hifreq = hicov > 0 ? hicnt / (double) hicov : 0;
                tvref.extrafreq = cnt.extracnt != 0 ? cnt.extracnt / (double) ttcov : 0;
                tvref.shift3 = 0;
                tvref.msi = 0;
                tvref.nm = Utils.roundHalfEven("0.0", cnt.nm / (double) cnt.cnt);
                tvref.hicnt = hicnt;
                tvref.hicov = hicov;
                tvref.duprate = duprate;

                //append variant record
                var.add(tvref);
                if (conf.debug) {
                    tvref.debugVariantsContentSimple(tmp, n);
                }
            }
            //Handle insertions separately
            Map<String, Variation> iv = iHash.get(p);
            if (iv != null) {
                List<String> ikeys = new ArrayList<>(iv.keySet());
                Collections.sort(ikeys);
                //Loop over insertion variants
                for (String n : ikeys) {
                    // String n = entV.getKey();
                    Variation cnt = iv.get(n);
                    //count of variants in forward strand
                    int fwd = cnt.getDir(false);
                    //count of variants in reverse strand
                    int rev = cnt.getDir(true);
                    //strand bias flag (0, 1 or 2)
                    int bias = strandBias(fwd, rev, conf.bias, conf.minb);
                    //mean base quality for variant
                    double vqual = roundHalfEven("0.0", cnt.qmean / cnt.cnt); // base quality
                    //mean mapping quality for variant
                    double mq = roundHalfEven("0.0", cnt.Qmean / (double)cnt.cnt); // mapping quality
                    //number of high-quality reads for variant
                    int hicnt = cnt.hicnt;
                    //number of low-quality reads for variant
                    int locnt = cnt.locnt;

                    // Also commented in Perl
                    // hicov += hicnt;

                    //adjust position coverage if variant count is more than position coverage and no more than
                    // position coverage + extracnt
                    int ttcov = tcov;
                    if (cnt.cnt > tcov && cnt.extracnt != 0 &&cnt.cnt - tcov < cnt.extracnt) {
                        ttcov = cnt.cnt;
                    }
                    if (ttcov < cnt.cnt) {
                        ttcov = cnt.cnt;
                        if (cov.containsKey(p + 1) && ttcov < cov.get(p + 1) - cnt.cnt) {
                            ttcov = cov.get(p + 1);
                            // Adjust the reference
                            getVariationMaybe(hash, p + 1, ref.get(p + 1)).dirPlus -= fwd;
                            getVariationMaybe(hash, p + 1, ref.get(p + 1)).dirMinus -= rev;
                        }
                        tcov = ttcov;
                    }

                    Variant tvref = new Variant();
                    tvref.n = n;
                    tvref.cov = cnt.cnt;
                    tvref.fwd = fwd;
                    tvref.rev = rev;
                    tvref.bias = String.valueOf(bias);
                    tvref.freq = cnt.cnt / (double) ttcov;
                    tvref.pmean = cnt.pmean / (double) cnt.cnt;
                    tvref.pstd = cnt.pstd;
                    tvref.qual = vqual;
                    tvref.qstd = cnt.qstd;
                    tvref.mapq = mq;
                    tvref.qratio = hicnt / (locnt != 0 ? locnt : 0.5d);
                    tvref.hifreq = hicov > 0 ? hicnt / (double)hicov : 0;
                    tvref.extrafreq = cnt.extracnt != 0 ? cnt.extracnt / (double)ttcov : 0;
                    tvref.shift3 = 0;
                    tvref.msi = 0;
                    tvref.nm = Utils.roundHalfEven("0.0", cnt.nm / (double)cnt.cnt);
                    tvref.hicnt = hicnt;
                    tvref.hicov = hicov;
                    tvref.duprate = duprate;

                    var.add(tvref);
                    if (conf.debug) {
                        tvref.debugVariantsContentInsertion(tmp, n);
                    }
                }
            }

            //sort variants by product of quality and coverage
            Collections.sort(var, new Comparator<Variant>() {
                @Override
                public int compare(Variant o1, Variant o2) {
                    return Double.compare(o2.qual * o2.cov, o1.qual * o1.cov);
                }
            });
            double maxfreq = 0;
            for (Variant tvar : var) {
                //If variant description string is 1-char base and it matches reference base at this position
                if (tvar.n.equals(String.valueOf(ref.get(p)))) {
                    //this is a reference variant
                    getOrPutVars(vars, p).ref = tvar;
                } else {
                    //append variant to VAR and put it to VARN with key tvar.n (variant description string)
                    getOrPutVars(vars, p).var.add(tvar);
                    getOrPutVars(vars, p).varn.put(tvar.n, tvar);
                    if (tvar.freq > maxfreq) {
                        maxfreq = tvar.freq;
                    }
                }
            }
            if (!conf.doPileup && maxfreq <= conf.freq && conf.ampliconBasedCalling == null) {
                if (!conf.bam.hasBam2()) {
                    vars.remove(p);
                    continue;
                }
            }
            // Make sure the first bias is always for the reference nucleotide

            //reference forward strand coverage
            int rfc = 0;
            //reference reverse strand coverage
            int rrc = 0;
            //description string for reference or best non-reference variant
            String genotype1 = "";

            //if reference variant has frequency more than $FREQ, set genotype1 to reference variant
            //else if non-reference variant exists, set genotype1 to first (largest quality*coverage) variant
            //otherwise set genotype1 to reference variant
            Vars varsAtp = getOrPutVars(vars, p);
            if (varsAtp.ref != null && varsAtp.ref.freq >= conf.freq) {
                genotype1 = varsAtp.ref.n;
            } else if (varsAtp.var.size() > 0) {
                genotype1 = varsAtp.var.get(0).n;
            } else {
                genotype1 = varsAtp.ref.n;
            }
            if (varsAtp.ref != null) {
                rfc = varsAtp.ref.fwd;
                rrc = varsAtp.ref.rev;
            }

            if (genotype1.startsWith("+")) {
                Matcher mm = DUP_NUM.matcher(genotype1);
                if (mm.find()) {
                    genotype1 = "+" + (Configuration.SVFLANK + toInt(mm.group(1)));
                }
                else {
                    genotype1 = "+" + (genotype1.length() - 1);
                }
            }
            //description string for any other variant
            String genotype2 = "";

            if (tcov > cov.get(p) && hash.containsKey(p + 1) && ref.containsKey(p + 1)
                    && hash.get(p + 1).containsKey(ref.get(p + 1).toString())) {
                Variation tpref = getVariationMaybe(hash, p + 1, ref.get(p + 1));
                rfc = tpref.dirPlus;
                rrc = tpref.dirMinus;
            }

            // only reference reads are observed.
            if (varsAtp.var.size() > 0) { //Condition: non-reference variants are found
                //Loop over non-reference variants
                for (Variant vref : varsAtp.var) {
                    //vref - variant reference

                    genotype2 = vref.n;
                    if (genotype2.startsWith("+")) {
                        genotype2 = "+" + (genotype2.length() - 1);
                    }
                    //variant description string
                    final String vn = vref.n;
                    //length of deletion in variant (0 if no deletion)
                    int dellen = 0;
                    Matcher matcher = BEGIN_MINUS_NUMBER.matcher(vn);
                    if (matcher.find()) {
                        dellen = toInt(matcher.group(1));
                    }
                    //effective position (??): p + dellen - 1 for deletion, p otherwise
                    int ep = p;
                    if (vn.startsWith("-")) {
                        ep = p + dellen - 1;
                    }
                    //reference sequence for variation (to be written to .vcf file)
                    String refallele = "";
                    //variant sequence (to be written to .vcf file)
                    String varallele = "";

                    // how many bp can a deletion be shifted to 3 prime
                    //3' shift (integer) for MSI adjustment
                    int shift3 = 0;
                    double msi = 0;
                    String msint = "";

                    int sp = p;

                    //if variant is an insertion
                    if (vn.startsWith("+")) {
                        //If no '&' and '#' symbols are found in variant string
                        //These symbols are in variant if a matched sequence follows insertion
                        if (!vn.contains("&") && !vn.contains("#") && !vn.contains("<dup")) {
                            //variant description string without first symbol '+'
                            String tseq1 = vn.substring(1);
                            //left 50 bases in reference sequence
                            String leftseq = joinRef(ref, p - 50 > 1 ? p - 50 : 1, p); // left 10 nt
                            int x = getOrElse(chrs, region.chr, 0);
                            //right 70 bases in reference sequence
                            String tseq2 = joinRef(ref, p + 1, (p + 70 > x ? x : p + 70));

                            Tuple.Tuple3<Double, Integer, String> tpl = findMSI(tseq1, tseq2, leftseq);
                            msi = tpl._1;
                            shift3 = tpl._2;
                            msint = tpl._3;

                            //Try to adjust for microsatellite instability
                            tpl = findMSI(leftseq, tseq2, null);
                            double tmsi = tpl._1;
                            String tmsint = tpl._3;
                            if (msi < tmsi) {
                                msi = tmsi;
                                // Don't change shift3
                                msint = tmsint;
                            }
                            if (msi <= shift3 / (double)tseq1.length()) {
                                msi = shift3 / (double)tseq1.length();
                            }
                        }

                        //Shift position to 3' if -3 option is set
                        if (conf.moveIndelsTo3) {
                            sp += shift3;
                            ep += shift3;
                        }
                        //reference allele is 1 base
                        refallele = ref.containsKey(p) ? ref.get(p).toString() : "";
                        //variant allele is reference base concatenated with insertion
                        varallele = refallele + vn.substring(1);
                        if (varallele.length() > conf.SVMINLEN) {
                            ep += varallele.length();
                            varallele = "<DUP>";
                        }
                        Matcher mm = DUP_NUM.matcher(varallele);
                        if (mm.find()) {
                            int dupCount = toInt(mm.group(1));
                            ep = sp + (2 * Configuration.SVFLANK + dupCount) - 1;
                            genotype2 = "+" + (2 * Configuration.SVFLANK + dupCount);
                            varallele = "<DUP>";
                        }
                    } else if (vn.startsWith("-")) { //deletion variant
                        Matcher matcherINV = INV_NUM.matcher(vn);
                        Matcher matcherStartMinusNum = BEGIN_MINUS_NUMBER_CARET.matcher(vn);

                        if (dellen < conf.SVMINLEN) {
                            //variant allele is in the record
                            //remove '-' and number from beginning of variant string
                            varallele = vn.replaceFirst("^-\\d+", "");

                            //left 70 bases in reference sequence
                            String leftseq = joinRef(ref, (p - 70 > 1 ? p - 70 : 1), p - 1); // left 10 nt
                            int chr0 = getOrElse(chrs, region.chr, 0);
                            //right 70 + dellen bases in reference sequence
                            String tseq = joinRef(ref, p, p + dellen + 70 > chr0 ? chr0 : p + dellen + 70);

                            //Try to adjust for microsatellite instability
                            Tuple.Tuple3<Double, Integer, String> tpl = findMSI(substr(tseq, 0, dellen),
                                    substr(tseq, dellen), leftseq);
                            msi = tpl._1;
                            shift3 = tpl._2;
                            msint = tpl._3;

                            tpl = findMSI(leftseq, substr(tseq, dellen), null);
                            double tmsi = tpl._1;
                            String tmsint = tpl._3;
                            if (msi < tmsi) {
                                msi = tmsi;
                                // Don't change shift3
                                msint = tmsint;
                            }
                            if (msi <= shift3 / (double) dellen) {
                                msi = shift3 / (double) dellen;
                            }
                            if (matcherINV.find()) {
                                varallele = "<INV>";
                                genotype2 = "<INV" + dellen + ">";
                            }
                        } else if (matcherStartMinusNum.find()) {
                            varallele = "<INV>";
                            genotype2 = "<INV" + dellen + ">";
                        } else {
                            varallele = "<DEL>";
                        }
                        //If no matched sequence or indel follows the variant
                        if (!vn.contains("&") && !vn.contains("#") && !vn.contains("^")) {
                            //Shift position to 3' if -3 option is set
                            if (conf.moveIndelsTo3) {
                                sp += shift3;
                            }
                            //variant allele is 1 base from reference string preceding p
                            if (!varallele.equals("<DEL>")) {
                                varallele = ref.containsKey(p - 1) ? ref.get(p - 1).toString() : "";
                            }
                            //prepend same base to reference allele
                            refallele = ref.containsKey(p - 1) ? ref.get(p - 1).toString() : "";
                            sp--;
                        }
                        Matcher mm = SOME_SV_NUMBERS.matcher(vn);
                        if (mm.find()) {
                            refallele = ref.containsKey(p) ? ref.get(p).toString() : "";
                        }
                        else if (dellen < conf.SVMINLEN) {
                            //append dellen bases from reference string to reference allele
                            refallele += joinRef(ref, p, p + dellen - 1);
                        }
                    } else { //Not insertion/deletion variant. SNP or MNP
                        //Find MSI adjustment
                        String tseq1 = joinRef(ref, p - 30 > 1 ? p - 30 : 1, p + 1);
                        int chr0 = getOrElse(chrs, region.chr, 0);
                        String tseq2 = joinRef(ref, p + 2, p + 70 > chr0 ? chr0 : p + 70);

                        Tuple.Tuple3<Double, Integer, String> tpl = findMSI(tseq1, tseq2, null);
                        msi = tpl._1;
                        shift3 = tpl._2;
                        msint = tpl._3;
                        //reference allele is 1 base from reference sequence
                        refallele = ref.containsKey(p) ? ref.get(p).toString() : "";
                        //variant allele is same as description string
                        varallele = vn;
                    }

                    Matcher mtch = AMP_ATGC.matcher(vn);
                    if (mtch.find()) { //If variant is followed by matched sequence
                        //following matching sequence
                        String extra = mtch.group(1);
                        //remove '&' symbol from variant allele
                        varallele = varallele.replaceFirst("&", "");
                        //append length(extra) bases from reference sequence to reference allele and genotype1
                        String tch = joinRef(ref, ep + 1, ep + extra.length());
                        refallele += tch;
                        genotype1 += tch;

                        //Adjust position
                        ep += extra.length();

                        mtch = AMP_ATGC.matcher(varallele);
                        if (mtch.find()) {
                            String vextra = mtch.group(1);
                            varallele = varallele.replaceFirst("&", "");
                            tch = joinRef(ref, ep + 1, ep + vextra.length());
                            refallele += tch;
                            genotype1 += tch;
                            ep += vextra.length();
                        }

                        //If description string starts with '+' sign, remove it from reference and variant alleles
                        if (vn.startsWith("+")) {
                            refallele = refallele.substring(1);
                            varallele = varallele.substring(1);
                            sp++;
                        }

                        if (varallele.equals("<DEL>") && refallele.length() > 1) {
                            refallele = ref.containsKey(sp) ? ref.get(sp).toString() : "";
                            if (cov.containsKey(sp - 1)) {
                                tcov = cov.get(sp - 1);
                            }
                            if (vref.cov > tcov ){
                                tcov = vref.cov;
                            }
                        }
                    }

                    //If variant is followed by short matched sequence and insertion/deletion
                    mtch = HASH_GROUP_CARET_GROUP.matcher(vn);
                    if (mtch.find()) {
                        //matched sequence
                        String mseq = mtch.group(1);
                        //insertion/deletion tail
                        String tail = mtch.group(2);

                        //adjust position by length of matched sequence
                        ep += mseq.length();

                        //append bases from reference sequence to reference allele
                        refallele += joinRef(ref, ep - mseq.length() + 1, ep);

                        //If tail is a deletion
                        mtch = BEGIN_DIGITS.matcher(tail);
                        if (mtch.find()) {
                            //append (deletion length) bases from reference sequence to reference allele
                            int d = toInt(mtch.group(1));
                            refallele += joinRef(ref, ep + 1, ep + d);

                            //shift position by deletion length
                            ep += d;
                        }

                        //clean special symbols from alleles
                        varallele = varallele.replaceFirst("#", "").replaceFirst("\\^(\\d+)?", "");

                        //replace '#' with 'm' and '^' with 'i' in genotypes
                        genotype1 = genotype1.replaceFirst("#", "m").replaceFirst("\\^", "i");
                        genotype2 = genotype2.replaceFirst("#", "m").replaceFirst("\\^", "i");
                    }
                    mtch = CARET_ATGNC.matcher(vn); // for deletion followed directly by insertion in novolign
                    if (mtch.find()) {
                        //remove '^' sign from varallele
                        varallele = varallele.replaceFirst("\\^", "");

                        //replace '^' sign with 'i' in genotypes
                        genotype1 = genotype1.replaceFirst("\\^", "i");
                        genotype2 = genotype2.replaceFirst("\\^", "i");
                    }

                    //preceding reference sequence
                    vref.leftseq = joinRef(ref, sp - 20 < 1 ? 1 : sp - 20, sp - 1); // left 20 nt
                    int chr0 = getOrElse(chrs, region.chr, 0);
                    //following reference sequence
                    vref.rightseq = joinRef(ref, ep + 1, ep + 20 > chr0 ? chr0 : ep + 20); // right 20 nt
                    //genotype description string
                    String genotype = genotype1 + "/" + genotype2;
                    //remove '&' and '#' symbols from genotype string
                    //replace '^' symbol with 'i' in genotype string
                    genotype = genotype
                            .replace("&", "")
                            .replace("#", "")
                            .replace("^", "i");
                    //convert extrafreq, freq, hifreq, msi fields to strings
                    vref.extrafreq = roundHalfEven("0.0000", vref.extrafreq);
                    vref.freq = roundHalfEven("0.0000", vref.freq);
                    vref.hifreq = roundHalfEven("0.0000", vref.hifreq);
                    vref.msi = roundHalfEven("0.000", msi);
                    vref.msint = msint.length();
                    vref.shift3 = shift3;
                    vref.sp = sp;
                    vref.ep = ep;
                    vref.refallele = refallele;
                    vref.varallele = varallele;
                    vref.genotype = genotype;
                    vref.tcov = tcov;
                    vref.rfc = rfc;
                    vref.rrc = rrc;

                    //bias is [0-2];[0-2] where first flag is for reference, second for variant
                    //if reference variant is not found, first flag is 0
                    if (varsAtp.ref != null) {
                        vref.bias = varsAtp.ref.bias + ";" + vref.bias;
                    } else {
                        vref.bias = "0;" + vref.bias;
                    }

                    adjustVariantCounts(p, vref);

                    if (conf.debug) {
                        StringBuilder sb = new StringBuilder();
                        for (String str : tmp) {
                            if (sb.length() > 0) {
                                sb.append(" & ");
                            }
                            sb.append(str);
                        }
                        vref.DEBUG = sb.toString();
                    }
                }
                //TODO: It is a "lazy" solution because current logic in realignment methods can't be changed simply for --nosv option
                if (conf.disableSV) {
                    varsAtp.var.removeIf(vref -> ANY_SV.matcher(vref.varallele).find());
                }
            } else if (varsAtp.ref != null) {
                Variant vref = varsAtp.ref; //no variant reads are detected.
                vref.tcov = tcov;
                vref.cov = 0;
                vref.freq = 0;
                vref.rfc = rfc;
                vref.rrc = rrc;
                vref.fwd = 0;
                vref.rev = 0;
                vref.msi = 0;
                vref.msint = 0;
                vref.bias += ";0";
                vref.shift3 = 0;
                vref.sp = p;
                vref.ep = p;
                vref.hifreq = roundHalfEven("0.0000", vref.hifreq);
                String r = ref.containsKey(p) ? ref.get(p).toString() : "";
                //both refallele and varallele are 1 base from reference string
                vref.refallele = r;
                vref.varallele = r;
                vref.genotype = r + "/" + r;
                vref.leftseq = "";
                vref.rightseq = "";
                vref.duprate = duprate;
                if (conf.debug) {
                    StringBuilder sb = new StringBuilder();
                    for (String str : tmp) {
                        if (sb.length() > 0) {
                            sb.append(" & ");
                        }
                        sb.append(str);
                    }
                    vref.DEBUG = sb.toString();
                }
            } else {
                varsAtp.ref = new Variant();
            }
        }
        if (conf.y) {
            System.err.println("TIME: Finish preparing vars:" + LocalDateTime.now());
        }
        return tuple(Rlen, vars);
    }

    /**
     * Adjust variant negative counts of fields FWD, REV, RFC, RRC to zeros and print the information message to console
     * @param p start position of variant
     * @param vref variant to adjust
     */
    private static void adjustVariantCounts(int p, Variant vref) {
        String message = "column in variant on position: " + p + " " + vref.refallele + "->" +
                vref.varallele + " was negative, adjusted to zero.";

        if (vref.rfc < 0 ) {
            vref.rfc = 0;
            System.err.println("Reference forward count " + message);
        }
        if (vref.rrc < 0) {
            vref.rrc = 0;
            System.err.println("Reference reverse count " + message);
        }
        if (vref.fwd < 0) {
            vref.fwd = 0;
            System.err.println("Variant forward count " + message);
        }
        if (vref.rev < 0 ) {
            vref.rev = 0;
            System.err.println("Variant reverse count " + message);
        }
    }

    private static int calcHicov(VariationMap<String, Variation> iv,
                                 VariationMap<String, Variation> v) {
        int hicov = 0;
        for (Map.Entry<String, Variation> vr : v.entrySet()) {
            if (vr.getKey().equals("SV")) {
                continue;
            }
            hicov += vr.getValue().hicnt;
        }
        if (iv != null) {
            for (Variation vr : iv.values()) {
                hicov += vr.hicnt;
            }
        }
        return hicov;
    }

    /**
     * Find microsatellite instability
     * Tandemly repeated short sequence motifs ranging from 1– 6(8 in our case) base pairs are called microsatellites.
     * Other frequently used terms for these DNA regions are simple sequences or short tandem repeats (STRs)
     * @param tseq1 variant description string
     * @param tseq2 right 70 bases in reference sequence
     * @param left left 50 bases in reference sequence
     * @return Tuple of (MSI count, No. of bases to be shifted to 3 prime for deletions due to alternative alignment,
     * MicroSattelite unit length in base pairs)
     */
    private static Tuple.Tuple3<Double, Integer, String> findMSI(String tseq1, String tseq2, String left) {

        //Number of nucleotides in microsattelite
        int nmsi = 1;
        //Number of bases to be shifted to 3 prime
        int shift3 = 0;
        String maxmsi = "";
        double msicnt = 0;
        while (nmsi <= tseq1.length() && nmsi <= 6) {
            //Microsattelite nucleotide sequence; trim nucleotide(s) from the end
            String msint = substr(tseq1, -nmsi, nmsi);
            Pattern pattern = Pattern.compile("((" + msint + ")+)$");
            Matcher mtch = pattern.matcher(tseq1);
            String msimatch = "";
            if (mtch.find()) {
                msimatch = mtch.group(1);
            }
            if (left != null && !left.isEmpty()) {
                mtch = pattern.matcher(left + tseq1);
                if (mtch.find()) {
                    msimatch = mtch.group(1);
                }
            }
            double curmsi = msimatch.length() / (double)nmsi;
            mtch = Pattern.compile("^((" + msint + ")+)").matcher(tseq2);
            if (mtch.find()) {
                curmsi += mtch.group(1).length() / (double)nmsi;
            }
            if (curmsi > msicnt) {
                maxmsi = msint;
                msicnt = curmsi;
            }
            nmsi++;
        }

        String tseq = tseq1 + tseq2;
        while (shift3 < tseq2.length() && tseq.charAt(shift3) == tseq2.charAt(shift3)) {
            shift3++;
        }
        return tuple(msicnt, shift3, maxmsi);
    }

}
