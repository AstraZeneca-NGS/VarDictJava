package com.astrazeneca.vardict;

import static com.astrazeneca.vardict.Tuple.tuple;
import static com.astrazeneca.vardict.Utils.*;
import static com.astrazeneca.vardict.VarDict.VarsType.ref;
import static com.astrazeneca.vardict.VarDict.VarsType.var;
import static com.astrazeneca.vardict.VarDict.VarsType.varn;
import static java.lang.Math.abs;
import static java.lang.String.format;
import static java.util.Collections.singletonList;
import static java.util.Collections.singletonMap;
import htsjdk.samtools.*;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;

import java.io.*;
import java.text.DecimalFormat;
import java.util.*;
import java.util.Map.Entry;
import java.util.concurrent.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import jregex.Replacer;

import com.astrazeneca.vardict.Tuple.Tuple2;
import com.astrazeneca.vardict.Tuple.Tuple3;
import com.astrazeneca.vardict.Tuple.Tuple4;

public class VarDict {

    final static Random RND = new Random(System.currentTimeMillis());

    final static Pattern SN = Pattern.compile("\\s+SN:(\\S+)");
    final static Pattern LN = Pattern.compile("\\sLN:(\\d+)");

    public static void start(Configuration conf) throws IOException {
        if (conf.printHeader) {
            System.out.println(join("\t",
                    "Sample", "Gene", "Chr", "Start", "End", "Ref", "Alt", "Depth", "AltDepth", "RefFwdReads",
                    "RefRevReads", "AltFwdReads", "AltRevReads", "Genotype", "AF", "Bias", "PMean", "PStd",
                    "QMean", "QStd", "5pFlankSeq", "3pFlankSeq"));
        }

        Tuple2<String, String> stpl = getSampleNames(conf.bam.getBamRaw(), conf.sampleName, conf.sampleNameRegexp);
        String sample = stpl._1;
        String samplem = stpl._2;

        Map<String, Integer> chrs = readChr(conf.bam.getBamX());

        if (conf.regionOfInterest != null) {
            Region region = buildRegion(conf.regionOfInterest, conf.numberNucleotideToExtend, chrs, conf.isZeroBasedDefined() ? conf.zeroBased : false);
            nonAmpVardict(singletonList(singletonList(region)), chrs, conf.ampliconBasedCalling, sample, samplem, conf);
        } else {
            Tuple3<String, Boolean, List<String>> tpl = readBedFile(conf);
            String ampliconBasedCalling = tpl._1;
            Boolean zeroBased = tpl._2;
            List<String> segraw = tpl._3;

            if (ampliconBasedCalling != null) {
                List<List<Region>> segs = toRegions(segraw, chrs, zeroBased != null ? zeroBased : false, conf.delimiter);
                if (conf.threads == 1)
                    ampVardictNotParallel(segs, chrs, ampliconBasedCalling, conf.bam.getBam1(), sample, conf);
                else
                    ampVardictParallel(segs, chrs, ampliconBasedCalling, conf.bam.getBam1(), sample, conf);
            } else {
                List<List<Region>> regions = toRegions(segraw, chrs, zeroBased, conf);
                nonAmpVardict(regions, chrs, null, sample, samplem, conf);
            }
        }
    }

    /**
     * Read map of chromosome lengths
     * @param bam BAM file name
     * @return Map of chromosome lengths. Key - chromosome name, value - length
     * @throws IOException
     */
    static Map<String, Integer> readChr(String bam) throws IOException {
        try (SamReader reader = SamReaderFactory.makeDefault().open(new File(bam))) {
            SAMFileHeader header = reader.getFileHeader();
            Map<String, Integer> chrs = new HashMap<>();
            for (SAMSequenceRecord record : header.getSequenceDictionary().getSequences()) {
                record.getSequenceLength();
                String sn = record.getSequenceName();
                int ln = record.getSequenceLength();
                chrs.put(sn, ln);
            }
            return chrs;
        }
    }

    private final static Comparator<Region> ISTART_COMPARATOR = new Comparator<Region>() {
        @Override
        public int compare(Region o1, Region o2) {
            return Integer.compare(o1.istart, o2.istart);
        }
    };

    private static List<List<Region>> toRegions(List<String> segraw, Map<String, Integer> chrs, Boolean zeroBased, Configuration conf) throws IOException {
        boolean zb = conf.isZeroBasedDefined() ? conf.zeroBased : false;
        List<List<Region>> segs = new LinkedList<>();
        BedRowFormat format = conf.bedRowFormat;
        for (String seg : segraw) {
            String[] splitA = seg.split(conf.delimiter);
            if (!conf.isColumnForChromosomeSet() && splitA.length == 4) {
                try {
                    int a1 = toInt(splitA[1]);
                    int a2 = toInt(splitA[2]);
                    if (a1 <= a2) {
                        format = CUSTOM_BED_ROW_FORMAT;
                        if (zeroBased == null) {
                            zb = true;
                        }
                    }
                } catch (NumberFormatException e) {
                }
            }
            String chr = splitA[format.chrColumn];
            chr = correctChr(chrs, chr);
            int cdss = toInt(splitA[format.startColumn]);
            int cdse = toInt(splitA[format.endColumn]);
            String gene = format.geneColumn < splitA.length ? splitA[format.geneColumn] : chr;

            String[] starts = splitA[format.thickStartColumn].split(",");
            String[] ends = splitA[format.thickEndColumn].split(",");
            List<Region> cds = new LinkedList<>();
            for (int i = 0; i < starts.length; i++) {
                int s = toInt(starts[i]);
                int e = toInt(ends[i]);
                if (cdss > e) {
                    continue;
                }
                if (cdse > e) {
                    break;
                }
                if (s < cdss)
                    s = cdss;
                if (e > cdse)
                    e = cdse;
                s -= conf.numberNucleotideToExtend;
                e += conf.numberNucleotideToExtend;
                if (zb && s < e) {
                    s++;
                }
                cds.add(new Region(chr, s, e, gene));
            }
            segs.add(cds);
        }
        return segs;
    }

    private static List<List<Region>> toRegions(List<String> segraw, Map<String, Integer> chrs, boolean zeroBased, String delimiter) throws IOException {
        List<List<Region>> segs = new LinkedList<>();
        Map<String, List<Region>> tsegs = new HashMap<>();
        for (String string : segraw) {
            String[] split = string.split(delimiter);
            String chr = split[0];
            chr = correctChr(chrs, chr);
            int start = toInt(split[1]);
            int end = toInt(split[2]);
            String gene = split[3];
            int istart = toInt(split[6]);
            int iend = toInt(split[7]);
            if (zeroBased && start < end) {
                start++;
                istart++;
            }
            List<Region> list = tsegs.get(chr);
            if (list == null) {
                list = new ArrayList<>();
                tsegs.put(chr, list);
            }
            list.add(new Region(chr, start, end, gene, istart, iend));
        }
        List<Region> list = new LinkedList<>();
        segs.add(list);
        int pend = -1;
        for (Map.Entry<String, List<Region>> entry : tsegs.entrySet()) {
            List<Region> regions = entry.getValue();
            Collections.sort(regions, ISTART_COMPARATOR);
            String pchr = null;
            for (Region region : regions) {
                if (pend != -1 && (!region.chr.equals(pchr) || region.istart > pend)) {
                    list = new LinkedList<>();
                    segs.add(list);
                }
                list.add(region);
                pchr = region.chr;
                pend = region.iend;
            }
        }

        return segs;
    }

    private static String correctChr(Map<String, Integer> chrs, String chr) {
        if (!chrs.containsKey(chr)) {
            if (chr.startsWith("chr")) {
                chr = chr.substring("chr".length());
            } else {
                chr = "chr" + chr;
            }
        }
        return chr;
    }

    private static Tuple3<String, Boolean, List<String>> readBedFile(Configuration conf) throws IOException, FileNotFoundException {
        String a = conf.ampliconBasedCalling;
        Boolean zeroBased = conf.zeroBased;
        List<String> segraw = new ArrayList<>();
        try (BufferedReader bedFileReader = new BufferedReader(new FileReader(conf.bed))) {
            String line;
            while ((line = bedFileReader.readLine()) != null) {
                if (line.startsWith("#")
                        || line.startsWith("browser")
                        || line.startsWith("track")) {
                    continue;
                }
                if (a == null) {
                    String[] ampl = line.split(conf.delimiter);
                    if (ampl.length > 7) {
                        try {
                            int a1 = toInt(ampl[1]);
                            int a2 = toInt(ampl[2]);
                            int a6 = toInt(ampl[6]);
                            int a7 = toInt(ampl[7]);
                            if (a6 >= a1 && a7 <= a2) {
                                a = "10:0.95";
                                if (!conf.isZeroBasedDefined()) {
                                    zeroBased = true;
                                }
                            }
                        } catch (NumberFormatException e) {
                            continue;
                        }
                    }
                }
                segraw.add(line);
            }
        }
        return tuple(a, zeroBased, segraw);
    }

    private static void nonAmpVardict(List<List<Region>> segs, Map<String, Integer> chrs, String ampliconBasedCalling, String sample, String samplem, Configuration conf) throws IOException {

        if (conf.bam.hasBam2()) {
            Tuple2<String, String> stpl = getSampleNames(sample, samplem, conf.bam.getBam1(), conf.bam.getBam2(), conf.sampleName, conf.sampleNameRegexp);
            sample = stpl._1;
            samplem = stpl._2;
            if (conf.threads == 1)
                somaticNotParallel(segs, chrs, ampliconBasedCalling, sample, samplem, conf);
            else
                somaticParallel(segs, chrs, ampliconBasedCalling, sample, samplem, conf);
        } else {
            if (conf.threads == 1)
                vardictNotParallel(segs, chrs, ampliconBasedCalling, sample, conf);
            else
                vardictParallel(segs, chrs, ampliconBasedCalling, sample, conf);
        }
    }

    private static void vardictParallel(final List<List<Region>> segs, final Map<String, Integer> chrs, final String ampliconBasedCalling, final String sample, final Configuration conf) throws IOException {
        final ExecutorService executor = Executors.newFixedThreadPool(conf.threads);
        final BlockingQueue<Future<OutputStream>> toPrint = new LinkedBlockingQueue<>(10);
        executor.submit(new Runnable() {

            @Override
            public void run() {
                try {
                    for (List<Region> list : segs) {
                        for (Region region : list) {
                            toPrint.put(executor.submit(new VardictWorker(region, chrs, new HashSet<String>(), ampliconBasedCalling, sample, conf)));
                        }
                    }
                    toPrint.put(NULL_FUTURE);
                } catch (InterruptedException e) {
                    e.printStackTrace();
                }
            }
        });

        try {
            while (true) {
                Future<OutputStream> wrk = toPrint.take();
                if (wrk == NULL_FUTURE) {
                    break;
                }
                System.out.print(wrk.get());
            }
        } catch (InterruptedException | ExecutionException e) {
            e.printStackTrace();
        }

        executor.shutdown();

    }

    private static void vardictNotParallel(final List<List<Region>> segs, final Map<String, Integer> chrs, final String ampliconBasedCalling, final String sample, final Configuration conf) throws IOException {
        for (List<Region> list : segs) {
            for (Region region : list) {
                Map<Integer, Character> ref = getREF(region, chrs, conf.fasta, conf.numberNucleotideToExtend);
                final Set<String> splice = new HashSet<>();
                Tuple2<Integer, Map<Integer, Vars>> tpl = toVars(region, conf.bam.getBam1(), ref, chrs, sample, splice, ampliconBasedCalling, 0, conf);
                vardict(region, tpl._2, sample, splice, conf, System.out);
            }
        }
    }

    private static void somaticParallel(final List<List<Region>> segs, final Map<String, Integer> chrs, final String ampliconBasedCalling, final String sample, String samplem, final Configuration conf) throws IOException {
        final ExecutorService executor = Executors.newFixedThreadPool(conf.threads);
        final BlockingQueue<Future<OutputStream>> toSamdict = new LinkedBlockingQueue<>(10);

        executor.submit(new Runnable() {

            @Override
            public void run() {
                try {
                    for (List<Region> list : segs) {
                        for (Region region : list) {
                            final Set<String> splice = new ConcurrentHashSet<>();
                            Map<Integer, Character> ref = getREF(region, chrs, conf.fasta, conf.numberNucleotideToExtend);
                            Future<Tuple2<Integer, Map<Integer, Vars>>> f1 = executor.submit(new ToVarsWorker(region, conf.bam.getBam1(), chrs, sample, splice, ampliconBasedCalling, ref, conf));
                            Future<OutputStream> f2 = executor.submit(new SomdictWorker(region, conf.bam.getBam2(), chrs, splice, ampliconBasedCalling, ref, conf, f1, sample));
                            toSamdict.put(f2);
                        }
                    }
                    toSamdict.put(NULL_FUTURE);
                } catch (IOException | InterruptedException e) {
                    e.printStackTrace();
                }
            }
        });

        try {
            while (true) {
                Future<OutputStream> wrk = toSamdict.take();
                if (wrk == NULL_FUTURE) {
                    break;
                }
                System.out.print(wrk.get());
            }
        } catch (InterruptedException | ExecutionException e) {
            e.printStackTrace();
        }
        executor.shutdown();
    }

    private static void somaticNotParallel(final List<List<Region>> segs,
            final Map<String, Integer> chrs, final String ampliconBasedCalling,
            final String sample, String samplem, final Configuration conf) throws IOException {
        for (List<Region> list : segs) {
            for (Region region : list) {
                final Set<String> splice = new ConcurrentHashSet<>();
                Map<Integer, Character> ref = getREF(region, chrs, conf.fasta, conf.numberNucleotideToExtend);
                Tuple2<Integer, Map<Integer, Vars>> t1 = toVars(region, conf.bam.getBam1(), ref, chrs, sample, splice, ampliconBasedCalling, 0, conf);
                Tuple2<Integer, Map<Integer, Vars>> t2 = toVars(region, conf.bam.getBam2(), ref, chrs, sample, splice, ampliconBasedCalling, t1._1, conf);
                somdict(region, t1._2, t2._2, sample, chrs, splice, ampliconBasedCalling, Math.max(t1._1, t2._1), conf, System.out);
            }
        }

    }

    final static Pattern SAMPLE_PATTERN = Pattern.compile("([^\\/\\._]+).sorted[^\\/]*.bam");
    final static Pattern SAMPLE_PATTERN2 = Pattern.compile("([^\\/]+)[_\\.][^\\/]*bam");

    private static Tuple2<String, String> getSampleNames(String bam, String sampleName, String regexp) {
        String sample = null;
        String samplem = "";

        if (sampleName != null) {
            sample = sampleName;
        } else {
            Pattern rn;
            if (regexp == null || regexp.isEmpty()) {
                rn = SAMPLE_PATTERN;
            } else {
                rn = Pattern.compile(regexp);
            }

            Matcher matcher = rn.matcher(bam);
            if (matcher.find()) {
                sample = matcher.group(1);
            }
        }

        if (sample == null) {
            Matcher matcher = SAMPLE_PATTERN2.matcher(bam);
            if (matcher.find()) {
                sample = matcher.group(1);
            }
        }

        return tuple(sample, samplem);
    }

    private static Tuple2<String, String> getSampleNames(String sample, String samplem, String bam1, String bam2, String sampleName, String regexp) {
        if (regexp != null) {
            Pattern rn = Pattern.compile(regexp);
            Matcher m = rn.matcher(bam1);
            if (m.find()) {
                sample = m.group(1);
            }
            m = rn.matcher(bam2);
            if (m.find()) {
                samplem = m.group(1);
            }
        } else {
            if (sampleName != null) {
                String[] split = sampleName.split("\\|");
                sample = split[0];
                if (split.length > 1) {
                    samplem = split[1];
                } else {
                    samplem = split[0] + "_match";
                }
            }
        }

        return tuple(sample, samplem);
    }

    /**
     * Paired sample variant calling
     * @param segs region
     * @param vars1 variants from BAM1
     * @param vars2 variants from BAM2
     * @param sample sample name
     * @param chrs map of chromosome lengths
     * @param splice set of strings representing introns in splice
     * @param ampliconBasedCalling string of maximum_distance:minimum_overlap for amplicon based calling
     * @param rlen max read length
     * @param conf Configuration
     * @param out output stream
     * @return maximum read length
     * @throws IOException
     */
    static int somdict(Region segs, Map<Integer, Vars> vars1, Map<Integer, Vars> vars2,
            String sample,
            Map<String, Integer> chrs,
            Set<String> splice,
            String ampliconBasedCalling,
            int rlen,
            Configuration conf, PrintStream out) throws IOException {

        Set<Integer> ps = new HashSet<>(vars1.keySet());
        ps.addAll(vars2.keySet());
        List<Integer> pp = new ArrayList<>(ps);
        Collections.sort(pp);

        for (Integer p : pp) {
            Vars v1 = vars1.get(p);
            Vars v2 = vars2.get(p);
            if (v1 == null && v2 == null) { // both samples have no coverage
                continue;
            }

            String vartype = "";

            if (v1 == null) { // no coverage for sample 1
                if (v2.var.isEmpty()) {
                    continue;
                }
                Variant var = v2.var.get(0);

                vartype = varType(var.refallele, var.varallele);
                if (!isGoodVar(var, v2.ref, vartype, splice, conf)) {
                    continue;
                }
                if (vartype.equals("Complex")) {
                    adjComplex(var);
                }
                out.println(join("\t", sample, segs.gene, segs.chr,
                        var.sp,
                        var.ep,
                        var.refallele,
                        var.varallele,

                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,

                        joinVar2(var, "\t"),
                        format("%.1f", var.nm),

                        var.shift3,
                        var.msi == 0 ? 0 : format("%.3f", var.msi),
                        var.msint,
                        var.leftseq,
                        var.rightseq,

                        segs.chr + ":" + segs.start + "-" + segs.end,
                        "Deletion", vartype));
            } else if (v2 == null) { // no coverage for sample 2
                if (v1.var.isEmpty()) {
                    continue;
                }
                Variant var = v1.var.get(0);
                vartype = varType(var.refallele, var.varallele);
                if (!isGoodVar(var, v1.ref, vartype, splice, conf)) {
                    continue;
                }
                if (vartype.equals("Complex")) {
                    adjComplex(var);
                }
                out.println(join("\t", sample, segs.gene, segs.chr,

                        var.sp,
                        var.ep,
                        var.refallele,
                        var.varallele,

                        joinVar2(var, "\t"),
                        format("%.1f", var.nm),

                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,

                        var.shift3,
                        var.msi == 0 ? 0 : format("%.3f", var.msi),
                        var.msint,
                        var.leftseq,
                        var.rightseq,

                        segs.chr + ":" + segs.start + "-" + segs.end,
                        "SampleSpecific", vartype));

            } else { // both samples have coverage
                if (v1.var.isEmpty() && v2.var.isEmpty()) {
                    continue;
                }
                if (v1.var.size() > 0) {
                    int n = 0;
                    while (n < v1.var.size()
                            && isGoodVar(v1.var.get(n), v1.ref, varType(v1.var.get(n).refallele, v1.var.get(n).varallele), splice, conf)) {
                        final Variant vref = v1.var.get(n);
                        final String nt = vref.n;
//                        if (nt.length() > 1
//                                && vref.refallele.length() == vref.varallele.length()
//                                && (!isGoodVar(getVarMaybe(vars2, p, varn, nt), null, null, splice, conf))
//                                && !vref.genotype.contains("-")
//                                && !vref.genotype.contains("m")
//                                && !vref.genotype.contains("i")) {
//
//                            String fnt = substr(nt, 0, -1).replaceFirst("&$", "");
//                            String lnt = substr(nt, 1).replaceFirst("^&", "");
//                            if (lnt.length() > 1) {
//                                lnt = lnt.charAt(0) + "&" + lnt.substring(1);
//                            }
//                            Variant vf = getVarMaybe(v2, varn, fnt);
//                            Variant vl = getVarMaybe(vars2, p + nt.length() - 2, varn, lnt);
//                            if (vf != null && isGoodVar(vf, v2.ref, null, splice, conf)) {
//                                vref.sp += vref.refallele.length() - 1;
//                                vref.refallele = substr(vref.refallele, -1);
//                                vref.varallele = substr(vref.varallele, -1);
//                            } else if (vl != null && isGoodVar(vl, getVarMaybe(vars2, p + nt.length() - 2, VarsType.ref), null, splice, conf)) {
//                                vref.ep += vref.refallele.length() - 1;
//                                vref.refallele = substr(vref.refallele, 0, -1);
//                                vref.varallele = substr(vref.varallele, 0, -1);
//                            }
//
//                        }
                        vartype = varType(vref.refallele, vref.varallele);
                        if (vartype.equals("Complex")) {
                            adjComplex(vref);
                        }
                        Variant v2nt = getVarMaybe(v2, varn, nt);
                        if (v2nt != null) {
                            String type;
                            if (isGoodVar(v2nt, v2.ref, vartype, splice, conf)) {
                                if (vref.freq > (1 - conf.lofreq) && v2nt.freq < 0.8d && v2nt.freq > 0.2d) {
                                    type = "LikelyLOH";
                                } else {
                                    type = "Germline";
                                }
                            } else {
                                if (v2nt.freq < conf.lofreq || v2nt.cov <= 1) {
                                    type = "LikelySomatic";
                                } else {
                                    type = "AFDiff";
                                }
                            }
                            if (isNoise(v2nt, conf.goodq, conf.lofreq) && vartype.equals("SNV")) {
                                type = "StrongSomatic";
                            }
                            out.println(join("\t", sample, segs.gene, segs.chr,
                                    vref.sp,
                                    vref.ep,
                                    vref.refallele,
                                    vref.varallele,

                                    joinVar2(vref, "\t"),
                                    format("%.1f", vref.nm),

                                    joinVar2(v2nt, "\t"),
                                    format("%.1f", v2nt.nm),

                                    v2nt.shift3,
                                    v2nt.msi == 0 ? 0 : format("%.3f", v2nt.msi),
                                    v2nt.msint,
                                    v2nt.leftseq,
                                    v2nt.rightseq,

                                    segs.chr + ":" + segs.start + "-" + segs.end,
                                    type, vartype));

                        } else { // sample 1 only, should be strong somatic
                            String type = "StrongSomatic";
                            v2nt = new Variant();
                            v2.varn.put(nt, v2nt); // Ensure it's initialized before passing to combineAnalysis
                            Tuple2<Integer, String> tpl = combineAnalysis(vref, v2nt, segs.chr, p, nt, chrs, sample, splice, ampliconBasedCalling, rlen, conf);
                            rlen = tpl._1;
                            String newtype = tpl._2;
                            if ("FALSE".equals(newtype)) {
                                n++;
                                continue;
                            }
                            if (newtype.length() > 0) {
                                type = newtype;
                            }
                            if (type.equals("StrongSomatic")) {
                                String tvf;
                                if (v2.ref == null) {
                                    Variant v2m = getVarMaybe(v2, var, 0);
                                    int tcov = v2m != null && v2m.tcov != 0 ? v2m.tcov : 0;
                                    tvf = join("\t", tcov, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
                                } else {
                                    Variant v2ref = vars2.get(p).ref;
                                    tvf = join("\t", joinVar2(v2ref, "\t"), format("%.1f", v2ref.nm));
                                }
                                out.println(join("\t", sample, segs.gene, segs.chr,
                                        vref.sp,
                                        vref.ep,
                                        vref.refallele,
                                        vref.varallele,

                                        joinVar2(vref, "\t"),
                                        format("%.1f", vref.nm),

                                        tvf,
                                        vref.shift3,
                                        vref.msi == 0 ? 0 : format("%.3f", vref.msi),
                                        vref.msint,
                                        vref.leftseq,
                                        vref.rightseq,

                                        segs.chr + ":" + segs.start + "-" + segs.end,
                                        "StrongSomatic", vartype));
                            } else {
                                out.println(join("\t", sample, segs.gene, segs.chr,

                                        vref.sp,
                                        vref.ep,
                                        vref.refallele,
                                        vref.varallele,

                                        joinVar2(vref, "\t"),
                                        format("%.1f", vref.nm),

                                        joinVar2(v2nt, "\t"),
                                        format("%.1f", v2nt.nm),

                                        vref.shift3,
                                        vref.msi == 0 ? 0 : format("%.3f", vref.msi),
                                        vref.msint,
                                        vref.leftseq,
                                        vref.rightseq,
                                        segs.chr + ":" + segs.start + "-" + segs.end,

                                        type, vartype));
                            }
                        }
                        n++;
                    }
                    if (n == 0) {
                        if (v2.var.isEmpty()) {
                            continue;
                        }
                        Variant v2var = v2.var.get(0);
                        vartype = varType(v2var.refallele, v2var.varallele);
                        if (!isGoodVar(v2var, v2.ref, vartype, splice, conf)) {
                            continue;
                        }
                        // potentail LOH
                        String nt = v2var.n;
                        Variant v1nt = getVarMaybe(v1, varn, nt);
                        if (v1nt != null) {
                            String type = v1nt.freq < conf.lofreq ? "LikelyLOH" : "Germline";
                            if ("Complex".equals(vartype)) {
                                adjComplex(v1nt);
                            }
                            out.println(join("\t", sample, segs.gene, segs.chr,
                                    v1nt.sp,
                                    v1nt.ep,
                                    v1nt.refallele,
                                    v1nt.varallele,

                                    joinVar2(v1nt, "\t"),
                                    format("%.1f", v1nt.nm),

                                    joinVar2(v2var, "\t"),
                                    format("%.1f", v2var.nm),

                                    v2var.shift3,
                                    v2var.msi == 0 ? 0 : format("%.3f", v2var.msi),
                                    v2var.msint,
                                    v2var.leftseq,
                                    v2var.rightseq,

                                    segs.chr + ":" + segs.start + "-" + segs.end,
                                    type, varType(v1nt.refallele, v1nt.varallele)
                                    ));
                        } else {
                            String th1;
                            Variant v1ref = v1.ref;
                            if (v1ref != null) {
                                th1 = join("\t",
                                        joinVar2(v1ref, "\t"),
                                        format("%.1f", v1ref.nm)
                                        );
                            } else {
                                th1 = join("\t", v1.var.get(0).tcov, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
                            }
                            if ("Complex".equals(vartype)) {
                                adjComplex(v2var);
                            }
                            out.println(join("\t", sample, segs.gene, segs.chr,
                                    v2var.sp,
                                    v2var.ep,
                                    v2var.refallele,
                                    v2var.varallele,
                                    th1,

                                    joinVar2(v2var, "\t"),
                                    format("%.1f", v2var.nm),

                                    v2var.shift3,
                                    v2var.msi == 0 ? 0 : format("%.3f", v2var.msi),
                                    v2var.msint,
                                    v2var.leftseq,
                                    v2var.rightseq,

                                    segs.chr + ":" + segs.start + "-" + segs.end,
                                    "StrongLOH", vartype
                                    ));
                        }
                    }
                } else if (v2.var.size() > 0) { // sample 1 has only reference
                    Variant v2var = v2.var.get(0);
                    vartype = varType(v2var.refallele, v2var.varallele);
                    if (!isGoodVar(v2var, v2.ref, vartype, splice, conf)) {
                        continue;
                    }
                    // potential LOH
                    String nt = v2var.n;
                    String type = "StrongLOH";
                    Variant v1nt = v1.varn.get(nt);
                    if (v1nt == null) {
                        v1nt = new Variant();
                        v1.varn.put(nt, v1nt);
                    }
                    v1nt.cov = 0;
                    Tuple2<Integer, String> tpl = combineAnalysis(v2.varn.get(nt), v1nt, segs.chr, p, nt, chrs, sample, splice, ampliconBasedCalling, rlen, conf);
                    rlen = tpl._1;
                    String newtype = tpl._2;
                    if ("FALSE".equals(newtype)) {
                        continue;
                    }
                    String th1;
                    if (newtype.length() > 0) {
                        type = newtype;
                        th1 = join("\t",
                                joinVar2(v1nt, "\t"),
                                format("%.1f", v1nt.nm)
                                );
                    } else {
                        Variant v1ref = v1.ref;
                        th1 = join("\t",
                                joinVar2(v1ref, "\t"),
                                format("%.1f", v1ref.nm)
                                );
                    }

                    if ("Complex".equals(vartype)) {
                        adjComplex(v2var);
                    }
                    out.println(join("\t", sample, segs.gene, segs.chr,
                            v2var.sp,
                            v2var.ep,
                            v2var.refallele,
                            v2var.varallele,

                            th1,

                            joinVar2(v2var, "\t"),
                            format("%.1f", v2var.nm),

                            v2var.shift3,
                            v2var.msi == 0 ? 0 : format("%.3f", v2var.msi),
                            v2var.msint,
                            v2var.leftseq,
                            v2var.rightseq,

                            segs.chr + ":" + segs.start + "-" + segs.end,
                            type, vartype

                            ));
                }
            }
        }

        return rlen;
    }

    /**
     * Taken a likely somatic indels and see whether combine two bam files still support somatic status. This is mainly for Indels
     * that softclipping overhang is too short to positively being called in one bam file, but can be called in the other bam file,
     * thus creating false positives
     *
     * @param var1 variant 1
     * @param var2 variant 2
     * @param chr chromosome
     * @param p position
     * @param nt nucleotide
     *
     * @param chrs map of chromosome lengths
     * @param splice set of strings representing introns in splice
     * @param ampliconBasedCalling string of maximum_distance:minimum_overlap for amplicon based calling
     * @param rlen max read length
     * @param conf Configuration
     * @return (new <code>rlen</code>, "FALSE" | "")
     * @throws IOException
     */
    static Tuple2<Integer, String> combineAnalysis(Variant var1, Variant var2, String chr, int p, String nt,
            Map<String, Integer> chrs, String sample, Set<String> splice, String ampliconBasedCalling, int rlen, Configuration conf) throws IOException {
        if (conf.y) {
            System.err.printf("Start Combine %s %s\n", p, nt);
        }
        Region region = new Region(chr, var1.sp - rlen, var1.ep + rlen, "");
        Map<Integer, Character> ref = getREF(region, chrs, conf.fasta, conf.numberNucleotideToExtend);
        Tuple2<Integer, Map<Integer, Vars>> tpl = toVars(region, conf.bam.getBam1() + ":" + conf.bam.getBam2(), ref,
                chrs, sample, splice, ampliconBasedCalling, rlen, conf);
        rlen = tpl._1;
        Map<Integer, Vars> vars = tpl._2;
        Variant vref = getVarMaybe(vars, p, varn, nt);
        if (vref != null) {
            if (conf.y) {
                System.err.printf("Combine: 1: %s comb: %s\n", var1.cov, vref.cov);
            }
            if (vref.cov - var1.cov >= conf.minr) {
                var2.tcov = vref.tcov - var1.tcov;
                if (var2.tcov < 0)
                    var2.tcov = 0;

                var2.cov = vref.cov - var1.cov;
                if (var2.cov < 0)
                    var2.cov = 0;

                var2.rfc = vref.rfc - var1.rfc;
                if (var2.rfc < 0)
                    var2.rfc = 0;

                var2.rrc = vref.rrc - var1.rrc;
                if (var2.rrc < 0)
                    var2.rrc = 0;

                var2.fwd = vref.fwd - var1.fwd;
                if (var2.fwd < 0)
                    var2.fwd = 0;

                var2.rev = vref.rev - var1.rev;
                if (var2.rev < 0)
                    var2.rev = 0;

                var2.pmean = (vref.pmean * vref.cov - var1.pmean * var1.cov) / var2.cov;
                var2.qual = (vref.qual * vref.cov - var1.qual * var1.cov) / var2.cov;
                var2.mapq = (vref.mapq * vref.cov - var1.mapq * var1.cov) / var2.cov;
                var2.hifreq = (vref.hifreq * vref.cov - var1.hifreq * var1.cov) / var2.cov;
                var2.extrafreq = (vref.extrafreq * vref.cov - var1.extrafreq * var1.cov) / var2.cov;
                var2.nm = (vref.nm * vref.cov - var1.nm * var1.cov) / var2.cov;

                var2.pstd = true;
                var2.qstd = true;

                if (var2.tcov <= 0) {
                    return tuple(rlen, "FALSE");
                }

                var2.freq = round(var2.cov / (double)var2.tcov, 4);
                var2.qratio = var1.qratio; // Can't back calculate and should be inaccurate
                var2.genotype = vref.genotype;
                var2.bias = strandBias(var2.rfc, var2.rrc, conf.bias, conf.minb) + ";" + strandBias(var2.fwd, var2.rev, conf.bias, conf.minb);
                return tuple(rlen, "Germline");
            } else if (vref.cov < var1.cov - 2) {
                if (conf.y) {
                    System.err.printf("Combine produce less: %s %s %s %s %s\n", chr, p, nt, vref.cov, var1.cov);
                }
                return tuple(rlen, "FALSE");
            } else {
                return tuple(rlen, "");
            }

        }
        return tuple(rlen, "FALSE");
    }

    /**
     * A variance is considered noise if the quality is below <code>goodq</code> and there're no more than 3 reads
     *
     * @param vref Variant
     * @param goodq quality threshold
     * @param lofreq The minimun alelle frequency allowed in normal for a somatic mutation
     * @return Returns <tt>true</tt> if variance is considered noise if the quality is below <code>goodq</code> and there're no more
     *         than 3 reads
     */
    static boolean isNoise(Variant vref, double goodq, double lofreq) {
        final double qual = vref.qual;
        if (((qual < 4.5d || (qual < 12 && !vref.qstd)) && vref.cov <= 3)
                || (qual < goodq && vref.freq < 2 * lofreq && vref.cov <= 1)) {

            vref.tcov -= vref.cov;
            vref.cov = 0;
            vref.fwd = 0;
            vref.rev = 0;
            vref.freq = 0;
            vref.hifreq = 0;

            return true;

        }
        return false;
    }

    /**
     * Single sample mode variant calling
     * @param region region of interest
     * @param vars map of variations
     * @param sample sample name
     * @param splice set of strings representing spliced regions
     * @param conf configuration
     * @param out output stream
     */
    static void vardict(Region region, Map<Integer, Vars> vars, String sample, Set<String> splice, Configuration conf, PrintStream out) {
        for (int p = region.start; p <= region.end; p++) {
            List<String> vts = new ArrayList<>();
            List<Variant> vrefs = new ArrayList<>();
            if (!vars.containsKey(p) || vars.get(p).var.isEmpty()) {
                if (!conf.doPileup) {
                    continue;
                }
                Variant vref = getVarMaybe(vars, p, ref);
                if (vref == null) {
                    out.println(join("\t", sample, region.gene, region.chr, p, p,
                            "", "", 0, 0, 0, 0, 0, 0, "", 0, "0;0", 0, 0, 0, 0, 0, "", 0, 0, 0, 0, 0, 0, "", "", 0, 0,
                            region.chr + ":" + region.start + "-" + region.end, ""
                            ));
                    continue;
                }
                vts.add("");
                vrefs.add(vref);
            } else {
                List<Variant> vvar = vars.get(p).var;
                Variant rref = getVarMaybe(vars, p, ref);
                for (int i = 0; i < vvar.size(); i++) {
                    Variant vref = vvar.get(i);
                    if (vref.refallele.contains("N")) {
                        continue;
                    }
                    String vartype = varType(vref.refallele, vref.varallele);
                    if (!isGoodVar(vref, rref, vartype, splice, conf)) {
                        if (!conf.doPileup) {
                            continue;
                        }
                    }

                    vts.add(vartype);
                    vrefs.add(vref);
                }
            }
            for (int vi = 0; vi < vts.size(); vi++) {
                String vartype = vts.get(vi);
                Variant vref = vrefs.get(vi);
                if ("Complex".equals(vartype)) {
                    adjComplex(vref);
                }
                out.println(join("\t", sample, region.gene, region.chr,
                        vref.sp, vref.ep, vref.refallele, vref.varallele,
                        vref.tcov,
                        vref.cov,
                        vref.rfc,
                        vref.rrc,
                        vref.fwd,
                        vref.rev,
                        vref.genotype,
                        format("%.4f", vref.freq),
                        vref.bias,
                        format("%.1f", vref.pmean),
                        vref.pstd ? 1 : 0,
                        format("%.1f", vref.qual),
                        vref.qstd ? 1 : 0,
                        format("%.1f", vref.mapq),
                        format("%.3f", vref.qratio),
                        vref.hifreq == 0 ? 0 : format("%.4f", vref.hifreq),
                        vref.extrafreq == 0 ? 0 : format("%.4f", vref.extrafreq),
                        vref.shift3,
                        vref.msi == 0 ? 0 : format("%.3f", vref.msi),
                        vref.msint,
                        format("%.1f", vref.nm),
                        vref.hicnt,
                        vref.hicov,
                        vref.leftseq.isEmpty() ? "0" : vref.leftseq,
                        vref.rightseq.isEmpty() ? "0" : vref.rightseq,
                        region.chr + ":" + region.start + "-" + region.end, vartype
                        ));
                if (conf.debug) {
                    out.println("\t" + vref.DEBUG);
                }
            }

        }
    }

    private static String[] retriveSubSeq(String fasta, String chr, int start, int end) throws IOException {
        try (IndexedFastaSequenceFile idx = new IndexedFastaSequenceFile(new File(fasta))) {
            ReferenceSequence seq = idx.getSubsequenceAt(chr, start, end);
            byte[] bases = seq.getBases();
            return new String[] { ">" + chr + ":" + start + "-" + end, bases != null ? new String(bases) : "" };
        }
    }

    private static final jregex.Pattern SOFT_CLIPPED = new jregex.Pattern("(\\d+)[MIS]");
    /**
     * regexp finds numbers followed by M (matched) or D (deleted) in CIGAR string
     */
    private static final jregex.Pattern ALIGNED_LENGTH = new jregex.Pattern("(\\d+)[MD]");

    private static Variation getVariationFromSeq(Sclip sclip, int idx, Character ch) {
        Map<Character, Variation> map = sclip.seq.get(idx);
        if (map == null) {
            map = new HashMap<>();
            sclip.seq.put(idx, map);
        }
        Variation variation = map.get(ch);
        if (variation == null) {
            variation = new Variation();
            map.put(ch, variation);
        }
        return variation;
    }

    /**
     * Increase count for given key
     * @param cnts map of counts
     * @param key key to add count
     * @param add amount to add
     */
    @SuppressWarnings({ "unchecked", "rawtypes" })
    private static void incCnt(Map cnts, Object key, int add) {
        Integer integer = (Integer)cnts.get(key);
        if (integer == null) {
            cnts.put(key, add);
        } else {
            cnts.put(key, integer + add);
        }
    }

    private static Variation getVariation(Map<Integer, Map<String, Variation>> hash, int start, String ref) {
        Map<String, Variation> map = hash.get(start);
        if (map == null) {
            map = new LinkedHashMap<>();
            hash.put(start, map);
        }
        Variation variation = map.get(ref);
        if (variation == null) {
            variation = new Variation();
            map.put(ref, variation);
        }
        return variation;
    }

    private static Variation getVariationMaybe(Map<Integer, Map<String, Variation>> hash, int start, Character ref) {
        if (ref == null)
            return null;

        Map<String, Variation> map = hash.get(start);
        if (map == null) {
            return null;
        }

        return map.get(ref.toString());
    }

    private static void subCnt(Variation vref, boolean dir, int rp, double q, int Q, int nm, Configuration conf) {
        // ref dir read_position quality
        vref.cnt--;
        vref.decDir(dir);
        vref.pmean -= rp;
        vref.qmean -= q;
        vref.Qmean -= Q;
        vref.nm -= nm;
        if (q >= conf.goodq) {
            vref.hicnt--;
        } else {
            vref.locnt--;
        }
    }

    private static void addCnt(Variation vref, boolean dir, int rp, double q, int Q, int nm, int goodq) {
        vref.cnt++;
        vref.incDir(dir);
        vref.pmean += rp;
        vref.qmean += q;
        vref.Qmean += Q;
        vref.nm += nm;
        if (q >= goodq) {
            vref.hicnt++;
        } else {
            vref.locnt++;
        }
    }

    public static class Sclip extends Variation {
        TreeMap<Integer, Map<Character, Integer>> nt = new TreeMap<>();
        TreeMap<Integer, Map<Character, Variation>> seq = new TreeMap<>();
        String sequence;
        boolean used;
    }

    /**
     * Intermediate variant structure
     */
    public static class Variation {
        /**
         * Variant count
         */
        int cnt;

        /**
         * Variant count on forward strand
         */
        int dirPlus;

        /**
         * Variant count on reverse strand
         */
        int dirMinus;

        /**
         * Sum of variant positions in read
         */
        int pmean;

        /**
         * Sum of base qualities for variant
         */
        double qmean;

        /**
         * Sum of mapping qualities for variant
         */
        int Qmean;

        /**
         * Sum of number of mismatches for variant
         */
        int nm;

        /**
         * Number of low-quality reads with the variant
         */
        int locnt;

        /**
         * Number of high-quality reads with the variant
         */
        int hicnt;

        /**
         * Flags that is true when variant is found in at least 2 different positions
         */
        boolean pstd;

        /**
         * Flags that is 1 when variant is read with at least 2 different qualities
         */
        boolean qstd;

        /**
         * Position in read for previous instance of this variant (used for pstd)
         */
        int pp;

        /**
         * Base quality for previous instance of this variant (used for qstd)
         */
        double pq;

        /**
         * Adjusted count for indels due to local realignment
         */
        int extracnt;

        /**
         * Increment count for direction
         * @param dir true for forward strand, false for reverse strand
         */
        public void incDir(boolean dir) {
            if (dir)
                this.dirMinus++;
            else
                this.dirPlus++;
        }

        /**
         * Decrement count for direction
         * @param dir true for forward strand, false for reverse strand
         */
        public void decDir(boolean dir) {
            if (dir)
                this.dirMinus--;
            else
                this.dirPlus--;
        }

        /**
         * Get variant count for direction
         * @param dir true for forward strand, false for reverse strand
         * @return variant count
         */
        public int getDir(boolean dir) {
            if (dir)
                return this.dirMinus;
            return this.dirPlus;
        }

        /**
         * Add count for direction
         * @param dir true for forward strand, false for reverse strand
         */
        public void addDir(boolean dir, int add) {
            if (dir)
                this.dirMinus += add;
            else
                this.dirPlus += add;
        }

        /**
         * Subtract count for direction
         * @param dir true for forward strand, false for reverse strand
         */
        public void subDir(boolean dir, int sub) {
            if (dir)
                this.dirMinus -= sub;
            else
                this.dirPlus -= sub;
        }

    }

    /**
     * Regexp finds number followed by S followed by number and I or D at the start of string
     */
    private static final jregex.Pattern BEGIN_NUMBER_S_NUMBER_IorD = new jregex.Pattern("^(\\d+)S(\\d+)([ID])");
    /**
     * Regexp finds number followed by I or D followed by number followed by S followed by end of string
     */
    private static final jregex.Pattern NUMBER_IorD_NUMBER_S_END = new jregex.Pattern("(\\d+)([ID])(\\d+)S$");
    /**
     * Regexp finds number followed by S followed by number followed by M followed by number followed by I or D at the start of string
     */
    private static final jregex.Pattern BEGIN_NUMBER_S_NUMBER_M_NUMBER_IorD = new jregex.Pattern("^(\\d+)S(\\d+)M(\\d+)([ID])");
    /**
     * Regexp finds number followed by I or D followed by number followed by M followed by number followed by S followed by end of string
     */
    private static final jregex.Pattern NUMBER_IorD_NUMBER_M_NUMBER_S_END = new jregex.Pattern("(\\d+)([ID])(\\d+)M(\\d+)S$");
    /**
     * Regexp finds digit-M-number-(I or D)-number-M
     */
    private static final jregex.Pattern BEGIN_DIGIT_M_NUMBER_IorD_NUMBER_M = new jregex.Pattern("^(\\d)M(\\d+)([ID])(\\d+)M");
    private static final Pattern BEGIN_DIGIT_M_NUMBER_IorD_NUMBER_M_ = Pattern.compile("^\\dM\\d+[ID]\\d+M");
    /**
     * Regexp replaces number-(I or D)-number-M
     */
    private static final jregex.Pattern NUMBER_IorD_DIGIT_M_END = new jregex.Pattern("(\\d+)([ID])(\\d)M$");

    /**
     * Regexp finds number-M-number-D-digit-M-number-I-digit-M-number-D.
     * ^(.*?) captures everything before the M-D-M-I-M-D complex
     */
    private static final jregex.Pattern D_M_D_DD_M_D_I_D_M_D_DD = new jregex.Pattern("^(.*?)(\\d+)M(\\d+)D(\\d)M(\\d+)I(\\d)M(\\d+)D");
    private static final Pattern D_M_D_DD_M_D_I_D_M_D_DD_prim = Pattern.compile("(\\d+)M(\\d+)D(\\d)M(\\d+)I(\\d)M(\\d+)D");
    /**
     * Regexp finds number-D-digit-M-number-D and optionally number-I
     */
    private static final Pattern DIG_D_DIG_M_DIG_DI_DIGI = Pattern.compile("(\\d+)D(\\d)M(\\d+)([DI])(\\d+I)?");
    /**
     * Regexp findx number-I-digit-M-number and optionally number-I
     */
    private static final Pattern DIG_I_dig_M_DIG_DI_DIGI = Pattern.compile("(\\d+)I(\\d)M(\\d+)([DI])(\\d+I)?");

    private static final Pattern BEGIN_DIGITS = Pattern.compile("^(\\d+)");
    private static final Pattern HASH_GROUP_CARET_GROUP = Pattern.compile("#(.+)\\^(.+)");
    /**
     * Regexp finds number-S-number-M at start of CIGAR string
     */
    private static final Pattern D_S_D_M = Pattern.compile("^(\\d+)S(\\d+)M");
    private static final Pattern D_M_D_S_END = Pattern.compile("\\d+M\\d+S$");
    /**
     * Regexp finds number-M-number-S at end of CIGAR string
     * ^(.*?) captures everything before M-S complex
     */
    private static final Pattern ANY_NUMBER_M_NUMBER_S_END = Pattern.compile("^(.*?)(\\d+)M(\\d+)S$");
    private static final Pattern START_DIG = Pattern.compile("^-(\\d+)");

    private static class SamView implements AutoCloseable {

        private SAMRecordIterator iterator;
        private SamReader reader;
        private int filter = 0;

        public SamView(String file, String samfilter, Region region, ValidationStringency stringency) {
            reader = SamReaderFactory.makeDefault().validationStringency(stringency).open(new File(file));
            iterator = reader.queryOverlapping(region.chr, region.start, region.end);
            if (!"".equals(samfilter)) {
                filter = Integer.decode(samfilter);
            }

        }

        public SAMRecord read() throws IOException {
            while(iterator.hasNext()) {
                SAMRecord record = iterator.next();
                if (filter != 0 && (record.getFlags() & filter) != 0) {
                    continue;
                }
                return record;
            }
            return null;
        }


        @Override
        public void close() throws IOException {
            reader.close();
        }
    }


    /**
     * Construct a variant structure given a region and BAM files.
     * @param region region
     * @param bam BAM file name
     * @param chrs map of chromosome lengths
     * @param splice set of strings representing spliced regions
     * @param ampliconBasedCalling string of maximum_distance:minimum_overlap for amplicon based calling
     * @param rlen max read length
     * @param ref reference in a given region
     * @param conf Configuration
     * @return Tuple of (noninsertion variant  structure, insertion variant structure, coverage, maxmimum read length)
     * @throws IOException
     */
    static Tuple4<Map<Integer, Map<String, Variation>>, Map<Integer, Map<String, Variation>>, Map<Integer, Integer>, Integer> parseSAM(Region region, String bam,
            Map<String, Integer> chrs, String sample, Set<String> splice, String ampliconBasedCalling, int rlen, Map<Integer, Character> ref, Configuration conf) throws IOException {

        String[] bams = bam.split(":");

        Map<Integer, Map<String, Variation>> hash = new HashMap<>();
        Map<Integer, Map<String, Variation>> iHash = new HashMap<>();
        Map<Integer, Integer> cov = new HashMap<>();
        Map<Integer, Sclip> sclip3 = new HashMap<>(); // soft clipped at 3'
        Map<Integer, Sclip> sclip5 = new HashMap<>(); // soft clipped at 5'
        Map<Integer, Map<String, Integer>> ins = new HashMap<>();
        Map<Integer, Map<String, Integer>> mnp = new HashMap<>(); // Keep track of MNPs
        Map<Integer, Map<String, Integer>> dels5 = new HashMap<>();
        Map<String, int[]> spliceCnt = new HashMap<String, int[]>();

        String chr = region.chr;
        if (conf.chromosomeNameIsNumber && chr.startsWith("chr")) { //remove prefix 'chr' if option -C is set
            chr = region.chr.substring("chr".length());
        }

        for (String bami : bams) {
            String samfilter = conf.samfilter == null || conf.samfilter.isEmpty() ? "" : conf.samfilter;
            try (SamView reader =  new SamView(bami, samfilter, region, conf.validationStringency)) {
                //dup contains already seen reads. For each seen read dup contains either POS-RNEXT-PNEXT or POS-CIGAR (if next segment in template is unmapped).
                Set<String> dup = new HashSet<>();
                //position of first matching base (POS in SAM)
                int dupp = -1;
                SAMRecord record;
                while ((record = reader.read()) != null) {
                    if (conf.isDownsampling() && RND.nextDouble() <= conf.downsampling) {
                        continue;
                    }

                    final String querySequence = record.getReadString();
                    final Flags flag = new Flags(record.getFlags());
                    final int mappingQuality = record.getMappingQuality();

                    if (conf.hasMappingQuality() && mappingQuality < conf.mappingQuality) { // ignore low mapping quality reads
                        continue;
                    }

                    if (flag.isNotPrimaryAlignment() && conf.samfilter != null) {
                        continue;
                    }

                    if (querySequence.length() == 1 && querySequence.charAt(0) == '*') {
                        continue;
                    }


                    final String mrnm = getMrnm(record);


                    // filter duplicated reads if option -t is set
                    if (conf.removeDuplicatedReads) {
                        if (record.getAlignmentStart() != dupp) {
                            dup.clear();
                        }
                        if (record.getMateAlignmentStart() < 10) {
                            String dupKey = record.getAlignmentStart() + "-" + mrnm + "-" + record.getMateAlignmentStart();
                            if (dup.contains(dupKey)) {
                                continue;
                            }
                            dup.add(dupKey);
                            dupp = record.getAlignmentStart();
                        } else if (flag.isUnmappedMate()) {
                            String dupKey = record.getAlignmentStart() + "-" + record.getCigarString();
                            if (dup.contains(dupKey)) {
                                continue;
                            }
                            dup.add(dupKey);
                            dupp = record.getAlignmentStart();
                        }
                    }

                    final Cigar readCigar = record.getCigar();
                    final int indel = getInsertionDeletionLenght(readCigar);

                    int tnm = 0;
                    Integer nmi = record.getIntegerAttribute(SAMTag.NM.name());
                    if (nmi != null) { // number of mismatches. Don't use NM since it includes gaps, which can be from indels
                        tnm = nmi - indel;
                        if (tnm > conf.mismatch) { // edit distance - indels is the # of mismatches
                            continue;
                        }
                    } else { //Skip the read if number of mismatches is not available
                        if (conf.y && !record.getCigarString().equals("*")) {
                            System.err.println("No XM tag for mismatches. " + record.getSAMString());
                        }
                        continue;
                    }


                    final String queryQuality = record.getBaseQualityString();
                    final boolean isMrnmEqual = record.getReferenceName().equals(record.getMateReferenceName());

                    //number of mismatches
                    final int nm = tnm;
                    int n = 0; // keep track the read position, including softclipped
                    int p = 0; // keep track the position in the alignment, excluding softclipped
                    boolean dir = flag.isReverseStrand();
                    if (ampliconBasedCalling != null) {
                        String[] split = ampliconBasedCalling.split(":");
                        //distance to amplicon (specified in -a option)
                        int dis;
                        //overlap fraction (in -a option)
                        double ovlp;
                        try {
                            dis = toInt(split[0]);
                            ovlp = Double.parseDouble(split.length > 1 ? split[1] : "");
                        } catch (NumberFormatException e) {
                            dis = 10;
                            ovlp = 0.95;
                        }
                        //rlen3 holds sum of lengths of matched and deleted segments
                        int rlen3 = getAlignedLenght(readCigar); // The total aligned length, excluding soft-clipped
                                                                                // bases and insertions
                        int segstart = record.getAlignmentStart();
                        int segend = segstart + rlen3 - 1;


                        if (readCigar.getCigarElement(0).getOperator() == CigarOperator.S) { //If read starts with soft-clipped sequence
                            //Ignore reads that overlap with region of interest by fraction less than ovlp
                            int ts1 = segstart > region.start ? segstart : region.start;
                            int te1 = segend < region.end ? segend : region.end;
                            if (Math.abs(ts1 - te1) / (double)(segend - segstart) > ovlp == false) {
                                continue;
                            }
                        } else if (readCigar.getCigarElement(readCigar.numCigarElements() - 1).getOperator() == CigarOperator.S) { //If read ends with contains soft-clipped sequence
                            //Ignore reads that overlap with region of interest by fraction less than ovlp
                            int ts1 = segstart > region.start ? segstart : region.start;
                            int te1 = segend < region.end ? segend : region.end;
                            if (Math.abs(te1 - ts1) / (double)(segend - segstart) > ovlp == false) {
                                continue;
                            }

                        } else { //no soft-clipping
                            //if RNEXT is identical to RNAME and TLEN is defined
                            if ( isMrnmEqual && record.getInferredInsertSize() != 0) {
                                if (record.getInferredInsertSize() > 0) {
                                    segend = segstart + record.getInferredInsertSize() - 1;
                                } else {
                                    segstart = record.getMateAlignmentStart();
                                    segend = record.getMateAlignmentStart() - record.getInferredInsertSize() - 1;
                                }
                            }
                            // No segment overlapping test since samtools should take care of it
                            int ts1 = segstart > region.start ? segstart : region.start;
                            int te1 = segend < region.end ? segend : region.end;
                            //ignore reads that are more than dis from region of interest and overlap is less than ovlp
                            if ((abs(segstart - region.start) > dis || abs(segend - region.end) > dis)
                                    || abs((ts1 - te1) / (double)(segend - segstart)) <= ovlp) {
                                continue;
                            }
                        }
                    }
                    if (flag.isUnmappedMate()) {
                        // to be implemented
                    } else {
                        if (isMrnmEqual) {
                            //if 'SA' tag (supplementary alignment) is present
                            if (record.getStringAttribute(SAMTag.SA.name()) != null) {
                                if (flag.isSupplementaryAlignment()) { // the supplementary alignment
                                    continue; // Ignore the supplementary for now so that it won't skew the coverage
                                }
                            }
                        }

                    }
                    final int position;
                    final Cigar cigar;

                    if (conf.performLocalRealignment) {
                        // Modify the CIGAR for potential mis-alignment for indels at the end of reads to softclipping and let VarDict's
                        // algorithm to figure out indels
                        Tuple2<Integer, String> mc = modifyCigar(indel, ref, record.getAlignmentStart(), record.getCigarString(), querySequence, queryQuality, conf.lowqual);
                        position = mc._1;
                        cigar = TextCigarCodec.decode(mc._2);
                    } else {
                        position = record.getAlignmentStart();
                        cigar = record.getCigar();
                    }

                    //adjusted start position
                    int start = position;
                    int offset = 0;
                     // Only match and insertion counts toward read length
                     // For total length, including soft-clipped bases
                    int rlen1 = getMatchInsertionLenght(cigar); // The read length for matched bases
                    if ( conf.minmatch != 0 && rlen1 < conf.minmatch) {
                        continue;
                    }
                    int rlen2 = getSoftClippedLenght(cigar); // The total length, including soft-clipped bases
                    if (rlen2 > rlen) { // Determine the read length
                        rlen = rlen2;
                    }

                    //Loop over CIGAR records
                    for (int ci = 0; ci < cigar.numCigarElements(); ci++) {

                        //length of segment in CIGAR
                        int m =  cigar.getCigarElement(ci).getLength();

                        //letter from CIGAR
                        final CigarOperator operator = getCigarOperator(cigar, ci);

                        switch (operator) {
                            case N: //N in CIGAR - skipped region from reference
                                //Skip the region and add string start-end to %SPLICE
                                String key = (start - 1) + "-" + (start + m - 1);
                                splice.add(key);
                                int[] cnt = spliceCnt.get(key);
                                if (cnt == null) {
                                    cnt = new int[] { 0 };
                                    spliceCnt.put(key, cnt);
                                }
                                cnt[0]++;

                                start += m;
                                offset = 0;
                                continue;

                            case S:
                                //First record in CIGAR
                                if (ci == 0) { // 5' soft clipped
                                    // align softclipped but matched sequences due to mis-softclipping
                                    /*
                                    Conditions:
                                    1). segment length > 1
                                    2). start between 0 and chromosome length,
                                    3). reference genome is known at this position
                                    4). reference and read bases match
                                    5). read quality is more than 10
                                     */
                                    while (m - 1 >= 0 && start - 1 > 0 && start - 1 <= chrs.get(chr)
                                            && isHasAndEquals(querySequence.charAt(m - 1), ref, start - 1)
                                            && queryQuality.charAt(m - 1) - 33 > 10) {
                                        //create variant if it is not present
                                        Variation variation = getVariation(hash, start - 1, ref.get(start - 1).toString());
                                        //add count
                                        addCnt(variation, dir, m, queryQuality.charAt(m - 1) - 33, mappingQuality, nm, conf.goodq);
                                        //increase coverage
                                        incCnt(cov, start - 1, 1);
                                        start--;
                                        m--;
                                    }
                                    if (m > 0) { //If there remains a soft-clipped sequence at the beginning (not everything was matched)
                                        int q = 0; //sum of read qualities (to get mean quality)
                                        int qn = 0; //number of quality figures
                                        int lowqcnt = 0; //number of low-quality bases in the sequence

                                        //loop over remaining soft-clipped sequence
                                        for (int si = m - 1; si >= 0; si--) {
                                            //stop if unknown base (N - any of ATGC) is found
                                            if (querySequence.charAt(si) == 'N') {
                                                break;
                                            }
                                            //tq - base quality
                                            int tq = queryQuality.charAt(si) - 33;
                                            if (tq <= 12)
                                                lowqcnt++;
                                            //Stop if a low-quality base is found
                                            if (lowqcnt > 1)
                                                break;

                                            q += tq;
                                            qn++;
                                        }
                                        //If we have at least 1 high-quality soft-clipped base within conf.buffer of region of interest
                                        if (qn >= 1 && qn > lowqcnt && start >= region.start - conf.buffer && start <= region.end + conf.buffer) {
                                            //add record to $sclip5
                                            Sclip sclip = sclip5.get(start);
                                            if (sclip == null) {
                                                sclip = new Sclip();
                                                sclip5.put(start, sclip);
                                            }
                                            for (int si = m - 1; m - si <= qn; si--) {
                                                Character ch = querySequence.charAt(si);
                                                int idx = m - 1 - si;
                                                Map<Character, Integer> cnts = sclip.nt.get(idx);
                                                if (cnts == null) {
                                                    cnts = new LinkedHashMap<>();
                                                    sclip.nt.put(idx, cnts);
                                                }
                                                incCnt(cnts, ch, 1);
                                                Variation seqVariation = getVariationFromSeq(sclip, idx, ch);
                                                addCnt(seqVariation, dir, si - (m - qn), queryQuality.charAt(si) - 33, mappingQuality, nm, conf.goodq);
                                            }
                                            addCnt(sclip, dir, m, q / (double)qn, mappingQuality, nm, conf.goodq);
                                        }

                                    }
                                    m = cigar.getCigarElement(ci).getLength();
                                } else if (ci == cigar.numCigarElements() - 1) { // 3' soft clipped
                                    /*
                                    Conditions:
                                    1). read position is less than sequence length
                                    2). reference base is defined for start
                                    3). reference base at start matches read base at n
                                    4). read quality is more than 10
                                     */
                                    while (n < querySequence.length()
                                            && isHasAndEquals(querySequence.charAt(n), ref, start)
                                            && queryQuality.charAt(n) - 33 > 10) {
                                        //initialize entry in $hash if not present
                                        Variation variation = getVariation(hash, start, ref.get(start).toString());
                                        //add count
                                        addCnt(variation, dir, rlen2 - p, queryQuality.charAt(n) - 33, mappingQuality, nm, conf.goodq);
                                        //add coverage
                                        incCnt(cov, start, 1);
                                        n++;
                                        start++;
                                        m--;
                                        p++;
                                    }
                                    if (querySequence.length() - n > 0) { //If there remains a soft-clipped sequence at the end (not everything was matched)
                                        int q = 0; //sum of read qualities (to get mean quality)
                                        int qn = 0; //number of quality figures
                                        int lowqcnt = 0; //number of low-quality bases in the sequence
                                        for (int si = 0; si < m; si++) { //loop over remaining soft-clipped sequence
                                            if (querySequence.charAt(n + si) == 'N') { //stop if unknown base (N - any of ATGC) is found
                                                break;
                                            }
                                            int tq = queryQuality.charAt(n + si) - 33; //base quality
                                            if (tq <= 12) {
                                                lowqcnt++;
                                            }
                                            //Stop if a low-quality base is found
                                            if (lowqcnt > 1) {
                                                break;
                                            }
                                            q += tq;
                                            qn++;
                                        }
                                        //If we have at least 1 high-quality soft-clipped base within conf.buffer of region of interest
                                        if (qn >= 1 && qn > lowqcnt && start >= region.start - conf.buffer && start <= region.end + conf.buffer) {
                                            //add record to $sclip3
                                            Sclip sclip = sclip3.get(start);
                                            if (sclip == null) {
                                                sclip = new Sclip();
                                                sclip3.put(start, sclip);
                                            }
                                            for (int si = 0; si < qn; si++) {
                                                Character ch = querySequence.charAt(n + si);
                                                int idx = si;
                                                Map<Character, Integer> cnts = sclip.nt.get(idx);
                                                if (cnts == null) {
                                                    cnts = new HashMap<>();
                                                    sclip.nt.put(idx, cnts);
                                                }
                                                incCnt(cnts, ch, 1);
                                                Variation variation = getVariationFromSeq(sclip, idx, ch);
                                                addCnt(variation, dir, qn - si, queryQuality.charAt(n + si) - 33, mappingQuality, nm, conf.goodq);
                                            }
                                            addCnt(sclip, dir, m, q / (double)qn, mappingQuality, nm, conf.goodq);
                                        }

                                    }

                                }
                                //move read position by m (length of segment in CIGAR)
                                n += m;
                                offset = 0;
                                start = position; // had to reset the start due to softclipping adjustment
                                continue;
                            case H: //Hard clipping - skip
                                offset = 0;
                                continue;
                            case I: { //Insertion
                                offset = 0;
                                //inserted segment of read sequence
                                StringBuilder s = new StringBuilder(substr(querySequence, n, m));
                                //quality of this segment
                                StringBuilder q = new StringBuilder(substr(queryQuality, n, m));
                                //sequence to be appended if next segment is matched
                                String ss = "";

                                // For multiple indels within 10bp

                                //offset for read position if next segment is matched
                                int multoffs = 0;
                                //offset for reference position if next segment is matched
                                int multoffp = 0;
                                int nmoff = 0;

                                /*
                                Condition:
                                1). CIGAR string has next entry
                                2). length of next CIGAR segment is less than conf.vext
                                3). next segment is matched
                                4). CIGAR string has one more entry after next one
                                5). this entry is insertion or deletion
                                 */
                                if (conf.performLocalRealignment && cigar.numCigarElements() > ci + 2
                                        && cigar.getCigarElement(ci + 1).getLength() <= conf.vext
                                        && cigar.getCigarElement(ci + 1).getOperator() == CigarOperator.M
                                        && (cigar.getCigarElement(ci + 2).getOperator()  == CigarOperator.I
                                            || cigar.getCigarElement(ci + 2).getOperator()  == CigarOperator.D)) {

                                    int mLen = cigar.getCigarElement(ci + 1).getLength();
                                    int indelLen = cigar.getCigarElement(ci + 2).getLength();
                                    CigarOperator indelOperator = cigar.getCigarElement(ci + 2).getOperator();
                                    //append to s '#' and part of read sequence corresponding to next CIGAR segment (matched one)
                                    s.append("#").append(substr(querySequence, n + m, mLen));
                                    //append next segment quality to q
                                    q.append(substr(queryQuality, n + m, mLen));

                                    //if an insertion is two segments ahead, append '^' + part of sequence corresponding to next-next segment
                                    //otherwise (deletion) append '^' + length of a next-next segment
                                    s.append('^').append(indelOperator == CigarOperator.I ? substr(querySequence, n + m + mLen, indelLen) : indelLen);
                                    //if an insertion is two segments ahead, append part of quality string sequence corresponding to next-next segment
                                    //otherwise (deletion) append first quality score of next segment
                                    q.append(indelOperator == CigarOperator.I ? substr(queryQuality, n + m + mLen, indelLen) : queryQuality.charAt(n + m + mLen));

                                    //add length of next segment to both multoffs and multoffp
                                    //add length of next-next segment to multoffp (for insertion) or to multoffs (for deletion)
                                    multoffs += mLen + (indelOperator == CigarOperator.D ? indelLen : 0);
                                    multoffp += mLen + (indelOperator == CigarOperator.I ? indelLen : 0);

                                    int ci6 = cigar.numCigarElements() > ci + 3 ? cigar.getCigarElement(ci + 3).getLength() : 0;
                                    if (ci6 != 0 && cigar.getCigarElement(ci + 3).getOperator()  == CigarOperator.M) {
                                        Tuple4<Integer, String, String, Integer> tpl = finndOffset(start + multoffs,
                                                n + m + multoffp, ci6, querySequence, queryQuality, ref, cov, conf.vext, conf.goodq);
                                        offset = tpl._1;
                                        ss = tpl._2;
                                        q.append(tpl._3);
                                    }
                                    //skip 2 CIGAR segments
                                    ci += 2;
                                } else {
                                    /*
                                    Condition:
                                    1). CIGAR string has next entry
                                    2). next CIGAR segment is matched
                                     */
                                    if (conf.performLocalRealignment && cigar.numCigarElements() > ci + 1 && cigar.getCigarElement(ci + 1).getOperator() == CigarOperator.M) {
                                        int vsn = 0;
                                        //Loop over next CIGAR segment (no more than conf.vext bases ahead)
                                        for (int vi = 0; vsn <= conf.vext && vi < cigar.getCigarElement(ci + 1).getLength(); vi++) {
                                            //If base is unknown, exit loop
                                            if (querySequence.charAt(n + m + vi) == 'N') {
                                                break;
                                            }
                                            //If base quality is less than conf.goodq, exit loop
                                            if (queryQuality.charAt(n + m + vi) - 33 < conf.goodq) {
                                                break;
                                            }
                                            //If reference sequence has base at this position and it matches read base, update offset
                                            if (ref.containsKey(start + vi)) {
                                                if (isNotEquals(querySequence.charAt(n + m + vi), ref.get(start + vi))) {
                                                    offset = vi + 1;
                                                    vsn = 0;
                                                } else {
                                                    vsn++;
                                                }
                                            }
                                        }
                                        if (offset != 0) { //If next CIGAR segment has good matching base
                                            //Append first offset bases of next segment to ss and q
                                            ss += substr(querySequence, n + m, offset);
                                            q.append(substr(queryQuality, n + m, offset));
                                            //Increase coverage for positions corresponding to first offset bases of next segment
                                            for (int osi = 0; osi < offset; osi++) {
                                                incCnt(cov, start + osi, 1);
                                            }
                                        }
                                    }
                                }

                                //offset should be reset to 0 on every loop, so it is non-zero only if previous part of code was executed
                                //Append '&' and $ss to s if next segment has good matching base
                                if (offset > 0) {
                                    s.append("&").append(ss);
                                }

                                //If start of the segment is within region of interest and the segment does not have unknown bases
                                if (start - 1 >= region.start && start - 1 <= region.end && !s.toString().contains("N")) {
                                    //add '+' + s to insertions at start - 1
                                    incCnt(getOrElse(ins, start - 1, new HashMap<String, Integer>()), "+" + s, 1);
                                    //add insertion to table of variations
                                    Variation hv = getVariation(iHash, start - 1, "+" + s); //variant structure for this insertion
                                    hv.incDir(dir);
                                    //add count
                                    hv.cnt++;
                                    //minimum of positions from start of read and end of read
                                    int tp = p < rlen1 - p ? p + 1 : rlen1 - p;

                                    //mean read quality of the segment
                                    double tmpq = 0;
                                    for (int i = 0; i < q.length(); i++) {
                                        tmpq += q.charAt(i) - 33;
                                    }
                                    tmpq = tmpq / q.length();

                                    //pstd is a flag that is 1 if the variant is covered by at least 2 read segments with different positions
                                    if (hv.pstd == false && hv.pp != 0 && tp != hv.pp) {
                                        hv.pstd = true;
                                    }
                                    //qstd is a flag that is 1 if the variant is covered by at least 2 segment reads with different qualities
                                    if (hv.qstd == false && hv.pq != 0 && tmpq != hv.pq) {
                                        hv.qstd = true;
                                    }
                                    hv.pmean += tp;
                                    hv.qmean += tmpq;
                                    hv.Qmean += mappingQuality;
                                    hv.pp = tp;
                                    hv.pq = tmpq;
                                    if (tmpq >= conf.goodq) {
                                        hv.hicnt++;
                                    } else {
                                        hv.locnt++;
                                    }
                                    hv.nm += nm - nmoff;

                                    // Adjust the reference count for insertion reads

                                    /*
                                    Condition:
                                    1). reference sequence has base for the position
                                    2). hash contains variant structure for the position
                                    3). read base at position n-1 matches reference at start-1
                                     */
                                    if (getVariationMaybe(hash, start - 1, ref.get(start - 1)) != null
                                            && isHasAndEquals(querySequence.charAt(n - 1), ref, start - 1)) {

                                        // subCnt(getVariation(hash, start - 1, ref.get(start - 1 ).toString()), dir, tp, tmpq,
                                        // Qmean, nm, conf);
                                        Variation tv = getVariation(hash, start - 1, String.valueOf(querySequence.charAt(n - 1)));
                                        //Substract count.
                                        subCnt(tv, dir, tp, queryQuality.charAt(n - 1) - 33, mappingQuality, nm, conf);
                                    }
                                    // Adjust count if the insertion is at the edge so that the AF won't > 1
                                    /*
                                    Condition:
                                    1). looking at second segment in CIGAR string
                                    2). first segment is a soft-clipping or a hard-clipping
                                    */
                                    if (ci == 1 && (cigar.getCigarElement(0).getOperator() == CigarOperator.S || cigar.getCigarElement(0).getOperator() == CigarOperator.H)) {
                                        //Add one more variant corresponding to base at start - 1 to hash
                                        Variation ttref = getVariation(hash, start - 1, ref.get(start - 1).toString());
                                        ttref.incDir(dir);
                                        ttref.cnt++;
                                        ttref.pstd = hv.pstd;
                                        ttref.qstd = hv.qstd;
                                        ttref.pmean += tp;
                                        ttref.qmean += tmpq;
                                        ttref.Qmean += mappingQuality;
                                        ttref.pp = tp;
                                        ttref.pq = tmpq;
                                        ttref.nm += nm - nmoff;
                                        incCnt(cov, start - 1, 1);
                                    }
                                }

                                //adjust read position by m (CIGAR segment length) + offset + multoffp
                                n += m + offset + multoffp;
                                p += m + offset + multoffp;
                                //adjust reference position by offset + multoffs
                                start += offset + multoffs;
                                continue;
                            }
                            case D: { //deletion
                                offset = 0;
                                //description string of deleted segment
                                StringBuilder s = new StringBuilder("-").append(m);
                                //sequence to be appended if next segment is matched
                                StringBuilder ss = new StringBuilder();
                                //quality of last base before deletion
                                char q1 = queryQuality.charAt(n - 1);
                                //quality of this segment
                                StringBuilder q = new StringBuilder();

                                // For multiple indels within $VEXT bp
                                //offset for read position if next segment is matched
                                int multoffs = 0;
                                //offset for reference position if next segment is matched
                                int multoffp = 0;
                                int nmoff = 0;

                                /*
                                Condition:
                                1). CIGAR string has next entry
                                2). length of next CIGAR segment is less than conf.vext
                                3). next segment is matched
                                4). CIGAR string has one more entry after next one
                                5). this entry is insertion or deletion
                                 */
                                if (conf.performLocalRealignment && cigar.numCigarElements() > ci + 2
                                        && cigar.getCigarElement(ci + 1).getLength() <= conf.vext
                                        && cigar.getCigarElement(ci + 1).getOperator() == CigarOperator.M
                                        && (cigar.getCigarElement(ci + 2).getOperator() == CigarOperator.I
                                         || cigar.getCigarElement(ci + 2).getOperator() == CigarOperator.D)) {

                                    int mLen = cigar.getCigarElement(ci + 1).getLength();
                                    int indelLen = cigar.getCigarElement(ci + 2).getLength();

                                    //append '#' + next matched segment from read
                                    s.append("#").append(substr(querySequence, n, mLen));
                                    //append quality string of next matched segment from read
                                    q.append(substr(queryQuality, n, mLen));

                                    //if an insertion is two segments ahead, append '^' + part of sequence corresponding to next-next segment
                                    //otherwise (deletion) append '^' + length of a next-next segment
                                    s.append('^').append(cigar.getCigarElement(ci + 2).getOperator() == CigarOperator.I ? substr(querySequence, n + mLen, indelLen) : indelLen);
                                    //same for quality string
                                    q.append(cigar.getCigarElement(ci + 2).getOperator() == CigarOperator.I ? substr(queryQuality, n + mLen, indelLen) : "");

                                    //add length of next segment to both read and reference offsets
                                    //add length of next-next segment to reference position (for insertion) or to read position(for deletion)
                                    multoffs += mLen + (cigar.getCigarElement(ci + 2).getOperator() == CigarOperator.D ? indelLen : 0);
                                    multoffp += mLen + (cigar.getCigarElement(ci + 2).getOperator() == CigarOperator.I ? indelLen : 0);
                                    if (cigar.numCigarElements() > ci + 3
                                            && cigar.getCigarElement(ci + 3).getOperator() == CigarOperator.M) {
                                        int vsn = 0;
                                        int tn = n + multoffp;
                                        int ts = start + multoffs + m;
                                        for (int vi = 0; vsn <= conf.vext && vi < cigar.getCigarElement(ci + 3).getLength(); vi++) {
                                            if (querySequence.charAt(tn + vi) == 'N') {
                                                break;
                                            }
                                            if (queryQuality.charAt(tn + vi) - 33 < conf.goodq) {
                                                break;
                                            }
                                            if (isHasAndEquals('N', ref, ts + vi)) {
                                                break;
                                            }
                                            Character refCh = ref.get(ts + vi);
                                            if (refCh != null) {
                                                if (isNotEquals(querySequence.charAt(tn + vi), refCh)) {
                                                    offset = vi + 1;
                                                    nmoff++;
                                                    vsn = 0;
                                                } else {
                                                    vsn++;
                                                }
                                            }
                                        }
                                        if (offset != 0) {
                                            ss.append(substr(querySequence, tn, offset));
                                            q.append(substr(queryQuality, tn, offset));
                                        }
                                    }
                                    // skip next 2 CIGAR segments
                                    ci += 2;
                                } else if (conf.performLocalRealignment && cigar.numCigarElements() > ci + 1 && cigar.getCigarElement(ci + 1).getOperator() == CigarOperator.I) {
                                    /*
                                    Condition:
                                    1). CIGAR string has next entry
                                    2). next CIGAR segment is an insertion
                                     */

                                    int insLen = cigar.getCigarElement(ci + 1).getLength();
                                    //Append '^' + next segment (inserted)
                                    s.append("^").append(substr(querySequence, n, insLen));
                                    //Append next segement to quality string
                                    q.append(substr(queryQuality, n, insLen));

                                    //Shift reference position by length of next segment
                                    //skip next CIGAR segment
                                    multoffp += insLen;
                                    if (cigar.numCigarElements() > ci + 2
                                            && cigar.getCigarElement(ci + 2).getOperator() == CigarOperator.M) {
                                        int mLen = cigar.getCigarElement(ci + 2).getLength();
                                        int vsn = 0;
                                        int tn = n + multoffp;
                                        int ts = start + m;
                                        for (int vi = 0; vsn <= conf.vext && vi < mLen; vi++) {
                                            char seqCh = querySequence.charAt(tn + vi);
                                            if (seqCh == 'N') {
                                                break;
                                            }
                                            if (queryQuality.charAt(tn + vi) - 33 < conf.goodq) {
                                                break;
                                            }
                                            Character refCh = ref.get(ts + vi);
                                            if (refCh != null) {
                                                if (isEquals('N', refCh)) {
                                                    break;
                                                }
                                                if (isNotEquals(seqCh, refCh)) {
                                                    offset = vi + 1;
                                                    nmoff++;
                                                    vsn = 0;
                                                } else {
                                                    vsn++;
                                                }
                                            }
                                        }
                                        if (offset != 0) {
                                            ss.append(substr(querySequence, tn, offset));
                                            q.append(substr(queryQuality, tn, offset));
                                        }
                                    }
                                    ci += 1;
                                } else {
                                    /*
                                    Condition:
                                    1). CIGAR string has next entry
                                    2). next CIGAR segment is matched
                                     */
                                    if (conf.performLocalRealignment && cigar.numCigarElements() > ci + 1 && cigar.getCigarElement(ci + 1).getOperator() == CigarOperator.M) {
                                        int mLen = cigar.getCigarElement(ci + 1).getLength();
                                        int vsn = 0;
                                        //Loop over next CIGAR segment (no more than conf.vext bases ahead)
                                        for (int vi = 0; vsn <= conf.vext && vi < mLen; vi++) {
                                            char seqCh = querySequence.charAt(n + vi);
                                            //If base is unknown, exit loop
                                            if (seqCh == 'N') {
                                                break;
                                            }
                                            //If base quality is less than $GOODQ, exit loop
                                            if (queryQuality.charAt(n + vi) - 33 < conf.goodq) {
                                                break;
                                            }
                                            //If reference sequence has base at this position and it matches read base, update offset
                                            Character refCh = ref.get(start + m + vi);
                                            if (refCh != null) {
                                                if (isEquals('N', refCh)) {
                                                    break;
                                                }
                                                if (isNotEquals(seqCh, refCh)) {
                                                    offset = vi + 1;
                                                    nmoff++;
                                                    vsn = 0;
                                                } else {
                                                    vsn++;
                                                }
                                            }

                                        }

                                        //If next CIGAR segment has good matching base
                                        if (offset != 0) {
                                            //Append first offset bases of next segment to ss and q
                                            ss.append(substr(querySequence, n, offset));
                                            q.append(substr(queryQuality, n, offset));
                                        }
                                    }
                                }
                                //offset should be reset to 0 on every loop, so it is non-zero only if previous next CIGAR segment is matched
                                // Append '&' and ss to s if next segment has good matching base
                                if (offset > 0) {
                                    s.append("&").append(ss);
                                }

                                //quality of first matched base after deletion
                                //append best of $q1 and $q2
                                if (n + offset >= queryQuality.length()) {
                                    q.append(q1);
                                } else {
                                    char q2 = queryQuality.charAt(n + offset);
                                    q.append(q1 > q2 ? q1 : q2);
                                }

                                //If reference position is inside region of interest
                                if (start >= region.start && start <= region.end) {
                                    //add variant structure for deletion at this position
                                    Variation hv = getVariation(hash, start, s.toString()); //variation structure
                                    //add record for deletion in deletions map
                                    increment(dels5, start, s.toString());
                                    hv.incDir(dir);
                                    //increase count
                                    hv.cnt++;

                                    //minimum of positions from start of read and end of read
                                    int tp = p < rlen1 - p ? p + 1 : rlen1 - p;

                                    //average quality of bases
                                    double tmpq = 0;

                                    for (int i = 0; i < q.length(); i++) {
                                        tmpq += q.charAt(i) - 33;
                                    }

                                    tmpq = tmpq / q.length();

                                    //pstd is a flag that is 1 if the variant is covered by at least 2 read segments with different positions
                                    if (hv.pstd == false && hv.pp != 0 && tp != hv.pp) {
                                        hv.pstd = true;
                                    }

                                    //qstd is a flag that is 1 if the variant is covered by at least 2 segment reads with different qualities
                                    if (hv.qstd == false && hv.pq != 0 && tmpq != hv.pq) {
                                        hv.qstd = true;
                                    }
                                    hv.pmean += tp;
                                    hv.qmean += tmpq;
                                    hv.Qmean += mappingQuality;
                                    hv.pp = tp;
                                    hv.pq = tmpq;
                                    hv.nm += nm - nmoff;
                                    if (tmpq >= conf.goodq) {
                                        hv.hicnt++;
                                    } else {
                                        hv.locnt++;
                                    }

                                    //increase coverage count for reference bases missing from the read
                                    for (int i = 0; i < m; i++) {
                                        incCnt(cov, start + i, 1);
                                    }
                                }

                                //adjust reference position by offset + multoffs
                                start += m + offset + multoffs;

                                //adjust read position by m (CIGAR segment length) + offset + multoffp
                                n += offset + multoffp;
                                p += offset + multoffp;
                                continue;
                            }
                            default:
                                break;
                        }
                        //End branch on CIGAR segment type

                        // Now dealing with matching part
                        int nmoff = 0;
                        int moffset = 0;
                        //Loop over bases of CIGAR segment
                        for (int i = offset; i < m; i++) {
                            //flag to trim reads at opt_T bases from start or end (depending on direction)
                            boolean trim = false;
                            if (conf.trimBasesAfter != 0) {
                                if (dir == false) {
                                    if (n > conf.trimBasesAfter) {
                                        trim = true;
                                    }
                                } else {
                                    if (rlen2 - n > conf.trimBasesAfter) {
                                        trim = true;
                                    }
                                }
                            }

                            //variation string. Initialize to first base of the read sequence
                            final char ch1 = querySequence.charAt(n);
                            String s = String.valueOf(ch1);
                            boolean startWithDelition = false;
                            //skip if base is unknown
                            if (ch1 == 'N') {
                                start++;
                                n++;
                                p++;
                                continue;
                            }

                            //sum of qualities for bases
                            double q = queryQuality.charAt(n) - 33;
                            //number of bases for quality calculation
                            int qbases = 1;
                            //number of bases in insertion for quality calculation
                            int qibases = 0;
                            // for more than one nucleotide mismatch
                            StringBuilder ss = new StringBuilder();
                            // More than one mismatches will only perform when all nucleotides have queryQuality > $GOODQ
                            // Update: Forgo the queryQuality check. Will recover later

                            /*
                            Condition:
                            1). start + 1 is in region of interest
                            2). base index is not more than segment length
                            3). reference sequence has base at start
                            4). base at n in read is not equal to base in reference sequence
                             */
                            while ((start + 1) >= region.start
                                    && (start + 1) <= region.end && (i + 1) < m
                                    && q >= conf.goodq
                                    && isHasAndNotEquals(ref, start, querySequence, n)
                                    && isNotEquals('N', ref.get(start))) {

                                //Break if base is unknown in the read
                                char nuc = querySequence.charAt(n + 1);
                                if (nuc == 'N') {
                                    break;
                                }
                                if (isHasAndEquals('N', ref, start + 1)) {
                                    break;
                                }

                                //Condition: base at n + 1 does not match reference base at start + 1
                                if (isNotEquals(ref.get(start + 1), nuc)) {

                                    //append the base from read
                                    ss.append(nuc);
                                    //add quality to total sum
                                    q += queryQuality.charAt(n + 1) - 33;
                                    //increase number of bases
                                    qbases++;
                                    //shift read position by 1
                                    n++;
                                    p++;
                                    i++;
                                    //shift reference position by 1
                                    start++;
                                    nmoff++;
                                } else { //if bases match, exit loop
                                    break;
                                }
                            }

                            //If multi-base mismatch is found, append it to s after '&'
                            if (ss.length() > 0) {
                                s += "&" + ss;
                            }
                            int ddlen = 0;

                            /*
                            Condition:
                            1). index $i is no farther than $VEXT from end of CIGAR segment
                            2). CIGAR string contains next entry
                            3). next CIGAR entry is a deletion
                            4). reference sequence contains a base at $start
                            5). either multi-nucleotide mismatch is found or read base at $n does not match reference base
                            6). read base at $n has good quality
                             */
                            if (conf.performLocalRealignment && m - i <= conf.vext
                                    && cigar.numCigarElements() > ci + 1 && cigar.getCigarElement(ci + 1).getOperator() == CigarOperator.D
                                    && ref.containsKey(start)
                                    && (ss.length() > 0 || isNotEquals(querySequence.charAt(n), ref.get(start)))
                                    && queryQuality.charAt(n) - 33 > conf.goodq) {

                                //loop until end of CIGAR segments
                                while (i + 1 < m) {
                                    //append next base to s and add its quality to q
                                    s += querySequence.charAt(n + 1);
                                    q += queryQuality.charAt(n + 1) - 33;
                                    //increase number of bases
                                    qbases++;

                                    //shift read and reference indices by 1
                                    i++;
                                    n++;
                                    p++;
                                    start++;
                                }

                                //remove '&' delimiter from s
                                s = s.replaceFirst("&", "");
                                //prepend s with deletion of length of next CIGAR segment + '&' delimiter
                                s = "-" + cigar.getCigarElement(ci + 1).getLength() + "&" + s;
                                startWithDelition = true;
                                ddlen = cigar.getCigarElement(ci + 1).getLength();
                                ci += 1;

                                //If CIGAR has insertion two segments ahead
                                if (cigar.numCigarElements() > ci + 1 && cigar.getCigarElement(ci + 1).getOperator() == CigarOperator.I) {
                                    //append '^' + next-next segment sequence
                                    s += "^" + substr(querySequence, n + 1, cigar.getCigarElement(ci + 1).getLength());

                                    //Loop over next-next segment
                                    int nextLen = cigar.getCigarElement(ci + 1).getLength();
                                    for (int qi = 1; qi <= nextLen; qi++) {
                                        //add base quality to total quality
                                        q += queryQuality.charAt(n + 1 + qi) - 33;
                                        //increase number of insertion bases
                                        qibases++;
                                    }

                                    //adjust read position by length of next-next segment
                                    n += nextLen;
                                    p += nextLen;
                                    ci += 1;
                                }
                                if (cigar.numCigarElements() > ci + 1 && cigar.getCigarElement(ci + 1).getOperator() == CigarOperator.M) {
                                    Tuple4<Integer, String, String, Integer> tpl =
                                            finndOffset(start + ddlen + 1, n + 1, cigar.getCigarElement(ci + 1).getLength(), querySequence, queryQuality, ref, cov, conf.vext, conf.goodq);
                                    int toffset = tpl._1;
                                    if (toffset != 0) {
                                        moffset = toffset;
                                        nmoff += tpl._4;
                                        s += "&" + tpl._2;
                                        String tq = tpl._3;
                                        for (int qi = 0; qi < tq.length(); qi++) {
                                            q += tq.charAt(qi) - 33;
                                            qibases++;
                                        }
                                    }
                                }
                            }
                            if (trim == false) {
                                //If start - qbases + 1 is in region of interest
                                final int pos = start - qbases + 1;
                                if (pos >= region.start && pos <= region.end) {
                                    //add variation record for $s
                                    Variation hv = getVariation(hash, pos, s); //reference to variant structure
                                    hv.incDir(dir);

                                    if(isBEGIN_ATGC_AMP_ATGCs_END(s)) {
                                        //if s is one base followed by '&' and one or more bases
                                        //add variant record for s to mnp
                                        increment(mnp, pos, s);
                                    }

                                    //increment count
                                    ++hv.cnt;

                                    //minimum of positions from start of read and end of read
                                    int tp = p < rlen1 - p ? p + 1 : rlen1 - p;

                                    //average quality of bases in the variation
                                    q = q / (qbases + qibases);

                                    //pstd is a flag that is 1 if the variant is covered by at least 2 read segments with different positions
                                    if (hv.pstd == false && hv.pp != 0 && tp != hv.pp) {
                                        hv.pstd = true;
                                    }

                                    //qstd is a flag that is 1 if the variant is covered by at least 2 segment reads with different qualities
                                    if (hv.qstd == false && hv.pq != 0 && q != hv.pq) {
                                        hv.qstd = true;
                                    }
                                    hv.pmean += tp;
                                    hv.qmean += q;
                                    hv.Qmean += mappingQuality;
                                    hv.pp = tp;
                                    hv.pq = q;
                                    hv.nm += nm - nmoff;
                                    if (q >= conf.goodq) {
                                        hv.hicnt++;
                                    } else {
                                        hv.locnt++;
                                    }

                                    //increase coverage for bases covered by the variation
                                    for (int qi = 1; qi <= qbases; qi++) {
                                        incCnt(cov, start - qi + 1, 1);
                                    }

                                    //If variation starts with a deletion ('-' character)
                                    if (startWithDelition) {
                                        //add variation to deletions map
                                        increment(dels5, pos, s);

                                       //increase coverage for next CIGAR segment
                                        for (int qi = 1; qi < ddlen; qi++) {
                                            incCnt(cov, start + qi, 1);
                                        }
                                    }
                                }
                            }

                            //If variation starts with a deletion ('-' character)
                            if (startWithDelition) {
                                start += ddlen;
                            }

                            //Shift reference position by 1 if CIGAR segment is not insertion
                            if (operator != CigarOperator.I) {
                                start++;
                            }
                            //Shift read position by 1 if CIGAR segment is not deletion
                            if (operator != CigarOperator.D) {
                                n++;
                                p++;
                            }
                        }
                        if (moffset != 0) {
                            offset = moffset;
                            n += moffset;
                            start += moffset;
                            p += moffset;
                        }
                        if (start > region.end) { //end if reference position is outside region of interest
                            break;
                        }
                    }
                }
            }
        }

        if (conf.outputSplicing) {
            for (Entry<String, int[]> entry : spliceCnt.entrySet()) {
                System.out.printf("%s\t%s\t%s\t%s\n", sample, region.chr, entry.getKey(), entry.getValue()[0]);
            }
            System.exit(0);
        }

        if (conf.performLocalRealignment) {
            if (conf.y)
                System.err.println("Start Realigndel");
            realigndel(hash, dels5, cov, sclip5, sclip3, ref, region.chr, chrs, rlen, bams, conf);
            if (conf.y)
                System.err.println("Start Realignins");
            realignins(hash, iHash, ins, cov, sclip5, sclip3, ref, region.chr, chrs, conf);
            if (conf.y)
                System.err.println("Start Realignlgdel");
            realignlgdel(hash, cov, sclip5, sclip3, ref, region.chr, chrs, rlen, bams, conf);
            if (conf.y)
                System.err.println("Start Realignlgins");
            realignlgins(hash, iHash, cov, sclip5, sclip3, ref, region.chr, chrs, rlen, bams, conf);
            if (conf.y)
                System.err.println("Start Realignlgins30");
            realignlgins30(hash, iHash, cov, sclip5, sclip3, ref, region.chr, chrs, rlen, bams, conf);
        }

        adjMNP(hash, mnp, cov, ref, sclip3, sclip5, conf);

        return tuple(hash, iHash, cov, rlen);
    }

    private static CigarOperator getCigarOperator(Cigar cigar, int ci) {
        CigarOperator operator =  cigar.getCigarElement(ci).getOperator();
        // Treat insertions at the edge as soft-clipping
        if ((ci == 0 || ci == cigar.numCigarElements() - 1) && operator == CigarOperator.I) {
            operator = CigarOperator.S;
        }
        return operator;
    }

    private static int getAlignedLenght(Cigar readCigar) {
        int lenght = 0;
        for (CigarElement element : readCigar.getCigarElements()) {
            if (element.getOperator() == CigarOperator.M || element.getOperator() == CigarOperator.D) {
                lenght += element.getLength();
            }
        }
        return lenght;
    }

    private static int getSoftClippedLenght(Cigar readCigar) {
        int lenght = 0;
        for (CigarElement element : readCigar.getCigarElements()) {
            CigarOperator operator = element.getOperator();
            if (operator == CigarOperator.M || operator == CigarOperator.I || operator == CigarOperator.S) {
                lenght += element.getLength();
            }
        }
        return lenght;
    }

    private static int getMatchInsertionLenght(Cigar readCigar) {
        int lenght = 0;
        for (CigarElement element : readCigar.getCigarElements()) {
            if (element.getOperator() == CigarOperator.M || element.getOperator() == CigarOperator.I) {
                lenght += element.getLength();
            }
        }
        return lenght;
    }

    private static int getInsertionDeletionLenght(Cigar readCigar) {
        int lenght = 0;
        for (CigarElement element : readCigar.getCigarElements()) {
            if (element.getOperator() == CigarOperator.I || element.getOperator() == CigarOperator.D) {
                lenght += element.getLength();
            }
        }
        return lenght;
    }

    private static String getMrnm(SAMRecord record) {
        if (record.getMateReferenceName() == null) {
            return "*";
        }

        if (record.getReferenceName().equals(record.getMateReferenceName())) {
            return "=";
        }
        return record.getMateReferenceName();
    }

    /**
     * Read BAM files and create variant structure
     * @param region region of interest
     * @param bam BAM file names (':' delimiter)
     * @param ref part of reference sequence (key - position, value - base)
     * @param chrs map of chromosome lengths
     * @param SPLICE set of strings representing spliced regions
     * @param ampliconBasedCalling string of maximum_distance:minimum_overlap for amplicon based calling
     * @param Rlen maximum read length
     * @param conf VarDict configuration
     * @return Tuple of (maxmimum read length, variant structure)
     * @throws IOException
     */
    static Tuple2<Integer, Map<Integer, Vars>> toVars(Region region, String bam, Map<Integer, Character> ref,
            Map<String, Integer> chrs, String sample, Set<String> SPLICE, String ampliconBasedCalling, int Rlen, Configuration conf) throws IOException {

        Tuple4<Map<Integer, Map<String, Variation>>, Map<Integer, Map<String, Variation>>, Map<Integer, Integer>, Integer> parseTpl =
                parseSAM(region, bam, chrs, sample, SPLICE, ampliconBasedCalling, Rlen, ref, conf);

        Map<Integer, Map<String, Variation>> hash = parseTpl._1;
        Map<Integer, Map<String, Variation>> iHash = parseTpl._2;
        Map<Integer, Integer> cov = parseTpl._3;
        Rlen = parseTpl._4;

        //the variant structure
        Map<Integer, Vars> vars = new HashMap<>();
        //Loop over positions
        for (Entry<Integer, Map<String, Variation>> entH : hash.entrySet()) {
            final int p = entH.getKey();
            Map<String, Variation> v = entH.getValue();

            if(v.isEmpty()) {
                continue;
            }

            //Skip if position is outside region of interest
            if (p < region.start || p > region.end) {
                continue;
            }

            //skip position if it has no coverage
            if (!cov.containsKey(p) || cov.get(p) == 0) {
                continue;
            }

            Set<String> vk = new HashSet<String>(v.keySet());
            if (iHash.containsKey(p)) {
                vk.add("I");
            }
            if (vk.size() == 1 && ref.containsKey(p) && vk.contains(ref.get(p).toString())) {
                if (!conf.doPileup && !conf.bam.hasBam2() && conf.ampliconBasedCalling == null) { // ignore if only reference were seen and no pileup to avoid computation
                    continue;
                }
            }

            //total position coverage
            final int tcov = cov.get(p);
            if (tcov == 0) { // ignore when there's no coverage
                continue;
            }

            //position coverage by high-quality reads
            final int hicov = calcHicov(v,iHash.get(p));

            //array of all variants for the position
            List<Variant> var = new ArrayList<>();
            //temporary array used for debugging
            List<String> tmp = new ArrayList<>();
            List<String> keys = new ArrayList<>(v.keySet());
            Collections.sort(keys);

            //Loop over all varints found for the position except insertions
            for (String n : keys) {
                // for (Entry<String, Variation> entV : v.entrySet()) {
                // String n = entV.getKey();
                // Variation cnt = entV.getValue();
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
                double vqual = round(cnt.qmean / cnt.cnt, 1); // base quality
                //mean mapping quality for variant
                double mq = cnt.Qmean / (double)cnt.cnt; // mapping quality
                //number of high-quality reads for variant
                int hicnt = cnt.hicnt;
                //number of low-quality reads for variant
                int locnt = cnt.locnt;
                /**
                 * Condition:
                 # 1). cnt.cnt > tcov                         - variant count is more than position coverage
                 # 2). cnt.cnt - tcov < cnt.extracnt          - variant count is no more than poistion coverage + extracnt
                 */
                int ttcov = tcov;
                if (cnt.cnt > tcov && cnt.cnt - tcov < cnt.extracnt) { //adjust position coverage if condition holds
                    ttcov = cnt.cnt;
                }

                //create variant record
                Variant tvref = new Variant();
                tvref.n = n;
                tvref.cov = cnt.cnt;
                tvref.fwd = fwd;
                tvref.rev = rev;
                tvref.bias = String.valueOf(bias);
                tvref.freq = cnt.cnt / (double)ttcov;
                tvref.pmean = round(cnt.pmean / (double)cnt.cnt, 1);
                tvref.pstd = cnt.pstd;
                tvref.qual = vqual;
                tvref.qstd = cnt.qstd;
                tvref.mapq = mq;
                tvref.qratio = hicnt / (locnt != 0 ? locnt : 0.5d);
                tvref.hifreq = hicov > 0 ? hicnt / (double)hicov : 0;
                tvref.extrafreq = cnt.extracnt != 0 ? cnt.extracnt / (double)ttcov : 0;
                tvref.shift3 = 0;
                tvref.msi = 0;
                tvref.nm = cnt.nm / (double)cnt.cnt;
                tvref.hicnt = hicnt;
                tvref.hicov = hicov;

                //append variant record
                var.add(tvref);
                if (conf.debug) { //debugging output
                    tmp.add(n
                            + ":" + (fwd + rev)
                            + ":F-" + fwd
                            + ":R-" + rev
                            + ":" + format("%.3f", tvref.freq)
                            + ":" + tvref.bias
                            + ":" + tvref.pmean
                            + ":" + tvref.pstd
                            + ":" + vqual
                            + ":" + tvref.qstd
                            + ":" + format("%.3f", tvref.hifreq)
                            + ":" + tvref.mapq
                            + ":" + tvref.qratio);
                }

            }

            //Handle insertions separately
            Map<String, Variation> iv = iHash.get(p);
            if (iv != null) {
                List<String> ikeys = new ArrayList<>(iv.keySet());
                Collections.sort(ikeys);
                // for (Entry<String, Variation> entV : iv.entrySet()) {
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
                    double vqual = round(cnt.qmean / cnt.cnt, 1); // base quality
                    //mean mapping quality for variant
                    double mq = cnt.Qmean / (double)cnt.cnt; // mapping quality
                    //number of high-quality reads for variant
                    int hicnt = cnt.hicnt;
                    //number of low-quality reads for variant
                    int locnt = cnt.locnt;

//                    hicov += hicnt;

                    //adjust position coverage if variant count is more than position coverage and no more than poistion coverage + extracnt
                    int ttcov = tcov;
                    if (cnt.cnt > tcov && cnt.cnt - tcov < cnt.extracnt) {
                        ttcov = cnt.cnt;
                    }

                    Variant tvref = new Variant();
                    tvref.n = n;
                    tvref.cov = cnt.cnt;
                    tvref.fwd = fwd;
                    tvref.rev = rev;
                    tvref.bias = String.valueOf(bias);
                    tvref.freq = cnt.cnt / (double)ttcov;
                    tvref.pmean = round(cnt.pmean / (double)cnt.cnt, 1);
                    tvref.pstd = cnt.pstd;
                    tvref.qual = vqual;
                    tvref.qstd = cnt.qstd;
                    tvref.mapq = mq;
                    tvref.qratio = hicnt / (locnt != 0 ? locnt : 0.5d);
                    tvref.hifreq = hicov > 0 ? hicnt / (double)hicov : 0;
                    tvref.extrafreq = cnt.extracnt != 0 ? cnt.extracnt / (double)ttcov : 0;
                    tvref.shift3 = 0;
                    tvref.msi = 0;
                    tvref.nm = cnt.nm / (double)cnt.cnt;
                    tvref.hicnt = hicnt;
                    tvref.hicov = hicov;

                    var.add(tvref);
                    if (conf.debug) {
                        tmp.add("I" + n
                                + ":" + (fwd + rev)
                                + ":F-" + fwd
                                + ":R-" + rev
                                + ":" + format("%.3f", tvref.freq)
                                + ":" + tvref.bias
                                + ":" + tvref.pmean
                                + ":" + tvref.pstd
                                + ":" + vqual
                                + ":" + tvref.qstd
                                + ":" + format("%.3f", tvref.hifreq)
                                + ":" + tvref.mapq
                                + ":" + tvref.qratio);
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
            if (varsAtp.ref != null) {
                if (varsAtp.ref.freq >= conf.freq) {
                    genotype1 = varsAtp.ref.n;
                } else if (varsAtp.var.size() > 0) {
                    genotype1 = varsAtp.var.get(0).n;
                }
                rfc = varsAtp.ref.fwd;
                rrc = varsAtp.ref.rev;
            } else if (varsAtp.var.size() > 0) {
                genotype1 = varsAtp.var.get(0).n;
            }
            if (genotype1.startsWith("+")) {
                genotype1 = "+" + (genotype1.length() - 1);
            }
            //description string for any other variant
            String genotype2 = "";

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
                    Matcher matcher = START_DIG.matcher(vn);
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
                        if (!vn.contains("&") && !vn.contains("#")) {
                            //variant description string without first symbol '+'
                            String tseq1 = vn.substring(1);
                            //left 50 bases in reference sequence
                            String leftseq = joinRef(ref, p - 50 > 1 ? p - 50 : 1, p); // left 10 nt
                            int x = getOrElse(chrs, region.chr, 0);
                            //right 70 bases in reference sequence
                            String tseq2 = joinRef(ref, p + 1, (p + 70 > x ? x : p + 70));

                            Tuple3<Double, Integer, String> tpl = findMSI(tseq1, tseq2, leftseq);
                            msi = tpl._1;
                            shift3 = tpl._2;
                            msint = tpl._3;

                            //Try to adjust for microsatellite instability
                            tpl = findMSI(leftseq, tseq2, null);
                            double tmsi = tpl._1;
                            int tshift3 = tpl._2;
                            String tmsint = tpl._3;
                            if (msi < tmsi) {
                                msi = tmsi;
                                shift3 = tshift3;
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
                    } else if (vn.startsWith("-")) { //deletion variant
                        //variant allele is in the record
                        //remove '-' and number from beginning of variant string
                        varallele = vn.replaceFirst("^-\\d+", "");

                        //left 70 bases in reference sequence
                        String leftseq = joinRef(ref, (p - 70 > 1 ? p - 70 : 1), p - 1); // left 10 nt
                        int chr0 = getOrElse(chrs, region.chr, 0);
                        //right 70 + dellen bases in reference sequence
                        String tseq = joinRef(ref, p, p + dellen + 70 > chr0 ? chr0 : p + dellen + 70);

                        //Try to adjust for microsatellite instability
                        Tuple3<Double, Integer, String> tpl = findMSI(substr(tseq, 0, dellen), substr(tseq, dellen), leftseq);
                        msi = tpl._1;
                        shift3 = tpl._2;
                        msint = tpl._3;

                        tpl = findMSI(leftseq, substr(tseq, dellen), leftseq);
                        double tmsi = tpl._1;
                        int tshift3 = tpl._2;
                        String tmsint = tpl._3;
                        if (msi < tmsi) {
                            msi = tmsi;
                            shift3 = tshift3;
                            msint = tmsint;
                        }
                        if (msi <= shift3 / (double)dellen) {
                            msi = shift3 / (double)dellen;
                        }

                        //If no matched sequence or indel follows the variant
                        if (!vn.contains("&") && !vn.contains("#") && !vn.contains("^")) {
                            //Shift position to 3' if -3 option is set
                            if (conf.moveIndelsTo3) {
                                sp += shift3;
                            }
                            //variant allel is 1 base from reference string preceding p
                            varallele = ref.containsKey(p - 1) ? ref.get(p - 1).toString() : "";
                            //prepend same base to reference allele
                            refallele = varallele;
                            sp--;
                        }
                        //append dellen bases from reference string to reference allele
                        refallele += joinRef(ref, p, p + dellen - 1);
                    } else { //Not insertion/deletion variant. SNP or MNP

                        //Find MSI adjustment
                        String tseq1 = joinRef(ref, p - 30 > 1 ? p - 30 : 1, p + 1);
                        int chr0 = getOrElse(chrs, region.chr, 0);
                        String tseq2 = joinRef(ref, p + 2, p + 70 > chr0 ? chr0 : p + 70);

                        Tuple3<Double, Integer, String> tpl = findMSI(tseq1, tseq2, null);
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
                    genotype = genotype.replace("&", "").replace("#", "").replace("^", "i");
                    //convert extrafreq, freq, hifreq, msi fields to strings
                    vref.extrafreq = vref.extrafreq;
                    vref.freq = round(vref.freq, 4);
                    vref.hifreq = round(vref.hifreq, 4);
                    vref.msi = msi;
                    vref.msint = msint.length();
                    vref.shift3 = shift3;
                    vref.sp = sp;
                    final Matcher pnMtch = PLUS_NUMBER.matcher(genotype);
                    if (pnMtch.find()) {
                        vref.ep = sp + toInt(pnMtch.group(1));
                    } else {
                        vref.ep = ep;
                    }
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
            } else {
                if (varsAtp.ref == null)
                    varsAtp.ref = new Variant();
            }
            if (varsAtp.ref != null) {
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
                vref.hifreq = vref.hifreq;
                String r = ref.containsKey(p) ? ref.get(p).toString() : "";
                //both refallele and varallele are 1 base from reference string
                vref.refallele = r;
                vref.varallele = r;
                vref.genotype = r + "/" + r;
                vref.leftseq = "";
                vref.rightseq = "";
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

        }

        return tuple(Rlen, vars);
    }

    private static int calcHicov(Map<String, Variation> v, Map<String, Variation> iv) {
        int hicov = 0;
        for (Variation vr : v.values()) {
            hicov += vr.hicnt;
        }
        if (iv != null) {
            for (Variation vr : iv.values()) {
                hicov += vr.hicnt;
            }
        }
        return hicov;
    }

    /**
     * Get part of reference sequence
     * @param region region of interest
     * @param chrs map of chromosome lengths
     * @param fasta file namo of reference genome in FASTA format
     * @param numberNucleotideToExtend number of base pairs to extend around region of interest
     * @return reference sequence map: key - position, value - base
     * @throws IOException
     */
    static Map<Integer, Character> getREF(Region region, Map<String, Integer> chrs, String fasta, int numberNucleotideToExtend) throws IOException {
        Map<Integer, Character> ref = new HashMap<Integer, Character>();

        int s_start = region.start - numberNucleotideToExtend - 700 < 1 ? 1 : region.start - numberNucleotideToExtend - 700;
        int len = chrs.containsKey(region.chr) ? chrs.get(region.chr) : 0;
        int s_end = region.end + numberNucleotideToExtend + 700 > len ?
                len : region.end + numberNucleotideToExtend + 700;

        String[] subSeq = retriveSubSeq(fasta, region.chr, s_start, s_end);
//        String header = subSeq[0];
        String exon = subSeq[1];
        for (int i = s_start; i < s_start + exon.length(); i++) { // TODO why '<=' ?
            ref.put(i, Character.toUpperCase(exon.charAt(i - s_start)));
        }

        return ref;
    }

    /**
     * Find microsatellite instability
     * Tandemly repeated short sequence motifs ranging from 1– 6(8 in our case) base pairs are called microsatellites.
     * Other frequently used terms for these DNA regions are simple sequences or short tandem repeats (STRs)
     * @param tseq1
     * @param tseq2
     * @param left
     * @return Tuple of (MSI count, No. of bases to be shifted to 3 prime for deletions due to alternative alignment, MicroSattelite unit length in base pairs)
     */
    static Tuple3<Double, Integer, String> findMSI(String tseq1, String tseq2, String left) {

        //number of nucleotides in microsattelite
        int nmsi = 1;
        //number of bases to be shifted to 3 prime
        int shift3 = 0;
        String maxmsi = "";
        double msicnt = 0;
        while (nmsi <= tseq1.length() && nmsi <= 8) {
            //microsattelite nucleotide sequence; trim nucleotide(s) from the end
            String msint = substr(tseq1, -nmsi);
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

    /**
     * Variants for position
     */
    static class Vars {
        /**
         * Reference variant
         */
        Variant ref;

        /**
         * List of all variants except reference variant
         */
        List<Variant> var = new ArrayList<>();

        /**
         * Map of all variants except reference variant.
         * Key - variant description string, value - variant
         */
        Map<String, Variant> varn = new HashMap<>();
    }

    /**
     * Calculate strand bias flag
     * @param fwd Variant count for forward strand
     * @param rev Variant count for reverse strand
     * @param bias
     * @param minb minimum reads for bias
     * @return 0 - small total count, only one of strands, 1 - strand bias, 2 - no strand bias
     */
    static int strandBias(int fwd, int rev, double bias, int minb) {

        if (fwd + rev <= 12) { // using p=0.01, because prop.test(1,12) = 0.01
            return fwd * rev > 0 ? 2 : 0;
        }

        return (fwd / (double)(fwd + rev) >= bias && rev / (double)(fwd + rev) >= bias && fwd >= minb && rev >= minb) ? 2 : 1;
    }

    /**
     * this will try to realign large insertions (typically larger than 30bp)
     * @param hash map of non-insertion variants produced by parseSAM method
     * @param iHash map of insertion variants produced by parseSAM method
     * @param cov coverage
     * @param sclip5 5' softclips
     * @param sclip3 3' softclips
     * @param ref map of reference bases
     * @param chr chromosome name
     * @param chrs map of chromosome lengths
     * @param rlen max read length
     * @param bams BAM file list
     * @param conf configuration
     * @throws IOException
     */
    static void realignlgins30(Map<Integer, Map<String, Variation>> hash,
            Map<Integer, Map<String, Variation>> iHash,
            Map<Integer, Integer> cov,
            Map<Integer, Sclip> sclip5,
            Map<Integer, Sclip> sclip3,
            Map<Integer, Character> ref,
            String chr,
            Map<String, Integer> chrs,
            int rlen,
            String[] bams,
            Configuration conf) throws IOException {

        List<Tuple3<Integer, Sclip, Integer>> tmp5 = new ArrayList<>();
        for (Entry<Integer, Sclip> ent5 : sclip5.entrySet()) {
            //ent5: (position, soft clipping 5' variant)
            tmp5.add(tuple(ent5.getKey(), ent5.getValue(), ent5.getValue().cnt));
        }
        Collections.sort(tmp5, COMP3);

        List<Tuple3<Integer, Sclip, Integer>> tmp3 = new ArrayList<>();
        for (Entry<Integer, Sclip> ent3 : sclip3.entrySet()) {
            tmp3.add(tuple(ent3.getKey(), ent3.getValue(), ent3.getValue().cnt));
        }
        Collections.sort(tmp3, COMP3);

        for (Tuple3<Integer, Sclip, Integer> t5 : tmp5) {
            final int p5 = t5._1;
            final Sclip sc5v = t5._2;
            final int cnt5 = t5._3;
            if (sc5v.used) {
                continue;
            }

            for (Tuple3<Integer, Sclip, Integer> t3 : tmp3) {
                final int p3 = t3._1;
                final Sclip sc3v = t3._2;
                final int cnt3 = t3._3;
                if (sc3v.used) {
                    continue;
                }
                if (p5 - p3 > rlen / 1.5) {
                    continue;
                }
                if (p3 - p5 > rlen - 10) { // if they're too far away, don't even try
                    continue;
                }
                final String seq5 = findconseq(sc5v, conf.y);
                final String seq3 = findconseq(sc3v, conf.y);
                //next until at least one of consensus sequences has length > 10
                if (seq5.length() <= 10 || seq3.length() <= 10) {
                    continue;
                }
                if (conf.y) {
                    System.err.printf("  Working lgins30: %s %s 3: %s %s 5: %s %s\n",
                            p3, p5, seq3, cnt3, new StringBuilder(seq5).reverse(), cnt5);
                }
                Tuple3<Integer, Integer, Integer> tpl = find35match(seq5, seq3, p5, p3, ref);
                int bp5 = tpl._1;
                int bp3 = tpl._2;
                //length of match
                int score = tpl._3;

                if (score == 0) {
                    continue;
                }
                //add matched part of seq3 (left clip)
                String ins = bp3 > 1 ? substr(seq3, 0, -bp3 + 1) : seq3;
                if (bp5 > 0) { //add not matched part of seq5 (right clip)
                    ins += new StringBuilder(substr(seq5, 0, bp5)).reverse();
                }
                if (islowcomplexseq(ins)) {
                    if (conf.y) {
                        System.err.println("  Discard low complex insertion found " + ins + ".");
                    }
                    continue;
                }
                int bi = 0;
                Variation vref;
                if (conf.y) {
                    System.err.printf("  Found candidate lgins30: %s %s %s\n", p3, p5, ins);
                }
                if (p5 > p3) {
                    if (seq3.length() > ins.length()
                            && !ismatch(substr(seq3, ins.length()), joinRef(ref, p5, p5 + seq3.length() - ins.length() + 2), 1, conf.y)) {
                        continue;
                    }
                    if (seq5.length() > ins.length()
                            && !ismatch(substr(seq5, ins.length()), joinRef(ref, p3 - seq5.length() - ins.length() - 2, p3 - 1), -1, conf.y)) {
                        continue;
                    }
                    if (conf.y) {
                        System.err.printf("  Found lgins30 complex: %s %s %s %s\n", p3, p5, ins.length(), ins);
                    }
                    String tmp = joinRef(ref, p3, p5 - 1);
                    if (tmp.length() > ins.length()) { // deletion is longer
                        ins = (p3 - p5) + "^" + ins;
                        bi = p3;
                        vref = getVariation(hash, p3, ins);
                    } else if (tmp.length() < ins.length()) {
//                        int p3s = p3 + tmp.length();
//                        int p3e = p3s + seq3.length() - ins.length() + 2;
                        ins = substr(ins, 0, ins.length() - tmp.length()) + "&" + substr(ins, p3 - p5);
                        ins = "+" + ins;
                        bi = p3 - 1;
                        vref = getVariation(iHash, bi, ins);
                    } else { // long MNP
                        ins = "-" + ins.length() + "^" + ins;
                        bi = p3;
                        vref = getVariation(hash, p3, ins);
                    }
                } else {
                    if (seq3.length() > ins.length()
                            && !ismatch(substr(seq3, ins.length()), joinRef(ref, p5, p5 + seq3.length() - ins.length() + 2), 1, conf.y)) {
                        continue;
                    }
                    if (seq5.length() > ins.length()
                            && !ismatch(substr(seq5, ins.length()), joinRef(ref, p3 - (seq5.length() - ins.length()) - 2, p3 - 1), -1, conf.y)) {
                        continue;
                    }
                    String tmp = ins.length() > p3 - p5 ? joinRef(ref, p5, p3)
                            : joinRef(ref, p5, p5 + (p3 - p5 - ins.length()) / 2); // Tandem duplication
                    if (conf.y) {
                        System.err.printf("Found lgins30: %s %s %s %s + %s\n", p3, p5, ins.length(), tmp, ins);
                    }
                    ins = "+" + tmp + ins;
                    bi = p5 - 1;
                    vref = getVariation(iHash, bi, ins);
                }
                sc3v.used = true;
                sc5v.used = true;
                vref.pstd = true;
                vref.qstd = true;
                incCnt(cov, bi, sc5v.cnt);
                if (conf.y) {
                    System.err.printf(" lgins30 Found: '%s' %s %s %s\n", ins, bi, bp3, bp5);
                }

                if (ins.startsWith("+")) {
                    Variation mvref = getVariationMaybe(hash, bi, ref.get(bi));
                    adjCnt(vref, sc3v, mvref, conf);
                    adjCnt(vref, sc5v, conf);
                    if (bams != null && bams.length > 0
                            && p3 - p5 >= 5 && p3 - p5 > rlen - 10
                            && mvref != null && mvref.cnt != 0
                            && vref.cnt > 2 * mvref.cnt
                            && noPassingReads(chr, p5, p3, bams, conf)) {
                        adjCnt(vref, mvref, mvref, conf);
                    }
                    Map<Integer, Map<String, Integer>> tins = new HashMap<>();
                    Map<String, Integer> map = new HashMap<>();
                    map.put(ins, vref.cnt);
                    tins.put(bi, map);
                    realignins(hash, iHash, tins, cov, sclip5, sclip3, ref, chr, chrs, conf);
                } else if (ins.startsWith("-")) {
                    adjCnt(vref, sc3v, getVariationMaybe(hash, bi, ref.get(bi)), conf);
                    adjCnt(vref, sc5v, conf);
                    Map<Integer, Map<String, Integer>> tdel = new HashMap<>();
                    Map<String, Integer> map = new HashMap<>();
                    map.put(ins, vref.cnt);
                    tdel.put(bi, map);
                    realigndel(hash, tdel, cov, sclip5, sclip3, ref, chr, chrs, rlen, bams, conf);
                } else {
                    adjCnt(vref, sc3v, conf);
                    adjCnt(vref, sc5v, conf);
                }
                break;
            }

        }
        if (conf.y) {
            System.err.println("Done: lgins30\n");
        }
    }

    /**
     * test whether the two soft-clipped reads match
     * Returns the breakpoints in 5 and 3 prime soft-clipped reads
     * @param seq5 consensus sequence 5' strand
     * @param seq3 consensus sequence 3' strand
     * @param p5 position 5'
     * @param p3 position 3'
     * @param ref map of reference bases (not used)
     * @return Tuple of (start position of 5' strand, start position of 3' strand, max match length)
     */
    static Tuple3<Integer, Integer, Integer> find35match(String seq5, String seq3, int p5, int p3, Map<Integer, Character> ref) {
        final int longmm = 3;
        int max = 0;
        int b3 = 0;
        int b5 = 0;
        for (int i = 0; i < seq5.length() - 8; i++) {
            for (int j = 1; j < seq3.length() - 8; j++) {
                int nm = 0;
                int n = 0;
                while (n + j <= seq3.length() && i + n <= seq5.length()) {
                    if (!substr(seq3, -j - n, 1).equals(substr(seq5, i + n, 1))) {
                        nm++;
                    }
                    if (nm > longmm) {
                        break;
                    }
                    n++;
                }
                if (n - nm > max
                        && n - nm > 8
                        && nm / (double)n < 0.1d
                        && (n + j >= seq3.length() || i + n >= seq5.length())) {

                    max = n - nm;
                    b3 = j;
                    b5 = i;
                    return tuple(b5, b3, max);
                }

            }
        }
        return tuple(b5, b3, max);
    }

    /**
     * Realign large insertions that are not present in alignment
     * @param hash map of non-insertion variants produced by parseSAM method
     * @param iHash map of insertion variants produced by parseSAM method
     * @param cov coverage
     * @param sclip5 5' softclips
     * @param sclip3 3' softclips
     * @param ref map of reference bases
     * @param chr chromosome name
     * @param chrs map of chromosome lengths
     * @param rlen read length
     * @param bams BAM file list
     * @param conf configuration
     * @throws IOException
     */
    static void realignlgins(Map<Integer, Map<String, Variation>> hash,
            Map<Integer, Map<String, Variation>> iHash,
            Map<Integer, Integer> cov,
            Map<Integer, Sclip> sclip5,
            Map<Integer, Sclip> sclip3,
            Map<Integer, Character> ref,
            String chr,
            Map<String, Integer> chrs,
            int rlen,
            String[] bams,
            Configuration conf) throws IOException {

        List<Tuple2<Integer, Sclip>> tmp = new ArrayList<>();
        for (Entry<Integer, Sclip> ent5 : sclip5.entrySet()) {
            tmp.add(tuple(ent5.getKey(), ent5.getValue()));
        }
        //sort by descending cnt
        Collections.sort(tmp, COMP2);

        for (Tuple2<Integer, Sclip> t : tmp) {
            int p = t._1;
            Sclip sc5v = t._2;
            //already been used in
            if (sc5v.used) {
                continue;
            }
            String seq = findconseq(sc5v, conf.y);
            if (seq.isEmpty()) {
                continue;
            }
            if (conf.y) {
                System.err.println("  Working lgins: 5: " + p + " " + seq);
            }
            if (B_A8.matcher(seq).find() || B_T8.matcher(seq).find()) {
                continue;
            }
            if (seq.length() < 12) {
                continue;
            }
            if (islowcomplexseq(seq)) {
                continue;
            }
            Tuple3<Integer, String, Integer> tpl = findbi(seq, p, ref, -1, chr, chrs);
            final int bi = tpl._1;
            final String ins = tpl._2;
            if (bi == 0) {
                continue;
            }
            if (conf.y) {
                System.err.printf("  Found candidate lgins from 5: %s +%s %s %s\n", bi, ins, p, seq);
            }
            final Variation iref = getVariation(iHash, bi, "+" + ins);
            iref.pstd = true;
            iref.qstd = true;
            adjCnt(iref, sc5v, conf);
            boolean rpflag = true; // A flag to indicate whether an insertion is a repeat
            for (int i = 0; i < ins.length(); i++) {
                if (!isEquals(ref.get(bi + 1 + i), ins.charAt(i))) {
                    rpflag = false;
                    break;
                }
            }
            incCnt(cov, bi, sc5v.cnt);
            int len = ins.length();
            if (ins.indexOf('&') != -1) {
                len--;
            }
            int seqLen = sc5v.seq.lastKey() + 1;
            for (int ii = len + 1; ii < seqLen; ii++) {
                int pii = bi - ii + len;
                if (!sc5v.seq.containsKey(ii)) {
                    continue;
                }
                for (Entry<Character, Variation> ent : sc5v.seq.get(ii).entrySet()) {
                    Character tnt = ent.getKey();
                    Variation tv = ent.getValue();
                    Variation tvr = getVariation(hash, pii, tnt.toString());
                    adjCnt(tvr, tv, conf);
                    tvr.pstd = true;
                    tvr.qstd = true;
                    incCnt(cov, pii, tv.cnt);
                }
            }
            sc5v.used = true;
            Map<Integer, Map<String, Integer>> tins = singletonMap(bi, singletonMap("+" + ins, iref.cnt));
            realignins(hash, iHash, tins, cov, sclip5, sclip3, ref, chr, chrs, conf);
            Variation mref = getVariationMaybe(hash, bi, ref.get(bi));
            if (rpflag && bams.length > 0 && ins.length() >= 5
                    && ins.length() < rlen - 10
                    && mref != null && mref.cnt != 0
                    && noPassingReads(chr, bi, bi + ins.length(), bams, conf)
                    && iref.cnt > 2 * mref.cnt) {
                adjCnt(iref, mref, mref, conf);
            }
        }

        tmp = new ArrayList<>();
        for (Entry<Integer, Sclip> ent3 : sclip3.entrySet()) {
            tmp.add(tuple(ent3.getKey(), ent3.getValue()));
        }
        Collections.sort(tmp, COMP2);
        for (Tuple2<Integer, Sclip> t : tmp) {
            final int p = t._1;
            final Sclip sc3v = t._2;
            if (sc3v.used) {
                continue;
            }
            String seq = findconseq(sc3v, conf.y);
            if (seq.isEmpty()) {
                continue;
            }
            if (conf.y) {
                System.err.println("  Working lgins 3: " + seq);
            }
            if (B_A7.matcher(seq).find() || B_T7.matcher(seq).find()) {
                continue;
            }
            if (seq.length() < 12) {
                continue;
            }
            if (islowcomplexseq(seq)) {
                continue;
            }
            Tuple3<Integer, String, Integer> tpl = findbi(seq, p, ref, 1, chr, chrs);
            final int bi = tpl._1;
            final String ins = tpl._2;
//            final int be = tpl._3;
            if (bi == 0) {
                continue;
            }
            final Variation iref = getVariation(iHash, bi, "+" + ins);
            iref.pstd = true;
            iref.qstd = true;
            final Variation lref = getVariationMaybe(hash, bi, ref.get(bi));
            adjCnt(iref, sc3v, lref, conf);
            boolean rpflag = true;
            for (int i = 0; i < ins.length(); i++) {
                if (!isEquals(ref.get(bi + 1 + i), ins.charAt(i))) {
                    rpflag = false;
                    break;
                }
            }
//            final int offset = bi == be ? (p - bi - 1) : -(p + be - bi);
            int len = ins.length();
            if (ins.indexOf('&') != -1) {
                len--;
            }
            int lenSeq = sc3v.seq.lastKey() + 1;
            for (int ii = len; ii < lenSeq; ii++) {
                int pii = p + ii - len;
                Map<Character, Variation> map = sc3v.seq.get(ii);
                if (map == null) {
                    continue;
                }
                for (Entry<Character, Variation> ent : map.entrySet()) {
                    Character tnt = ent.getKey();
                    Variation tv = ent.getValue();
                    Variation vref = getVariation(hash, pii, tnt.toString());
                    adjCnt(vref, tv, conf);
                    vref.pstd = true;
                    vref.qstd = true;
                    incCnt(cov, pii, tv.cnt);
                }
            }
            sc3v.used = true;
            Map<Integer, Map<String, Integer>> tins = singletonMap(bi, singletonMap("+" + ins, iref.cnt));
            realignins(hash, iHash, tins, cov, sclip5, sclip3, ref, chr, chrs, conf);
            Variation mref = getVariationMaybe(hash, bi, ref.get(bi));
            if (rpflag && bams.length > 0 && ins.length() >= 5 && ins.length() < rlen - 10
                    && mref != null && mref.cnt != 0
                    && noPassingReads(chr, bi, bi + ins.length(), bams, conf)
                    && iref.cnt > 2 * mref.cnt) {

                adjCnt(iref, mref, mref, conf);
            }
        }
    }

    /**
     * Find the insertion
     * @param seq sequence
     * @param p position
     * @param ref map of reference bases
     * @param dir direction
     * @param chr chromosome name
     * @param chrs map of chromosome lengths
     * @return Tuple of (BI (insert starting position), INS (insert sequence), BI2 ( = BI))
     */
    static Tuple3<Integer, String, Integer> findbi(String seq, int p, Map<Integer, Character> ref, final int dir, String chr, Map<String, Integer> chrs) {
        final int maxmm = 3; // maximum mismatches allowed
        final int dirExt = dir == -1 ? 1 : 0;
        int score = 0;
        int bi = 0;
        String ins = "";
        int bi2 = 0;

        for (int n = 6; n < seq.length(); n++) {
            if (p + 6 >= chrs.get(chr)) {
                break;
            }
            int mm = 0;
            int i = 0;
            Set<Character> m = new HashSet<>();
            for (i = 0; i + n < seq.length(); i++) {
                if (p + dir * i - dirExt < 1) {
                    break;
                }
                if (p + dir * i - dirExt > chrs.get(chr)) {
                    break;
                }
                if (isNotEquals(seq.charAt(i + n), ref.get(p + dir * i - dirExt))) {
                    mm++;
                } else {
                    m.add(seq.charAt(i + n));
                }
                if (mm > maxmm) {
                    break;
                }
            }
            int mnt = m.size();
            if (mnt < 2) { // at least three different NT for overhang sequences, weeding out low complexity seq
                continue;
            }
            if ((mnt >= 3 && i + n >= seq.length() - 1 && i >= 8 && mm / (double)i < 0.15)
                    || (mnt >= 2 && mm == 0 && i + n == seq.length() && n >= 20 && i >= 8)) {

                StringBuilder insert = new StringBuilder(substr(seq, 0, n));
                StringBuilder extra = new StringBuilder();
                int ept = 0;
                while (n + ept + 1 < seq.length() && (!isEquals(seq.charAt(n + ept), ref.get(p + ept * dir - dirExt))
                        || !isEquals(seq.charAt(n + ept + 1), ref.get(p + (ept + 1) * dir - dirExt)))) {
                    extra.append(seq.charAt(n + ept));
                    ept++;
                }
                if (dir == -1) {
                    insert.append(extra);
                    insert.reverse();
                    if (extra.length() > 0) {
                        insert.insert(insert.length() - extra.length(), "&");
                    }
                    if (mm == 0 && i + n == seq.length()) {
                        bi = p - 1 - extra.length();
                        ins = insert.toString();
                        bi2 = p - 1;
                        if (extra.length() == 0) {
                            Tuple3<Integer, String, Integer> tpl = adjInsPos(bi, ins, ref);
                            bi = tpl._1;
                            ins = tpl._2;
                            bi2 = tpl._3;
                        }
                        return tuple(bi, ins, bi2);
                    } else if (i - mm > score) {
                        bi = p - 1 - extra.length();
                        ins = insert.toString();
                        bi2 = p - 1;
                        score = i - mm;
                    }
                } else {
                    int s = -1;
                    if (extra.length() > 0) {
                        insert.append("&").append(extra);
                    } else {
                        while (s >= -n && isEquals(charAt(insert, s), ref.get(p + s))) {
                            s--;
                        }
                        if (s < -1) {
                            String tins = substr(insert.toString(), s + 1, 1 - s);
                            insert.delete(insert.length() + s + 1, insert.length());
                            insert.insert(0, tins);
                        }

                    }
                    if (mm == 0 && i + n == seq.length()) {
                        bi = p + s;
                        ins = insert.toString();
                        bi2 = p + s + extra.length();
                        if (extra.length() == 0) {
                            Tuple3<Integer, String, Integer> tpl = adjInsPos(bi, ins, ref);
                            bi = tpl._1;
                            ins = tpl._2;
                            bi2 = tpl._3;
                        }
                        return tuple(bi, ins, bi2);
                    } else if (i - mm > score) {
                        bi = p + s;
                        ins = insert.toString();
                        bi2 = p + s + extra.length();
                        score = i - mm;
                    }
                }
            }
        }
        if (bi2 == bi && ins.length() > 0 && bi != 0) {
            Tuple3<Integer, String, Integer> tpl = adjInsPos(bi, ins, ref);
            bi = tpl._1;
            ins = tpl._2;
        }
        return tuple(bi, ins, bi2);
    }

    /**
     * Adjust the insertion position if necessary
     * @param bi starting position of insert
     * @param ins insert sequence
     * @param ref map of reference bases
     * @return Tuple of (int bi, String ins, int bi)
     */
    static Tuple3<Integer, String, Integer> adjInsPos(int bi, String ins, Map<Integer, Character> ref) {
        int n = 1;
        int len = ins.length();
        while (isEquals(ref.get(bi), ins.charAt(ins.length() - n))) {
            n++;
            if (n > len) {
                n = 1;
            }
            bi--;
        }
        if (n > 1) {
            ins = substr(ins, 1 - n) + substr(ins, 0, 1 - n);
        }
        return tuple(bi, ins, bi);
    }

    private static final Comparator<Tuple2<Integer, Sclip>> COMP2 = new Comparator<Tuple2<Integer, Sclip>>() {
        @Override
        public int compare(Tuple2<Integer, Sclip> o1, Tuple2<Integer, Sclip> o2) {
            int f = Integer.compare(o2._2.cnt, o1._2.cnt);
            if (f != 0)
                return f;
            return Integer.compare(o1._1, o2._1);
        }
    };

    private static final Comparator<Tuple3<Integer, Sclip, Integer>> COMP3 = new Comparator<Tuple3<Integer, Sclip, Integer>>() {
        @Override
        public int compare(Tuple3<Integer, Sclip, Integer> o1, Tuple3<Integer, Sclip, Integer> o2) {
            return Integer.compare(o2._3, o1._3);
        }
    };

    private final static Pattern B_A7 = Pattern.compile("^.AAAAAAA");
    private final static Pattern B_T7 = Pattern.compile("^.TTTTTTT");
    private final static Pattern B_A8 = Pattern.compile("^.AAAAAAAA");
    private final static Pattern B_T8 = Pattern.compile("^.TTTTTTTT");

    /**
     * Realign large deletions that are not present in alignment
     * @param hash variant structure produced by parseSam method
     * @param cov coverage
     * @param sclip5 5' softclips
     * @param sclip3 3' softclips
     * @param ref map of reference bases
     * @param chr chromosome name
     * @param chrs map of chromosome lengths
     * @param rlen read length
     * @param bams BAM file list
     * @param conf configuration
     * @throws IOException
     */
    static void realignlgdel(Map<Integer, Map<String, Variation>> hash,
            Map<Integer, Integer> cov,
            Map<Integer, Sclip> sclip5,
            Map<Integer, Sclip> sclip3,
            Map<Integer, Character> ref,
            String chr,
            Map<String, Integer> chrs,
            final int rlen,
            String[] bams,
            Configuration conf) throws IOException {

        final int longmm = 3;

        List<Tuple2<Integer, Sclip>> tmp = new ArrayList<>();
        for (Entry<Integer, Sclip> ent5 : sclip5.entrySet()) {
            tmp.add(tuple(ent5.getKey(), ent5.getValue()));
        }
        Collections.sort(tmp, COMP2);
        for (Tuple2<Integer, Sclip> t : tmp) {
            int p = t._1;
            Sclip sc5v = t._2;
            int cnt = t._2.cnt;
            if (sc5v.used) {
                continue;
            }
            String seq = findconseq(sc5v, conf.y);
            if (seq.isEmpty()) {
                continue;
            }
            if (B_A7.matcher(seq).find() || B_T7.matcher(seq).find()) {
                continue;
            }
            if (seq.length() < 7) {
                continue;
            }
            if (islowcomplexseq(seq)) {
                continue;
            }
            int bp = findbp(seq, p - 5, ref, conf.indelsize, -1, chr, chrs, conf.y);
            final int dellen = p - bp;
            if (bp == 0) {
                continue;
            }
            StringBuilder extra = new StringBuilder();
            int en = 0;
            while (!Character.valueOf(seq.charAt(en)).equals(ref.get(bp - en - 1)) && en < seq.length()) {
                extra.append(seq.charAt(en));
                en++;
            }
            final String gt;
            if (extra.length() > 0) {
                extra = extra.reverse();
                gt = "-" + dellen + "&" + extra;
                bp -= extra.length();
            } else {
                gt = String.valueOf(-dellen);
            }
            if (conf.y) {
                System.err.printf("  Found Realignlgdel: %s %s 5' %s %s %s\n", bp, gt, p, seq, cnt);
            }
            final Variation tv = getVariation(hash, bp, gt);
            tv.qstd = true; // more accurate implementation lat
            tv.pstd = true; // more accurate implementation lat
            for (int tp = bp; tp < bp + dellen; tp++) {
                incCnt(cov, tp, sc5v.cnt);
            }
            adjCnt(tv, sc5v, conf);
            sc5v.used = true;

            // Work on the softclipped read at 3'
            int n = 0;
            while (ref.containsKey(bp + n)
                    && ref.containsKey(bp + dellen + n)
                    && isEquals(ref.get(bp + n), ref.get(bp + dellen + n))) {
                n++;
            }
            int sc3p = bp + n;
            StringBuilder str = new StringBuilder();
            int mcnt = 0;
            while (mcnt <= longmm
                    && ref.containsKey(bp + n)
                    && ref.containsKey(bp + dellen + n)
                    && isNotEquals(ref.get(bp + n), ref.get(bp + dellen + n))) {
                str.append(ref.get(bp + dellen + n));
                n++;
                mcnt++;
            }
            if (str.length() == 1) {
                while (ref.containsKey(bp + n)
                        && ref.containsKey(bp + dellen + n)
                        && isEquals(ref.get(bp + n), ref.get(bp + dellen + n))) {
                    n++;
                }
                sc3p = bp + n;
            }
            if (sclip3.containsKey(sc3p) && !sclip3.get(sc3p).used) {
                Sclip sclip = sclip3.get(sc3p);
                if (sc3p > bp) {
                    adjCnt(tv, sclip, getVariationMaybe(hash, bp, ref.get(bp)), conf);
                } else {
                    adjCnt(tv, sclip, conf);
                }

                if (sc3p == bp) {
                    for (int tp = bp; tp < bp + dellen; tp++) {
                        incCnt(cov, tp, sclip3.get(sc3p).cnt);
                    }
                }
                for (int ip = bp + 1; ip < sc3p; ip++) {
                    Variation vv = getVariation(hash, ip, ref.get(dellen + ip).toString());
                    rmCnt(vv, sclip);
                    if (vv.cnt == 0) {
                        hash.get(ip).remove(ref.get(dellen + ip).toString());
                    }
                    if (hash.get(ip).size() == 0) {
                        hash.remove(ip);
                    }
                }
                sclip.used = true;
            }
            Map<Integer, Map<String, Integer>> dels5 = singletonMap(bp, singletonMap(gt, tv.cnt));
            realigndel(hash, dels5, cov, sclip5, sclip3, ref, chr, chrs, rlen, bams, conf);
            if (conf.y) {
                System.err.printf("  Found lgdel done: %s %s %s 5' %s %s\n\n", bp, gt, p, seq, tv.cnt);
            }
        }

        tmp = new ArrayList<>();
        for (Entry<Integer, Sclip> ent3 : sclip3.entrySet()) {
            tmp.add(tuple(ent3.getKey(), ent3.getValue()));
        }
        Collections.sort(tmp, COMP2);
        for (Tuple2<Integer, Sclip> t : tmp) {
            int p = t._1;
            final Sclip sc3v = t._2;
            final int cnt = sc3v.cnt;
            if (sc3v.used) {
                continue;
            }
            String seq = findconseq(sc3v, conf.y);
            if (seq.isEmpty()) {
                continue;
            }
            if (B_A7.matcher(seq).find() || B_T7.matcher(seq).find()) {
                continue;
            }
            if (seq.length() < 7) {
                continue;
            }
            if (islowcomplexseq(seq)) {
                continue;
            }
            int bp = findbp(seq, p + 5, ref, conf.indelsize, 1, chr, chrs, conf.y);
            final int dellen = bp - p;
            if (bp == 0) {
                continue;
            }
            StringBuilder extra = new StringBuilder();
            int en = 0;
            while (en < seq.length() && isNotEquals(seq.charAt(en), ref.get(bp + en))) {
                extra.append(seq.charAt(en));
                en++;
            }
            bp = p; // Set it to 5'
            final String gt;
            if (extra.length() > 0) {
                gt = "-" + dellen + "&" + extra;
            } else {
                gt = "-" + dellen;
                while (isEquals(ref.get(bp - 1), ref.get(bp + dellen - 1))) {
                    bp--;
                }
            }
            if (conf.y) {
                System.err.printf("  Found Realignlgdel: %s %s 3' %s %s %s\n", bp, gt, p, seq, cnt);
            }
            Variation tv = getVariation(hash, bp, gt);
            tv.qstd = true; // more accurate implementation later
            tv.pstd = true; // more accurate implementation later
            for (int tp = bp; tp < bp + dellen; tp++) {
                incCnt(cov, tp, sc3v.cnt);
            }
            sc3v.pmean += dellen * sc3v.cnt;
            adjCnt(tv, sc3v, conf);
            sc3v.used = true;

            Map<Integer, Map<String, Integer>> dels5 = new HashMap<>();
            HashMap<String, Integer> map = new HashMap<>();
            map.put(gt, tv.cnt);
            dels5.put(bp, map);
            realigndel(hash, dels5, cov, sclip5, sclip3, ref, chr, chrs, rlen, bams, conf);
            if (conf.y) {
                System.err.printf("  Found lgdel: %s %s $p 3' %s %s\n\n", bp, gt, p, tv.cnt);
            }
        }
        if (conf.y) {
            System.err.println("  Done: Realignlgdel\n");
        }
    }

    private static boolean isHasAndEquals(Character ch1, Map<Integer, Character> ref, int index) {
        Character refc = ref.get(index);
        if (refc == null)
            return false;
        return refc.equals(ch1);
    }

    private static boolean isHasAndEquals(Map<Integer, Character> ref, int index1, String str, int index2) {
        Character refc = ref.get(index1);
        if (refc == null)
            return false;
        return refc.equals(str.charAt(index2));
    }

    private static boolean isHasAndNotEquals(Character ch1, Map<Integer, Character> ref, int index) {
        Character refc = ref.get(index);
        if (refc == null)
            return false;
        return !refc.equals(ch1);
    }

    private static boolean isHasAndNotEquals(Map<Integer, Character> ref, int index1, String str, int index2) {
        Character refc = ref.get(index1);
        if (refc == null)
            return false;
        return !refc.equals(str.charAt(index2));
    }

    private static boolean isEquals(Character ch1, Character ch2) {
        if (ch1 == null && ch2 == null)
            return true;
        if (ch1 == null || ch2 == null)
            return false;
        return ch1.equals(ch2);
    }

    private static boolean isNotEquals(Character ch1, Character ch2) {
        return !isEquals(ch1, ch2);
    }

    private static void rmCnt(Variation vref, Variation tv) {
        vref.cnt -= tv.cnt;
        vref.hicnt -= tv.hicnt;
        vref.locnt -= tv.locnt;
        vref.pmean -= tv.pmean;
        vref.qmean -= tv.qmean;
        vref.Qmean -= tv.Qmean;
        vref.subDir(true, tv.getDir(true));
        vref.subDir(false, tv.getDir(false));
        correctCnt(vref);
    }

    /**
     * Find breakpoint
     * @param seq sequence
     * @param sp start position of sequence
     * @param ref map of reference bases
     * @param dis distance (always = $INDELSIZE)
     * @param dir direction
     * @param chr chromosome name
     * @param chrs map of chromosome lengths
     * @param debugLog print debug message if true
     * @return breakpoint position
     */
    static int findbp(String seq, int sp,
            Map<Integer, Character> ref,
            int dis, int dir, String chr,
            Map<String, Integer> chrs,
            boolean debugLog) {

        final int maxmm = 3; // maximum mismatches allowed
        int bp = 0;
        int score = 0;
        int idx = chrs.containsKey(chr) ? chrs.get(chr) : 0;
        for (int n = 0; n < dis; n++) {
            int mm = 0;
            int i = 0;
            Set<Character> m = new HashSet<>();
            for (i = 0; i < seq.length(); i++) {
                if (sp + dir * n + dir * i < 1) {
                    break;
                }
                if (sp + dir * n + dir * i > idx) {
                    break;
                }
                if (isEquals(seq.charAt(i), ref.get(sp + dir * n + dir * i))) {
                    m.add(seq.charAt(i));
                } else {
                    mm++;
                }
                if (mm > maxmm - n / 100) {
                    break;
                }
            }
            if (m.size() < 3) {
                continue;
            }
            if (mm <= maxmm - n / 100 && i >= seq.length() - 2 && i >= 8 + n / 10 && mm / (double)i < 0.12) {
                int lbp = sp + dir * n - (dir < 0 ? dir : 0);
                if (mm == 0 && i == seq.length()) {
                    if (debugLog) {
                        System.err.printf("  Findbp: %s %s %s %s %s\n", seq, sp, lbp, mm, i);
                    }
                    return lbp;
                } else if (i - mm > score) {
                    bp = lbp;
                    score = i - mm;
                }
            }
        }
        if (debugLog && bp != 0) {
            System.err.printf("  Findbp with mismatches: %s %s %s %s %s\n", seq, sp, bp, dir, score);
        }
        return bp;
    }

    private static int count(String str, char chr) {
        int cnt = 0;
        for (int i = 0; i < str.length(); i++) {
            if (str.charAt(i) == chr) {
                cnt++;
            }
        }
        return cnt;

    }

    /**
     * If any of the base is represented by 75%, it's a low complexity seq
     *
     * Low-complexity sequences are simple repeats such as ATATATATAT
     * or regions that are highly enriched for just one letter, e.g. AAACAAAAAAAAGAAAAAAC.
     * @param seq sequence
     * @return true if sequence has low complexity
     */
    static boolean islowcomplexseq(String seq) {
        int len = seq.length();
        if (len == 0)
            return true;
        if (count(seq, 'A') / (double)len > 0.75)
            return true;
        if (count(seq, 'T') / (double)len > 0.75)
            return true;
        if (count(seq, 'G') / (double)len > 0.75)
            return true;
        return count(seq, 'C') / (double)len > 0.75;
    }

    private static List<Object[]> fillTmp(Map<Integer, Map<String, Integer>> changes) {
        List<Object[]> tmp = new ArrayList<>();
        for (Entry<Integer, Map<String, Integer>> ent : changes.entrySet()) {
            int p = ent.getKey();
            Map<String, Integer> v = ent.getValue();
            for (Entry<String, Integer> entV : v.entrySet()) {
                String vn = entV.getKey();
                int cnt = entV.getValue();
                // int ecnt = 0;
                // Matcher mtch = ATGC_E.matcher(vn);
                // if (mtch.find()) {
                // ecnt = mtch.group(1).length();
                // }
                tmp.add(new Object[] { p, vn, cnt, /* ecnt */});
            }
        }
        Collections.sort(tmp, REALIGNDEL_COMPARATOR);
        return tmp;
    }

    private static final Pattern BEGIN_PLUS_ATGC = Pattern.compile("^\\+([ATGC]+)");
    private static final Pattern AMP_ATGC = Pattern.compile("&([ATGC]+)");
    private static final Pattern HASH_ATGC = Pattern.compile("#([ATGC]+)");
    private static final Pattern CARET_ATGNC = Pattern.compile("\\^([ATGNC]+)");

    private static final Pattern BEGIN_MINUS_NUMBER = Pattern.compile("^-(\\d+)");
    private static final Pattern BEGIN_MINUS_NUMBER_ANY = Pattern.compile("^-\\d+(.*)");
    private static final Pattern UP_NUMBER_END = Pattern.compile("\\^(\\d+)$");
    private static final Pattern ATGSs_AMP_ATGSs_END = Pattern.compile("(\\+[ATGC]+)&[ATGC]+$");
    private static final Pattern PLUS_NUMBER = Pattern.compile("\\+(\\d+)");

    private static final Comparator<Object[]> REALIGNDEL_COMPARATOR = new Comparator<Object[]>() {
        @Override
        public int compare(Object[] o1, Object[] o2) {
            int x1 = (Integer)o1[2];
            int x2 = (Integer)o2[2];
            int f = Integer.compare(x2, x1);
            if (f != 0)
                return f;
            x1 = (Integer)o1[0];
            x2 = (Integer)o2[0];
            f =  Integer.compare(x1, x2);
            if (f != 0)
                return f;
            String s1 = (String)o1[1];
            String s2 = (String)o2[1];
            return s2.compareTo(s1);

        }
    };

    /**
     * Realign insertions
     * @param hash map of non-insertion variants produced by parseSAM method
     * @param iHash map of insertion variants produced by parseSAM method
     * @param ins insertion variants
     * @param cov coverage
     * @param sclip5 5' softclips
     * @param sclip3 3' softclips
     * @param ref map of reference bases
     * @param chr chromosome name
     * @param chrs map of chromosome lengths
     * @param conf configuration
     */
    static void realignins(Map<Integer, Map<String, Variation>> hash,
            Map<Integer, Map<String, Variation>> iHash,
            Map<Integer, Map<String, Integer>> ins,
            Map<Integer, Integer> cov,
            Map<Integer, Sclip> sclip5,
            Map<Integer, Sclip> sclip3,
            Map<Integer, Character> ref,
            String chr,
            Map<String, Integer> chrs,
            Configuration conf) {

        List<Object[]> tmp = fillTmp(ins);
        for (Object[] objects : tmp) {
            Integer p = (Integer)objects[0];
            String vn = (String)objects[1];
            Integer icnt = (Integer)objects[2];
            if (conf.y) {
                System.err.println(format("  Realign Ins: %s %s %s", p, vn, icnt));
            }
            String insert;
            Matcher mtch = BEGIN_PLUS_ATGC.matcher(vn);
            if (mtch.find()) {
                insert = mtch.group(1);
            } else {
                continue;
            }
            String extra = "";
            mtch = AMP_ATGC.matcher(vn);
            if (mtch.find()) {
                extra = mtch.group(1);
            }
            String compm = ""; // the match part for a complex variant
            mtch = HASH_ATGC.matcher(vn);
            if (mtch.find()) {
                compm = mtch.group(1);
            }
//            String newins = ""; // the adjacent insertion
//            mtch = CARET_ATGC_E.matcher(vn);
//            if (mtch.find()) {
//                newins = mtch.group(1);
//            }

            int newdel = 0; // the adjacent deletion
            try {
                newdel = toInt(vn);
            } catch (NumberFormatException ignore) {
            }
            String tn = vn.replaceFirst("^\\+", "")
                    .replaceFirst("&", "")
                    .replaceFirst("#", "")
                    .replaceFirst("\\^\\d+$", "")
                    .replaceFirst("\\^", "");

            int wustart = p - 100 - vn.length() + 1;
            if (wustart <= 1) {
                wustart = 1;
            }
            String wupseq = joinRef(ref, wustart, p) + tn; // 5prime flanking seq
            /*
            5' flanking region is a region of DNA that is adjacent to the 5' end of the gene.
            The 5' flanking region contains the promoter, and may contain enhancers or other protein binding sites.
            It is the region of DNA that is not transcribed into RNA.
             */
            Integer tend = chrs.get(chr);
            int sanend = p + vn.length() + 100;
            if (tend != null && tend < sanend) {
                sanend = tend;
            }
            String sanpseq = tn + joinRef(ref, p + extra.length() + 1 + compm.length() + newdel, sanend); // 3prime flanking seq
            /*
            3' flanking region is a region of DNA which is NOT copied into the mature mRNA, but which is present adjacent to 3' end of the gene.
            It was originally thought that the 3' flanking DNA was not transcribed at all, but it was discovered to be transcribed into RNA, but
            quickly removed during processing of the primary transcript to form the mature mRNA. The 3' flanking region often contains sequences
            which affect the formation of the 3' end of the message. It may also contain enhancers or other sites to which proteins may bind.
             */
            MMResult findMM3 = findMM3(ref, p + 1, sanpseq, insert.length() + compm.length(), sclip3); // mismatches, mismatch
                                                                                                       // positions, 5 or 3 ends
            MMResult findMM5 = findMM5(ref, p + extra.length() + compm.length() + newdel, wupseq, insert.length() + compm.length(), sclip5);

            List<Tuple3<String, Integer, Integer>> mm3 = findMM3.mm;
            List<Integer> sc3p = findMM3.scp;
            int nm3 = findMM3.nm;
            int misp3 = findMM3.misp;
            String misnt3 = findMM3.misnt;

            List<Tuple3<String, Integer, Integer>> mm5 = findMM5.mm;
            List<Integer> sc5p = findMM5.scp;
            int nm5 = findMM5.nm;
            int misp5 = findMM5.misp;
            String misnt5 = findMM5.misnt;

            List<Tuple3<String, Integer, Integer>> mmm = new ArrayList<>(mm3);
            mmm.addAll(mm5);
            Variation vref = getVariation(iHash, p, vn);
            for (Tuple3<String, Integer, Integer> tuple3 : mmm) {
                //mismatch nucleotide
                String mm = tuple3._1;
                //start position of clip that constains mm
                Integer mp = tuple3._2;
                //end (3 or 5)
                Integer me = tuple3._3;

                if (mm.length() > 1) {
                    mm = mm.charAt(0) + "&" + mm.substring(1);
                }
                if (!hash.containsKey(mp)) {
                    continue;
                }

                Variation tv = hash.get(mp).get(mm);
                if (tv == null) {
                    continue;
                }
                if (tv.cnt == 0) {
                    continue;
                }
                if (tv.qmean / tv.cnt < conf.goodq) {
                    continue;
                }
                if (tv.pmean / tv.cnt > (me == 3 ? nm3 + 4 : nm5 + 4)) { // opt_k;
                    continue;
                }
                if (tv.cnt >= icnt + insert.length() || tv.cnt / icnt >= 8) {
                    continue;
                }
                if (conf.y) {
                    System.err.printf("    insMM: %s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", mm, mp, me, nm3, nm5, vn, icnt, tv.cnt, tv.qmean, tv.pmean, cov.get(p));
                }
                // Adjust ref cnt so that AF won't > 1
                if (mp > p && me == 5) {
                    incCnt(cov, p, tv.cnt);
                }

                Variation lref = null;
                if (mp > p && me == 3 &&
                        hash.containsKey(p) &&
                        ref.containsKey(p) &&
                        hash.get(p).containsKey(ref.get(p).toString())) {

                    lref = hash.get(p).get(ref.get(p).toString());
                }
                adjCnt(vref, tv, lref, conf);
                hash.get(mp).remove(mm);
                if (hash.get(mp).isEmpty()) {
                    hash.remove(mp);
                }
            }
            if (misp3 != 0 && mm3.size() == 1 && hash.containsKey(misp3)) {
                hash.get(misp3).remove(misnt3);
            }
            if (misp5 != 0 && mm5.size() == 1 && hash.containsKey(misp5)) {
                hash.get(misp5).remove(misnt5);
            }
            for (Integer sc5pp : sc5p) {
                Sclip tv = sclip5.get(sc5pp);
                if (tv != null && !tv.used) {
                    String seq = findconseq(tv, conf.y);
                    if (conf.y) {
                        System.err.printf("    ins5: %s %s %s %s %s %s\n", p, sc5pp, seq, wupseq, icnt, tv.cnt);
                    }
                    if (!seq.isEmpty() && ismatch(seq, wupseq, -1, conf.y)) {
                        if (conf.y) {
                            System.err.printf("    ins5: %s %s $s %s %s %s used\n", p, sc5pp, seq, wupseq, icnt, tv.cnt);
                        }
                        if (sc5pp > p) {
                            incCnt(cov, p, tv.cnt);
                        }
                        adjCnt(vref, tv, conf);
                        tv.used = true;
                    }
                }
            }
            for (Integer sc3pp : sc3p) {
                Sclip tv = sclip3.get(sc3pp);
                if (conf.y) {
                    System.err.printf("    33: %s %s %s %s\n", p, sc3pp, vn, sanpseq);
                }
                if (tv != null && !tv.used) {
                    String seq = findconseq(tv, conf.y);
                    if (conf.y) {
                        System.err.printf("    ins3: %s %s %s %s %s %s %s\n", p, sc3pp, seq, sanpseq, vn, icnt, tv.cnt);
                    }
                    if (!seq.isEmpty() && ismatch(seq, substr(sanpseq, sc3pp - p - 1), 1, conf.y)) {
                        if (conf.y) {
                            System.err.printf("    ins3: %s %s %s %s %s %s used\n", p, sc3pp, seq, vn, icnt, tv.cnt);
                        }
                        if (sc3pp <= p) {
                            incCnt(cov, p, tv.cnt);
                        }
                        Variation lref = null;
                        if (sc3pp > p &&
                                hash.containsKey(p) &&
                                ref.containsKey(p) &&
                                hash.get(p).containsKey(ref.get(p).toString())) {

                            lref = hash.get(p).get(ref.get(p).toString());
                        }
                        adjCnt(vref, tv, lref, conf);
                        tv.used = true;
                    }

                }
            }
        }

        for (Object[] objects : tmp) {
            Integer p = (Integer)objects[0];
            String vn = (String)objects[1];
            if (!iHash.containsKey(p)) {
                continue;
            }
            Variation vref = iHash.get(p).get(vn);
            if (vref == null) {
                continue;
            }
            Matcher mtch = ATGSs_AMP_ATGSs_END.matcher(vn);
            if (mtch.find()) {
                String tn = mtch.group(1);
                Variation tref = iHash.get(p).get(tn);
                if (tref != null) {
                    if (vref.cnt < tref.cnt) {
                        adjCnt(tref, vref, getVariationMaybe(hash, p, ref.get(p)), conf);
                        iHash.get(p).remove(vn);
                    }
                }
            }
        }

    }

    /**
     * Realign deletions if already present in alignment
     * @param hash variant structure produced by parseSam method
     * @param dels5 deletion variants
     * @param cov coverage
     * @param sclip5 5' softclips
     * @param sclip3 3' softclips
     * @param ref map of reference bases
     * @param chr chromosome name
     * @param chrs map of chromosome lengths
     * @param rlen read length
     * @param bams BAM file list
     * @param conf configuration
     * @throws IOException
     */
    static void realigndel(Map<Integer, Map<String, Variation>> hash,
            Map<Integer, Map<String, Integer>> dels5,
            Map<Integer, Integer> cov,
            Map<Integer, Sclip> sclip5,
            Map<Integer, Sclip> sclip3,
            Map<Integer, Character> ref,
            String chr,
            Map<String, Integer> chrs,
            final int rlen,
            String[] bams,
            Configuration conf) throws IOException {

        // int longmm = 3; //Longest continued mismatches typical aligned at the end
        List<Object[]> tmp = fillTmp(dels5);

        for (Object[] objects : tmp) {
            Integer p = (Integer)objects[0];
            String vn = (String)objects[1];
            Integer dcnt = (Integer)objects[2];
            if (conf.y) {
                System.err.printf("  Realigndel for: %s %s %s cov: %s\n", p, vn, dcnt, cov.get(p));
            }
            final Variation vref = getVariation(hash, p, vn);
            int dellen = 0;
            Matcher mtch = BEGIN_MINUS_NUMBER.matcher(vn);
            if (mtch.find()) {
                dellen = toInt(mtch.group(1));
            }
            mtch = UP_NUMBER_END.matcher(vn);
            if (mtch.find()) {
                dellen += toInt(mtch.group(1));
            }
            String extrains = "";
            mtch = CARET_ATGNC.matcher(vn);
            if (mtch.find()) {
                extrains = mtch.group(1);
            }
            String extra = "";
            mtch = BEGIN_MINUS_NUMBER_ANY.matcher(vn);
            if (mtch.find()) {
                extra = mtch.group(1).replaceAll("\\^|&|#", "");
            }

            int wustart = p - dellen - 100;
            if (wustart <= 1) {
                wustart = 1;
            }

            String wupseq = joinRef(ref, wustart, p - 1) + extra; // 5' flanking seq
            int sanend = p + 2 * dellen + 100;
            if (sanend > chrs.get(chr)) {
                sanend = chrs.get(chr);
            }
            String sanpseq = extra + joinRef(ref, p + dellen + extra.length() - extrains.length(), sanend); // 3' flanking seq
            MMResult r3 = findMM3(ref, p, sanpseq, dellen, sclip3); // mismatches, mismatch positions, 5 or 3 ends
            MMResult r5 = findMM5(ref, p + dellen + extra.length() - extrains.length() - 1, wupseq, dellen, sclip5);

            List<Tuple3<String, Integer, Integer>> mm3 = r3.mm;
            List<Integer> sc3p = r3.scp;
            int nm3 = r3.nm;
            int misp3 = r3.misp;
            String misnt3 = r3.misnt;

            List<Tuple3<String, Integer, Integer>> mm5 = r5.mm;
            List<Integer> sc5p = r5.scp;
            int nm5 = r5.nm;
            int misp5 = r5.misp;
            String misnt5 = r5.misnt;
            if (conf.y) {
                System.err.printf("  Mismatches: misp3: %s-%s misp5: %s-%s sclip3: %s sclip5: %s\n",
                        misp3, misnt3, misp5, misnt5, Utils.toString(sc3p), Utils.toString(sc5p));
            }

            List<Tuple3<String, Integer, Integer>> mmm = new ArrayList<>(mm3);
            mmm.addAll(mm5);
            for (Tuple3<String, Integer, Integer> tuple : mmm) {
                String mm = tuple._1;
                Integer mp = tuple._2;
                Integer me = tuple._3;
                if (mm.length() > 1) {
                    mm = mm.charAt(0) + "&" + mm.substring(1);
                }
                if (hash.get(mp) == null) {
                    continue;
                }
                Variation tv = hash.get(mp).get(mm);
                if (tv == null) {
                    continue;
                }
                if (tv.cnt == 0) {
                    continue;
                }
                if (tv.qmean / tv.cnt < conf.goodq) {
                    continue;
                }

                if (tv.pmean / tv.cnt > (me == 3 ? nm3 + 4 : nm5 + 4)) {
                    continue;
                }
                if (tv.cnt >= dcnt + dellen || tv.cnt / dcnt >= 8) {
                    continue;
                }
                if (conf.y) {
                    System.err.printf("  Realigndel Adj: %s %s %s %s %s %s %s %s cov: %s\n",
                            mm, mp, me, nm3, nm5, p, tv.cnt, tv.qmean, cov.get(p));
                }
                // Adjust ref cnt so that AF won't > 1
                if (mp > p && me == 5) {
                    double f = tv.pmean != 0 ? (mp - p) / (tv.pmean / (double)tv.cnt) : 1;
                    if (f > 1) {
                        f = 1;
                    }
                    incCnt(cov, p, (int)(tv.cnt * f));
                    adjRefCnt(tv, getVariationMaybe(hash, p, ref.get(p)), dellen);
                }
                Variation lref = (mp > p && me == 3) ? (hash.containsKey(p) && hash.get(p).containsKey(ref.get(p).toString()) ? hash.get(p).get(ref.get(p).toString()) : null) : null;
                adjCnt(vref, tv, lref, conf);
                hash.get(mp).remove(mm);
                if (hash.get(mp).isEmpty()) {
                    hash.remove(mp);
                }
                if (conf.y) {
                    System.err.printf("  Realigndel AdjA: %s %s %s %s %s %s %s %s cov: %s\n",
                            mm, mp, me, nm3, nm5, p, tv.cnt, tv.qmean, cov.get(p));
                }
            }
            if (misp3 != 0 && mm3.size() == 1 && hash.containsKey(misp3)) {
                hash.get(misp3).remove(misnt3);
            }
            if (misp5 != 0 && mm5.size() == 1 && hash.containsKey(misp5)) {
                hash.get(misp5).remove(misnt5);
            }

            for (Integer sc5pp : sc5p) {
                if (sclip5.containsKey(sc5pp) && !sclip5.get(sc5pp).used) {
                    Sclip tv = sclip5.get(sc5pp);
                    String seq = findconseq(tv, conf.y);
                    if (conf.y) {
                        System.err.printf("  Realigndel 5: %s %s Seq: '%s' %s %s %s %s %s cov: %s\n",
                                p, sc5pp, seq, new StringBuilder(wupseq).reverse(), tv.cnt, dcnt, vn, p, cov.get(p));
                    }
                    if (!seq.isEmpty() && ismatch(seq, wupseq, -1, conf.y)) {
                        if (sc5pp > p) {
                            incCnt(cov, p, tv.cnt);
                        }
                        adjCnt(vref, tv, conf);
                        sclip5.get(sc5pp).used = true;
                        if (conf.y) {
                            System.err.printf("  Realigndel 5: %s %s %s %s %s %s %s %s used cov: %s\n",
                                    p, sc5pp, seq, new StringBuilder(wupseq).reverse(), tv.cnt, dcnt, vn, p, cov.get(p));
                        }
                    }
                }
            }

            for (Integer sc3pp : sc3p) {
                if (sclip3.containsKey(sc3pp) && !sclip3.get(sc3pp).used) {
                    Sclip tv = sclip3.get(sc3pp);
                    String seq = findconseq(tv, conf.y);
                    if (conf.y) {
                        System.err.printf("  Realigndel 3: %s %s seq '%s' %s %s %s %s %s %s %s\n",
                                p, sc3pp, seq, sanpseq, tv.cnt, dcnt, vn, p, dellen, substr(sanpseq, sc3pp - p));
                    }
                    if (!seq.isEmpty() && ismatch(seq, substr(sanpseq, sc3pp - p), 1, conf.y)) {
                        if (conf.y) {
                            System.err.printf("  Realigndel 3: %s %s %s %s %s %s %s %s used\n", p, sc3pp, seq, sanpseq, tv.cnt, dcnt, vn, p);
                        }
                        if (sc3pp <= p) {
                            incCnt(cov, p, tv.cnt);
                        }
                        Variation lref = sc3pp <= p ? null : getVariationMaybe(hash, p, ref.get(p));
                        adjCnt(vref, tv, lref, conf);
                        sclip3.get(sc3pp).used = true;
                    }
                }
            }
            // int pe = p + dellen + extra.length() + compm.length();
            int pe = p + dellen + extra.length() - extrains.length();
            Variation h = getVariationMaybe(hash, p, ref.get(p));
            if (bams != null && bams.length > 0
                    && pe - p >= 5
                    && pe - p < rlen - 10
                    && h != null && h.cnt != 0
                    && noPassingReads(chr, p, pe, bams, conf)
                    && vref.cnt > 2 * h.cnt * (1 - (pe - p) / (double)rlen)) {

                adjCnt(vref, h, h, conf);
            }

        }

        for (int i = tmp.size() - 1; i >= 0; i--) {
            Object[] os = tmp.get(i);
            int p = (Integer)os[0];
            String vn = (String)os[1];
            if (!hash.containsKey(p)) {
                continue;
            }
            Variation vref = hash.get(p).get(vn);
            if (vref == null) {
                continue;
            }
            Matcher matcher = MINUS_NUMBER_AMP_ATGCs_END.matcher(vn);
            if (matcher.find()) {
                String tn = matcher.group(1);
                Variation tref = hash.get(p).get(tn);
                if (tref != null) {
                    if (vref.cnt < tref.cnt) {
                        adjCnt(tref, vref, conf);
                        hash.get(p).remove(vn);
                    }
                }
            }
        }
    }

    private static final Pattern MINUS_NUMBER_AMP_ATGCs_END = Pattern.compile("(-\\d+)&[ATGC]+$");

    /**
     * check whether there're reads supporting wild type in deletions
     * Only for indels that have micro-homology
     * @param chr chromosome name
     * @param s start position
     * @param e end position
     * @param bams BAM file list
     * @param conf configuration
     * @return true if any read was found in chr:s-e
     * @throws IOException
     */
    static boolean noPassingReads(String chr, int s, int e, String[] bams, Configuration conf) throws IOException {
        int cnt = 0;
        int midcnt = 0; // Reads end in the middle
        int dlen = e - s;
        String dlenqr = dlen + "D";
        Region region = new Region(chr, s, e, "");
        for (String bam : bams) {
            try (SamView reader = new SamView(bam, "", region, conf.validationStringency)) {
                SAMRecord record;
                while ((record = reader.read()) != null) {
                    if (record.getCigarString().contains(dlenqr)) {
                        continue;
                    }
                    int rs = record.getAlignmentStart();
                    int rlen = getAlignedLenght(record.getCigar()); // The total aligned length, excluding soft-clipped bases and
                                                                      // insertions
                    int re = rs + rlen;
                    if (re > e + 2 && rs < s - 2) {
                        cnt++;
                    }
                    if (rs < s - 2 && re > s && re < e) {
                        midcnt++;
                    }
                }

            }
        }
        if (conf.y) {
            System.err.printf("    Passing Read CNT: %s %s %s %s %s\n", cnt, chr, s, e, midcnt);
        }
        return cnt <= 0;
    }

    /**
     * Find if sequences match
     * @param seq1 first sequence
     * @param seq2 second sequence
     * @param dir direction of seq2 (1 or -1)
     * @param debugLog print debug message if true
     * @return true if seq1 matches seq2 with no more than 2 mismathes
     */
    static boolean ismatch(String seq1, String seq2, int dir, boolean debugLog) {
        if (debugLog) {
            System.err.printf("    Matching %s %s %s\n", seq1, seq2, dir);
        }
        seq2 = seq2.replaceAll("#|\\^", "");
        int mm = 0;
        Map<Character, Boolean> nts = new HashMap<>();
        for (int n = 0; n < seq1.length() && n < seq2.length(); n++) {
            nts.put(seq1.charAt(n), true);
            if (seq1.charAt(n) != substr(seq2, dir * n - (dir == -1 ? 1 : 0), 1).charAt(0)) {
                mm++;
            }
        }

        return (mm <= 2 && mm / (double)seq1.length() < 0.15);
    }

    /**
     * Find the consensus sequence in soft-clipped reads. Consensus is called if
     * the matched nucleotides are &gt;90% of all softly clipped nucleotides.
     * @param scv soft-clipped sequences
     * @param debugLog print debug message if true
     * @return consensus sequence
     */
    static String findconseq(Sclip scv, boolean debugLog) {
        if (scv.sequence != null) {
            return scv.sequence;
        }

        int total = 0;
        int match = 0;
        StringBuilder seq = new StringBuilder();
        boolean flag = false;
        for (Map.Entry<Integer, Map<Character, Integer>> nve : scv.nt.entrySet()) {
            Integer i = nve.getKey();
            Map<Character, Integer> nv = nve.getValue();
            int max = 0;
            double maxq = 0;
            Character mnt = null;
            int tt = 0;
            for (Entry<Character, Integer> ent : nv.entrySet()) {
                Character nt = ent.getKey();
                int ncnt = ent.getValue();
                tt += ncnt;
                if (scv.seq.containsKey(i) && scv.seq.get(i).containsKey(nt) && scv.seq.get(i).get(nt).qmean > maxq) {
                    max = ncnt;
                    mnt = nt;
                    maxq = scv.seq.get(i).get(nt).qmean;
                }
            }
            if ((tt - max > 2 || max <= tt - max) && max / (double)tt < 0.8) {
                if (flag)
                    break;
                flag = true;
            }
            total += tt;
            match += max;
            if (mnt != null) {
                seq.append(mnt);
            }
        }

        Integer ntSize = scv.nt.lastKey();
        if (total != 0
                && match / (double)total > 0.9
                && seq.length() / 1.5 > ntSize - seq.length()
                && (seq.length() / (double)ntSize > 0.8
                        || ntSize - seq.length() < 10
                        || seq.length() > 25)) {
            scv.sequence = seq.toString();
        } else {
            scv.sequence = "";
        }
        if (debugLog) {
            System.err.printf("  candidate consensus: %s M: %s T: %s Final: %s\n", seq, match, total, scv.sequence);
        }
        return scv.sequence;

    }

    /**
     * Given a variant sequence, find the mismatches and potential softclipping positions
     * @param ref map of reference bases
     * @param p position
     * @param wupseq sequence
     * @param len not used
     * @param sclip5 map of 5' softclips
     * @return
     */
    static MMResult findMM5(Map<Integer, Character> ref, int p, String wupseq, int len, Map<Integer, Sclip> sclip5) {
        String seq = wupseq.replaceAll("#|\\^", "");
        int longmm = 3;
        List<Tuple3<String, Integer, Integer>> mm = new ArrayList<>(); // mismatches, mismatch positions, 5 or 3 ends
        int n = 0;
        int mn = 0;
        int mcnt = 0;
        StringBuilder str = new StringBuilder();
        List<Integer> sc5p = new ArrayList<>();
        while (isHasAndNotEquals(charAt(seq, -1 - n), ref, p - n) && mcnt < longmm) {
            str.insert(0, charAt(seq, -1 - n));
            mm.add(tuple(str.toString(), p - n, 5));
            n++;
            mcnt++;
        }
        sc5p.add(p + 1);
        // Adjust clipping position if only one mismatch
        int misp = 0;
        Character misnt = null;
        if (str.length() == 1) {
            while (isHasAndEquals(charAt(seq, -1 - n), ref, p - n)) {
                n++;
                mn++;
            }
            if (mn > 1) {
                int n2 = 0;
                while (-1 - n - 1 - n2 >= 0 && isHasAndEquals(charAt(seq, -1 - n - 1 - n2), ref, p - n - 1 - n2)) {
                    n2++;
                }
                if (n2 > 2) {
                    sc5p.add(p - n - n2);
                    misp = p - n;
                    misnt = charAt(seq, -1 - n);
                    if (sclip5.containsKey(p - n - n2)) {
                        sclip5.get(p - n - n2).used = true;
                    }
                    mn += n2;
                } else {
                    sc5p.add(p - n);
                    if (sclip5.containsKey(p - n)) {
                        sclip5.get(p - n).used = true;
                    }
                }

            }
        }
        return new MMResult(mm, sc5p, mn, misp, misnt == null ? "" : misnt.toString());
    }

    static class MMResult {
        private final List<Tuple3<String, Integer, Integer>> mm;
        private final List<Integer> scp;
        private final int nm;
        private final int misp;
        private final String misnt;

        public MMResult(List<Tuple3<String, Integer, Integer>> mm, List<Integer> scp, int nm, int misp, String misnt) {
            this.mm = mm;
            this.scp = scp;
            this.nm = nm;
            this.misp = misp;
            this.misnt = misnt;
        }

    }

    /**
     * Given a variant sequence, find the mismatches and potential softclipping positions
     * @param ref map of reference bases
     * @param p position
     * @param sanpseq sequence
     * @param len not used
     * @param sclip3 map of 3' softclips
     * @return array of [str (mismatches), Tbp (soft clip start position), 3], sc3p - array of soft clip positions, mn - number of nucleotides to adjust
     */
    static MMResult findMM3(Map<Integer, Character> ref, int p, String sanpseq, int len, Map<Integer, Sclip> sclip3) {
        String seq = sanpseq.replaceAll("#|\\^", ""); // ~ s/#|\^//g;
        final int longmm = 3;
        List<Tuple3<String, Integer, Integer>> mm = new ArrayList<>(); // mismatches, mismatch positions, 5 or 3 ends
        int n = 0;
        int mn = 0;
        int mcnt = 0;
        List<Integer> sc3p = new ArrayList<>();
        StringBuilder str = new StringBuilder();
        while (n < seq.length() && isEquals(ref.get(p + n), seq.charAt(n))) {
            n++;
        }
        sc3p.add(p + n);
        int Tbp = p + n;
        while (mcnt <= longmm && n < seq.length() && isNotEquals(ref.get(p + n), seq.charAt(n))) {
            str.append(seq.charAt(n));
            mm.add(tuple(str.toString(), Tbp, 3));
            n++;
            mcnt++;

        }
        // Adjust clipping position if only one mismatch
        int misp = 0;
        Character misnt = null;
        if (str.length() == 1) {
            while (n < seq.length() && isHasAndEquals(seq.charAt(n), ref, p + n)) {
                n++;
                mn++;
            }
            if (mn > 1) {
                int n2 = 0;
                while (n + n2 + 1 < seq.length() && isHasAndEquals(seq.charAt(n + n2 + 1), ref, p + n + 1 + n2)) {
                    n2++;
                }
                if (n2 > 2 && n + n2 + 1 < seq.length()) {
                    sc3p.add(p + n + n2);
                    misp = p + n;
                    misnt = Character.valueOf(seq.charAt(n));
                    if (sclip3.containsKey(p + n + n2)) {
                        sclip3.get(p + n + n2).used = true;
                    }
                    mn += n2;
                } else {
                    sc3p.add(p + n);
                    if (sclip3.containsKey(p + n)) {
                        sclip3.get(p + n).used = true;
                    }
                }
            }
        }
        // print STDERR "MM3: $seq $len $p '@sc3p'\n";

        return new MMResult(mm, sc3p, mn, misp, misnt == null ? "" : misnt.toString());
    }

    /**
     * Construct reference sequence
     * @param ref map of reference bases
     * @param from start position
     * @param to end position
     * @return reference sequence starting at from and ending at to
     */
    private static String joinRef(Map<Integer, Character> ref, int from, int to) {
        StringBuilder sb = new StringBuilder();
        for (int i = from; i <= to; i++) {
            Character ch = ref.get(i);
            if (ch != null) {
                sb.append(ch);
            }
        }
        return sb.toString();
    }

    // Find closest mismatches to combine with indels
    private static Tuple4<Integer, String, String, Integer> finndOffset(int refp, int readp, int mlen, String rdseq, String qstr,
            Map<Integer, Character> ref,
            Map<Integer, Integer> cov,
            int vext, int goodq) {
        int offset = 0;
        String ss = "";
        String q = "";
        int tnm = 0;
        int vsn = 0;
        for (int vi = 0; vsn <= vext && vi < mlen; vi++) {
            if (rdseq.charAt(readp + vi) == 'N') {
                break;
            }
            if (qstr.charAt(readp + vi) - 33 < goodq) {
                break;
            }
            Character refCh = ref.get(refp + vi);
            if (refCh != null) {
                char ch = rdseq.charAt(readp + vi);
                if (isNotEquals(ch, refCh)) {
                    offset = vi + 1;
                    tnm++;
                    vsn = 0;
                } else {
                    vsn++;
                }
            }
        }
        if (offset > 0) {
            ss = substr(rdseq, readp, offset);
            q = substr(qstr, readp, offset);
            for (int osi = 0; osi < offset; osi++) {
                incCnt(cov, refp + osi, 1);
            }
        }

        return tuple(offset, ss, q, tnm);

    }

    /**
     * Adjust MNP when there're breakpoints within MNP (multi-nucleotide polymorphism)
     * @param hash map of variants produced by parseSAM method
     * @param mnp map of MNP variants
     * @param cov coverage
     * @param ref
     * @param sclip5
     * @param sclip3
     * @param conf configuration
     */
    static void adjMNP(Map<Integer, Map<String, Variation>> hash,
            Map<Integer, Map<String, Integer>> mnp,
            Map<Integer, Integer> cov, Map<Integer, Character> ref, Map<Integer, Sclip> sclip3, Map<Integer, Sclip> sclip5, Configuration conf) {

        for (Map.Entry<Integer, Map<String, Integer>> entry : mnp.entrySet()) {
            final Integer p = entry.getKey();
            Map<String, Integer> v = entry.getValue();

            for (Map.Entry<String, Integer> en : v.entrySet()) {
                final String vn = en.getKey();
                final Map<String, Variation> hashP = hash.get(p);
                if (hashP == null) {
                    continue;
                }
                final Variation vref = hashP.get(vn);
                if (vref == null ) { // The variant is likely already been used by indel realignment
                    continue;
                }
                final String mnt = vn.replaceFirst("&", "");
                for (int i = 0; i < mnt.length() - 1; i++) {
                    String left = substr(mnt, 0, i + 1);
                    String right = substr(mnt, -(mnt.length() - i - 1));
                    {
                        Variation tref = hashP.get(left);
                        if (tref != null) {
                            if (tref.cnt < vref.cnt && tref.pmean / tref.cnt <= i + 1) {
                                if (conf.y) {
                                    System.err.printf(" AdjMnt Left: %s %s %s\n", p, vn, tref.cnt);
                                }
                                adjCnt(vref, tref, conf);
                                hashP.remove(left);
                            }
                        }
                    }
                    if (hash.containsKey(p + i + 1)) {
                        Variation tref = hash.get(p + i + 1).get(right);
                        if (tref != null) {
                            if (tref.cnt < vref.cnt && tref.pmean / tref.cnt <= mnt.length() - i - 1) {
                                if (conf.y) {
                                    System.err.printf(" AdjMnt Right: %s %s %s\n", p, vn, tref.cnt);
                                }
                                adjCnt(vref, tref, conf);
                                incCnt(cov, p, tref.cnt);
                                hash.get(p + i + 1).remove(right);
                            }
                        }
                    }

                }
                if (sclip3.containsKey(p)) {
                    final Sclip sc3v = sclip3.get(p);
                    if (!sc3v.used) {
                        final String seq = findconseq(sc3v, conf.y);
                        if (seq.startsWith(mnt)) {
                            if(seq.length() == mnt.length() || ismatchref(seq.substring(mnt.length()), ref, p + mnt.length(), 1, conf.y)) {
                                adjCnt(hash.get(p).get(vn), sc3v, conf);
                                incCnt(cov, p, sc3v.cnt);
                                sc3v.used = true;
                            }
                        }
                    }
                }
                if (sclip5.containsKey(p + mnt.length())) {
                    final Sclip sc5v = sclip5.get(p + mnt.length());
                    if (!sc5v.used) {
                        String seq =  findconseq(sc5v, conf.y);
                        if (!seq.isEmpty() && seq.length() >= mnt.length()) {
                            seq =  new StringBuffer(seq).reverse().toString();
                            if (seq.endsWith(mnt)) {
                                if (seq.length() == mnt.length() || ismatchref(seq.substring(0, seq.length() - mnt.length()), ref, p - 1, -1, conf.y)) {
                                    adjCnt(hash.get(p).get(vn), sc5v, conf);
                                    incCnt(cov, p, sc5v.cnt);
                                    sc5v.used = true;
                                }
                            }
                        }

                    }
                }
            }
        }

    }

    private static boolean ismatchref(String seq, Map<Integer, Character> ref, int p, int dir, boolean debugLog) {
        if (debugLog) {
            System.err.println(format("      Matching REF %s %s %s", seq, p, dir));
        }
        int mm = 0;
        for (int n = 0; n < seq.length(); n++) {
            final Character refCh = ref.get(p + dir * n);
            if (refCh == null || charAt(seq, dir == 1 ? n : dir * n - 1) != refCh) {
                mm++;
            }
        }
        return mm <= 3 && mm / (double)seq.length() < 0.15;
    }

    /**
     * Adjust count. Don't adjust reference
     * @param vref variant to add counts
     * @param tv variant
     * @param conf configuration
     */
    private static void adjCnt(Variation vref, Variation tv, Configuration conf) {
        adjCnt(vref, tv, null, conf);
    }

    /**
     * Adjust the count,  If ref is not null, the count for reference is also adjusted.
     * Variant counts of tv are added to vref and removed from ref
     * @param vref variant to add counts
     * @param tv variant
     * @param ref reference variant
     * @param conf configuration
     */
    private static void adjCnt(Variation vref, Variation tv, Variation ref, Configuration conf) {
        vref.cnt += tv.cnt;
        vref.extracnt += tv.cnt;
        vref.hicnt += tv.hicnt;
        vref.locnt += tv.locnt;
        vref.pmean += tv.pmean;
        vref.qmean += tv.qmean;
        vref.Qmean += tv.Qmean;
        vref.nm += tv.nm;
        vref.pstd = true;
        vref.qstd = true;
        vref.addDir(true, tv.getDir(true));
        vref.addDir(false, tv.getDir(false));

        if (conf.y) {
            String refCnt = ref != null ? String.valueOf(ref.cnt) : "NA";
            System.err.printf("    AdjCnt: '+' %s %s %s %s Ref: %s\n", vref.cnt, tv.cnt, vref.getDir(false), tv.getDir(false), refCnt);
            System.err.printf("    AdjCnt: '-' %s %s %s %s Ref: %s\n", vref.cnt, tv.cnt, vref.getDir(true), tv.getDir(true), refCnt);
        }

        if (ref == null)
            return;

        ref.cnt -= tv.cnt;
        ref.hicnt -= tv.hicnt;
        ref.locnt -= tv.locnt;
        ref.pmean -= tv.pmean;
        ref.qmean -= tv.qmean;
        ref.Qmean -= tv.Qmean;
        ref.nm -= tv.nm;
        ref.subDir(true, tv.getDir(true));
        ref.subDir(false, tv.getDir(false));
        correctCnt(ref);
    }

    /**
     * correct counts for negative values
     * @param ref variant
     */
    private static void correctCnt(Variation ref) {
        if (ref.cnt < 0)
            ref.cnt = 0;
        if (ref.hicnt < 0)
            ref.hicnt = 0;
        if (ref.locnt < 0)
            ref.locnt = 0;
        if (ref.pmean < 0)
            ref.pmean = 0;
        if (ref.qmean < 0)
            ref.qmean = 0;
        if (ref.Qmean < 0)
            ref.Qmean = 0;
        if (ref.getDir(true) < 0)
            ref.addDir(true, -ref.getDir(true));
        if (ref.getDir(false) < 0)
            ref.addDir(false, -ref.getDir(false));
    }

    /**
     * Adjust the reference count.
     * @param tv non-reference variant
     * @param ref reference variant
     * @param len length for adhustment factor
     */
    static void adjRefCnt(Variation tv, Variation ref, int len) {
        if (ref == null) {
            return;
        }
        double f = tv.pmean != 0 ? (tv.pmean / (double)tv.cnt - len + 1) / (tv.pmean / (double)tv.cnt) : 0; // the adjustment factor
        if (f < 0) {
            return;
        }

        if (f > 1) {
            f = 1;
        }

        ref.cnt -= f * tv.cnt;
        ref.hicnt -= f * tv.hicnt;
        ref.locnt -= f * tv.locnt;
        ref.pmean -= f * tv.pmean;
        ref.qmean -= f * tv.qmean;
        ref.Qmean -= f * tv.Qmean;
        ref.nm -= f * tv.nm;
        ref.subDir(true, (int)(f * tv.getDir(true)));
        ref.subDir(false, (int)(f * tv.getDir(false)));
        correctCnt(ref);
    }

    private static void increment(Map<Integer, Map<String, Integer>> counters, int idx, String s) {
        Map<String, Integer> map = counters.get(idx);
        if (map == null) {
            map = new HashMap<>();
            counters.put(idx, map);
        }
        incCnt(map, s, 1);
    }

    public static final BedRowFormat DEFAULT_BED_ROW_FORMAT = new BedRowFormat(2, 6, 7, 9, 10, 12);
    private static final BedRowFormat CUSTOM_BED_ROW_FORMAT = new BedRowFormat(0, 1, 2, 3, 1, 2);

    static class BedRowFormat {
        public final int chrColumn;
        public final int startColumn;
        public final int endColumn;
        public final int thickStartColumn;
        public final int thickEndColumn;
        public final int geneColumn;

        public BedRowFormat(int chrColumn, int startColumn, int endColumn, int thickStartColumn, int thickEndColumn, int geneColumn) {
            this.chrColumn = chrColumn;
            this.startColumn = startColumn;
            this.endColumn = endColumn;
            this.thickStartColumn = thickStartColumn;
            this.thickEndColumn = thickEndColumn;
            this.geneColumn = geneColumn;
        }

    }

    private static class ToVarsWorker implements Callable<Tuple2<Integer, Map<Integer, Vars>>> {
        final Region region;
        final String bam;
        final Map<String, Integer> chrs;
        final String sample;
        final Set<String> splice;
        final String ampliconBasedCalling;
        final Configuration conf;
        Map<Integer, Character> ref;

        public ToVarsWorker(Region region, String bam, Map<String, Integer> chrs, String sample, Set<String> splice, String ampliconBasedCalling, Map<Integer, Character> ref, Configuration conf) {
            super();
            this.region = region;
            this.bam = bam;
            this.chrs = chrs;
            this.sample = sample;
            this.splice = splice;
            this.ampliconBasedCalling = ampliconBasedCalling;
            this.conf = conf;
            this.ref = ref;
        }

        @Override
        public Tuple2<Integer, Map<Integer, Vars>> call() throws Exception {
            if (ref == null)
                ref = getREF(region, chrs, conf.fasta, conf.numberNucleotideToExtend);
            return toVars(region, bam, ref, chrs, sample, splice, ampliconBasedCalling, 0, conf);
        }

    }

    private static class SomdictWorker implements Callable<OutputStream> {

        final Future<Tuple2<Integer, Map<Integer, Vars>>> first;
        final ToVarsWorker second;
        final String sample;

        public SomdictWorker(Region region, String bam, Map<String, Integer> chrs, Set<String> splice, String ampliconBasedCalling, Map<Integer, Character> ref, Configuration conf,
                Future<Tuple2<Integer, Map<Integer, Vars>>> first,
                String sample) {
            this.first = first;
            this.second = new ToVarsWorker(region, bam, chrs, sample, splice, ampliconBasedCalling, ref, conf);
            this.sample = sample;
        }

        @Override
        public OutputStream call() throws Exception {
            Tuple2<Integer, Map<Integer, Vars>> t2 = second.call();
            Tuple2<Integer, Map<Integer, Vars>> t1 = first.get();
            ByteArrayOutputStream baos = new ByteArrayOutputStream();
            PrintStream out = new PrintStream(baos);
            somdict(second.region, t1._2, t2._2, sample, second.chrs, second.splice, second.ampliconBasedCalling, Math.max(t1._1, t2._1), second.conf, out);
            out.close();
            return baos;
        }

    }

    private static class VardictWorker implements Callable<OutputStream> {

        final String sample;
        private Region region;
        private Configuration conf;
        private Map<String, Integer> chrs;
        private Set<String> splice;
        private String ampliconBasedCalling;

        public VardictWorker(Region region, Map<String, Integer> chrs, Set<String> splice, String ampliconBasedCalling, String sample, Configuration conf) {
            super();
            this.region = region;
            this.chrs = chrs;
            this.splice = splice;
            this.ampliconBasedCalling = ampliconBasedCalling;
            this.sample = sample;
            this.conf = conf;
        }

        @Override
        public OutputStream call() throws Exception {
            Map<Integer, Character> ref = getREF(region, chrs, conf.fasta, conf.numberNucleotideToExtend);
            Tuple2<Integer, Map<Integer, Vars>> tpl = toVars(region, conf.bam.getBam1(), ref, chrs, sample, splice, ampliconBasedCalling, 0, conf);
            ByteArrayOutputStream baos = new ByteArrayOutputStream();
            PrintStream out = new PrintStream(baos);
            vardict(region, tpl._2, sample, splice, conf, out);
            out.close();
            return baos;
        }

    }

    private static class AmpVardictWorker implements Callable<OutputStream> {
        final Map<Integer, List<Tuple2<Integer, Region>>> pos;
        final Region rg;
        final List<Future<Tuple2<Integer, Map<Integer, Vars>>>> workers;
        final ToVarsWorker worker;
        final String sample;

        public AmpVardictWorker(Map<Integer, List<Tuple2<Integer, Region>>> pos, Region rg, String sample, List<Future<Tuple2<Integer,
                Map<Integer, Vars>>>> workers,
                ToVarsWorker worker) {
            this.pos = pos;
            this.rg = rg;
            this.workers = workers;
            this.worker = worker;
            this.sample = sample;
        }

        @Override
        public OutputStream call() throws Exception {
            Tuple2<Integer, Map<Integer, Vars>> last = worker.call();
            List<Map<Integer, Vars>> vars = new ArrayList<>();
            for (Future<Tuple2<Integer, Map<Integer, Vars>>> future : workers) {
                vars.add(future.get()._2);
            }
            vars.add(last._2);
            ByteArrayOutputStream baos = new ByteArrayOutputStream();
            PrintStream out = new PrintStream(baos);
            ampVardict(rg, vars, pos, sample, worker.splice, worker.conf, out);
            out.close();
            return baos;
        }
    }

    final static Future<OutputStream> NULL_FUTURE = new FutureTask<>(new Callable<OutputStream>() {
        @Override
        public OutputStream call() throws Exception {
            return null;
        }
    });

    private static void ampVardictParallel(final List<List<Region>> segs, final Map<String, Integer> chrs, final String ampliconBasedCalling,
            final String bam1, final String sample, final Configuration conf) throws IOException {

        final ExecutorService executor = new ForkJoinPool(); //Executors.newFixedThreadPool(conf.threads);
        final BlockingQueue<Future<OutputStream>> toPrint = new LinkedBlockingQueue<>(21);

        executor.submit(new Runnable() {
            @Override
            public void run() {
                try {
                    for (List<Region> regions : segs) {
                        Map<Integer, List<Tuple2<Integer, Region>>> pos = new HashMap<>();
                        int j = 0;
                        Region rg = null;
                        List<Future<Tuple2<Integer, Map<Integer, Vars>>>> workers = new ArrayList<>(regions.size() - 1);
                        final Set<String> splice = new ConcurrentHashSet<>();
                        for (Region region : regions) {
                            rg = region; // ??
                            for (int p = region.istart; p <= region.iend; p++) {
                                List<Tuple2<Integer, Region>> list = pos.get(p);
                                if (list == null) {
                                    list = new ArrayList<>();
                                    pos.put(p, list);
                                }
                                list.add(tuple(j, region));
                            }
                            ToVarsWorker toVars = new ToVarsWorker(region, bam1, chrs, sample, splice, ampliconBasedCalling, null, conf);
                            if (workers.size() == regions.size() - 1) {
                                toPrint.put(executor.submit(new AmpVardictWorker(pos, rg, sample, workers, toVars)));
                            } else {
                                workers.add(executor.submit(toVars));
                            }
                            j++;

                        }
                    }
                    toPrint.put(NULL_FUTURE);
                } catch (InterruptedException e) {
                    e.printStackTrace();
                }
            }
        });

        try {
            while (true) {
                Future<OutputStream> seg = toPrint.take();
                if (seg == NULL_FUTURE) {
                    break;
                }
                System.out.print(seg.get());
            }
        } catch (InterruptedException | ExecutionException e) {
            e.printStackTrace();
        }
        executor.shutdown();
    }

    private static void ampVardictNotParallel(final List<List<Region>> segs, final Map<String, Integer> chrs, final String ampliconBasedCalling,
            final String bam1, final String sample, final Configuration conf) throws IOException {

        for (List<Region> regions : segs) {
            Map<Integer, List<Tuple2<Integer, Region>>> pos = new HashMap<>();
            int j = 0;
            Region rg = null;
            final Set<String> splice = new HashSet<>();
            List<Map<Integer, Vars>> vars = new ArrayList<>();
            for (Region region : regions) {
                rg = region; // ??
                for (int p = region.istart; p <= region.iend; p++) {
                    List<Tuple2<Integer, Region>> list = pos.get(p);
                    if (list == null) {
                        list = new ArrayList<>();
                        pos.put(p, list);
                    }
                    list.add(tuple(j, region));
                }
                vars.add(toVars(region, bam1, getREF(region, chrs, conf.fasta, conf.numberNucleotideToExtend), chrs, sample, splice, ampliconBasedCalling, 0, conf)._2);
                j++;
            }
            ampVardict(rg, vars, pos, sample, splice, conf, System.out);
        }
    }

    /**
     * amplicon variant calling
     *
     * @param rg
     *            region
     * @param vars
     *            result of {@link #toVars} calling
     * @param positions map of positon =&gt; (list of (region number, region))
     * @param sample sample name
     * @param splice set of strings representing spliced regions
     * @param conf configuration
     * @param out output stream
     */
    static void ampVardict(Region rg, List<Map<Integer, Vars>> vars, Map<Integer, List<Tuple2<Integer, Region>>> positions,
            final String sample, final Set<String> splice, final Configuration conf,
            PrintStream out) {

        List<Integer> pp = new ArrayList<>(positions.keySet());
        Collections.sort(pp);
        for (Integer p : pp) {

            final List<Tuple2<Integer, Region>> v = positions.get(p);

            // good variants
            List<Tuple2<Variant, String>> gvs = new ArrayList<>();
            //reference variants
            List<Variant> ref = new ArrayList<>();
            String nt = null;
            double maxaf = 0;
            //vartype may take values SNV (Single Nucleotide Variant), Complex (or MNV (Multiple Nucleotide Variant)), Insertion, Deletion
            String vartype = "SNV";
            boolean flag = false;
            Variant vref;
            //DNA sequencing coverage
            int nocov = 0;
            //max DNA sequencing coverage (max depth)
            int maxcov = 0;
            //good amplicon
            Set<String> goodmap = new HashSet<>();
            List<Integer> vcovs = new ArrayList<>();
            //amps map of amplicons.
            //An amplicon is a piece of DNA or RNA that is the source and/or product of
            //natural or artificial amplification or replication events.
            for (Tuple2<Integer, Region> amps : v) {
                final int amp = amps._1;
                //chromosome name
                final String chr = amps._2.chr;
                //start index
                final int S = amps._2.start;
                //end index
                final int E = amps._2.end;

                Vars vtmp = vars.get(amp).get(p);
                List<Variant> l = vtmp == null ? null : vtmp.var;
                Variant refAmpP = vtmp == null ? null : vtmp.ref;
                if (l != null && !l.isEmpty()) {
                    Variant tv = l.get(0);
                    vcovs.add(tv.tcov);
                    if (tv.tcov > maxcov) {
                        maxcov = tv.tcov;
                    }
                    vartype = varType(tv.refallele, tv.varallele);
                    if (isGoodVar(tv, refAmpP, vartype, splice, conf)) {
                        gvs.add(tuple(tv, chr + ":" + S + "-" + E));
                        if (nt != null && !tv.n.equals(nt)) {
                            flag = true;
                        }
                        if (tv.freq > maxaf) {
                            maxaf = tv.freq;
                            nt = tv.n;
                            vref = tv;
                        }
                        goodmap.add(format("%s-%s-%s", amp, tv.refallele, tv.varallele));
                    }
                } else if (refAmpP != null) {
                    vcovs.add(refAmpP.tcov);
                } else {
                    vcovs.add(0);
                }
                if (refAmpP != null) {
                    ref.add(refAmpP);
                }
            }

            //Depth (coverage) in DNA sequencing refers to the number of times a nucleotide is read during the sequencing process.
            // Coverage is the average number of reads representing a given nucleotide in the reconstructed sequence.
            for (int t : vcovs) {
                if (t < maxcov / 50) { //The amplicon that has depth less than 1/50 of the max depth will be considered not working and thus not used.
                    nocov++;
                }
            }

            if (gvs.size() > 1) {
                Collections.sort(gvs, GVS_COMPARATOR);
            }
            if (ref.size() > 1) {
                Collections.sort(ref, VAR_TCOV_COMPARATOR);
            }

            if (gvs.isEmpty()) { // Only referenece
                if (conf.doPileup) {
                    if (!ref.isEmpty()) {
                        vref = ref.get(0);
                    } else {
                        out.println(join("\t", sample, rg.gene, rg.chr, p, p, "", "", 0, 0, 0, 0, 0, 0, "", 0,
                                "0;0", 0, 0, 0, 0, 0, "", 0, 0, 0, 0, 0, 0, "", "", 0, 0,
                                rg.chr + ":", +p + "-" + p, "", 0, 0, 0, 0));
                        continue;
                    }
                } else {
                    continue;
                }
            } else {
                vref = gvs.get(0)._1;
            }
            if (flag) { // different good variants detected in different amplicons
                String gdnt = gvs.get(0)._1.n;
                List<Tuple2<Variant, String>> gcnt = new ArrayList<>();
                for (Tuple2<Integer, Region> amps : v) {
                    final Vars vtmp = vars.get(amps._1).get(p);
                    final Variant variant = vtmp == null ? null : vtmp.varn.get(gdnt);
                    if (variant != null && isGoodVar(variant, vtmp.ref, null, splice, conf)) {
                        gcnt.add(tuple(variant, amps._2.chr + ":" + amps._2.start + "-" + amps._2.end));
                    }
                }
                if (gcnt.size() == gvs.size()) {
                    flag = false;
                }
                Collections.sort(gcnt, GVS_COMPARATOR);
                gvs = gcnt;
            }

            //bad variants
            List<Tuple2<Variant, String>> badv = new ArrayList<>();
            int gvscnt = gvs.size();
            if (gvscnt != v.size() || flag) {
                for (Tuple2<Integer, Region> amps : v) {
                    int amp = amps._1;
                    Region reg = amps._2;
                    if (goodmap.contains(format("%s-%s-%s", amp, vref.refallele, vref.varallele))) {
                        continue;
                    }
                    // my $tref = $vars[$amp]->{ $p }->{ VAR }->[0]; ???
                    if (vref.sp >= reg.istart && vref.ep <= reg.iend) {

                        String regStr = reg.chr + ":" + reg.start + "-" + reg.end;

                        if (vars.get(amp).containsKey(p) && vars.get(amp).get(p).var.size() > 0) {
                            badv.add(tuple(vars.get(amp).get(p).var.get(0), regStr));
                        } else if (vars.get(amp).containsKey(p) && vars.get(amp).get(p).ref != null) {
                            badv.add(tuple(vars.get(amp).get(p).ref, regStr));
                        } else {
                            badv.add(tuple((Variant)null, regStr));
                        }
                    } else if ((vref.sp < reg.iend && reg.iend < vref.ep)
                            || (vref.sp < reg.istart && reg.istart < vref.ep)) { // the variant overlap with amplicon's primer
                        if (gvscnt > 1)
                            gvscnt--;
                    }
                }
            }
            if (flag && gvscnt < gvs.size()) {
                flag = false;
            }
            if (vartype.equals("Complex")) {
                adjComplex(vref);
            }
            out.print(join("\t", sample, rg.gene, rg.chr,
                    joinVar1(vref, "\t"), gvs.get(0)._2, vartype, gvscnt, gvscnt + badv.size(), nocov, flag ? 1 : 0));
            if (conf.debug) {
                out.print("\t" + vref.DEBUG);
            }
            if (conf.debug) {
                for (int gvi = 0; gvi < gvs.size(); gvi++) {
                    Tuple2<Variant, String> tp = gvs.get(gvi);
                    out.print("\tGood" + gvi + " " + join(" ", joinVar2(tp._1, " "), tp._2));
                }
                for (int bvi = 0; bvi < badv.size(); bvi++) {
                    Tuple2<Variant, String> tp = badv.get(bvi);
                    out.print("\tBad" + bvi + " " + join(" ", joinVar2(tp._1, " "), tp._2));
                }
            }
            out.println();
        }
    }

    final static Comparator<Variant> VAR_TCOV_COMPARATOR = new Comparator<Variant>() {
        @Override
        public int compare(Variant o1, Variant o2) {
            return Integer.compare(o2.tcov, o1.tcov);
        }
    };


    final static Comparator<Tuple2<Variant, String>> GVS_COMPARATOR = new Comparator<Tuple2<Variant, String>>() {
        @Override
        public int compare(Tuple2<Variant, String> a, Tuple2<Variant, String> b) {
            return Double.compare(b._1.freq, a._1.freq);
        }

    };


    private static String joinVar2(Variant var, String delm) {
        StringBuilder sb = new StringBuilder();

        sb.append(var.tcov).append(delm);
        sb.append(var.cov).append(delm);
        sb.append(var.rfc).append(delm);
        sb.append(var.rrc).append(delm);
        sb.append(var.fwd).append(delm);
        sb.append(var.rev).append(delm);
        sb.append(var.genotype).append(delm);
        sb.append(var.freq == 0 ? 0 : format("%.4f", var.freq)).append(delm);
        sb.append(var.bias).append(delm);
        sb.append(format("%.1f", var.pmean)).append(delm);
        sb.append(var.pstd ? 1 : 0).append(delm);
        sb.append(format("%.1f", var.qual)).append(delm);
        sb.append(var.qstd ? 1 : 0).append(delm);
        sb.append(format("%.1f", var.mapq)).append(delm);
        sb.append(format("%.3f", var.qratio)).append(delm);
        sb.append(var.hifreq == 0 ? 0 : format("%.4f", var.hifreq)).append(delm);
        sb.append(var.extrafreq == 0 ? 0 : format("%.4f", var.extrafreq));

        return sb.toString();
    }

    private static String joinVar1(Variant var, String delm) {
        StringBuilder sb = new StringBuilder();

        sb.append(var.sp).append(delm);
        sb.append(var.ep).append(delm);
        sb.append(var.refallele).append(delm);
        sb.append(var.varallele).append(delm);

        sb.append(joinVar2(var, delm)).append(delm);

        sb.append(var.shift3).append(delm);
        sb.append(var.msi == 0 ? 0 : format("%.3f", var.msi)).append(delm);
        sb.append(var.msint).append(delm);
        sb.append(format("%.1f", var.nm)).append(delm);
        sb.append(var.hicnt).append(delm);
        sb.append(var.hicov).append(delm);
        sb.append(var.leftseq).append(delm);
        sb.append(var.rightseq);

        return sb.toString();
    }

    /**
     * Adjust the complex variant
     * @param vref variant
     */
    static void adjComplex(Variant vref) {
        String refnt = vref.refallele;
        String varnt = vref.varallele;
        int n = 0;
        while (refnt.length() - n > 1 && varnt.length() - n > 1 && refnt.charAt(n) == varnt.charAt(n)) {
            n++;
        }
        if (n > 0) {
            vref.sp += n;
            vref.refallele = substr(refnt, n);
            vref.varallele = substr(varnt, n);
            vref.leftseq += substr(refnt, 0, n);
            vref.leftseq = substr(vref.leftseq, n);
        }
        refnt = vref.refallele;
        varnt = vref.varallele;
        n = 1;
        while (refnt.length() - n > 0 && varnt.length() - n > 0 && substr(refnt, -n, 1).equals(substr(varnt, -n, 1))) {
            n++;
        }
        if (n > 1) {
            vref.ep -= n - 1;
            vref.refallele = substr(refnt, 0, 1 - n);
            vref.varallele = substr(varnt, 0, 1 - n);
            vref.rightseq = substr(refnt, 1 - n, n - 1) + substr(vref.rightseq, 0, 1 - n);
        }

    }

    /**
     * VARN - array of variants, REF - reference variant, VAR - map of (variant description string =&gt; variant)
     */
    static enum VarsType {
        varn, ref, var,
    }

    private static Vars getOrPutVars(Map<Integer, Vars> map, int position) {
        Vars vars = map.get(position);
        if (vars == null) {
            vars = new Vars();
            map.put(position, vars);
        }
        return vars;
    }

    private static Variant getVarMaybe(Map<Integer, Vars> vars, int key, VarsType type, Object... keys) {
        if (vars.containsKey(key)) {
            switch (type) {
                case var:
                    if (vars.get(key).var.size() > (Integer)keys[0]) {
                        return vars.get(key).var.get((Integer)keys[0]);
                    }
                case varn:
                    return vars.get(key).varn.get(keys[0]);
                case ref:
                    return vars.get(key).ref;
            }
        }
        return null;
    }

    private static Variant getVarMaybe(Vars vars, VarsType type, Object... keys) {
        switch (type) {
            case var:
                if (vars.var.size() > (Integer)keys[0]) {
                    return vars.var.get((Integer)keys[0]);
                }
            case varn:
                return vars.varn.get(keys[0]);
            case ref:
                return vars.ref;
        }
        return null;
    }

    /**
     * Returns <tt>true</tt> whether a variant meet specified criteria
     * @param vref variant
     * @param rref reference variant
     * @param type Type of variant
     * @param splice set of strings representing introns in splice
     * @param conf Configuration (contains preferences for min, freg, filter and etc)
     * @return <tt>true</tt> if variant meet specified criteria
     */
    static boolean isGoodVar(Variant vref, Variant rref, String type,
            Set<String> splice,
            Configuration conf) {
        if (vref == null || vref.refallele == null || vref.refallele.isEmpty())
            return false;

        if (type == null || type.isEmpty()) {
            type = varType(vref.refallele, vref.varallele);
        }
        if (vref.freq < conf.freq
                || vref.hicnt < conf.minr
                || vref.pmean < conf.readPosFilter
                || vref.qual < conf.goodq) {
            return false;
        }

        if (rref != null && rref.hicnt > conf.minr && vref.freq < 0.25d) {
            //The reference allele has much better mean mapq than var allele, thus likely false positives
            double d = vref.mapq + vref.refallele.length() + vref.varallele.length();
            if ((d - 2 < 5 && rref.mapq > 20)
                    || (1 + d) / (rref.mapq + 1) < 0.25d) {
                return false;
            }
        }

        if (type.equals("Deletion") && splice.contains(vref.sp + "-" + vref.ep)) {
            return false;
        }
        if (vref.qratio < conf.qratio) {
            return false;
        }
        if (vref.freq > 0.35d) {
            return true;
        }
        if (vref.mapq < conf.mapq) {
            return false;
        }
        if (vref.msi >= 13 && vref.freq <= 0.275d && vref.msint == 1) {
            return false;
        }
        if (vref.msi >= 8 && vref.freq <= 0.2d && vref.msint > 1) {
            return false;
        }
        if (vref.bias.equals("2;1") && vref.freq < 0.25d) {
            if (type == null || type.equals("SNV") || (vref.refallele.length() <= 3 && vref.varallele.length() <= 3)) {
                return false;
            }
        }

        return true;
    }

    /**
     * Find variant type based on variant sequence
     * @param ref reference sequence
     * @param var variant sequence
     * @return variant type
     */
    private static String varType(String ref, String var) {
        if (ref.length() == 1 && var.length() == 1) {
            return "SNV";
        } else if (ref.length() == 0 || var.length() == 0) { //issue #19
            return "Complex";
        } else if (ref.charAt(0) != var.charAt(0)) {
            return "Complex";
        } else if (ref.length() == 1 && var.length() > 1 && var.startsWith(ref)) {
            return "Insertion";
        } else if (ref.length() > 1 && var.length() == 1 && ref.startsWith(var)) {
            return "Deletion";
        }
        return "Complex";
    }

    /**
     * Create region from command-line option
     * @param region Region string from command line
     * @param numberNucleotideToExtend number of nucleotides to extend region
     * @param chrs
     * @param zeroBased Are positions zero-based or 1-based
     * @return Region
     */
    private static Region buildRegion(String region, final int numberNucleotideToExtend, Map<String, Integer> chrs, final boolean zeroBased) {
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

    /**
     * Modify the CIGAR for potential mis-alignment for indels at the end of reads to softclipping and let VarDict's algorithm to figure out indels
     * @param indel lenght of insert/delete
     * @param ref Map of reference sequence (key - position, value - base)
     * @param oPosition Position of first matched base in sequence
     * @param oCigar Original CIGAR string
     * @param querySeq Base sequence
     * @param queryQual
     * @param lowqual
     * @return Tuple of (adjusted position of first matched base, modified CIGAR string)
     */
    static Tuple2<Integer, String> modifyCigar(int indel, Map<Integer, Character> ref,
            final int oPosition, final String oCigar, final String querySeq,
            final String queryQual, final int lowqual) {

        int position = oPosition;
        String cigarStr = oCigar;
        /**
         * flag is set to true if CIGAR string is modified and should be looked at again
         */
        boolean flag = true;
        while (flag && indel > 0) {
            flag = false;
            jregex.Matcher mm = BEGIN_NUMBER_S_NUMBER_IorD.matcher(cigarStr);
            if (mm.find()) { // If CIGAR starts with soft-clipping followed by insertion or deletion
                /*
                If insertion follows soft-clipping, add the inserted sequence to soft-clipped start
                Otherwise increase position by number of deleted bases
                 */
                String tslen = toInt(mm.group(1)) + (mm.group(3).equals("I") ? toInt(mm.group(2)) : 0) + "S";
                position += mm.group(3).equals("D") ? 2 : 0;
                //Regexp replaces found CIGAR sequence with tslen (number + S)
                Replacer r = BEGIN_NUMBER_S_NUMBER_IorD.replacer(tslen);
                cigarStr = r.replace(cigarStr);
                flag = true;
            }
            mm = NUMBER_IorD_NUMBER_S_END.matcher(cigarStr);
            if (mm.find()) { // If CIGAR ends with insertion or deletion followed by soft-clipping
                //Replace insertion or deletion with soft-clipping
                String tslen = toInt(mm.group(3)) + (mm.group(2).equals("I") ? toInt(mm.group(1)) : 0) + "S";
                //Regexp replaces found CIGAR sequence with $tslen (number + S)
                Replacer r = NUMBER_IorD_NUMBER_S_END.replacer(tslen);
                cigarStr = r.replace(cigarStr);
                flag = true;
            }
            mm = BEGIN_NUMBER_S_NUMBER_M_NUMBER_IorD.matcher(cigarStr);
            if (mm.find()) { // If CIGAR starts with soft-clipping followed by matched sequence and insertion or deletion
                int tmid = toInt(mm.group(2));
                if (tmid <= 10) { // If matched sequence length is no more than 10, replace everything with soft-clipping
                    String tslen = toInt(mm.group(1)) + tmid + (mm.group(4).equals("I") ? toInt(mm.group(3)) : 0) + "S";
                    position += tmid + (mm.group(4).equals("D") ? toInt(mm.group(3)) : 0);
                    Replacer r = BEGIN_NUMBER_S_NUMBER_M_NUMBER_IorD.replacer(tslen);
                    cigarStr = r.replace(cigarStr);
                    flag = true;
                }
            }
            mm = NUMBER_IorD_NUMBER_M_NUMBER_S_END.matcher(cigarStr);
            if (mm.find()) { // If CIGAR ends with insertion or deletion, matched sequence and soft-clipping
                int tmid = toInt(mm.group(3));
                if (tmid <= 10) { // If matched sequence length is no more than 10, replace everything with soft-clipping
                    String tslen = toInt(mm.group(4)) + tmid + (mm.group(2).equals("I") ? toInt(mm.group(1)) : 0) + "S";
                    Replacer r = NUMBER_IorD_NUMBER_M_NUMBER_S_END.replacer(tslen);
                    cigarStr = r.replace(cigarStr);
                    flag = true;
                }
            }

            // The following two clauses to make indels at the end of reads as softly
            // clipped reads and let VarDict's algorithm identify indels
            mm = BEGIN_DIGIT_M_NUMBER_IorD_NUMBER_M.matcher(cigarStr); //If CIGAR starts with 1-9 bases long matched sequence, insertion or deletion and matched sequence
            if (mm.find()) {
                int tmid = toInt(mm.group(1));
                int mlen = toInt(mm.group(4));
                if (tmid <= 8) {
                    /*
                    If first matched sequence length is no more than 8, this sequence and insertion/deletion are
                    replaced with soft-clipping up to 1st matching base in last matched sequence
                    For deletion position is adjusted by deletion length
                     */
                    int tslen = tmid + (mm.group(3).equals("I") ? toInt(mm.group(2)) : 0);
                    position += tmid + (mm.group(3).equals("D") ? toInt(mm.group(2)) : 0);
                    int n0 = 0;
                    while (n0 < mlen
                            && isHasAndNotEquals(querySeq.charAt(tslen + n0), ref, position + n0)) {
                        n0++;
                    }
                    tslen += n0;
                    mlen -= n0;
                    position += n0;
                    //Regexp replaces digit-M-number-(I or D)-number-M with tslen + S + mlen + M
                    cigarStr = BEGIN_DIGIT_M_NUMBER_IorD_NUMBER_M_.matcher(cigarStr).replaceFirst(tslen + "S" + mlen + "M");
                    flag = true;
                }
            }
            mm = NUMBER_IorD_DIGIT_M_END.matcher(cigarStr);
            if (mm.find()) { //If CIGAR ends with insertion or deletion and 1-9 bases long matched sequence
                int tmid = toInt(mm.group(3));
                if (tmid <= 8) { //If matched sequence length is no more than 8, insertion/deletion and matched sequence are replaced with soft-clipping
                    String tslen = tmid + (mm.group(2).equals("I") ? toInt(mm.group(1)) : 0) + "S";
                    Replacer r = NUMBER_IorD_DIGIT_M_END.replacer(tslen);
                    cigarStr = r.replace(cigarStr);
                    flag = true;
                }
            }

            // Combine two deletions and insertion into one complex if they are close
            mm = D_M_D_DD_M_D_I_D_M_D_DD.matcher(cigarStr);
            if (mm.find()) { //If CIGAR string contains matched sequence, deletion, short (<=9 bases) matched sequence, insertion, short (<=9 bases) matched sequence and deletion
                //length of internal matched sequences
                int mid = toInt(mm.group(4)) + toInt(mm.group(6));
                if (mid <= 10) {
                    //length of both matched sequences and insertion
                    int tslen = mid + toInt(mm.group(5));
                    //length of deletions and internal matched sequences
                    int dlen = toInt(mm.group(3)) + mid + toInt(mm.group(7));
                    //prefix of CIGAR string before the M-D-M-I-M-D complex
                    String ov5 = mm.group(1);
                    //offset of first deletion in the read
                    int rdoff = toInt(mm.group(2));
                    //offset of first deletion in the reference sequence
                    int refoff = position + rdoff;
                    //offset of first deletion in the read corrected by possibly matching bases
                    int RDOFF = rdoff;
                    if (!ov5.isEmpty()) { //If the complex is not at start of CIGAR string
                        rdoff += sum(globalFind(SOFT_CLIPPED, ov5)); // read position
                        refoff += sum(globalFind(ALIGNED_LENGTH, ov5)); // reference position
                    }
                    //number of bases after refoff/rdoff that match in reference and read
                    int rn = 0;
                    while (rdoff + rn < querySeq.length()
                            && isHasAndEquals(querySeq.charAt(rdoff + rn), ref, refoff + rn)) {
                        rn++;
                    }
                    RDOFF += rn;
                    dlen -= rn;
                    tslen -= rn;
                    //If length of internal matched sequences is no more than 10, replace M-D-M-I-M-D complex with M-D-I
                    cigarStr = D_M_D_DD_M_D_I_D_M_D_DD_prim.matcher(cigarStr).replaceFirst(RDOFF + "M" + dlen + "D" + tslen + "I");
                    flag = true;
                }
            }
            // Combine two close deletions (<10bp) into one
            Matcher cm = DIG_D_DIG_M_DIG_DI_DIGI.matcher(cigarStr);
            if (cm.find()) { //If CIGAR string contains deletion, short (<= 9 bases) matched sequence, deletion and possibly insertion
                int g2 = toInt(cm.group(2));
                int g3 = toInt(cm.group(3));
                String op = cm.group(4);

                //length of both deletions and matched sequence
                int dlen = toInt(cm.group(1)) + g2;
                //matched sequence length
                int ilen = g2;
                if (op.equals("I")) {
                    ilen += g3;
                } else { // op == "D"
                    dlen += g3;
                    //insertion string
                    String istr = cm.group(5);
                    if (istr != null) { //If insertion is present after 2nd deletion, add its length to $ilen
                        ilen += toInt(istr.substring(0, istr.length() - 1));
                    }
                }
                //Replace D-M-D-I? complex with deletion and insertion
                cigarStr = cm.replaceFirst(dlen + "D" + ilen + "I");
                flag = true;
            }

            // Combine two close indels (<10bp) into one
            cm = DIG_I_dig_M_DIG_DI_DIGI.matcher(cigarStr);
            if (cm.find()) { //If CIGAR string contains insertion, short (<=9 bases) matched sequence, deletion and possibly insertion
                String op = cm.group(4);
                int g2 = toInt(cm.group(2));
                int g3 = toInt(cm.group(3));

                //length of matched sequence and deletion
                int dlen = g2;
                //length of first insertion and matched sequence
                int ilen = toInt(cm.group(1)) + g2;
                if (op.equals("I")) {
                    ilen += g3;
                } else { // op == "D"
                    dlen += g3;
                    //last insertion string
                    String istr = cm.group(5);
                    if (istr != null) { //If insertion is present after deletion, add its length to ilen
                        ilen += toInt(istr.substring(0, istr.length() - 1));
                    }
                }
                //Replace I-M-D-I? complex with deletion and insertion
                cigarStr = cm.replaceFirst(dlen + "D" + ilen + "I");
                flag = true;
            }
        }

        //The following two clauses to capture sometimes mis-softly clipped reads by aligner
        Matcher mtch = ANY_NUMBER_M_NUMBER_S_END.matcher(cigarStr);
        if (mtch.find()) {
            //prefix of CIGAR string before last matched sequence
            String ov5 = mtch.group(1);
            //length of matched sequence
            int mch = toInt(mtch.group(2));
            //length of soft-clipping
            int soft = toInt(mtch.group(3));
            //offset of soft-clipped sequence in the reference string (position + length of matched)
            int refoff = position + mch;
            //offset of soft-clipped sequence in the read
            int rdoff = mch;
            if (!ov5.isEmpty()) { //If prefix is present
                //Add all matched, insertion and soft-clipped lengths to read position
                rdoff += sum(globalFind(SOFT_CLIPPED, ov5)); // read position
                //Add all matched and deletion lengths to reference position
                refoff += sum(globalFind(ALIGNED_LENGTH, ov5)); // reference position
            }
            //number of bases after refoff/rdoff that match in reference and read sequences
            int rn = 0;
            Set<Character> RN = new HashSet<Character>();
            while( rn + 1 < soft
                    && isHasAndEquals(ref, refoff + rn + 1, querySeq, rdoff + rn + 1)
                    && queryQual.charAt(rdoff + rn + 1) - 33 > lowqual ) {
                rn++;
                RN.add(ref.get(refoff + rn + 1));
            }
            int rn_nt = RN.size(); // don't adjust if homopolymer
            if ( (rn > 3 && rn_nt > 1) || (isHasAndEquals(ref, refoff, querySeq, rdoff))) { //If more than 3 bases match after refoff/rdoff or base at refoff/rdoff match
                mch += rn + 1;
                soft -= rn + 1;
                if (soft > 0) {
                    cigarStr = D_M_D_S_END.matcher(cigarStr).replaceFirst(mch + "M" + soft + "S");
                } else {
                    cigarStr = D_M_D_S_END.matcher(cigarStr).replaceFirst(mch + "M");
                }
            }
            if (rn == 0) {
                while (rn < mch && isHasAndNotEquals(ref, refoff - rn - 1, querySeq, rdoff - rn - 1)) {
                    rn++;
                }
                if (rn > 0 && rn < mch) {
                    soft += rn;
                    mch -= rn;
                    cigarStr = D_M_D_S_END.matcher(cigarStr).replaceFirst(mch + "M" + soft + "S");
                }
            }
        }

        mtch = D_S_D_M.matcher(cigarStr);
        if (mtch.find()) {
            //length of matched sequence
            int mch = toInt(mtch.group(2));
            //length of soft-clipping
            int soft = toInt(mtch.group(1));
            //number of bases before matched sequence that match in reference and read sequences
            int rn = 0;
            Set<Character> RN = new HashSet<Character>();
            while (rn + 1 < soft && isHasAndEquals(ref, position - rn - 2, querySeq, soft - rn - 2)
                    && queryQual.charAt(soft - rn - 2) - 33 > lowqual) {
                rn++;
                RN.add(ref.get(position - rn - 2));
            }
            int rn_nt = RN.size();
            if ((rn > 3 && rn_nt > 1) || (isHasAndEquals(ref, position - 1, querySeq, soft -1 ))) {//If more than 3 bases match before matched sequence or last base of clipped sequence matches
                //Replace the S-M complex with either match or match-soft clip
                mch += rn + 1;
                soft -= rn + 1;
                if (soft > 0) {
                    cigarStr = mtch.replaceFirst(soft + "S" + mch + "M");
                } else {
                    cigarStr = mtch.replaceFirst(mch + "M");
                }
                position -= rn + 1;
                rn++;
            }
            if (rn == 0) {
                while (rn < mch && isHasAndNotEquals(ref, position + rn, querySeq, soft + rn)) {
                    rn++;
                }
                if (rn > 0 && rn < mch) {
                    soft += rn;
                    mch -= rn;
                    cigarStr = mtch.replaceFirst(soft + "S" + mch + "M");
                    position += rn;
                }
            }
        }
        return tuple(position, cigarStr);
    }

    private static boolean isBEGIN_ATGC_AMP_ATGCs_END(String s) {
        if (s.length() > 2) {
            char ch1 = s.charAt(0);
            char ch2 = s.charAt(1);
            if (ch2 == '&' && isATGC(ch1)) {
                for (int i = 2; i < s.length(); i++) {
                    if(!isATGC(s.charAt(i))) {
                        return false;
                    }
                }
                return true;
            }
        }
        return false;
    }

    private static boolean isATGC(char ch) {
        switch (ch) {
            case 'A':
            case 'T':
            case 'G':
            case 'C':
                return true;

            default:
                return false;
        }

    }

    public static void main(String[] args) {
        System.err.println(new DecimalFormat("#.##").format(12.357));
    }

}
