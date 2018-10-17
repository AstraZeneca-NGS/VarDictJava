package com.astrazeneca.vardict;

import com.astrazeneca.vardict.collection.ConcurrentHashSet;
import com.astrazeneca.vardict.collection.Tuple.Tuple2;
import com.astrazeneca.vardict.collection.Tuple.Tuple3;
import com.astrazeneca.vardict.data.*;
import com.astrazeneca.vardict.modules.*;
import com.astrazeneca.vardict.variations.*;
import htsjdk.samtools.*;

import java.io.*;
import java.text.DecimalFormat;
import java.util.*;
import java.util.concurrent.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static com.astrazeneca.vardict.ReferenceResource.getREF;
import static com.astrazeneca.vardict.modules.ToVarsBuilder.toVars;
import static com.astrazeneca.vardict.variations.Variant.callingEmptySimpleVariant;
import static com.astrazeneca.vardict.variations.Variant.joinEmptyVariantWithTcov;
import static com.astrazeneca.vardict.variations.Variant.joinVariantWithNM;
import static com.astrazeneca.vardict.variations.VariationUtils.*;
import static com.astrazeneca.vardict.collection.Tuple.tuple;
import static com.astrazeneca.vardict.Utils.*;
import static com.astrazeneca.vardict.data.Patterns.*;
import static com.astrazeneca.vardict.variations.VariationUtils.VarsType.varn;
import static com.astrazeneca.vardict.variations.VariationUtils.VarsType.var;
import static java.lang.String.format;
import static java.util.Collections.singletonList;

public class VarDict {

    /**
     * The main method to start Vardict in any mode.
     * @param conf prepared configuration object from arguments
     * @throws IOException
     */
    public static void start(Configuration conf) throws IOException {
        if (conf.printHeader) {
            System.out.println(join("\t",
                    "Sample", "Gene", "Chr", "Start", "End", "Ref", "Alt", "Depth", "AltDepth", "RefFwdReads",
                    "RefRevReads", "AltFwdReads", "AltRevReads", "Genotype", "AF", "Bias", "PMean", "PStd",
                    "QMean", "QStd", "MQ", "Sig_Noise", "HiAF", "ExtraAF", "shift3", "MSI", "MSI_NT", "NM",
                    "HiCnt", "HiCov", "5pFlankSeq", "3pFlankSeq", "Seg", "VarType", "Duprate", "SV_info"));
        }

        Tuple2<String, String> stpl = getSampleNames(conf.bam.getBamRaw(), conf.sampleName, conf.sampleNameRegexp);
        String sample = stpl._1;
        String samplem = stpl._2;

        Map<String, Integer> chrs = readChr(conf.bam.getBamX());

        if (conf.regionOfInterest != null) {
            Region region = Region.buildRegion(conf.regionOfInterest, conf.numberNucleotideToExtend, chrs,
                    conf.isZeroBasedDefined() ? conf.zeroBased : false);
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

    /**
     * Method splits list of lines from BED file to list of Regions in non-amplicon mode.
     * @param segraw list of lines from BED file
     * @param chrs map of chromosome lengths
     * @param zeroBased true if coordinates in BED file start from 0
     * @param conf Vardict Configuration (contains bed row format and the number of nucleotides to extend for each segment)
     * @return list of segments split by Regions
     */
    private static List<List<Region>> toRegions(List<String> segraw,
                                                Map<String, Integer> chrs,
                                                Boolean zeroBased,
                                                Configuration conf) {
        boolean zb = conf.isZeroBasedDefined() ? conf.zeroBased : false;
        List<List<Region>> segs = new LinkedList<>();
        BedRowFormat format = conf.bedRowFormat;
        for (String seg : segraw) {
            String[] columnValues = seg.split(conf.delimiter);
            // For Custom Bed Row format columns number must be 4 and parameter -c doesn't set
            if (!conf.isColumnForChromosomeSet() && columnValues.length == 4) {
                try {
                    int startRegion_a1 = toInt(columnValues[1]);
                    int endRegion_a2 = toInt(columnValues[2]);
                    if (startRegion_a1 <= endRegion_a2) {
                        format = CUSTOM_BED_ROW_FORMAT;
                        if (zeroBased == null) {
                            zb = true;
                        }
                    }
                } catch (NumberFormatException e) {
                    System.err.println("Incorrect format of BED file. 2 and 3 columns must contain region start and end.");
                    throw e;
                }
            }
            String chr = columnValues[format.chrColumn];
            chr = correctChr(chrs, chr);
            int cdsStart = toInt(columnValues[format.startColumn]);
            int cdsEnd = toInt(columnValues[format.endColumn]);
            String gene = format.geneColumn < columnValues.length ? columnValues[format.geneColumn] : chr;

            String[] starts = columnValues[format.thickStartColumn].split(",");
            String[] ends = columnValues[format.thickEndColumn].split(",");
            List<Region> cds = new LinkedList<>();
            for (int i = 0; i < starts.length; i++) {
                int regionStart = toInt(starts[i]);
                int regionEnd = toInt(ends[i]);
                if (cdsStart > regionEnd) {
                    continue;
                }
                if (cdsEnd > regionEnd) {
                    break;
                }
                if (regionStart < cdsStart)
                    regionStart = cdsStart;
                if (regionEnd > cdsEnd)
                    regionEnd = cdsEnd;
                regionStart -= conf.numberNucleotideToExtend;
                regionEnd += conf.numberNucleotideToExtend;
                // increment start if zero based parameter true (coordinates start from 0)
                if (zb && regionStart < regionEnd) {
                    regionStart++;
                }
                cds.add(new Region(chr, regionStart, regionEnd, gene));
            }
            segs.add(cds);
        }
        return segs;
    }

    private static List<List<Region>> toRegions(List<String> segraw,
                                                Map<String, Integer> chrs,
                                                boolean zeroBased,
                                                String delimiter) {
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

    /**
     * Method reads BED file line by line and checks if an amplicon mode sets in Configuration or by BED file.
     *
     * @param conf Vardict Configuration (contains amplicon based calling parameter)
     * @return tuple of amplicon parameters (distance to the edges and overlap fraction), zero based parameter
     * and list of lines from BED file
     * @throws IOException
     */
    private static Tuple3<String, Boolean, List<String>> readBedFile(Configuration conf) throws IOException {
        String ampliconParameters = conf.ampliconBasedCalling;
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
                // Check that the amplicon parameters are not set in Configuration
                if (ampliconParameters == null) {
                    String[] columnValues = line.split(conf.delimiter);
                    // Amplicon mode is on only if bed file contains exactly 8 columns and columns #7 and #8 are numbers
                    if (columnValues.length == 8) {
                        Matcher column6Matcher = INTEGER_ONLY.matcher(columnValues[6]);
                        Matcher column7Matcher = INTEGER_ONLY.matcher(columnValues[7]);
                        if (column6Matcher.find() && column7Matcher.find()) {
                            try {
                                int startRegion_a1 = toInt(columnValues[1]);
                                int endRegion_a2 = toInt(columnValues[2]);
                                int startAmplicon_a6 = toInt(columnValues[6]);
                                int endAmplicon_a7 = toInt(columnValues[7]);
                                if (startAmplicon_a6 >= startRegion_a1 && endAmplicon_a7 <= endRegion_a2) {
                                    // Read pair is considered belonging the amplicon if the edges are less than 10 bp
                                    // to the amplicon and overlap fraction is at least 0.95 by default
                                    ampliconParameters = Configuration.DEFAULT_AMPLICON_PARAMETERS;
                                    if (!conf.isZeroBasedDefined()) {
                                        zeroBased = true;
                                    }
                                }
                            } catch (NumberFormatException e) {
                                System.err.println("Incorrect format of BED file for amplicon mode. It must be 8 columns " +
                                        "and 2, 3, 7 and 8 columns must contain region and amplicon starts and ends.");
                                throw e;
                            }
                        }
                    }
                }
                segraw.add(line);
            }
        }
        return tuple(ampliconParameters, zeroBased, segraw);
    }

    private static void nonAmpVardict(List<List<Region>> segs,
                                      Map<String, Integer> chrs,
                                      String ampliconBasedCalling,
                                      String sample,
                                      String samplem,
                                      Configuration conf) throws IOException {

        if (conf.bam.hasBam2()) {
            Tuple2<String, String> stpl = getSampleNames(sample, samplem, conf.bam.getBam1(), conf.bam.getBam2(),
                    conf.sampleName, conf.sampleNameRegexp);
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

    private static void vardictParallel(final List<List<Region>> segs,
                                        final Map<String, Integer> chrs,
                                        final String ampliconBasedCalling,
                                        final String sample,
                                        final Configuration conf) {
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

    private static void vardictNotParallel(final List<List<Region>> segs,
                                           final Map<String, Integer> chrs,
                                           final String ampliconBasedCalling,
                                           final String sample,
                                           final Configuration conf) throws IOException {
        for (List<Region> list : segs) {
            for (Region region : list) {
                Reference reference = getREF(region, chrs, conf);
                final Set<String> splice = new HashSet<>();
                Tuple2<Integer, Map<Integer, Vars>> tpl = toVars(region, conf.bam.getBam1(), reference, chrs, sample,
                        splice, ampliconBasedCalling, 0, conf);
                vardict(region, tpl._2, sample, splice, conf, System.out);
            }
        }
    }

    private static void somaticParallel(final List<List<Region>> segs,
                                        final Map<String, Integer> chrs,
                                        final String ampliconBasedCalling,
                                        final String sample,
                                        String samplem,
                                        final Configuration conf) {
        final ExecutorService executor = Executors.newFixedThreadPool(conf.threads);
        final BlockingQueue<Future<OutputStream>> toSamdict = new LinkedBlockingQueue<>(10);

        executor.submit(new Runnable() {

            @Override
            public void run() {
                try {
                    for (List<Region> list : segs) {
                        for (Region region : list) {
                            final Set<String> splice = new ConcurrentHashSet<>();
                            Reference reference = getREF(region, chrs, conf);
                            Future<Tuple2<Integer, Map<Integer, Vars>>> f1 = executor.submit(new ToVarsWorker(region,
                                    conf.bam.getBam1(), chrs, sample, splice, ampliconBasedCalling, reference, conf));
                            Future<OutputStream> f2 = executor.submit(new SomdictWorker(region,
                                    conf.bam.getBam2(), chrs, splice, ampliconBasedCalling, reference, conf, f1, sample));
                            toSamdict.put(f2);
                        }
                    }
                    toSamdict.put(NULL_FUTURE);
                } catch (InterruptedException e) {
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
                                           final Map<String, Integer> chrs,
                                           final String ampliconBasedCalling,
                                           final String sample,
                                           String samplem,
                                           final Configuration conf) throws IOException {
        for (List<Region> list : segs) {
            for (Region region : list) {
                final Set<String> splice = new ConcurrentHashSet<>();
                Reference reference = getREF(region, chrs, conf);
                Tuple2<Integer, Map<Integer, Vars>> t1 = toVars(region, conf.bam.getBam1(), reference, chrs, sample,
                        splice, ampliconBasedCalling, 0, conf);
                Tuple2<Integer, Map<Integer, Vars>> t2 = toVars(region, conf.bam.getBam2(), reference, chrs, sample,
                        splice, ampliconBasedCalling, t1._1, conf);
                somdict(region, t1._2, t2._2, sample, chrs, splice, ampliconBasedCalling, Math.max(t1._1, t2._1), conf, System.out);
            }
        }
    }

    private static Tuple2<String, String> getSampleNames(String bam,
                                                         String sampleName,
                                                         String regexp) {
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

    private static Tuple2<String, String> getSampleNames(String sample,
                                                         String samplem,
                                                         String bam1,
                                                         String bam2,
                                                         String sampleName,
                                                         String regexp) {
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
                       Map<String, Integer> chrs, Set<String> splice,
                       String ampliconBasedCalling, int rlen,
                       Configuration conf, PrintStream out) throws IOException {

        Set<Integer> ps = new HashSet<>(vars1.keySet());
        ps.addAll(vars2.keySet());
        List<Integer> pp = new ArrayList<>(ps);
        Collections.sort(pp);

        for (Integer p : pp) {
            if (p < segs.start || p > segs.end) {
                continue;
            }
            Vars v1 = vars1.get(p);
            Vars v2 = vars2.get(p);
            if (v1 == null && v2 == null) { // both samples have no coverage
                continue;
            }

            if (v1 == null) { // no coverage for sample 1
                if (v2.var.isEmpty()) {
                    continue;
                }
                for (Variant variant : v2.var) {
                    variant.vartype = variant.varType();
                    if (!isGoodVar(variant, v2.ref, variant.vartype, splice, conf)) {
                        continue;
                    }
                    if (variant.vartype.equals("Complex")) {
                        variant.adjComplex();
                    }

                    String tvf = joinEmptyVariantWithTcov(0);
                    variant.callingOneSample(segs, sample, out, tvf, "Deletion", v2.sv, true);
                }

            } else if (v2 == null) { // no coverage for sample 2
                if (v1.var.isEmpty()) {
                    continue;
                }
                for (Variant variant : v1.var) {
                    variant.vartype = variant.varType();
                    if (!isGoodVar(variant, v1.ref, variant.vartype, splice, conf)) {
                        continue;
                    }
                    if (variant.vartype.equals("Complex")) {
                        variant.adjComplex();
                    }
                    String tvf = joinEmptyVariantWithTcov(0);
                    variant.callingOneSample(segs, sample, out, tvf, "SampleSpecific", v1.sv, false);
                }

            } else { // both samples have coverage
                if (v1.var.isEmpty() && v2.var.isEmpty()) {
                    continue;
                }
                if (v1.var.size() > 0) {
                    int n = 0;
                    while (n < v1.var.size()
                            && isGoodVar(v1.var.get(n), v1.ref, v1.var.get(n).varType(), splice, conf)) {
                        final Variant vref = v1.var.get(n);
                        final String nt = vref.n;

                        vref.vartype = vref.varType();
                        if (vref.vartype.equals("Complex")) {
                            vref.adjComplex();
                        }
                        Variant v2nt = getVarMaybe(v2, varn, nt);
                        if (v2nt != null) {
                            String type;
                            if (isGoodVar(v2nt, v2.ref, vref.vartype, splice, conf)) {
                                if (vref.freq > (1 - conf.lofreq) && v2nt.freq < 0.8d && v2nt.freq > 0.2d) {
                                    type = "LikelyLOH";
                                } else {
                                    if (v2nt.freq < conf.lofreq || v2nt.cov <= 1) {
                                        type = "LikelySomatic";
                                    } else {
                                        type = "Germline";
                                    }
                                }
                            } else {
                                if (v2nt.freq < conf.lofreq || v2nt.cov <= 1) {
                                    type = "LikelySomatic";
                                } else {
                                    type = "AFDiff";
                                }
                            }
                            if (v2nt.isNoise(conf.goodq, conf.lofreq) && vref.vartype.equals("SNV")) {
                                type = "StrongSomatic";
                            }

                            String info = join("\t",
                                    joinVariantWithNM(vref),
                                    joinVariantWithNM(v2nt));
                            vref.constructBothSamplesWithSecondVariant(segs, sample, out, info, v2nt, type, v1.sv, v2.sv);

                        } else { // sample 1 only, should be strong somatic
                            String tvf = "";
                            if (!v2.var.isEmpty()) {
                                Variant v2r = getVarMaybe(v2, var, 0);
                                int tcov = v2r != null && v2r.tcov != 0 ? v2r.tcov : 0;
                                int rfc = v2r != null && v2r.rfc != 0 ? v2r.rfc : 0;
                                int rrc = v2r != null && v2r.rrc != 0 ? v2r.rrc : 0;
                                tvf = join("\t", tcov, 0, rfc, rrc, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
                            } else if (v2.ref != null) {
                                if (v2.ref == null) {
                                    Variant v2m = getVarMaybe(v2, var, 0);
                                    int tcov = v2m != null && v2m.tcov != 0 ? v2m.tcov : 0;
                                    tvf = joinEmptyVariantWithTcov(tcov);
                                } else {
                                    Variant v2ref = v2.ref;
                                    tvf = joinVariantWithNM(v2ref);
                                }
                            }
                            String type = "StrongSomatic";
                            jregex.Matcher mm = MINUS_NUM_NUM.matcher(nt);
                            if (!vref.vartype.equals("SNV") && nt.length() > 10 || mm.find()) {
                                v2nt = new Variant();
                                v2.varn.put(nt, v2nt); // Ensure it's initialized before passing to combineAnalysis
                                if (vref.cov < conf.minr + 3 && !nt.contains("<")) {
                                    Tuple2<Integer, String> tpl = combineAnalysis(vref, v2nt, segs.chr, p, nt, chrs, sample,
                                            splice, ampliconBasedCalling, rlen, conf);
                                    rlen = tpl._1;
                                    String newtype = tpl._2;
                                    if ("FALSE".equals(newtype)) {
                                        n++;
                                        continue;
                                    }
                                    if (newtype.length() > 0) {
                                        type = newtype;
                                    }
                                }
                            }
                            if (type.equals("StrongSomatic")) {
                                String info = join("\t",
                                        joinVariantWithNM(vref),
                                        tvf);
                                vref.constructBothSamplesWithSecondVariant(segs, sample, out, info, vref,
                                        "StrongSomatic", v1.sv, v2.sv);
                            } else {
                                String info = join("\t",
                                        joinVariantWithNM(vref),
                                        joinVariantWithNM(v2nt));
                                vref.constructBothSamplesWithoutSecondVariant(segs, sample, out, info, type, v2nt,
                                        v1.sv, v2.sv);
                            }
                        }
                        n++;
                    }
                    if (n == 0) {
                        if (v2.var.isEmpty()) {
                            continue;
                        }
                        for (Variant v2var : v2.var) {
                            v2var.vartype = v2var.varType();
                            if (!isGoodVar(v2var, v2.ref, v2var.vartype, splice, conf)) {
                                continue;
                            }
                            // potential LOH
                            String nt = v2var.n;
                            Variant v1nt = getVarMaybe(v1, varn, nt);
                            if (v1nt != null) {
                                String type = v1nt.freq < conf.lofreq ? "LikelyLOH" : "Germline";
                                if ("Complex".equals(v2var.vartype)) {
                                    v1nt.adjComplex();
                                }

                                v1nt.vartype = v1nt.varType();
                                String info = join("\t",
                                        joinVariantWithNM(v1nt),
                                        joinVariantWithNM(v2var));
                                v1nt.constructBothSamplesWithSecondVariant(segs, sample, out, info, v2var, type, v1.sv, v2.sv);
                            } else {
                                String th1;
                                Variant v1var = getVarMaybe(v1, var, 0);
                                int tcov = v1var != null && v1var.tcov != 0 ? v1var.tcov : 0;

                                Variant v1ref = v1.ref;
                                int fwd = v1ref != null ? v1ref.fwd : 0;
                                int rev = v1ref != null ? v1ref.rev : 0;
                                th1 = join("\t", tcov, 0, fwd, rev, 0, 0);

                                String genotype = v1var != null ? v1var.genotype :
                                      (v1ref != null ? v1ref.n + "/" + v1ref.n : "N/N");
                                th1 = join("\t", th1, genotype, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);

                                if ("Complex".equals(v2var.vartype)) {
                                    v2var.adjComplex();
                                }
                                String info = join("\t",
                                        th1,
                                        joinVariantWithNM(v2var));
                                v2var.constructBothSamplesWithZeroDuprateSecondVariant(segs, sample, out, info, "StrongLOH", v2.sv);
                            }
                        }
                    }
                } else if (v2.var.size() > 0) { // sample 1 has only reference
                    for (Variant v2var : v2.var) {
                        v2var.vartype = v2var.varType();
                        if (!isGoodVar(v2var, v2.ref, v2var.vartype, splice, conf)) {
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
                        String newtype = "";

                        jregex.Matcher mm = MINUS_NUM_NUM.matcher(nt);
                        if (v2.varn.get(nt).cov < conf.minr + 3 && !nt.contains("<")
                                && nt.length() > 10 && mm.find()) {
                            Tuple2<Integer, String> tpl = combineAnalysis(v2.varn.get(nt), v1nt, segs.chr, p, nt,
                                    chrs, sample, splice, ampliconBasedCalling, rlen, conf);
                            rlen = tpl._1;
                            newtype = tpl._2;
                            if ("FALSE".equals(newtype)) {
                                continue;
                            }
                        }
                        String th1;
                        if (newtype.length() > 0) {
                            type = newtype;
                            th1 = joinVariantWithNM(v1nt);
                        } else {
                            Variant v1ref = v1.ref;
                            th1 = joinVariantWithNM(v1ref);
                        }

                        if ("Complex".equals(v2var.vartype)) {
                            v2var.adjComplex();
                        }
                        String info = join("\t",
                                th1,
                                joinVariantWithNM(v2var));
                        v2var.constructBothSamplesWithZeroDuprateSecondVariant(segs, sample, out, info, type, v2.sv);
                    }
                }
            }
        }

        return rlen;
    }

    /**
     * Taken a likely somatic indels and see whether combine two bam files still support somatic status.
     * This is mainly for Indels that softclipping overhang is too short to positively being called
     * in one bam file, but can be called in the other bam file, thus creating false positives
     *
     * @param var1 variant 1
     * @param var2 variant 2
     * @param chr chromosome
     * @param p position
     * @param nt nucleotide
     *
     * @param chrs map of chromosome lengths
     * @param sample sample name
     * @param splice set of strings representing introns in splice
     * @param ampliconBasedCalling string of maximum_distance:minimum_overlap for amplicon based calling
     * @param rlen max read length
     * @param conf Configuration
     * @return (new <code>rlen</code>, "FALSE" | "")
     * @throws IOException
     */
    static Tuple2<Integer, String> combineAnalysis(Variant var1, Variant var2,
                                                   String chr, int p, String nt,
                                                   Map<String, Integer> chrs,
                                                   String sample, Set<String> splice,
                                                   String ampliconBasedCalling, int rlen,
                                                   Configuration conf) throws IOException {
        if (conf.y) {
            System.err.printf("Start Combine %s %s\n", p, nt);
        }
        // Don't do it for structural variants
        if (var1.ep - var1.sp > conf.SVMINLEN) {
            return tuple(rlen, "");
        }

        Region region = new Region(chr, var1.sp - rlen, var1.ep + rlen, "");
        Reference reference = getREF(region, chrs, conf);
        Tuple2<Integer, Map<Integer, Vars>> tpl = toVars(region, conf.bam.getBam1() + ":" + conf.bam.getBam2(),
                reference, chrs, sample, splice, ampliconBasedCalling, rlen, conf);
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

                if (var2.cov != 0) {
                    var2.pmean = (vref.pmean * vref.cov - var1.pmean * var1.cov) / var2.cov;
                    var2.qual = (vref.qual * vref.cov - var1.qual * var1.cov) / var2.cov;
                    var2.mapq = (vref.mapq * vref.cov - var1.mapq * var1.cov) / var2.cov;
                    var2.hifreq = (vref.hifreq * vref.cov - var1.hifreq * var1.cov) / var2.cov;
                    var2.extrafreq = (vref.extrafreq * vref.cov - var1.extrafreq * var1.cov) / var2.cov;
                    var2.nm = (vref.nm * vref.cov - var1.nm * var1.cov) / var2.cov;
                } else {
                    var2.pmean = 0;
                    var2.qual = 0;
                    var2.mapq = 0;
                    var2.hifreq = 0;
                    var2.extrafreq = 0;
                    var2.nm = 0;
                }

                var2.pstd = true;
                var2.qstd = true;

                if (var2.tcov <= 0) {
                    return tuple(rlen, "FALSE");
                }

                var2.freq = var2.cov / (double)var2.tcov;
                var2.qratio = var1.qratio; // Can't back calculate and should be inaccurate
                var2.genotype = vref.genotype;
                var2.bias = strandBias(var2.rfc, var2.rrc, conf.bias, conf.minb) + ";" +
                            strandBias(var2.fwd, var2.rev, conf.bias, conf.minb);
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
     * Single sample mode variant calling
     * @param region region of interest
     * @param vars map of variations
     * @param sample sample name
     * @param splice set of strings representing spliced regions
     * @param conf configuration
     * @param out output stream
     */
    static void vardict(Region region, Map<Integer, Vars> vars,
                        String sample,
                        Set<String> splice,
                        Configuration conf,
                        PrintStream out) {
        for (Map.Entry<Integer, Vars> ent : vars.entrySet()) {
        //for (int p = region.start; p <= region.end; p++) {
            int p = ent.getKey();
            Vars pv = ent.getValue();
            List<String> vts = new ArrayList<>();
            List<Variant> vrefs = new ArrayList<>();
            if (pv.sv.isEmpty()) {
                if (p < region.start || p > region.end){
                    continue;
                }
            }
            if (pv != null && pv.var.isEmpty()) {
                if (!conf.doPileup) {
                    continue;
                }
                //Variant vref = getVarMaybe(vars, p, ref);
                Variant vref = pv.ref;
                if (vref == null) {
                    callingEmptySimpleVariant(region, sample, out, p);
                    continue;
                }
                vref.vartype = "";
                vrefs.add(vref);
            } else {
                List<Variant> vvar = pv.var;
                //Variant rref = getVarMaybe(vars, p, ref);
                for (Variant vref : vvar) {
                    if (vref.refallele.contains("N")) {
                        continue;
                    }
                    vref.vartype = vref.varType();
                    if (!isGoodVar(vref, pv.ref, vref.vartype, splice, conf)) {
                        if (!conf.doPileup) {
                            continue;
                        }
                    }
                    vrefs.add(vref);
                }
            }
            for (int vi = 0; vi < vrefs.size(); vi++) {
                Variant vref = vrefs.get(vi);
                if ("Complex".equals(vref.vartype)) {
                    vref.adjComplex();
                }
                vref.callingSimpleVariant(sample, region, out, pv.sv);

                if (conf.debug) {
                    out.println("\t" + vref.DEBUG);
                }
            }

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
        Reference reference;

        public ToVarsWorker(Region region, String bam, Map<String, Integer> chrs, String sample, Set<String> splice,
                            String ampliconBasedCalling, Reference reference, Configuration conf) {
            super();
            this.region = region;
            this.bam = bam;
            this.chrs = chrs;
            this.sample = sample;
            this.splice = splice;
            this.ampliconBasedCalling = ampliconBasedCalling;
            this.conf = conf;
            this.reference = reference;
        }

        @Override
        public Tuple2<Integer, Map<Integer, Vars>> call() throws Exception {
            if (reference == null) {
                reference = getREF(region, chrs, conf);
            }
            return toVars(region, bam, reference, chrs, sample, splice, ampliconBasedCalling, 0, conf);
        }
    }

    private static class SomdictWorker implements Callable<OutputStream> {

        final Future<Tuple2<Integer, Map<Integer, Vars>>> first;
        final ToVarsWorker second;
        final String sample;

        public SomdictWorker(Region region, String bam, Map<String, Integer> chrs, Set<String> splice,
                             String ampliconBasedCalling, Reference reference, Configuration conf,
                Future<Tuple2<Integer, Map<Integer, Vars>>> first,
                String sample) {
            this.first = first;
            this.second = new ToVarsWorker(region, bam, chrs, sample, splice, ampliconBasedCalling, reference, conf);
            this.sample = sample;
        }

        @Override
        public OutputStream call() throws Exception {
            Tuple2<Integer, Map<Integer, Vars>> t2 = second.call();
            Tuple2<Integer, Map<Integer, Vars>> t1 = first.get();
            ByteArrayOutputStream baos = new ByteArrayOutputStream();
            PrintStream out = new PrintStream(baos);
            somdict(second.region, t1._2, t2._2, sample, second.chrs, second.splice, second.ampliconBasedCalling,
                    Math.max(t1._1, t2._1), second.conf, out);
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

        public VardictWorker(Region region, Map<String, Integer> chrs, Set<String> splice, String ampliconBasedCalling,
                             String sample, Configuration conf) {
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
            Reference reference = getREF(region, chrs, conf);
            Tuple2<Integer, Map<Integer, Vars>> tpl = toVars(region, conf.bam.getBam1(), reference,
                    chrs, sample, splice, ampliconBasedCalling, 0, conf);
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

        public AmpVardictWorker(Map<Integer, List<Tuple2<Integer, Region>>> pos, Region rg, String sample,
                List<Future<Tuple2<Integer,
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
        public OutputStream call() {
            return null;
        }
    });

    private static void ampVardictParallel(final List<List<Region>> segs, final Map<String, Integer> chrs,
            final String ampliconBasedCalling,
            final String bam1, final String sample, final Configuration conf) {

        final ExecutorService executor = Executors.newFixedThreadPool(conf.threads);
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

    private static void ampVardictNotParallel(final List<List<Region>> segs,
                                              final Map<String, Integer> chrs,
                                              final String ampliconBasedCalling,
                                              final String bam1,
                                              final String sample,
                                              final Configuration conf) throws IOException {

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
                vars.add(toVars(region, bam1, getREF(region, chrs, conf), chrs, sample, splice, ampliconBasedCalling, 0, conf)._2);
                j++;
            }
            ampVardict(rg, vars, pos, sample, splice, conf, System.out);
        }
    }

    /**
     * Amplicon variant calling
     *
     * @param rg region
     * @param vars result of {@link ToVarsBuilder ToVarsBuilder#toVars} calling
     * @param positions map of position =&gt; (list of (region number, region))
     * @param sample sample name
     * @param splice set of strings representing spliced regions
     * @param conf configuration
     * @param out output stream
     */
    static void ampVardict(Region rg, List<Map<Integer, Vars>> vars,
                           Map<Integer, List<Tuple2<Integer, Region>>> positions,
                           final String sample,
                           final Set<String> splice,
                           final Configuration conf,
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
            //vartype may take values SNV (Single Nucleotide Variant), Complex (or MNV (Multiple Nucleotide Variant)),
            // Insertion, Deletion
            String vartype = "SNV";
            boolean flag = false;
            Variant vref;
            List<Variant> vrefList = new ArrayList<>();
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
                    for (Variant tv : l) {
                        vcovs.add(tv.tcov);
                        if (tv.tcov > maxcov) {
                            maxcov = tv.tcov;
                        }
                        vartype = tv.varType();
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
                //The amplicon that has depth less than 1/50 of the max depth will be considered
                // not working and thus not used.
                if (t < maxcov / 50) {
                    nocov++;
                }
            }

            if (gvs.size() > 1) {
                Collections.sort(gvs, GVS_COMPARATOR);
            }
            if (ref.size() > 1) {
                Collections.sort(ref, VAR_TCOV_COMPARATOR);
            }

            if (gvs.isEmpty()) { // Only reference
                if (conf.doPileup) {
                    if (!ref.isEmpty()) {
                        vrefList.add(ref.get(0));
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
                for (Tuple2<Variant, String> goodVariant : gvs) {
                    vrefList.add(goodVariant._1);
                }
            }
            List<Tuple2<Variant, String>> goodVariants = gvs;
            for (int i = 0; i < vrefList.size(); i++) {
                vref = vrefList.get(i);
                if (flag) { // different good variants detected in different amplicons
                    String gdnt = gvs.get(i)._1.n;
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
                    goodVariants = gcnt;
                }

                //bad variants
                List<Tuple2<Variant, String>> badv = new ArrayList<>();
                int gvscnt = goodVariants.size();
                if (gvscnt != v.size() || flag) {
                    for (Tuple2<Integer, Region> amps : v) {
                        int amp = amps._1;
                        Region reg = amps._2;
                        if (goodmap.contains(format("%s-%s-%s", amp, vref.refallele, vref.varallele))) {
                            continue;
                        }
                        // In perl it doesn't used but not commented
                        // my $tref = $vars[$amp]->{ $p }->{ VAR }->[0];
                        if (vref.sp >= reg.istart && vref.ep <= reg.iend) {

                            String regStr = reg.chr + ":" + reg.start + "-" + reg.end;

                            if (vars.get(amp).containsKey(p) && vars.get(amp).get(p).var.size() > 0) {
                                badv.add(tuple(vars.get(amp).get(p).var.get(0), regStr));
                            } else if (vars.get(amp).containsKey(p) && vars.get(amp).get(p).ref != null) {
                                badv.add(tuple(vars.get(amp).get(p).ref, regStr));
                            } else {
                                badv.add(tuple((Variant) null, regStr));
                            }
                        } else if ((vref.sp < reg.iend && reg.iend < vref.ep)
                                || (vref.sp < reg.istart && reg.istart < vref.ep)) { // the variant overlap with amplicon's primer
                            if (gvscnt > 1)
                                gvscnt--;
                        }
                    }
                }
                if (flag && gvscnt < goodVariants.size()) {
                    flag = false;
                }
                vartype = vref.varType();
                if (vartype.equals("Complex")) {
                    vref.adjComplex();
                }
                out.print(join("\t",
                        vref.joinSimpleVariant(sample, rg), goodVariants.get(0)._2, vartype, gvscnt, gvscnt + badv.size(), nocov, flag ? 1 : 0));
                if (conf.debug) {
                    out.print("\t" + vref.DEBUG);
                }
                if (conf.debug) {
                    for (int gvi = 0; gvi < goodVariants.size(); gvi++) {
                        Tuple2<Variant, String> tp = goodVariants.get(gvi);
                        out.print("\tGood" + gvi + " " + join(" ", tp._1.joinVar2(" "), tp._2));
                    }
                    for (int bvi = 0; bvi < badv.size(); bvi++) {
                        Tuple2<Variant, String> tp = badv.get(bvi);
                        out.print("\tBad" + bvi + " " + join(" ", tp._1.joinVar2(" "), tp._2));
                    }
                }
                out.println();
            }
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

    public static void main(String[] args) {
        System.err.println(new DecimalFormat("#.##").format(12.357));
    }
}
