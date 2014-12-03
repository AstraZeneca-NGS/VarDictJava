package com.epam;

import static java.lang.String.format;
import static com.epam.VarDict.VarsType.*;

import java.io.*;
import java.util.*;
import java.util.Map.Entry;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class VarDict {

    final static Pattern SN = Pattern.compile("\\s+SN:(\\S+)");
    final static Pattern LN = Pattern.compile("\\sLN:(\\d+)");

    public static void main(String[] args) throws Exception {

        String xx = "200SIxx4544IDyyy555";
        List<String> l = new LinkedList<>();
        Matcher xxm = Pattern.compile("(\\d+)([A-Z])").matcher(xx);
        while (xxm.find()) {
            l.add(xxm.group(1));
            l.add(xxm.group(2));
        }

        System.err.println(l);


        if (true)
            return;

        Matcher m = D_S_D_ID.matcher("123S77DIxxxxxxxxxxxxxx4456SIyyyy");
        if (m.find()) {
            String tslen = Integer.parseInt(m.group(1)) + Integer.parseInt(m.group(2)) + "S";
            System.err.println(m.replaceFirst(tslen));
        }


        System.err.println("asdfasdf\tSA:Z:ererer  asdasd".matches(".*\\tSA:Z:(\\S+).*"));



        System.err.println("123Ssdfgsdfg".matches("^(\\d+)S.*"));
        System.err.println("saffasf123S".matches("^.*(\\d+)S$"));

//      if ( $a[5] =~ /^(\d+)S/ ) {



        Matcher nmMatcher = NUMBER_MISMATCHES.matcher("nM:i:2345asdd");
        if (nmMatcher.find()) {
            System.err.println(nmMatcher.group(1));
        }



        Matcher matcher = INDEL.matcher("123I45D666ID");
        while(matcher.find()) {
            System.err.println(matcher.group(1));
        }

        Matcher matcher2 = INDEL.matcher("123I45D666ID");
        System.err.println(matcher2.matches());
//        if (matcher2.matches()) {
//            System.err.println("yyy");
//            for (int i = 1; i <= matcher2.groupCount(); i++) {
//                System.err.println(matcher2.group(i));
//            }
//        }

        if (true)
            return;

        Bam bam = new Bam("/Users/kirst/Work/bcbio.coverage/test/data/aligned-reads.bam");

        Map<String, Integer> chr = readCHR(bam.getBamX());

        for (Map.Entry<String, Integer> string : chr.entrySet()) {
            System.out.println(string.getKey() + " == " + string.getValue());
        }
    }


    public void start(Configuration conf) throws IOException {
        Tuple2<String, String> stpl = getSample(conf.bam.getBamRaw(), conf.sampleName, conf.sampleNameRegexp);
        String sample = stpl._1();
        String samplem = stpl._2();
        Map<String, Integer> chrs = readCHR(conf.bam.getBamX());
        if (conf.regionOfInterest != null) {
            List<List<Region>> regions = buildRegions(conf.regionOfInterest, conf.numberNucleotideToExtend, conf.zeroBased);
            finalCycle(regions, chrs, conf.ampliconBasedCalling, sample, samplem, conf);
        } else {
            Tuple2<String, List<String>> tpl = readBedFile(conf);
            String ampliconBasedCalling = tpl._1();
            List<String> segraw = tpl._2();

            if (ampliconBasedCalling != null) {
                aa(segraw, chrs, ampliconBasedCalling, sample, conf);
            } else {
                List<List<Region>> regions = buildRegionsFromFile(segraw, conf);
                finalCycle(regions, chrs, ampliconBasedCalling, sample, samplem, conf);
            }
        }
    }


    public static Map<String, Integer> readCHR(String bam) throws IOException {
        try (SamtoolsReader reader = new SamtoolsReader("view", "-H", bam)) {
            Map<String, Integer> chrs = new HashMap<>();
            String line;
            while ((line = reader.read()) != null) {
                if (line.startsWith("@SQ")) {
                    Matcher sn = SN.matcher(line);
                    Matcher ln = LN.matcher(line);
                    if (sn.find() && ln.find()) {
                        chrs.put(sn.group(1), toInt(ln.group(1)));
                    }
                }
            }
            return chrs;
        }
    }

    public static class Configuration {
        String rawBam;
        String delimiter;
        String bed;
        int numberNucleotideToExtend;
        boolean zeroBased; // -z,  default true if set -R
        String ampliconBasedCalling; //-a
        int columnForChromosome = -1;
        String sampleNameRegexp; // -n
        String sampleName; //-N
        String fasta;
        Bam bam;
        Double downsampling;
        boolean chromosomeNameIsNumber;
        Integer mappingQuality;//-Q
        boolean removeDuplicatedReads; //-t
        int mismatch; //-m, default = 8
        boolean y; //-y TODO ???
        int goodq; // -q, default = 23
        final int buffer = 200;
        int vext = 3; // -X, default 3
        int trimBasesAfter = 0; // -T, Trim bases after [INT] bases in the reads
        boolean performLocalRealignment; // -k, default false
        int indelsize = 120; // -I, default 120
        double bias = 0.05d; // The cutoff to decide whether a positin has read strand bias
        int minb = 2; // -B, default 2. The minimum reads for bias calculation
        int minr = 2; // -r, he minimum # of variance reads, default 2 //TODO -p
        boolean debug = false; // -D
        double freq = 0.5; // -f and -p
        boolean  moveIndelsTo3 = false;
        String samfilter = "0x500"; //-F
        String regionOfInterest; //-R chr:start[-end]
        int readPosFilter = 5; // - P The read position filter, default 5
        double qratio = 1.5; //-o
        double mapq = 0; // -O The minimun mean mapping quality to be considered, default 0
        boolean doPileup = false; // -p Do pileup regarless the frequency

        public boolean isColumnForChromosomeSet() {
            return columnForChromosome >= 0;
        }

        public boolean isDownsampling() {
            return downsampling != null;
        }

        public boolean hasMappingQuality() {
            return mappingQuality != null;
        }



    }

    public static class Region {
        final String chr;
        final int start;
        final int end;
        final String gene;
        int istart;
        int iend;

        public Region(String chr, int start, int end, String gene) {
            this.chr = chr;
            this.start = start;
            this.end = end;
            this.gene = gene;
        }

        public Region(String chr, int start, int end, String gene, int istart, int iend) {
            this(chr, start, end, gene);
            this.istart = istart;
            this.iend = iend;
        }

        public String getChr() {
            return chr;
        }

        public String getGene() {
            return gene;
        }

        public int getStart() {
            return start;
        }

        public int getEnd() {
            return end;
        }

        public int getIStart() {
            return istart;
        }

        public int getIEnd() {
            return iend;
        }
    }

    final static Comparator<Region> ISTART_COMPARATOR = new Comparator<Region>() {
        @Override
        public int compare(Region o1, Region o2) {
            return Integer.compare(o1.istart, o2.istart);
        }
    };

    public static List<List<Region>> buildRegionsFromFile(List<String> segraw, Configuration conf) throws IOException {
        boolean zeroBased = conf.zeroBased;
        List<List<Region>> segs = new LinkedList<>();
        BedRowFormat format = DEFAULT_BED_ROW_FORMAT;
        for (String seg : segraw) {
            String[] splitA = seg.split(conf.delimiter);
            if (!conf.isColumnForChromosomeSet() && splitA.length == 4) {
                try {
                    int a1 = toInt(splitA[1]);
                    int a2 = toInt(splitA[2]);
                    if (a1 <= a2) {
                        format = CUSTOM_BED_ROW_FORMAT;
                        zeroBased = true;
                    }
                } catch (NumberFormatException e) {
                }
            }
            String chr = splitA[format.chrColumn];
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
                if (zeroBased && s < e) {
                    s++;
                }
                cds.add(new Region(chr, s, e, gene));
            }
            segs.add(cds);
        }
        return segs;
    }


    private static void aa(List<String> segraw, Map<String, Integer> chrs, String ampliconBasedCalling, String sample, Configuration conf) throws IOException {
        List<List<Region>> segs = new LinkedList<>();
        Map<String, List<Region>> tsegs = new HashMap<>();
        for (String string : segraw) {
            String[] split = string.split(conf.delimiter);
            String chr = split[0];
            int start = toInt(split[1]);
            int end = toInt(split[2]);
            String gene = split[3];
            int istart = toInt(split[6]);
            int iend = toInt(split[7]);
            if (conf.zeroBased && start < end) {
                start++;
                end++;
            }
            List<Region> list = tsegs.get(chr);
            if (list == null) {
                list = new ArrayList<>();
                tsegs.put(chr, list);
            }
            list.add(new Region(chr, start, end, gene, istart, iend));
        }
        for (Map.Entry<String, List<Region>> entry : tsegs.entrySet()) {
            List<Region> regions = entry.getValue();
            Collections.sort(regions, ISTART_COMPARATOR);
            String pchr = null;
            int pend = -1;
            List<Region> list = new LinkedList<>();
            segs.add(list);
            for (Region region : regions) {
                if (pend != -1 && (!region.getChr().equals(pchr) || region.getIStart() > pend)) {
                    list = new LinkedList<>();
                    segs.add(list);
                }
                list.add(region);
                pchr = region.getChr();
                pend = region.getIEnd();
            }
        }
        ampVardict(segs, chrs, ampliconBasedCalling, sample, conf.bam.getBam1(), conf);
    }


    private static Tuple2<String, List<String>> readBedFile(Configuration conf) throws IOException, FileNotFoundException {
        String a = conf.ampliconBasedCalling;
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
                    if (ampl.length == 8) {
                        try {
                            int a1 = toInt(ampl[1]);
                            int a2 = toInt(ampl[2]);
                            int a6 = toInt(ampl[6]);
                            int a7 = toInt(ampl[7]);
                            if (a6 >= a1 && a7 <= a2) {
                                a = "10:0.95";
                            }
                        } catch (NumberFormatException e) {
                            continue;
                        }
                    }
                    segraw.add(line);
                }

            }
        }
        return Tuple2.newTuple(a, segraw);
    }

    public static class Bam {
        private final String[] bamNames;
        private final String[] bams;
        private final String bamRaw;

        public Bam(String value) {
            bamRaw = value;
            bamNames = value.split("\\|");
            bams = bamNames[0].split(":");
        }

        public String getBam1() {
            return bamNames[0];
        }

        public String getBam2() {
            return hasBam2() ? bamNames[1] : null;
        }

        public String getBamX() {
            return bams[0];
        }

        public boolean hasBam2() {
            return bamNames.length > 1;
        }

        public String getBamRaw() {
            return bamRaw;
        }

    }

    public void finalCycle(List<List<Region>> segs, Map<String, Integer> chrs, String ampliconBasedCalling, String sample, String samplem, Configuration conf) throws IOException {
        if (conf.bam.hasBam2()) {
            Tuple2<String, String> stpl = getSampleFromBam(sample, samplem, conf.bam.getBam1(), conf.bam.getBam2(), conf.sampleName, conf.sampleNameRegexp);
            sample = stpl._1();
            samplem = stpl._2();
        }
        Map<String, Integer> splice = new HashMap<>();
        int rlen = 0;
        for (List<Region> list : segs) {
            for (Region region : list) {
                if (conf.bam.hasBam2()) {
                    Tuple2<Integer, Vars> tpl1 = toVars(region, conf.bam.getBam1(), chrs, splice, ampliconBasedCalling, rlen, conf);
                    Tuple2<Integer, Vars> tpl2 = toVars(region, conf.bam.getBam2(), chrs, splice, ampliconBasedCalling, tpl1._1(), conf);
                    somdict(region, tpl1._2(), tpl2._2(), sample, splice, conf);
                } else {
                    Tuple2<Integer, Vars> tpl = toVars(region, conf.bam.getBam1(), chrs, splice, ampliconBasedCalling, rlen, conf);
                    vardict(region, tpl._2());
                }
            }
        }
    }

    final static Pattern SAMPLE_PATTERN = Pattern.compile("([^\\/\\._]+).sorted[^\\/]*.bam");
    final static Pattern SAMPLE_PATTERN2 = Pattern.compile("([^\\/]+)[_\\.][^\\/]*bam");

    public static Tuple2<String, String> getSample(String bam, String sampleName, String regexp) {
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

        return Tuple2.newTuple(sample, samplem);
    }

    public static Tuple2<String, String> getSampleFromBam(String sample, String samplem, String bam1, String bam2, String sampleName, String regexp) {
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

        return Tuple2.newTuple(sample, samplem);
    }

    private static boolean hasVar(Vars vars, int key) {
        return vars.var.containsKey(key)
                || vars.varn.containsKey(key)
                || vars.ref.containsKey(key);
    }


    private void somdict(Region segs, Vars vars1, Vars vars2, String sample, Map<String, Integer> splice, Configuration conf) {

        double fisherp = 0.01d;

        for (int p = 0; p <= segs.end; p++) {
            if (!hasVar(vars1, p) && !hasVar(vars2, p)) { // both samples have no coverag
                continue;
            }

            String vartype = "";

            if (!hasVar(vars1, p)) { // no coverage for sample 1
                if (vars2.var.containsKey(p)) {
                    Var var = vars2.var.get(p).get(0);
                    vartype = varType(var.refallele, var.varallele);
                    if (!isGoodVar(var, getVarMaybe(vars2, p, VarsType.ref), vartype, splice, conf)) {
                        continue;
                    }
                    if (vartype.equals("Complex")) {
                        adjComplex(var);
                    }
                    System.out.println(join("\t", sample, segs.gene, segs.chr, var.sp, var.ep, var.refallele, var.varallele,
                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                            var.sp, var.ep, var.refallele, var.varallele,
                            var.shift3, var.msi, var.msint, var.leftseq, var.rightseq,
                            segs.chr + ":" + segs.start + "-" + segs.end,
                            "Deletion", vartype));
                }
            } else if (!hasVar(vars2, p)) { // no coverage for sample 2
                if (vars1.var.containsKey(p)) {
                    Var var = vars1.var.get(p).get(0);
                    vartype = varType(var.refallele, var.varallele);
                    if (!isGoodVar(var, getVarMaybe(vars1, p, VarsType.ref), vartype, splice, conf)) {
                        continue;
                    }
                    if (vartype.equals("Complex")) {
                        adjComplex(var);
                    }
                    System.out.println(join("\t", sample, segs.gene, segs.chr, var.sp, var.ep, var.refallele, var.varallele,
                            var.tcov, var.cov, var.rfc, var.rrc, var.fwd, var.rev, var.genotype, var.freq, var.bias, var.pmean, var.pstd, var.qual, var.qstd,
                            var.mapq, var.qratio, var.hifreq, var.extrafreq, var.nm,
                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                            var.shift3, var.msi, var.msint, var.leftseq, var.rightseq,
                            segs.chr + ":" + segs.start + "-" + segs.end,
                            "SampleSpecific", vartype));
                }

            } else { // both samples have coverage
                if (!vars1.var.containsKey(p) && !vars2.var.containsKey(p)) {
                    continue;
                }
                if (vars1.var.containsKey(p)) {
                    int n = 0;
                    List<Var> v1 = vars1.var.get(p);
                    while(n < v1.size() && isGoodVar(v1.get(n), getVarMaybe(vars1, p, VarsType.ref), varType(v1.get(n).refallele, v1.get(n).varallele), splice, conf) ) {
                        Var vref = v1.get(n);
                        String nt = vref.n;
                        if ( nt.length() > 1
                                && vref.refallele.length() == vref.varallele.length()
                                && !vref.genotype.contains("-")
                                && !vref.genotype.contains("m")
                                && !vref.genotype.contains("i") ) {

                            String fnt = substr(nt, 0, -1).replaceFirst("&$", "");
                            String lnt = substr(nt, 1).replaceFirst("^&", "");
                            if (lnt.length() > 1) {
                                lnt = lnt.charAt(0) + "&" + lnt.substring(1);
                            }
                            Var vf = getVarMaybe(vars2, p, varn, fnt);
                            Var vl = getVarMaybe(vars2, p + nt.length() - 2, varn, lnt);
                            if (vf != null && isGoodVar(vf, getVarMaybe(vars2, p, VarsType.ref), null, splice, conf)) {
                                vref.sp += vref.refallele.length() - 1;
                                vref.refallele = substr(vref.refallele, -1);
                                vref.varallele = substr(vref.varallele, -1);
                            } else if (vl != null && isGoodVar(vl, getVarMaybe(vars2, p + nt.length() - 2, VarsType.ref), null, splice, conf)) {
                                vref.ep += vref.refallele.length() - 1;
                                vref.refallele = substr(vref.refallele, 0, -1);
                                vref.varallele = substr(vref.varallele, 0, -1);
                            }

                        }
                        vartype = varType(vref.refallele, vref.varallele);
                        if (vartype.equals("Complex")) {
                            adjComplex(vref);
                        }
                        //TODO !!!

                    }

                }
            }
        }


    }

    private void vardict(Region region, Vars vars) {
        // TODO Auto-generated method stub

    }

    public static String[] retriveSubSeq(String fasta, String chr, int start, int end) throws IOException {

        try (SamtoolsReader reader = new SamtoolsReader("faidx", fasta, chr + ":" + start + "-" + end)) {
            String header = reader.read();
            String exon = reader.read();
            return new String[] { header, exon.replaceAll("\\s+", "") };
        }

    }

    final static Random RND = new Random(System.currentTimeMillis());
    private static final Pattern INDEL = Pattern.compile("(\\d+)[ID]");
    private static final Pattern NUMBER_MISMATCHES = Pattern.compile("(?i)NM:i:(\\d+)");
    private static final Pattern MATCH_INSERTION = Pattern.compile("(\\d+)[MI]");
    private static final Pattern SOFT_CLIPPED = Pattern.compile("(\\d+)[MIS]");
    private static final Pattern ALIGNED_LENGTH = Pattern.compile("(\\d+)[MD]");
    private static final Pattern CIGAR_PAIR = Pattern.compile("(\\d+)([A-Z])");



    public static <K, V> V getOrElse(Map<K, V> map, K key, V or) {
        V v = map.get(key);
        if (v == null) {
            v = or;
            map.put(key, v);
        }
        return v;
    }

    private static Variation getVariation(Sclip sclip, int idx, Character ch) {
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
            map = new HashMap<>();
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

        return  map.get(ref.toString());
    }

    private static void subCnt(Variation vref, boolean dir, int rp, int q, int Q, int nm, Configuration conf) {
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


    private static void addCnt(Variation vref, boolean dir, int rp, int q, int Q, int nm, int goodq) {
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
        Map<Integer, Map<Character, Integer>> nt = new HashMap<>();
        Map<Integer, Map<Character, Variation>> seq = new HashMap<>();
        String sequence;
        boolean used;
    }


    public static class Variation {
        int cnt;
        Map<Boolean, int[]> dir = new HashMap<Boolean, int[]>() {{
            put(false, new int[1]);
            put(true, new int[1]);
        }};
        int pmean;
        int qmean;
        int Qmean;
        int nm;
        int locnt;
        int hicnt;

        boolean pstd;
        boolean qstd;
        int pp;
        int pq;
        int extracnt;

        public void incDir(boolean dir) {
            this.dir.get(dir)[0]++;
        }
        public void decDir(boolean dir) {
            this.dir.get(dir)[0]--;
        }
        public int getDir(boolean dir) {
            return this.dir.get(dir)[0];
        }
        public void addDir(boolean dir, int add) {
            this.dir.get(dir)[0] += add;
        }
        public void subDir(boolean dir, int sub) {
            this.dir.get(dir)[0] -= sub;
        }

    }

    public static class Flag {
        private int flag;

        public Flag(int flag) {
            this.flag = flag;
        }

        public boolean isSupplementaryAlignment() {
            return (flag & 0x800) != 0;
        }

        public boolean isUnmappedMate() {
            return (flag & 0x8) != 0;
        }

        public boolean isReverseStrand() {
            return (flag & 0x10) != 0;
        }

        public boolean isNotPrimaryAlignment() {
            return  (flag & 0x100) != 0;
        }

    }

    private static final Pattern D_S_D_ID = Pattern.compile("^(\\d+)S(\\d+)([ID])");
    private static final Pattern D_ID_D_S = Pattern.compile("(\\d+)([ID])(\\d+)S$");
    private static final Pattern D_S_D_M_ID = Pattern.compile("^(\\d+)S(\\d+)M(\\d+)([ID])");
    private static final Pattern D_ID_D_M_S = Pattern.compile("(\\d+)([ID])(\\d+)M(\\d+)S$");
    private static final Pattern D_M_D_ID_D_M = Pattern.compile("^(\\d)M(\\d+)([ID])(\\d+)M");
    private static final Pattern D_ID_DD_M = Pattern.compile("(\\d+)([ID])(\\d)M$");
    private static final Pattern D_M_D_DD_M_D_I_D_M_D_DD = Pattern.compile("^(.*?)(\\d+)M(\\d+)D(\\d)M(\\d+)I(\\d)M(\\d+)D");
    private static final Pattern D_M_D_DD_M_D_I_D_M_D_DD_prim = Pattern.compile("(\\d+)M(\\d+)D(\\d)M(\\d+)I(\\d)M(\\d+)D");
    private static final Pattern D_MIS = Pattern.compile("(\\d+)[MIS]");
    private static final Pattern D_MD = Pattern.compile("(\\d+)[MD]");
    private static final Pattern D_DD_D_M_D_DD_DI = Pattern.compile("(\\d+)D(\\d)M(\\d+)D(\\d+I)?");
    private static final Pattern D_I_DD_M_D_d_DI = Pattern.compile("(\\d+)I(\\d)M(\\d+)D(\\d+I)?");

    private static final Pattern BEGIN_DIGITS = Pattern.compile("^(\\d+)");
    private static final Pattern HASH_GROUP_CARET_GROUP = Pattern.compile("#(.+)\\^(.+)");


    public static class SamtoolsReader implements AutoCloseable {

        private Process proc;
        private BufferedReader reader;

        public SamtoolsReader(String... args) throws IOException {
            List<String> list = new ArrayList<String>(1 + args.length);
            list.add("samtools");
            for (String arg : args) {
                list.add(arg);
            }
            ProcessBuilder builder = new ProcessBuilder(list);
            builder.redirectErrorStream(true);
            proc = builder.start();
            reader = new BufferedReader(new InputStreamReader(proc.getInputStream()));
        }

        public String read() throws IOException {
            return reader.readLine();
        }

        @Override
        public void close() throws IOException {
            reader.close();
            proc.getInputStream().close();
            proc.getOutputStream().close();
            proc.destroy();
        }
    }

    public static class Samrecord {

        final String qname; //0
        final Flag flag;//1
        final String rname;//2
        final int position;//3
        final int mapq;//4
        final String cigar; //5
        final String mrnm; //6
        final int mpos; //7
        final int isize; //8
        final String querySeq; //9
        final String queryQual; //10

        int recordSize;

        public Samrecord(String record) {
            String[] a = record.split("\t");
            recordSize = a.length;

            qname = a[0];
            flag = new Flag(toInt(a[1]));
            rname = a.length > 2 ? a[2] : null;
            position = a.length > 3 ? toInt(a[3]) : 0;
            mapq = a.length > 4 ? toInt(a[4]) : 0;
            cigar = a.length > 5 ? a[5] : "";
            mrnm = a.length > 6 ? a[6] : "";
            mpos = a.length > 7 ? toInt(a[7]) : 0;
            isize = a.length > 8 ? toInt(a[8]) : 0;
            querySeq = a.length > 9 ? a[9] : "";
            queryQual = a.length > 10 ? a[10] : "";
        }

        public boolean isDefined(int num) {
            return num < recordSize;
        }

    }


    private static Tuple2<Integer, Vars> toVars(Region region, String bam, Map<String, Integer> chrs, Map<String, Integer> SPLICE, String ampliconBasedCalling, int Rlen, Configuration conf) throws IOException {
        String[] bams = bam.split(":");
        Map<Integer, Character> ref = new HashMap<>();
        Map<Integer, Map<String, Variation>> hash = new HashMap<>();
        Map<Integer, Map<String, Variation>> iHash = new HashMap<>();
        Map<Integer, Integer> cov = new HashMap<>();
        Map<Integer, Sclip> sclip3 = new HashMap<>();
        Map<Integer, Sclip> sclip5 = new HashMap<>();
        Map<Integer, Map<String, Integer>> ins = new HashMap<>();
        Map<Integer, Map<String, Integer>> mnp = new HashMap<>(); // Keep track of MNPs
        Map<Integer, Map<String, Integer>> dels5 = new HashMap<>(); // Keep track of MNPs
        for (String bami : bams) {
            int s_start = region.getStart() - conf.numberNucleotideToExtend - 700 < 1 ? 1 : region.getStart() - conf.numberNucleotideToExtend - 700;
            int len = chrs.containsKey(region.getChr()) ? chrs.get(region.getChr()) : 0;
            int s_end = region.getEnd() + conf.numberNucleotideToExtend + 700 > len ?
                    len : region.getEnd() + conf.numberNucleotideToExtend + 700;

            String[] subSeq = retriveSubSeq(conf.fasta, region.getChr(), s_start, s_end);
            String header = subSeq[0];
            String exon = subSeq[1];
            for(int i = s_start; i <= s_start + exon.length(); i++) { //TODO why '<=' ?
                ref.put(i, Character.toUpperCase(exon.charAt(i - s_start)));
            }

            String chr = region.getChr();
            if (conf.chromosomeNameIsNumber && chr.startsWith("chr")) {
                chr = region.getChr().substring("chr".length());
            }
            String samfilter = conf.samfilter != null ? "-F " + conf.samfilter : "";
            try (SamtoolsReader reader = new SamtoolsReader("faidx", "view", samfilter, bami, chr + ":" + s_start + "-" + s_end)) {
                Map<String, Boolean> dup = new HashMap<>();
                int dupp = -1;
                String line;
                while ((line = reader.read()) != null) {
                    if (conf.isDownsampling() && RND.nextDouble() <= conf.downsampling) {
                        continue;
                    }
                    Samrecord row = new Samrecord(line);


                    if (conf.hasMappingQuality() && row.mapq < conf.mappingQuality) { // ignore low mapping quality reads
                        continue;
                    }

                    if ("*".equals(row.querySeq)) {
                        continue;
                    }

                    // filter duplicated reads if option -t is set
                    if (conf.removeDuplicatedReads) {
                        if (row.position != dupp) {
                            dup.clear();
                        }
                        if (row.isDefined(7) && row.mpos < 10) {
                            if (dup.containsKey(row.position + "-" + row.mrnm + "-" + row.mpos)) {
                                continue;
                            }
                            dup.put(row.position + "-" + row.mrnm + "-" + row.mpos, Boolean.TRUE);
                            dupp = row.position;
                        } else if (row.flag.isUnmappedMate()) {
                            if (dup.containsKey(row.position + "-" + row.cigar)) {
                                continue;
                            }
                            dup.put(row.position + "-" + row.cigar, Boolean.TRUE);
                            dupp = row.position;
                        }
                    }

                    int nm = 0;
                    int idlen = 0;
                    for (String s : globalFind(INDEL, row.cigar)) {
                        idlen += toInt(s);
                    }
                    Matcher nmMatcher = NUMBER_MISMATCHES.matcher(line);
                    if (nmMatcher.find()) { // number of mismatches. Don't use NM since it includes gaps, which can be from indels
                        nm = toInt(nmMatcher.group(1)) - idlen;
                        if (nm > conf.mismatch) { // edit distance - indels is the # of mismatches
                            continue;
                        }
                    } else {
                        if (conf.y && !row.cigar.equals("*")) {
                            System.err.println("No XM tag for mismatches. " + line);
                            continue;
                        }
                    }
                    int n = 0; // keep track the read position, including softclipped
                    int p = 0; // keep track the position in the alignment, excluding softclipped
                    boolean dir = row.flag.isReverseStrand();
                    if ( ampliconBasedCalling != null ) {
                        String[] split = ampliconBasedCalling.split(":");
                        int dis;
                        double ovlp;
                        try {
                            dis = toInt(split[0]);
                            ovlp = Double.parseDouble(split.length > 1 ? split[1] : "");
                        } catch(NumberFormatException e) {
                            dis = 10;
                            ovlp = 0.95;
                        }
                        int rlen3 = sum(globalFind(ALIGNED_LENGTH, row.cigar)); // The total aligned length, excluding soft-clipped bases and insertions
                        int segstart = row.position;
                        int segend = segstart + rlen3 - 1;

                        if (row.cigar.matches("^(\\d+)S.*")) {
                            int ts1 = segstart > s_start ? segstart : s_start;
                            int te1 = segend < s_end ? segend : s_end;
                            if (Math.abs(ts1 - te1) / (double)(segend - segstart) > ovlp == false) {
                                continue;
                            }
                        } else if (row.cigar.matches("^.*(\\d+)S$")) {
                            int ts1 = segstart > s_start ? segstart : s_start;
                            int te1 = segend < s_end ? segend : s_end;
                            if (Math.abs(te1 - ts1) / (double)(segend - segstart) > ovlp == false) {
                                continue;
                            }

                        } else {
                          if (row.mrnm.equals("=") && row.isDefined(8)) {
                              if (row.isize > 0) {
                                  segend = segstart + row.isize -1;
                              } else {
                                  segstart = row.mpos;
                                  segend = segstart - row.isize - 1;
                              }
                          }
                          // No segment overlapping test since samtools should take care of it
                          int ts1 = segstart > s_start ? segstart : s_start;
                          int te1 = segend < s_end ? segend : s_end;
                          if ((((Math.abs(segstart - s_start) <= dis) && (Math.abs(segend - s_end) <= dis))
                                  && Math.abs(ts1 - te1) / (double)(segend - segstart) > ovlp) == false) {
                              continue;
                          }
                        }
                    }
                    if (row.flag.isUnmappedMate()) {
                        // to be implemented
                    } else {
                        if (row.mrnm.equals("=")) {
                            if (line.matches(".*\\tSA:Z:(\\S+).*")) {
                                if (row.flag.isSupplementaryAlignment()) { // the supplementary alignment
                                    continue; // Ignore the supplmentary for now so that it won't skew the coverage
                                }
                            }
                        }

                    }

                    //Modify the CIGAR for potential mis-alignment for indels at the end of reads to softclipping and let VarDict's algorithm to figure out indels
                    int position = row.position;
                    String cigarStr = row.cigar;
                    int offset = 0;
                    boolean flag = true;
                    while (flag) {
                        flag = false;
                        Matcher cm = D_S_D_ID.matcher(cigarStr);
                        if (cm.find()) {
                            String tslen = toInt(cm.group(1)) + (cm.group(3).equals("I") ? toInt(cm.group(2)) : 0) + "S";
                            position = cm.group(3).equals("D") ? 2 : 0;
                            cigarStr = cm.replaceFirst(tslen);
                            flag = true;
                        }
                        cm = D_ID_D_S.matcher(cigarStr);
                        if (cm.find()) {
                            String tslen = toInt(cm.group(3)) + (cm.group(2).equals("I") ? toInt(cm.group(1)) : 0) + "S";
                            cigarStr = cm.replaceFirst(tslen);
                            flag = true;
                        }
                        cm = D_S_D_M_ID.matcher(cigarStr);
                        if (cm.find()) {
                            int tmid = toInt(cm.group(2));
                            if (tmid <= 10 ) {
                                String tslen = toInt(cm.group(1)) + tmid + (cm.group(4).equals("I") ? toInt(cm.group(3)) : 0) + "S";
                                position = tmid + (cm.group(4).equals("D") ? toInt(cm.group(3)) : 0);
                                cigarStr = cm.replaceFirst(tslen);
                                flag = true;
                            }
                        }
                        cm = D_ID_D_M_S.matcher(cigarStr);
                        if (cm.find()) {
                            int tmid = toInt(cm.group(3));
                            if (tmid <= 10 ) {
                                String tslen = toInt(cm.group(4)) + tmid + (cm.group(2).equals("I") ? toInt(cm.group(1)) : 0) + "S";
                                cigarStr = cm.replaceFirst(tslen);
                                flag = true;
                            }
                        }



                        // The following two clauses to make indels at the end of reads as softly
                        // clipped reads and let VarDict's algorithm identify indels
                        cm = D_M_D_ID_D_M.matcher(cigarStr);
                        if (cm.find()) {
                            int tmid = toInt(cm.group(1));
                            int mlen = toInt(cm.group(4));
                            if (tmid <= 8 ) {
                                int tslen = tmid + (cm.group(3).equals("I") ? toInt(cm.group(2)) : 0);
                                position = tmid + (cm.group(3).equals("D") ? toInt(cm.group(2)) : 0);
                                int n0 = 0;
                                while (n0 < mlen && Character.valueOf(row.querySeq.charAt(tslen + n0)).equals(ref.get(position + n0))) {
                                    n0++;
                                }
                                tslen += n0;
                                mlen -= n0;
                                position += n0;
                                cigarStr.replaceFirst("^\\dM\\d+[ID]\\d+M", tslen + "S" + mlen);
                                flag = true;
                            }
                        }
                        cm = D_ID_DD_M.matcher(cigarStr);
                        if (cm.find()) {
                            int tmid = toInt(cm.group(3));
                            if (tmid <= 8 ) {
                                String tslen = tmid + (cm.group(2).equals("I") ? toInt(cm.group(1)) : 0) + "S";
                                cigarStr = cm.replaceFirst(tslen);
                                flag = true;
                            }
                        }

                        //Combine two deletions and insertion into one complex if they are close
                        cm = D_M_D_DD_M_D_I_D_M_D_DD.matcher(cigarStr);
                        if (cm.find()) {
                            int mid = toInt(cm.group(4)) + toInt(cm.group(6));
                            if (mid <= 10) {
                                int tslen = mid + toInt(cm.group(5));
                                int dlen = toInt(cm.group(3)) + mid + toInt(cm.group(7));
                                String ov5 = cm.group(1);
                                int rdoff = toInt(cm.group(2));
                                int refoff = position + rdoff;
                                int RDOFF = rdoff;
                                if (!ov5.isEmpty()) {
                                    Matcher matcher = D_MIS.matcher(ov5);
                                    while (matcher.find()) {
                                        rdoff += toInt(matcher.group(1));
                                    }
                                    matcher = D_MD.matcher(ov5);
                                    while (matcher.find()) {
                                        refoff += toInt(matcher.group(1));
                                    }
                                }
                                int rn = 0;
                                while (rdoff + rn < row.querySeq.length()
                                        && Character.valueOf(row.querySeq.charAt(rdoff + rn)).equals(ref.get(refoff + rn))) {
                                    rn++;
                                }
                                RDOFF += rn;
                                dlen -= rn;
                                tslen -= rn;
                                cigarStr = D_M_D_DD_M_D_I_D_M_D_DD_prim.matcher(cigarStr).replaceFirst(RDOFF + "M" + dlen + "D" + tslen + "I");
                                flag = true;
                            }
                        }
                        // Combine two close deletions (<10bp) into one
                        cm = D_DD_D_M_D_DD_DI.matcher(cigarStr);
                        if (cm.find()) {
                            int ilen = toInt(cm.group(2));
                            int dlen = ilen + toInt(cm.group(1)) + toInt(cm.group(3));
                            String istr = cm.group(4);
                            if (istr != null) {
                                ilen += toInt(istr.substring(0, istr.length() - 1));
                            }
                            cigarStr = cm.replaceFirst(dlen + "D" + ilen + "I");
                            flag = true;
                        }

                        // Combine two close indels (<10bp) into one
                        cm = D_I_DD_M_D_d_DI.matcher(cigarStr);
                        if (cm.find()) {
                            int dlen = toInt(cm.group(2)) + toInt(cm.group(3));
                            int ilen = toInt(cm.group(1)) + toInt(cm.group(2));
                            String istr = cm.group(4);
                            if (istr != null) {
                                ilen += toInt(istr.substring(0, istr.length() - 1));
                            }
                            cigarStr = cm.replaceFirst(dlen + "D" + ilen);
                            flag = true;
                        }

                    }

                    List<String> cigar = new ArrayList<>();
                    Matcher cigarM = CIGAR_PAIR.matcher(cigarStr);
                    while (cigarM.find()) {
                        cigar.add(cigarM.group(1));
                        cigar.add(cigarM.group(2));
                    }

                    int start = position;

                    List<String> segs = globalFind(MATCH_INSERTION, cigarStr); // Only match and insertion counts toward read length
                    List<String> segs2 = globalFind(SOFT_CLIPPED, cigarStr); // For total length, including soft-clipped bases
                    int rlen = sum(segs); //The read length for matched bases
                    int rlen2 = sum(segs2); //The total length, including soft-clipped bases
                    if (rlen2 > Rlen) { // Determine the read length
                        Rlen = rlen2;
                    }

                    for(int ci = 0; ci < cigar.size(); ci += 2) {
                        int m = toInt(cigar.get(ci));
                        String operation = cigar.get(ci + 1);
                        // Treat insertions at the edge as soft-clipping
                        if ( (ci == 0 || ci == cigar.size() - 2) && operation.equals("I")) {
                            operation = "S";
                        }

                        switch (operation) {
                            case "N":
                                String key = (start - 1) + "-" + (start + m - 1);
                                Integer splice = SPLICE.get(key);
                                if (splice == null) {
                                    SPLICE.put(key, 1);
                                } else {
                                    SPLICE.put(key, splice + 1);
                                }

                                start += m;
                                offset = 0;
                                continue;

                            case "S":
                                if (ci == 0) { // 5' soft clipped
                                    // align softclipped but matched sequences due to mis-softclipping
                                    while (m - 1 >= 0 && start - 1 > 0 && start - 1 <= chrs.get(chr)
                                            && ref.containsKey(start - 1)
                                            && ref.get(start - 1).equals(row.querySeq.charAt(m - 1))
                                            && row.queryQual.charAt(m - 1) - 33 > 10) {
                                        Variation variation = getVariation(hash, start - 1, ref.get(start - 1).toString());
                                        if (variation.cnt != 0) {
                                            variation.cnt = 0;
                                        }
                                        addCnt(variation, dir, m, row.queryQual.charAt(m - 1) - 33, row.mapq, nm, conf.goodq);
                                        incCnt(cov, start - 1, 1);
                                        start--;
                                        m--;
                                    }
                                    if (m > 0) {
                                        int q = 0;
                                        int qn = 0;
                                        int lowqcnt = 0;
                                        for (int si = m - 1; si >= 0; si--) {
                                            if (row.querySeq.charAt(si) == 'N') {
                                                break;
                                            }
                                            int tq = row.queryQual.charAt(si) - 33;
                                            if (tq <= 12)
                                                lowqcnt++;
                                            if (lowqcnt > 1)
                                                break;

                                            q += tq;
                                            qn++;
                                        }
                                        if (qn >= 1 && qn > lowqcnt && start >= s_start - conf.buffer && start <= s_end + conf.buffer) {
                                            Sclip sclip = sclip5.get(start);
                                            if (sclip == null) {
                                                sclip = new Sclip();
                                                sclip5.put(start, sclip);
                                            }
                                            for (int si = m - 1; m - si <= qn; si--) {
                                                Character ch = row.querySeq.charAt(si);
                                                int idx = m - 1 - si;
                                                Map<Character, Integer> cnts = sclip.nt.get(idx);
                                                if (cnts == null) {
                                                    cnts = new HashMap<>();
                                                    sclip.nt.put(idx, cnts);
                                                }
                                                incCnt(cnts, ch, 1);
                                                Variation variation = getVariation(sclip, idx, ch);
                                                if (variation.cnt != 0) {
                                                    variation.cnt = 0;
                                                }
                                                addCnt(variation, dir, si - (m - qn), row.queryQual.charAt(si) - 33, row.mapq, nm, conf.goodq);
                                            }
                                            if (sclip.cnt != 0)
                                                sclip.cnt = 0;
                                            addCnt(sclip, dir, m, q / qn, row.mapq, nm, conf.goodq);
                                        }

                                    }
                                    m = toInt(cigar.get(ci));
                                } else if (ci == cigar.size() - 2) { // 3' soft clipped
                                    int qmean = row.queryQual.charAt(n) - 33;
                                    while (n < row.querySeq.length()
                                            && ref.containsKey(start)
                                            && ref.get(start).equals(row.querySeq.charAt(n))
                                            && qmean > 10) {

                                        Variation variation = getVariation(hash, start, ref.get(start).toString());
                                        if (variation.cnt != 0) {
                                            variation.cnt = 0;
                                        }
                                        addCnt(variation, dir, rlen2 - p, qmean, row.mapq, nm, conf.goodq);
                                        incCnt(cov, start, 1);

                                        n++;
                                        start++;
                                        m--;
                                        p++;
                                    }
                                    if (row.querySeq.length() - n > 0) {
                                        int q = 0;
                                        int qn = 0;
                                        int lowqcnt = 0;
                                        for (int si = 0; si < m; si++) {
                                            if ( row.querySeq.charAt(n+si) == 'N' ) {
                                                break;
                                            }
                                            int tq = row.queryQual.charAt(n + si) - 33;
                                            if (tq <= 12) {
                                                lowqcnt++;
                                            }
                                            if ( lowqcnt > 1 ) {
                                                break;
                                            }
                                            q += tq;
                                            qn++;
                                        }
                                        if (qn >= 1 && qn > lowqcnt && start >= s_start - conf.buffer && start <= s_end + conf.buffer) {
                                            Sclip sclip = sclip3.get(start);
                                            if (sclip == null) {
                                                sclip = new Sclip();
                                                sclip3.put(start, sclip);
                                            }
                                            for (int si = 0; si < qn; si++) {
                                                Character ch = row.querySeq.charAt(n + si);
                                                int idx = si;
                                                Map<Character, Integer> cnts = sclip.nt.get(idx);
                                                if (cnts == null) {
                                                    cnts = new HashMap<>();
                                                    sclip.nt.put(idx, cnts);
                                                }
                                                incCnt(cnts, ch, 1);
                                                Variation variation = getVariation(sclip, idx, ch);
                                                if (variation.cnt != 0) {
                                                    variation.cnt = 0;
                                                }
                                                addCnt(variation, dir, qn -  si, row.queryQual.charAt(n + si) - 33, row.mapq, nm, conf.goodq);
                                            }
                                            if (sclip.cnt != 0)
                                                sclip.cnt = 0;
                                            addCnt(sclip, dir, m, q / qn, row.mapq, nm, conf.goodq);
                                        }

                                    }

                                }
                                n += m;
                                offset = 0;
                                start = row.position;  // had to reset the start due to softclipping adjustment
                                continue;
                            case "H":
                                continue;
                            case "I": {
                                StringBuilder s = new StringBuilder(substr(row.querySeq, n, m));
                                StringBuilder q = new StringBuilder(substr(row.queryQual, n, m));
                                StringBuilder ss = new StringBuilder();

                                // For multiple indels within 10bp
                                int multoffs = 0;
                                int multoffp = 0;

                                if (cigar.size() > ci + 5
                                        && toInt(cigar.get(ci + 2)) <= conf.vext
                                        && cigar.get(ci + 3).contains("M")
                                        && (cigar.get(ci + 5).contains("I") || cigar.get(ci + 5).contains("D"))) {

                                    int ci2 = toInt(cigar.get(ci + 2));
                                    int ci4 = toInt(cigar.get(ci + 4));
                                    s.append("#").append(substr(row.querySeq, n + m, ci2));
                                    q.append(substr(row.queryQual, n + m, ci2));
                                    s.append('^').append(cigar.get(ci + 5).equals("I") ? substr(row.querySeq, n + m + ci2, ci4) : ci4);
                                    q.append(cigar.get(ci + 5).equals("I") ? substr(row.queryQual, n + m + ci2, ci4) : row.queryQual.charAt(n + m + ci2));
                                    multoffs += ci2 + (cigar.get(ci + 5).equals("D") ? ci4 : 0);
                                    multoffp += ci2 + (cigar.get(ci + 5).equals("I") ? ci4 : 0);
                                    ci += 4;
                                } else {
                                    if (cigar.size() > ci + 3 && cigar.get(ci + 3).contains("M")) {
                                        int ci2 = toInt(cigar.get(ci + 2));
                                        for (int vi = 0; vi < conf.vext && vi < ci2; vi++) {
                                            if (row.querySeq.charAt(n + m + vi) == 'N') {
                                                break;
                                            }
                                            if (row.queryQual.charAt(n + m + vi) - 33 < conf.goodq) {
                                                break;
                                            }
                                            if (ref.containsKey(start + vi) && String.valueOf(row.querySeq.charAt(n + m + vi)).equals(ref.get(start + vi))) {
                                                offset = vi + 1;
                                            }
                                        }
                                        if (offset != 0) {
                                            ss.append(substr(row.querySeq, n+m, offset));
                                            q.append(substr(row.queryQual, n+m, offset));
                                            for(int osi = 0; osi < offset; osi++ ) {
                                                incCnt(cov, start + osi, 1);
                                            }
                                        }
                                    }
                                }
                                if ( offset > 0) {
                                    s.append("&").append(ss);
                                }

                                if (start - 1 >= s_start && start - 1 <= s_end && !s.toString().contains("N")) {
                                    incCnt(getOrElse(ins, start - 1, new HashMap<String, Integer>()), "+" + s, 1);
                                    Variation hv = getVariation(iHash, start - 1, "+" + s);
                                    hv.incDir(dir);
                                    hv.cnt++;
                                    int tp = p < rlen - p ? p + 1 : rlen - p;
                                    int tmpq = 0;
                                    for (int i = 0; i < q.length(); i++) {
                                        tmpq += q.charAt(i) - 33;
                                    }
                                    tmpq /= q.length();
                                    if (hv.pstd == false && hv.pp != 0 && tp != hv.pp) {
                                        hv.pstd = true;
                                    }
                                    if (hv.qstd == false && hv.pq != 0 && tmpq != hv.pq) {
                                        hv.qstd = true;
                                    }
                                    hv.pmean += tp;
                                    hv.qmean += tmpq;
                                    hv.Qmean += row.mapq;
                                    hv.pp = tp;
                                    hv.pq = tmpq;
                                    if (tmpq >= conf.goodq) {
                                        hv.hicnt++;
                                    } else {
                                        hv.locnt++;
                                    }
                                    hv.nm += nm;

                                    // Adjust the reference count for insertion reads
                                    if ( ref.containsKey(start - 1)
                                            && hash.containsKey(start - 1)
                                            && hash.get(start - 1).containsKey(ref.get(start - 1))
                                            && String.valueOf(row.querySeq.charAt(n-1)).equals(ref.get(start - 1))) {

//                                        subCnt(getVariation(hash, start - 1, ref.get(start - 1 ).toString()), dir, tp, tmpq, Qmean, nm, conf);
                                        subCnt(getVariation(hash, start - 1, String.valueOf(row.querySeq.charAt(n - 1))), dir, tp, row.queryQual.charAt(n - 1) - 33, row.mapq, nm, conf);
                                    }
                                    // Adjust count if the insertion is at the edge so that the AF won't > 1
                                    if (ci == 2 && (cigar.get(1).contains("S") || cigar.get(1).contains("H"))) {
                                        Variation ttref = getVariation(hash, start - 1, ref.get(start - 1).toString());
                                        ttref.incDir(dir);
                                        ttref.cnt++;
                                        ttref.pstd = hv.pstd;
                                        ttref.qstd = hv.qstd;
                                        ttref.pmean += tp;
                                        ttref.qmean += tmpq;
                                        ttref.Qmean += row.mapq;
                                        ttref.pp = tp;
                                        ttref.pq = tmpq;
                                        ttref.nm += nm;
                                        incCnt(cov, start - 1, 1);
                                    }
                                }
                                n += m + offset + multoffp;
                                p += m + offset + multoffp;
                                start += offset + multoffs;
                                }
                                continue;
                            case "D":
                                {
                                    StringBuilder s = new StringBuilder("-").append(m);
                                    StringBuilder ss = new StringBuilder();
                                    char q1 = row.queryQual.charAt(n - 1);
                                    StringBuilder q = new StringBuilder();

                                    // For multiple indels within $VEXT bp
                                    int multoffs = 0;
                                    int multoffp = 0;
                                    if (cigar.size() > ci + 5
                                            && toInt(cigar.get(ci + 2)) <= conf.vext
                                            && cigar.get(ci + 3).contains("M")
                                            && (cigar.get(ci + 5).contains("I") || cigar.get(ci + 5).contains("D"))) {

                                        int ci2 = toInt(cigar.get(ci + 2));
                                        int ci4 = toInt(cigar.get(ci + 4));

                                        s.append("#").append(substr(row.querySeq, n, ci2));
                                        q.append(substr(row.queryQual, n, ci2));
                                        s.append('^').append(cigar.get(ci + 5).equals("I") ? substr(row.querySeq, n + ci2, ci4) : ci4);
                                        q.append(cigar.get(ci + 5).equals("I") ? substr(row.queryQual, n + ci2, ci4) : "");
                                        multoffs += ci2 + (cigar.get(ci + 5).equals("D") ? ci4 : 0);
                                        multoffp += ci2 + (cigar.get(ci + 5).equals("I") ? ci4 : 0);
                                        ci += 4;
                                    } else if (cigar.size() > ci + 3 && cigar.get(ci + 3).equals("I")) {
                                       int ci2 = toInt(cigar.get(ci + 2));
                                       s.append("^").append(substr(row.querySeq, n, ci2));
                                       q.append(substr(row.queryQual, n, ci2));
                                       multoffp += ci2;
                                       ci += 2;
                                    } else {
                                        if (cigar.size() > ci + 3 && cigar.get(ci + 3).contains("M")) {
                                            int ci2 = toInt(cigar.get(ci + 2));
                                            for (int vi = 0; vi < conf.vext && vi < ci2; vi++) {
                                                if (row.querySeq.charAt(n + vi) == 'N') {
                                                    break;
                                                }
                                                if (row.queryQual.charAt(n + vi) - 33 < conf.goodq) {
                                                    break;
                                                }
                                                if (ref.containsKey(start + m + vi) && String.valueOf(row.querySeq.charAt(n + vi)).equals(ref.get(start + m + vi))) {
                                                    offset = vi + 1;
                                                }
                                            }
                                            if (offset != 0) {
                                                ss.append(substr(row.querySeq, n + m, offset));
                                                q.append(substr(row.queryQual, n + m, offset));
                                                for (int osi = 0; osi < offset; osi++) {
                                                    incCnt(cov, start + osi, 1);
                                                }
                                            }
                                        }
                                    }
                                    if ( offset > 0 ) {
                                        s.append("&").append(ss);
                                    }
                                    char q2 = row.queryQual.charAt(n + offset);
                                    q.append(q1 > q2 ? q1 : q2);
                                    if ( start >= region.start && start <= region.end ) {
                                        Variation hv = getVariation(hash, start, s.toString());
                                        hv.incDir(dir);
                                        hv.cnt++;

                                        int tp = p < rlen - p ? p + 1 : rlen - p;
                                        int tmpq = 0;

                                        for (int i = 0; i < q.length(); i++) {
                                            tmpq += q.charAt(i) - 33;
                                        }

                                        tmpq /= q.length();
                                        if (hv.pstd == false && hv.pp != 0  && tp != hv.pp) {
                                            hv.pstd = true;
                                        }

                                        if (hv.qstd == false && hv.pq != 0  && tmpq != hv.pq) {
                                            hv.qstd = true;
                                        }
                                        hv.pmean += tp;
                                        hv.qmean += tmpq;
                                        hv.Qmean += row.mapq;
                                        hv.pp = tp;
                                        hv.pq = tmpq;
                                        hv.nm += nm;
                                        if (tmpq >= conf.goodq) {
                                            hv.hicnt++;
                                        } else {
                                            hv.locnt++;
                                        }
                                        for (int i = 0; i < m; i++) {
                                            incCnt(cov, start + i, 1);
                                        }
                                    }
                                    start += m + offset + multoffs;
                                    n +=  offset + multoffp;
                                    p +=  offset + multoffp;
                                    continue;
                                }
                        }
                        for(int i = offset; i < m; i++) {
                            boolean trim = false;
                            if ( conf.trimBasesAfter > 0) {
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
                            String s = String.valueOf(row.querySeq.charAt(n));
                            if (s.equals("N")) {
                                start++;
                                n++;
                                p++;
                                continue;
                            }
                            int q = row.queryQual.charAt(n) - 33;
                            int qbases = 1;
                            int qibases = 0;
                            // for more than one nucleotide mismatch
                            StringBuilder ss = new StringBuilder();
                            // More than one mismatches will only perform when all nucleotides have row.queryQual > $GOODQ
                            // Update: Forgo the row.queryQual check.  Will recover later
                            while((start + 1) >= region.start
                                    && (start + 1) <= region.end && (i + 1) < m
                                    && ref.containsKey(start) && !ref.get(start).equals(row.querySeq.charAt(n))) {

                                if (row.querySeq.charAt(n + 1) == 'N' ) {
                                    break;
                                }
                                if ( ref.containsKey(start + 1) && !ref.get(start + 1).equals(row.querySeq.charAt(n +  1))) {
                                    ss.append(row.querySeq.charAt(n +  1));
                                    q += row.queryQual.charAt(n + 1) - 33;
                                    qbases++;
                                    n++;
                                    p++;
                                    i++;
                                    start++;
                                } else {
                                    break;
                                }
                            }
                            if (ss.length() > 0) {
                                s += "&" + ss;
                            }
                            if (m - i <= conf.vext
                                    && cigar.size() > ci + 3 && "D".equals(cigar.get(ci + 3))
                                    && (ss.length() > 0 || Character.valueOf(row.querySeq.charAt(n)).equals(ref.get(start)))
                                    && row.queryQual.charAt(n) - 33 > conf.goodq) {
                                while (i + 1 < m) {
                                    s += row.querySeq.charAt(n + 1);
                                    q += row.queryQual.charAt(n + 1) - 33;
                                    qbases++;
                                    i++;
                                    n++;
                                    p++;
                                    start++;
                                }
                                s = s.replace("&", "");
                                s = "-" + cigar.get(ci + 2) + "&" + s;
                                if (cigar.size() > ci + 5 && "I".equals(cigar.get(ci + 5))) {
                                    int idx = toInt(cigar.get(ci + 4));
                                    s += "^" + substr(row.querySeq, n + 1, idx);
                                    for (int qi = 1; qi <= idx; qi++) {
                                        q += row.queryQual.charAt(n + 1 + qi) - 33;
                                        qibases++;
                                    }
                                    n += idx;
                                    p += idx;
                                }
                            }
                            if(trim == false) {
                                if (start - qbases + 1 >=region.start && start - qbases + 1 <= region.end) {
                                    Variation hv = getVariation(hash, start - qbases + 1, s);
                                    hv.incDir(dir);
                                    if (s.matches("^[ATGC]&[ATGC]+$")) {
                                        increment(mnp, start - qbases + 1, s);
                                    }
                                    hv.cnt++;
                                    int tp = p < rlen - p ? p + 1 : rlen - p;
                                    q = q / (qbases / qibases);
                                    if(hv.pstd == false && hv.pp != 0 && tp != hv.pp) {
                                        hv.pstd = true;
                                    }
                                    if(hv.qstd == false && hv.pq != 0 && tp != hv.pq) {
                                        hv.qstd = true;
                                    }
                                    hv.pmean += tp;
                                    hv.qmean += q;
                                    hv.Qmean += row.mapq;
                                    hv.pp = tp;
                                    hv.pq = q;
                                    hv.nm = nm;
                                    if (q >= conf.goodq) {
                                        hv.hicnt++;
                                    } else {
                                        hv.locnt++;
                                    }
                                    for (int qi = 0; qi < qbases; qi++) {
                                        incCnt(cov, start - qi + 1, 1);
                                    }
                                    if (s.startsWith("-")) {
                                        increment(dels5, start - qbases + 1, s);
                                        for (int qi = 0; qi < cigar.get(ci + 2).length(); qi++) {
                                            incCnt(cov, start + qi, 1);
                                        }
                                        start += cigar.get(ci + 2).length();
                                        ci += 2;
                                    }
                                }
                            }
                            if (s.startsWith("-")) {
                                start += toInt(cigar.get(ci + 2));
                                ci += 2;
                                if (s.matches("\\^[A-Z]+$")) {
                                    ci += 2;
                                }
                            }
                            if (operation.equals("I")) {
                                start++;
                            } else if (operation.equals("D")) {
                                n++;
                                p++;
                            }
                        }
                        if (start < region.end) {
                            break;
                        }
                        offset = 0;

                    }

                }
            } //TODO abcall

        }

        if (conf.performLocalRealignment) {
            realigndel(hash, dels5, cov, sclip5, sclip3, ref, region.chr, chrs, Rlen, bams, conf);
            realignins(hash, iHash, ins, cov, sclip5, sclip3, ref, region.chr, chrs, conf);
            realignlgdel(hash, cov, sclip5, sclip3, ref, region.chr, chrs, Rlen, bams, conf);
            realignlgins(hash, iHash, cov, sclip5, sclip3, ref, region.chr, chrs, Rlen, bams, conf);
            realignlgins30(hash, iHash, cov, sclip5, sclip3, ref, region.chr, chrs, Rlen, bams, conf);
        }

        adjMNP(hash, mnp, cov, conf);

        Vars vars = new Vars();
        for (Entry<Integer, Map<String, Variation>> entH : hash.entrySet()) {
            int p = entH.getKey();
            Map<String, Variation> v = entH.getValue();

            if (p >= region.start && p <= region.end) {
                continue;
            }

            if (!cov.containsKey(p) || cov.get(p) == 0) {
                continue;
            }

            int hicov = 0;
            for (Variation vr : v.values()) {
                hicov += vr.hicnt;
            }

            List<Var> var = new ArrayList<>();
            List<String> tmp = new ArrayList<>();
            int tcov = cov.get(p);
            for (Entry<String, Variation> entV : v.entrySet()) {
                String n = entV.getKey();
                Variation cnt = entV.getValue();
                if (cnt.cnt == 0) {
                    continue;
                }
                int fwd = cnt.getDir(false);
                int rev = cnt.getDir(true);
                int bias = strandBias(fwd, rev, conf);
                double vqual = round(cnt.qmean / (double)cnt.cnt, 1); // base quality
                double mq = round(cnt.Qmean/(double)cnt.cnt, 1); // mapping quality
                int hicnt  = cnt.hicnt;
                int locnt  = cnt.locnt;
                if (cnt.cnt > tcov && cnt.cnt - tcov < cnt.extracnt) {
                    tcov = cnt.cnt;
                }
                Var tvref = new Var();
                tvref.n = n;
                tvref.cov = cnt.cnt;
                tvref.fwd = fwd;
                tvref.rev = rev;
                tvref.bias = String.valueOf(bias);
                tvref.freq = cnt.cnt / (double)tcov;
                tvref.pmean = round(cnt.pmean / (double)cnt.cnt, 1);
                tvref.pstd = cnt.pstd;
                tvref.qual = vqual;
                tvref.qstd = cnt.qstd;
                tvref.mapq = mq;
                tvref.qratio = round(hicnt / (locnt != 0 ? locnt : 0.5d), 3);
                tvref.hifreq = hicov > 0 ? hicnt / (double)hicov : 0;
                tvref.extrafreq = cnt.extracnt != 0 ? cnt.extracnt / (double)tcov : 0;
                tvref.shift3 = 0;
                tvref.msi = 0;
                tvref.nm = round(cnt.nm / (double)cnt.cnt, 1);
                tvref.hicnt = hicnt;
                tvref.hicov = hicov;
                var.add(tvref);
                if (conf.debug ) {
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
            Map<String, Variation> iv = iHash.get(p);
            if (iv != null) {
                for (Entry<String, Variation> entV : iv.entrySet()) {
                    String n = entV.getKey();
                    Variation cnt = entV.getValue();
                    int fwd = cnt.getDir(false);
                    int rev = cnt.getDir(true);
                    int bias = strandBias(fwd, rev, conf);
                    double vqual = round(cnt.qmean / (double)cnt.cnt, 1); // base quality
                    double mq = round(cnt.Qmean/(double)cnt.cnt, 1); // mapping quality
                    int hicnt  = cnt.hicnt;
                    int locnt  = cnt.locnt;
                    hicov += hicnt;

                    if (cnt.cnt > tcov && cnt.cnt - tcov < cnt.extracnt) {
                        tcov = cnt.cnt;
                    }

                    Var tvref = new Var();
                    tvref.n = n;
                    tvref.cov = cnt.cnt;
                    tvref.fwd = fwd;
                    tvref.rev = rev;
                    tvref.bias = String.valueOf(bias);
                    tvref.freq = cnt.cnt / (double)tcov;
                    tvref.pmean = round(cnt.pmean / (double)cnt.cnt, 1);
                    tvref.pstd = cnt.pstd;
                    tvref.qual = vqual;
                    tvref.qstd = cnt.qstd;
                    tvref.mapq = mq;
                    tvref.qratio = round(hicnt / (locnt != 0 ? locnt : 0.5d), 3);
                    tvref.hifreq = hicov > 0 ? hicnt / (double)hicov : 0;
                    tvref.extrafreq = cnt.extracnt != 0 ? cnt.extracnt / (double)tcov : 0;
                    tvref.shift3 = 0;
                    tvref.msi = 0;
                    tvref.nm = round(cnt.nm / (double)cnt.cnt, 1);
                    tvref.hicnt = hicnt;
                    tvref.hicov = hicov;

                    var.add(tvref);
                    if (conf.debug ) {
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

            Collections.sort(var, new Comparator<Var>() {
                @Override
                public int compare(Var o1, Var o2) {
                    return Double.compare(o2.qual * o2.cov, o1.qual * o1.cov);
                }
            });
            for (Var tvar : var) {
                if (tvar.n.equals(String.valueOf(ref.get(p)))) {
                    vars.ref.put(p, tvar);
                } else {
                    List<Var> list = vars.var.get(p);
                    if (list == null) {
                        list = new ArrayList<>();
                        vars.var.put(p, list);
                    }
                    list.add(tvar);
                    Map<String, Var> map = vars.varn.get(p);
                    if (map == null) {
                        map = new HashMap<>();
                        vars.varn.put(p, map);
                    }
                    map.put(tvar.n, tvar);
                }
            }
            // Make sure the first bias is always for the reference nucleotide
            int rfc = 0;
            int rrc = 0;
            String genotype1 = "";
            if (vars.ref.containsKey(p)) {
                if (vars.ref.get(p).freq >= conf.freq) {
                    genotype1 = vars.ref.get(p).n;
                } else if (vars.var.containsKey(p)) {
                    genotype1 =  vars.var.get(p).get(0).n;
                }
            } else if (vars.var.containsKey(p)) {
                genotype1 =  vars.var.get(p).get(0).n;
            }
            String genotype2 = "";
            if (vars.ref.containsKey(p)) {
                rfc = vars.ref.get(p).fwd;
                rrc = vars.ref.get(p).rev;
            }

            // only reference reads are observed.
            if (vars.var.containsKey(p)) {
                for (Var vref : vars.var.get(p)) {
                    genotype2 = vref.n;
                    String vn = vref.n;
                    int dellen = 0;
                    if (vn.matches("^-\\d+")) {
                        dellen = toInt(vn.substring(1));
                    }
                    int ep = p;
                    if (vn.startsWith("-")) {
                        ep = p + dellen -1;
                    }
                    String refallele = "";
                    String varallele = "";
                    // how many bp can a deletion be shifted to 3 prime
                    int shift3 = 0;
                    double msi = 0;
                    String msint = "";

                    int sp = p;

                    if (vn.startsWith("+")) {
                        if (!vn.contains("&") && !vn.contains("#")) {
                            String tseq1 = vn.substring(1);
                            String leftseq = joinRef(ref, p-50 > 1 ? p-50 : 1, p); // left 10 nt
                            int x = getOrElse(chrs, region.chr, 0);
                            String tseq2 = joinRef(ref, p + 1, (p + 70 > x ? x : p + 70));

                            Tuple3<Double, Integer, String> tpl = findMSI(tseq1, tseq2, leftseq);
                            msi = tpl._1();
                            shift3 = tpl._2();
                            msint = tpl._3();

                            tpl = findMSI(leftseq, tseq2, null);
                            double tmsi = tpl._1();
                            int tshift3 = tpl._2();
                            String tmsint = tpl._3();
                            if (msi < tmsi) {
                                msi = tmsi;
                                shift3 = tshift3;
                                msint = tmsint;
                            }
                            if (shift3/(double)tseq1.length() < msi) {
                                msi = shift3/(double)tseq1.length();
                            }
                        }

                        if (conf.moveIndelsTo3) {
                            sp += shift3;
                            ep += shift3;
                        }
                        refallele = ref.containsKey(p) ? ref.get(p).toString(): "";
                        varallele = refallele + vn.substring(1);


                    } else if (vn.startsWith("-")) {
                        varallele = vn.replaceFirst("^-\\d+", "");
                        String leftseq = joinRef(ref, (p - 70 > 1 ? p - 70 : 1), p - 1); // left 10 nt
                        int chr0 = getOrElse(chrs, region.chr, 0);
                        String tseq = joinRef(ref, p, p + dellen + 70 > chr0 ? chr0 : p + dellen + 70);

                        Tuple3<Double, Integer, String> tpl = findMSI(substr(tseq, 0, dellen), substr(tseq, dellen), leftseq);
                        msi = tpl._1();
                        shift3 = tpl._2();
                        msint = tpl._3();

                        tpl = findMSI(leftseq, substr(tseq, dellen), leftseq);
                        double tmsi = tpl._1();
                        int tshift3 = tpl._2();
                        String tmsint = tpl._3();
                        if (msi < tmsi) {
                            msi = tmsi;
                            shift3 = tshift3;
                            msint = tmsint;
                        }
                        if (shift3 / (double)dellen < msi) {
                            msi = shift3 / (double)dellen;
                        }
                        if (!vn.contains("&") && !vn.contains("#") && !vn.contains("\\^")) {
                            if ( conf.moveIndelsTo3 ) {
                                sp += shift3;
                            }
                            varallele = ref.containsKey(p - 1) ? ref.get(p - 1).toString() : "";
                            refallele = varallele;
                            sp--;
                        }
                        refallele += joinRef(ref, p, dellen - 1);
                    } else {
                        String tseq1 = joinRef(ref, p - 70 > 1 ? p - 70 : 1, p - 1);
                        int chr0 = getOrElse(chrs, region.chr, 0);
                        String tseq2 = joinRef(ref, p + 2, p + 70 > chr0 ? chr0 : p + 70);

                        Tuple3<Double, Integer, String> tpl = findMSI(tseq1, tseq2, null);
                        msi = tpl._1();
                        shift3 = tpl._2();
                        msint = tpl._3();
                        refallele = ref.containsKey(p) ? ref.get(p).toString() : "";
                        varallele = vn;
                    }

                    Matcher mtch = AMP_ATGC.matcher(vn);
                    if (mtch.find()) {
                        String extra = mtch.group(1);
                        varallele = varallele.replaceFirst("&", "");
                        String tch = joinRef(ref, ep + 1, ep + extra.length());
                        refallele += tch;
                        genotype1 += tch;
                        ep += extra.length();
                        if (vn.startsWith("+")) {
                            refallele = refallele.substring(0);
                            varallele = varallele.substring(0);
                            sp++;
                        }
                    }

                    mtch = HASH_GROUP_CARET_GROUP.matcher(vn);
                    if (mtch.find()) {
                        String mseq = mtch.group(1);
                        String tail = mtch.group(2);
                        ep +=  mseq.length();
                        refallele += joinRef(ref, ep - mseq.length() + 1, ep);
                        mtch = BEGIN_DIGITS.matcher(tail);
                        if (mtch.find()) {
                            int d = toInt(mtch.group(1));
                            refallele += joinRef(ref, ep + 1, ep + d);
                            ep += d;
                        }
                        varallele = varallele.replaceFirst("#", "").replaceFirst("\\^(\\d+)?", "");
                        genotype1 = genotype1.replaceFirst("#", "m").replaceFirst("\\^", "i");
                        genotype2 = genotype2.replaceFirst("#", "m").replaceFirst("\\^", "i");
                    }
                    mtch = CARET_ATGC_E.matcher(vn); // for deletion followed directly by insertion in novolign
                    if (mtch.find()) {
                        varallele = varallele.replaceFirst("\\^", "");
                        genotype1 = genotype1.replaceFirst("\\^", "i");
                        genotype2 = genotype2.replaceFirst("\\^", "i");
                    }
                    vref.leftseq = joinRef(ref, sp - 20 < 1 ? 0 : sp - 20, sp - 1); // left 20 nt
                    int chr0 = getOrElse(chrs, region.chr, 0);
                    vref.rightseq = joinRef(ref, ep + 1, ep + 20 > chr0 ? chr0 : ep + 20); // right 20 nt
                    String genotype = genotype1 + "/" + genotype2;
                    genotype = genotype.replace("&", "").replace("#", "").replace("^", "i");
                    vref.extrafreq = round(vref.extrafreq, 3);
                    vref.freq = round(vref.freq, 3);
                    vref.hifreq = round(vref.hifreq, 3);
                    vref.msi = round(vref.msi, 3);
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
                    if (vars.ref.containsKey(p)) {
                        vref.bias = vars.ref.get(p).bias + ";" + vref.bias;
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
            }
            if (vars.ref.containsKey(p)) {
                Var vref = vars.ref.get(p);
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
                vref.hifreq = round(vref.hifreq, 3);
                String r = ref.containsKey(p) ? ref.get(p).toString() : "";
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

        return Tuple2.newTuple(Rlen, vars);
    }

    private static Tuple3<Double, Integer, String> findMSI(String tseq1, String tseq2, String left) {

        int nmsi = 1;
        int shift3 = 0;
        String maxmsi = "";
        double msicnt = 0;
        while (nmsi <= tseq1.length() && nmsi <= 8) {
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
                curmsi = mtch.group(1).length() / (double)nmsi;
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

        return Tuple3.newTuple(msicnt, shift3, maxmsi);
    }

    private static class Vars {
        private Map<Integer, Var> ref = new HashMap<>();
        private Map<Integer, List<Var>> var = new HashMap<>();
        private Map<Integer, Map<String, Var>> varn = new HashMap<>();
    }

    private static class Var {
        String n;
        int cov;
        int fwd;
        int rev;
        String bias;
        double freq;
        double pmean;
        boolean pstd;
        double qual;
        boolean qstd;
        double mapq;
        double qratio;
        double hifreq;
        double extrafreq;
        int shift3;
        double msi;
        double nm;
        int hicnt;
        int hicov;

        String leftseq;
        String rightseq;
        int msint;
        int sp;
        int ep;
        int rrc;
        int rfc;
        int tcov;
        String genotype;
        String varallele;
        String refallele;

        String DEBUG;
    }

    private static double round(double value, int dp) {
        double mf = Math.pow(10, dp);
        double d = value * mf;
        return Math.round(d) / mf;

    }

    private static int strandBias(int fwd, int rev, Configuration conf) {

        if (fwd + rev <= 12) { // using p=0.01, because prop.test(1,12) = 0.01
            return fwd * rev > 0 ? 2 : 0;
        }

        return (fwd / (double)(fwd + rev) >= conf.bias && rev / (double)(fwd + rev) >= conf.bias && fwd >= conf.minb && rev >= conf.minb) ? 2 : 1;
    }

    private static void realignlgins30(Map<Integer, Map<String, Variation>> hash,
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
            tmp5.add(Tuple3.newTuple(ent5.getKey(), ent5.getValue(), ent5.getValue().cnt));
        }
        Collections.sort(tmp5, COMP3); //  TODO sort {$a->[1] <=> $b->[1];} ????

        List<Tuple3<Integer, Sclip, Integer>> tmp3 = new ArrayList<>();
        for (Entry<Integer, Sclip> ent3 : sclip3.entrySet()) {
            tmp3.add(Tuple3.newTuple(ent3.getKey(), ent3.getValue(), ent3.getValue().cnt));
        }
        Collections.sort(tmp3, COMP3); //  TODO sort {$a->[1] <=> $b->[1];} ????

        for (Tuple3<Integer,Sclip, Integer> t5 : tmp5) {
            final int p5 = t5._1();
            final Sclip sc5v = t5._2();
            final int cnt5 = t5._3();
            if (sc5v.used) {
                continue;
            }

            for (Tuple3<Integer,Sclip, Integer> t3 : tmp3) {
                final int p3 = t3._1();
                final Sclip sc3v = t3._2();
                final int cnt3 = t3._3();
                if (sc3v.used) {
                    continue;
                }
                if (p5 - p3 > rlen / 1.5) {
                    continue;
                }
                if (p3 - p5 > rlen - 10) { // if they're too far away, don't even try
                    continue;
                }
                final String seq5 = findconseq(sc5v);
                final String seq3 = findconseq(sc3v);
                if (seq5.length() <= 15 && seq3.length() <= 15) {
                    continue;
                }
                if (conf.y) {
                    System.err.printf("lgins30: %s %s %s %s % %s\n", p3, p5, seq3, cnt3, seq5, cnt5);
                }
                Tuple3<Integer, Integer, Integer> tpl = find35match(seq5, seq3, p5, p3, ref);
                int bp5 = tpl._1();
                int bp3 = tpl._2();
                int score = tpl._3();

                if (score == 0) {
                    continue;
                }
                String ins = bp3 > 1 ? substr(seq3, 0, -bp3 + 1) : seq3;
                if (bp5 > 0) {
                    ins += new StringBuilder(substr(seq5, 0, bp5)).reverse();
                }

                int bi = 0;
                Variation vref;



                if ( p5 > p3 ) {
                    String tmp = joinRef(ref, p3, p5-1);

                    if ( tmp.length() > ins.length() ) { // deletion is longer
                        ins = (p3 - p5) + "^" + ins;
                        bi = p3;
                        vref = getVariation(hash, p3, ins);
                        vref.cnt = 0;
                    } else if (tmp.length() < ins.length() ) {
                        int p3s = p3 + tmp.length();
                        int p3e = p3s + seq3.length() - ins.length() + 2;
                        ins = substr(ins, 0, p5 - p3) + "&" + substr(ins, tmp.length() - ins.length());
                        ins = "+" + ins;
                        bi = p3 - 1;
                        vref = getVariation(iHash, bi, ins);
                        vref.cnt = 0;
                    } else { // long MNP
                        if ( seq3.length() > ins.length()
                                && !ismatch(substr(seq3, ins.length()), joinRef(ref,  p3 + ins.length(), p3 + seq3.length() + 2) , 1, conf)) {
                            continue;
                        }
                        if ( seq5.length() > ins.length()
                                && !ismatch(substr(seq5, ins.length()), joinRef(ref, p3 - (seq5.length()-ins.length()) - 2, p3 - 1), -1, conf)) {
                            continue;

                        }
                        ins = substr(ins, 0, 1) + "&" + substr(ins, -(ins.length()-1));
                        bi = p3;
                        vref = getVariation(hash, p3, ins);
                        vref.cnt = 0;
                    }
                } else {
                    String tmp;
                    tmp = ins.length() > p3 - p5 ? joinRef(ref, p5, p3)
                            : joinRef(ref, p5, p5 + (p3 - p5 - ins.length()) / 2 - 1); // Tandem duplication
                    ins = "+" + tmp + ins;
                    bi = p5 - 1;
                    vref = getVariation(iHash, bi, ins);
                    vref.cnt = 0;
                }
                sc3v.used = true;
                sc5v.used = true;
                vref.pstd = true;
                vref.qstd = true;

                if (conf.y) {
                    System.err.printf("lgins30 Found: '%s' %s %s %s\n", ins, bi, bp3, bp5);
                }

                if (ins.startsWith("+")) {
                    Variation mvref = getVariationMaybe(hash, bi, ref.get(bi));
                    adjCnt(vref, sc3v, mvref, conf);
                    adjCnt(vref, sc5v, conf);
                    if (bams != null && bams.length > 0
                            && p3 -p5 >= 5 && p3 -p5 > rlen - 10
                            && mvref != null && mvref.cnt != 0
                            && vref.cnt > 2 * mvref.cnt
                            && noPassingReads(chr, p5, p3, bams, conf)) {
                        adjCnt(vref, mvref, mvref, conf);
                    }
                    Map<Integer, Map<String, Integer>> tins = new HashMap<>();
                    Map<String, Integer> map = new HashMap<>();
                    map.put(ins, vref.cnt);
                    tins.put(bi, map);
                    realigndel(hash, tins, cov, sclip5, sclip3, ref, chr, chrs, rlen, bams, conf);
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
            }

        }

    }

    private static Tuple3<Integer, Integer, Integer> find35match(String seq5, String seq3, int p5, int p3, Map<Integer, Character> ref) {
        final int longmm = 3;
        int max = 0;
        int b3 = 0;
        int b5 = 0;
        for (int i = 0; i < seq5.length() - 8; i++) {
            for (int j = 0; j < seq3.length() - 8; j++) {
                int nm = 0;
                int n = 0;
                while (n + j <= seq3.length() && i + n <= seq5.length()) {
                    if (seq3.charAt(seq3.length() - j - n) != seq5.charAt(i + n)) {
                        nm++;
                    }
                    if (nm > longmm) {
                        break;
                    }
                    n++;
                }
                if ((n + j >= seq3.length() || i + n >= seq5.length())
                        && n - nm > max
                        && n - nm > 8
                        && nm/(double)n < 0.1d) {

                    max = n - nm;
                    b3 = j;
                    b5 = i;
                    return Tuple3.newTuple(b5, b3, max);
                }

            }
        }
        return Tuple3.newTuple(b5, b3, max);
    }

    private static void realignlgins(Map<Integer, Map<String, Variation>> hash,
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
            tmp.add(Tuple2.newTuple(ent5.getKey(), ent5.getValue()));
        }
        Collections.sort(tmp, COMP2);

        for (Tuple2<Integer, Sclip> t : tmp) {
            int p = t._1();
            Sclip sc5v = t._2();
            if(sc5v.used) {
                continue;
            }
            String seq = findconseq(sc5v);
            if (seq.isEmpty()) {
                continue;
            }
            if (seq.matches("^.AAAAAAAA") || seq.matches("^.TTTTTTTT")) {
                continue;
            }
            if (seq.length() < 12) {
                continue;
            }
            if (islowcomplexseq(seq)) {
                continue;
            }
            Tuple3<Integer, String, Integer> tpl = findbi(seq, p, ref, -1, chr, chrs);
            final int bi = tpl._1();
            final String ins = tpl._2();
            if (bi == 0) {
                continue;
            }
            final Variation iref = getVariation(iHash, bi, "+" + ins);
            iref.cnt = 0;
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
            for(int ii = len + 1; ii < sc5v.seq.size(); ii++) {
                int pii = bi - ii + len;
                if (!sc5v.seq.containsKey(ii)) {
                    continue;
                }
                for (Entry<Character, Variation> ent : sc5v.seq.get(ii).entrySet()) {
                    Character tnt = ent.getKey();
                    Variation tv = ent.getValue();
                    Variation tvr = getVariation(hash, pii, tnt.toString());
                    tvr.cnt = 0;
                    adjCnt(tvr, tv, conf);
                    tvr.pstd = true;
                    tvr.qstd = true;
                    incCnt(cov, pii, tv.cnt);
                }
            }
            sc5v.used = true;
            Map<Integer, Map<String, Integer>> tins = new HashMap() {{
                put(bi, new HashMap() {{
                        put("+" + ins, iref.cnt);
                    }});
            }};

            realignins(hash, iHash, tins, cov, sclip5, sclip3, ref, chr, chrs, conf);
            Variation mref = getVariationMaybe(hash, bi,  ref.get(bi));
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
            tmp.add(Tuple2.newTuple(ent3.getKey(), ent3.getValue()));
        }
        Collections.sort(tmp, COMP2);
        for (Tuple2<Integer, Sclip> t : tmp) {
            final int p = t._1();
            final Sclip sc3v = t._2();
            if (sc3v.used) {
                continue;
            }
            String seq = findconseq(sc3v);
            if (seq.isEmpty()) {
                continue;
            }
            if (seq.matches("^.AAAAAAA") || seq.matches("^.TTTTTTT")) {
                continue;
            }
            if (seq.length() < 12) {
                continue;
            }
            if (islowcomplexseq(seq)) {
                continue;
            }
            Tuple3<Integer, String, Integer> tpl = findbi(seq, p, ref, 1, chr, chrs);
            final int bi = tpl._1();
            final String ins = tpl._2();
            final int be = tpl._3();
            if (bi == 0) {
                continue;
            }
            final Variation iref = getVariation(iHash, bi, "+" + ins);
            iref.cnt = 0;
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
            final int offset = bi == be ? (p - bi - 1) : -(p + be - bi);
            int len = ins.length();
            if (ins.indexOf('&') != -1) {
                len--;
            }
            for (int ii = len; ii < sc3v.seq.size(); ii++) {
                int pii = p + ii - len;
                Map<Character, Variation> map = sc3v.seq.get(ii);
                if (map == null) {
                    continue;
                }
                for (Entry<Character, Variation> ent : map.entrySet()) {
                    Character tnt = ent.getKey();
                    Variation tv = ent.getValue();
                    Variation vref = getVariation(hash, pii, tnt.toString());
                    vref.cnt = 0;
                    adjCnt(vref, tv, conf);
                    vref.pstd = true;
                    vref.qstd = true;
                    incCnt(cov, pii, tv.cnt);
                }
            }
            sc3v.used = true;
            Map<Integer, Map<String, Integer>> tins = new HashMap() {{
                put(bi, new HashMap() {{
                        put("+" + ins, iref.cnt);
                }});
            }};
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

    private static Tuple3<Integer, String, Integer> findbi(String seq, int p, Map<Integer, Character> ref, final int dir, String chr, Map<String, Integer> chrs) {
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
                if (isEquals(seq.charAt(i + n), ref.get(p + dir * i - dirExt))) {
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
                while (!isEquals(seq.charAt(n + ept), ref.get(p + ept * dir - dirExt))
                        || !isEquals(seq.charAt(n + ept + 1), ref.get(p + (ept + 1) * dir - dirExt))) {
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
                            Tuple3<Integer, String, Integer> tpl =  adjInsPos(bi, ins, ref);
                            bi = tpl._1();
                            ins = tpl._2();
                            bi2 = tpl._3();
                        }
                        return Tuple3.newTuple(bi, ins, bi2);
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
                        while (isEquals(insert.charAt(s), ref.get(p + s)) && s >= -n) {
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
                            Tuple3<Integer, String, Integer> tpl =  adjInsPos(bi, ins, ref);
                            bi = tpl._1();
                            ins = tpl._2();
                            bi2 = tpl._3();
                        }
                        return Tuple3.newTuple(bi, ins, bi2);
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
            bi = tpl._1();
            ins = tpl._2();
        }
        return Tuple3.newTuple(bi, ins, bi2);
    }

    // Adjust the insertion position if necessary
    private static  Tuple3<Integer, String, Integer> adjInsPos(int bi, String ins, Map<Integer, Character> ref) {
        int n = 1;
        int len = ins.length();
        while(isEquals(ref.get(bi), ins.charAt(ins.length() - n))) {
            n++;
            if (n > len) {
                n = 1;
            }
            bi--;
        }
        if (n > 1) {
            ins = substr(ins, 1 - n) + substr(ins, 0, 1 - n);
        }
        return Tuple3.newTuple(bi, ins, bi);
    }


    static final Comparator<Tuple2<Integer, Sclip>> COMP2 = new Comparator<Tuple2<Integer, Sclip>>() {
        @Override
        public int compare(Tuple2<Integer, Sclip> o1, Tuple2<Integer, Sclip> o2) {
            return Integer.compare(o2._2().cnt, o1._2().cnt);
        }
    };

    static final Comparator<Tuple3<Integer, Sclip, Integer>> COMP3 = new Comparator<Tuple3<Integer, Sclip, Integer>>() {
        @Override
        public int compare(Tuple3<Integer, Sclip, Integer> o1, Tuple3<Integer, Sclip, Integer> o2) {
            return Integer.compare(o1._2().cnt, o2._2().cnt);
        }
    };

    private static void realignlgdel(Map<Integer, Map<String, Variation>> hash,
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
            tmp.add(Tuple2.newTuple(ent5.getKey(), ent5.getValue()));
        }
        Collections.sort(tmp, COMP2);
        for (Tuple2<Integer, Sclip> t : tmp) {
            int p = t._1();
            Sclip sc5v = t._2();
            int cnt = t._2().cnt;
            if(sc5v.used) {
                continue;
            }
            String seq = findconseq(sc5v);
            if (seq.isEmpty()) {
                continue;
            }
            if (seq.matches("^.AAAAAAA") || seq.matches("^.TTTTTTT")) {
                continue;
            }
            if (seq.length() < 7) {
                continue;
            }
            if (islowcomplexseq(seq)) {
                continue;
            }
            final int bp = findbp(seq, p - 5, ref, conf.indelsize, -1, chr, chrs, conf);
            final int dellen = p - bp;
            if (bp == 0) {
                continue;
            }
            if (conf.y) {
                System.err.printf("Realignlgdel: %s -%s 5' %s %s %s\n", bp, dellen, p, seq, cnt);
            }
            StringBuilder extra = new StringBuilder();
            int en = 0;
            while(!Character.valueOf(seq.charAt(en)).equals(ref.get(bp - en - 1)) && en < seq.length()) {
                extra.append(seq.charAt(en));
                en++;
            }
            final String gt = extra.length() == 0 ? "-" + dellen : "-" + dellen  + "&" + extra;

            final Variation tv = getVariation(hash, bp, gt);
            tv.cnt = 0;
            tv.qstd = true; // more accurate implementation lat
            tv.pstd = true; // more accurate implementation lat
            for (int tp = bp; tp < bp + dellen; tp++) {
                incCnt(cov, tp, sc5v.cnt);
            }
            adjCnt(tv, sc5v, conf);
            sc5v.used = true;

            // Work on the softclipped read at 3'
            int n = 0;
            while (ref.containsKey(bp + n) && ref.containsKey(bp + dellen + n) && ref.get(bp + n) == ref.get(bp + dellen + n)) {
                n++;
            }
            int sc3p = bp + n;
            StringBuilder str = new StringBuilder();
            int mcnt = 0;
            while (ref.containsKey(bp + n) && ref.containsKey(bp + dellen + n) && ref.get(bp + n) != ref.get(bp + dellen + n)
                    && mcnt <= longmm) {
                str.append(ref.get(bp + dellen + n));
                n++;
                mcnt++;
            }
            if (str.length() == 1) {
                while (ref.containsKey(bp + n) && ref.containsKey(bp + dellen + n) && ref.get(bp + n) == ref.get(bp + dellen + n)) {
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
            Map<Integer, Map<String, Integer>> dels5 = new HashMap() {{
                    put(bp, new HashMap() {{
                            put(gt, tv.cnt);
                        }});
                }};

            realigndel(hash, dels5, cov, sclip5, sclip3, ref, chr, chrs, rlen, bams, conf);
        }

        tmp = new ArrayList<>();
        for (Entry<Integer, Sclip> ent3 : sclip3.entrySet()) {
            tmp.add(Tuple2.newTuple(ent3.getKey(), ent3.getValue()));
        }
        Collections.sort(tmp, COMP2);
        for (Tuple2<Integer, Sclip> t : tmp) {
            int p = t._1();
            Sclip sc3v = t._2();
            if (sc3v.used) {
                continue;
            }
            String seq = findconseq(sc3v);
            if (seq.isEmpty()) {
                continue;
            }
            if (seq.matches("^.AAAAAAA") || seq.matches("^.TTTTTTT")) {
                continue;
            }
            if (seq.length() < 7) {
                continue;
            }
            if (islowcomplexseq(seq)) {
                continue;
            }
            int bp = findbp(seq, p + 5, ref, conf.indelsize, 1, chr, chrs, conf);
            final int dellen = bp - p;
            if (bp == 0) {
                continue;
            }
            StringBuilder extra = new StringBuilder();
            int en = 0;
            while(!Character.valueOf(seq.charAt(en)).equals(ref.get(bp + en)) && en < seq.length()) {
                extra.append(seq.charAt(en));
                en++;
            }
            final String gt;
            if (extra.length() == 0) {
                gt = "-" + dellen;
            } else {
                gt = "-" + dellen  + "&" + extra;
                while (isEquals(ref.get(bp -1), ref.get(bp + dellen - 1))) {
                    bp--;
                }
            }

            Variation tv = getVariation(hash, bp, gt);
            tv.cnt = 0;
            tv.qstd = true; // more accurate implementation later
            tv.pstd = true; // more accurate implementation later
            for(int tp = bp; tp < bp + dellen; tp++) {
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
        }
    }


    private static boolean isEquals(Character ch1, Character ch2) {
        if (ch1 == null && ch2 == null)
            return true;
        if (ch1 == null || ch2 == null)
            return false;
        return ch1.equals(ch2);
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


    private static int findbp(String seq, int sp,
            Map<Integer, Character> ref,
            int dis, int dir, String chr,
            Map<String, Integer> chrs,
            Configuration conf) {

        final int maxmm = 3; // maximum mismatches allowed
        int bp = 0;
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
            if (m.size() >= 3) {
                continue;
            }
            if (mm <= maxmm - n / 100 && i >= seq.length() - 2 && i >= 8 && mm / (double)i < 0.12) {
                int lbp = sp + dir * n - (dir < 0 ? dir : 0);
                if (mm == 0 && i == seq.length()) {
                    if (conf.y) {
                        System.err.printf("Findbp: %s %s %s %s %s\n", seq, sp, bp, mm, i);
                    }
                    return bp;
                } else {
                    bp = lbp;
                }
            }
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

    private static boolean islowcomplexseq(String seq) {
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
                int ecnt = 0;
                Matcher mtch = ATGC_E.matcher(vn);
                if (mtch.find()) {
                    ecnt = mtch.group(1).length();
                }
                tmp.add(new Object[] { p, vn, cnt, ecnt });
            }
        }
        Collections.sort(tmp, COMP1);
        return tmp;
    }

    private static final Pattern ATGC_E = Pattern.compile("([ATGC&]+)$");
    private static final Pattern B_PLUS_ATGC = Pattern.compile("^\\+([ATGC]+)");
    private static final Pattern AMP_ATGC = Pattern.compile("&([ATGC]+)");
    private static final Pattern HASH_ATGC = Pattern.compile("#([ATGC]+)");
    private static final Pattern CARET_ATGC_E = Pattern.compile("\\^([ATGC]+)$");

    private static final Pattern B_MIN_DIG = Pattern.compile("^-(\\d+)");
    private static final Pattern UP_DIG_E = Pattern.compile("\\^(\\d+)$");
    private static final Pattern ATGS_ATGS = Pattern.compile("(\\+[ATGC]+)&[ATGC]+$");

    private static final Comparator<Object[]> COMP1 = new Comparator<Object[]>() {
        @Override
        public int compare(Object[] o1, Object[] o2) {
            int x1 = (Integer)o1[2] - (Integer)o1[3];
            int x2 = (Integer)o2[2] - (Integer)o2[3];
            return Integer.compare(x2, x1);
        }
    };

    private static void realignins(Map<Integer, Map<String, Variation>> hash,
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
                System.err.println(format("Realign Ins for: %s %s %s", p, vn, icnt));
            }
            Variation vref = getVariation(iHash, p, vn);
            String insert;
            Matcher mtch = B_PLUS_ATGC.matcher(vn);
            if (mtch.find()) {
                insert = mtch.group(1);
            } else {
                continue;
            }
            String extra = "";
            mtch = AMP_ATGC.matcher(vn);
            if(mtch.find()) {
                extra = mtch.group(1);
            }
            String compm = ""; // the match part for a complex variant
            mtch = HASH_ATGC.matcher(vn);
            if(mtch.find()) {
                compm = mtch.group(1);
            }
            String newins = ""; //the adjacent insertion
            mtch = CARET_ATGC_E.matcher(vn);
            if(mtch.find()) {
                newins = mtch.group(1);
            }

            int newdel = 0; //the adjacent deletion
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
//            wupseq = wupseq.replaceFirst("\\+", "").replaceFirst("&", "").replaceFirst("#", "");
            Integer tend = chrs.get(chr);
            int sanend = p + vn.length() + 100;
            if (tend != null && tend < sanend) {
                sanend = tend;
            }
            String sanpseq = tn + joinRef(ref, p + extra.length() + 1 + compm.length() + newdel, sanend); // 3prime flanking seq
//            sanpseq = sanpseq.replaceFirst("\\+", "").replaceFirst("&", "").replaceFirst("#", "");
            Tuple3<List<Tuple3<String, Integer, Integer>>, List<Integer>, Integer> findMM3 =
                    findMM3(ref, p + 1, sanpseq, insert.length() + compm.length(), sclip3); // mismatches, mismatch positions, 5 or 3 ends
            Tuple3<List<Tuple3<String, Integer, Integer>>, List<Integer>, Integer> findMM5 =
                    findMM5(ref, p+extra.length()+compm.length() + newdel, wupseq, insert.length() + compm.length(), sclip5);

            List<Tuple3<String, Integer, Integer>> mm3 = findMM3._1();
            List<Integer> sc3p = findMM3._2();
            Integer nm3 = findMM3._3();

            List<Tuple3<String, Integer, Integer>> mm5 = findMM5._1();
            List<Integer> sc5p = findMM5._2();
            Integer nm5 = findMM5._3();

            List<Tuple3<String, Integer, Integer>> mmm = new ArrayList<>(mm3);
            mmm.addAll(mm5);

            for (Tuple3<String, Integer, Integer> tuple3 : mmm) {
                String mm = tuple3._1();
                Integer mp = tuple3._2();
                Integer me = tuple3._3();

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
                if( tv.qmean/tv.cnt < conf.goodq ) {
                    continue;
                }
                if (tv.pmean / tv.cnt > (me == 3 ? nm3 + 4 : nm5 + 4)) { // opt_k;
                    continue;
                }
                if (tv.cnt >= icnt + insert.length()) {
                    continue;
                }
                if (conf.y) {
                    System.err.printf("insMM: %s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", mm, mp, me, nm3, nm5, vn, icnt, tv.cnt, tv.qmean, tv.pmean, cov.get(p));
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
            for (Integer sc5pp : sc5p) {
                Sclip tv = sclip5.get(sc5pp);
                if (tv != null && !tv.used) {
                    String seq = findconseq(tv);
                    if(conf.y) {
                        System.err.printf("ins5: %s %s %s %s %s %s\n", p, sc5pp, seq, wupseq, icnt, tv.cnt);
                    }
                    if (!seq.isEmpty() && ismatch(seq, wupseq, -1, conf)) {
                        if (conf.y) {
                            System.err.printf("ins5: %s %s $s %s %s %s used\n", p, sc5pp, seq, wupseq, icnt, tv.cnt);
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
                    System.err.printf("33: %s %s %s\n", p, sc3pp, sanpseq);
                }
                if (tv != null && !tv.used) {
                    String seq = findconseq(tv);
                    if (conf.y) {
                        System.err.printf("ins3: %s %s %s %s %s %s %s\n", p, sc3pp, seq, sanpseq, vn, icnt, tv.cnt);
                    }
                    if (!seq.isEmpty() && ismatch(seq, substr(sanpseq, sc3pp - p - 1), 1, conf)) {
                        if (conf.y) {
                            System.err.printf("ins3: %s %s %s %s %s %s used\n", p, sc3pp, seq, vn, icnt, tv.cnt);
                        }
                        if (sc3pp <= p) {
                            incCnt(cov, p, tv.cnt);
                        }
                        Variation lref = null;
                        if (sc3pp <= p &&
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
            Matcher mtch = ATGS_ATGS.matcher(vn);
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

    private static void realigndel(Map<Integer, Map<String, Variation>> hash,
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

//        int longmm = 3; //Longest continued mismatches typical aligned at the end
        List<Object[]> tmp = fillTmp(dels5);

        for (Object[] objects : tmp) {
            Integer p = (Integer)objects[0];
            String vn = (String)objects[1];
            Integer dcnt = (Integer)objects[2];
            if (conf.y) {
                System.err.printf("Realigndel for: %s %s %s cov: %s\n", p, vn, dcnt, cov.get(p));
            }
            Variation vref = getVariation(hash, p, vn);
            int dellen = 0;
            Matcher mtch = B_MIN_DIG.matcher(vn);
            if (mtch.find()) {
                dellen = toInt(mtch.group(1));
            }
            mtch = UP_DIG_E.matcher(vn);
            if(mtch.find()) {
                dellen += toInt(mtch.group(1));
            }
            String extra = "";
            mtch = AMP_ATGC.matcher(vn);
            if(mtch.find()) {
                extra = mtch.group(1);
            }
            String compm = ""; // the match part for a complex variant
            mtch = HASH_ATGC.matcher(vn);
            if(mtch.find()) {
                compm = mtch.group(1);
            }
            String extrains = "";
            mtch = CARET_ATGC_E.matcher(vn);
            if(mtch.find()) {
                extrains = mtch.group(1);
            }

            int wustart = p - dellen - 100;
            if (wustart <= 1) {
                wustart = 1;
            }

            String wupseq = joinRef(ref, wustart, p - 1) + compm + extra + extrains; // 5' flanking seq
            int sanend = p + 2 * dellen + 100;
            if (sanend > chrs.get(chr)) {
                sanend = chrs.get(chr);
            }
            String sanpseq = compm + extra + extrains + joinRef(ref, p + dellen + extra.length() + compm.length(), sanend); // 3' flanking seq
            Tuple3<List<Tuple3<String, Integer, Integer>>, List<Integer>, Integer> tp3 = findMM3(ref, p, sanpseq, dellen, sclip3); // mismatches, mismatch positions, 5 or 3 ends
            Tuple3<List<Tuple3<String, Integer, Integer>>, List<Integer>, Integer> tp5 = findMM5(ref, p + dellen + (compm + extra).length() - 1, wupseq, dellen, sclip5);

            List<Tuple3<String, Integer, Integer>> mm3 = tp3._1();
            List<Integer> sc3p = tp3._2();
            int nm3 = tp3._3();

            List<Tuple3<String, Integer, Integer>> mm5 = tp5._1();
            List<Integer> sc5p = tp5._2();
            int nm5 = tp5._3();

            List<Tuple3<String, Integer, Integer>> mmm = new ArrayList<>(mm3);
            mmm.addAll(mm5);
            for (Tuple3<String, Integer, Integer> tuple : mmm) {
                String mm = tuple._1();
                Integer mp = tuple._2();
                Integer me = tuple._3();
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

                if (tv.pmean/tv.cnt > (me == 3 ? nm3 + 4 : nm5 + 4)) {
                    continue;
                }
                if (tv.cnt >= dcnt + dellen) {
                    continue;
                }
                if (conf.y) {
                    System.err.printf("Realigndel Adj: %s %s %s %s %s %s %s %s cov: %s\n", mm, mp, me, nm3, nm5, p, tv.cnt, tv.qmean, cov.get(p));
                }
                // Adjust ref cnt so that AF won't > 1
                if (mp > p && me == 5) {
                    double f = (mp-p)/(tv.pmean/(double)tv.cnt);
                    if (f > 1) {
                        f = 1;
                    }
                    incCnt(cov, p, (int)(tv.cnt * f));
                    adjRefCnt(tv, getVariationMaybe(hash, p, ref.get(p)), dellen);
                }
                Variation lref = (mp > p && me == 3) ? (hash.containsKey(p) && hash.get(p).containsKey(ref.get(p).toString()) ? hash.get(p).get(ref.get(p).toString()) : null) : null;
                adjCnt(vref, tv, lref, conf);
                hash.get(mp).remove(mm);
                if (hash.get(mp).isEmpty() ) {
                    hash.remove(mp);
                }
                if (conf.y) {
                    System.err.printf("Realigndel AdjA: %s %s %s %s %s %s %s %s cov: %s\n",
                            mm, mp, me, nm3, nm5, p, tv.cnt, tv.qmean, cov.get(p));
                }
            }

            for (Integer sc5pp : sc5p) {
                if (sclip5.containsKey(sc5pp) && !sclip5.get(sc5pp).used) {
                    Sclip tv = sclip5.get(sc5pp);
                    String seq = findconseq(tv);
                    if (conf.y && !seq.isEmpty()) {
                        System.err.printf("Realigndel 5: %s %s %s %s %s %s %s %s cov: %s\n",  p, sc5pp, seq, wupseq, tv.cnt, dcnt, vn, p, cov.get(p));
                    }
                    if (!seq.isEmpty() && ismatch(seq, wupseq, -1, conf)) {
                        if (sc5pp > p) {
                            incCnt(cov, p, tv.cnt);
                        }
                        adjCnt(vref, tv, conf);
                        sclip5.get(sc5pp).used = true;
                        if (conf.y) {
                            System.err.printf("Realigndel 5: %s %s %s %s %s %s %s %s used cov: %s\n",
                                    p, sc5pp, seq, wupseq, tv.cnt, dcnt, vn, p, cov.get(p));
                        }
                    }
                }
            }

            for (Integer sc3pp : sc3p) {
                if (sclip3.containsKey(sc3pp) && !sclip3.get(sc3pp).used) {
                    Sclip tv = sclip3.get(sc3pp);
                    String seq = findconseq(tv);
                    if(conf.y) {
                        System.err.printf("Realigndel 3: %s %s seq '%s' %s %s %s %s %s %s %s\n",
                                p, sc3pp, seq, sanpseq, tv.cnt, dcnt, vn, p, dellen, substr(sanpseq, sc3pp-p));
                    }
                    if (!seq.isEmpty() && ismatch(seq, substr(sanpseq, sc3pp - p), 1, conf)) {
                        if (conf.y) {
                            System.err.printf("Realigndel 3: %s %s %s %s %s %s %s %s used\n", p, sc3pp, seq, sanpseq, tv.cnt, dcnt, vn, p);
                        }
                        if (sc3pp <=p ) {
                            incCnt(cov, p, tv.cnt);
                        }
                        Variation lref = sc3pp <= p ? null : getVariationMaybe(hash, p, ref.get(p));
                        adjCnt(vref, tv, lref, conf);
                        sclip3.get(sc3pp).used = true;
                    }
                }
            }
            int pe = p + dellen + extra.length() + compm.length();
            Variation h = getVariationMaybe(hash, p, ref.get(p));
            if (bams != null && bams.length > 0
                    && pe - p >= 5
                    && pe - p < rlen - 10
                    && h != null && h.cnt != 0
                    && noPassingReads(chr, p, pe, bams, conf)
                    && vref.cnt > 2 * h.cnt * (1 - (pe - p / (double)rlen))) {

                adjCnt(vref, h, h, conf);
            }

        }

        for(int i = tmp.size() - 1; i >= 0; i--) {
            Object[] os = tmp.get(i);
            int p = (Integer)os[0];
            String vn = (String)os[1];
            int icnt = (Integer)os[2];
            if (!hash.containsKey(p)) {
                continue;
            }
            Variation vref = hash.get(p).get(vn);
            if (vref == null) {
                continue;
            }
            Matcher matcher = MINUS_D_AMP_ATGC_E.matcher(vn);
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

    private static final Pattern MINUS_D_AMP_ATGC_E = Pattern.compile("(-\\d+)&[ATGC]+$");

    // check whether there're reads supporting wild type in deletions
    // Only for deletions that have micro-homology
    private static boolean noPassingReads(String chr, int s, int e, String[] bams, Configuration conf) throws IOException {
        int cnt = 0;
        for (String bam : bams) {
            try (SamtoolsReader reader = new SamtoolsReader("view", bam, chr + ":" + s + "-" + e)) {
                String line;
                while ((line = reader.read()) != null) {
                    String[] a = line.split("\t");
                    int rs = toInt(a[3]);
                    if (a[5].matches("^(\\d+)M$")) {
                        int re = rs + toInt(substr(a[5], 0, -1));
                        if (re > e + 2 && rs < s - 2) {
                            cnt++;
                        }
                    }

                }

            }
        }
        if (conf.y) {
            System.err.printf("Passing Read CNT: %s %s %s %s\n", cnt, chr, s, e);
        }
        return cnt <= 0;
    }


    private static boolean ismatch(String seq1, String seq2, int dir, Configuration conf) {
        if (conf.y) {
            System.err.printf("Matching %s %s %s\n", seq1, seq2, dir);
        }
        seq2 = seq2.replaceAll("#|\\^", "");
        int mm = 0;
        Map<Character, Boolean> nts = new HashMap<>();
        for(int n = 0; n < seq1.length() && n < seq2.length(); n++) {
            nts.put( seq1.charAt(n), true);
            if (seq1.charAt(n) != substr(seq2, dir * n - (dir == -1 ? 1 : 0), 1).charAt(0)) {
                mm++;
            }
        }

        return (mm <= 2 && mm/(double)seq1.length() < 0.15);
    }


    // Find the consensus sequence in soft-clipped reads.  Consensus is called if
    // the matched nucleotides are >90% of all softly clipped nucleotides.
    private static String findconseq(Sclip scv) {
        if (scv.sequence != null) {
            return scv.sequence;
        }

        int total = 0;
        int match = 0;
        StringBuilder seq = new StringBuilder();
        for (Map<Character, Integer> nv : scv.nt.values()) {
            int max = 0;
            Character mnt = null;
            int tt = 0;
            for (Entry<Character, Integer> ent : nv.entrySet()) {
                Character nt = ent.getKey();
                int ncnt = ent.getValue();
                tt += ncnt;
                if (ncnt > max) {
                    max = ncnt;
                    mnt = nt;
                }
            }
            if (tt - max > 1 && max / (double)tt < 0.8) {
                break;
            }
            if ((tt - max > 2 || max <= tt - max) && max / (double)tt < 0.8) {
                break;
            }
            total += tt;
            match += max;
            seq.append(mnt);
        }

        scv.sequence = seq.toString();

        if (total > 0
                && match / (double) total > 0.9
                && seq.length()/1.5d > scv.nt.size() - seq.length()
                && (seq.length() /(double) scv.nt.size() > 0.8 || scv.nt.size() - seq.length() < 10 || seq.length() > 25)) {
            return scv.sequence;

        }

        return "";

    }

    private static Tuple3<List<Tuple3<String, Integer, Integer>>, List<Integer>, Integer> findMM5(Map<Integer, Character> ref, int p, String wupseq, int len, Map<Integer, Sclip> sclip5) {
        String seq = wupseq.replaceAll("#|\\^", "");
        int longmm = 3;
        List<Tuple3<String, Integer, Integer>> mm = new ArrayList<>(); // mismatches, mismatch positions, 5 or 3 ends
        int n = 0;
        int mn = 0;
        int mcnt = 0;
        StringBuilder str = new StringBuilder();
        List<Integer> sc5p = new ArrayList<>();
        while (ref.containsKey(p - n) && !ref.get(p - n).toString().equals(substr(seq, -1 - n, 1)) && mcnt < longmm) {
            str.insert(0, substr(seq, -1 - n, 1));
            mm.add(Tuple3.newTuple(str.toString(), p - n, 5 ));
            n++;
            mcnt++;
        }
        sc5p.add(p + 1);
        // Adject clipping position if only one mismatch
        if ( str.length() == 1 ) {
            while( ref.containsKey(p-n) && ref.get(p-n).toString().equals(substr(seq, -1-n, 1) )) {
                n++;
                mn++;
            }
            if (mn > 1) {
                sc5p.add(p -n);
            }
            if (sclip5.containsKey(p-n) && mn > 1 ) {
                sclip5.get(p-n).used = true;
            }
        }
        return Tuple3.newTuple(mm, sc5p, mn);

    }

    // Given a variant sequence, find the mismatches and potential softclipping positions
    private static Tuple3<List<Tuple3<String, Integer, Integer>>, List<Integer>, Integer> findMM3(Map<Integer, Character> ref, int p, String sanpseq, int len, Map<Integer, Sclip> sclip3) {
        String seq = sanpseq.replaceAll("#|\\^", ""); // ~ s/#|\^//g;
        final int longmm = 3;
        List<Tuple3<String, Integer, Integer>> mm = new ArrayList<>(); //mismatches, mismatch positions, 5 or 3 ends
        int n = 0;
        int mn = 0;
        int mcnt = 0;
        List<Integer> sc3p = new ArrayList<>();
        StringBuilder str = new StringBuilder();
        while(ref.containsKey(p + n) && ref.get(p + n).equals(seq.charAt(n))) {
            n++;
        }
        sc3p.add(p + n);
        int Tbp = p + n;
        while (ref.containsKey(p + n) && !ref.get(p + n).equals(seq.charAt(n)) && mcnt <= longmm && n < seq.length()) {
            str.append(seq.charAt(n));
            mm.add(Tuple3.newTuple(str.toString(), Tbp, 3 ));
            n++;
            mcnt++;

        }
        //  Adject clipping position if only one mismatch
        if (str.length() == 1) {
            while(ref.containsKey(p + n) && ref.get(p + n).equals(seq.charAt(n))) {
                n++;
                mn++;
            }
            if (mn > 1) {
                sc3p.add(p + n);
            }
            if (sclip3.containsKey(p + n) && mn > 1) {
                sclip3.get(p + n).used = true;
            }
        }
        //print STDERR "MM3: $seq $len $p '@sc3p'\n";

        return Tuple3.newTuple(mm, sc3p, mn);
    }

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


    private static void adjMNP(Map<Integer, Map<String, Variation>> hash,
            Map<Integer, Map<String, Integer>> mnp,
            Map<Integer, Integer> cov, Configuration conf) {

        for (Map.Entry<Integer, Map<String, Integer>> entry: mnp.entrySet()) {
            Integer p = entry.getKey();
            Map<String, Integer> v = entry.getValue();

            for (Map.Entry<String, Integer> en: v.entrySet()) {
                String vn = en.getKey();
                String mnt = vn.replace("&", "");
                for (int i = 0; i < mnt.length() - 1; i++) {
                    String left = substr(mnt, 0, i + 1);
                    String right = substr(mnt, -(mnt.length() - i - 1));
                    Variation vref = getVariation(hash, p, vn);
                    if (hash.containsKey(p) && hash.get(p).containsKey(left)) {
                        Variation tref = getVariation(hash, p, left);
                        if (tref.cnt < vref.cnt && tref.pmean / tref.cnt <= i + 1) {
                            adjCnt(vref, tref, conf);
                            hash.get(p).remove(left);
                        }
                    }
                    if (hash.containsKey(p + i + 1) && hash.get(p + i + 1).containsKey(right)) {
                        Variation tref = getVariation(hash, p + i + 1, right);
                        if (tref.cnt < vref.cnt && tref.pmean / tref.cnt <= mnt.length() - i - 1) {
                            adjCnt(vref, tref, conf);
                            incCnt(cov, p, tref.cnt);
                            hash.get(p + i + 1).remove(right);
                        }
                    }

                }

            }
        }

    }

    private static void adjCnt(Variation vref, Variation tv, Configuration conf) {
        adjCnt(vref, tv, null, conf);
    }

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
            System.err.printf("AdjCnt: '+' %s %s %s %s Ref: %s\n", vref.cnt, tv.cnt, vref.getDir(false), tv.getDir(false), refCnt);
            System.err.printf("AdjCnt: '-' %s %s %s %s Ref: %s\n", vref.cnt, tv.cnt, vref.getDir(true), tv.getDir(true), refCnt);
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

    // Adjust the reference count.
    private static void adjRefCnt(Variation tv, Variation ref, int len) {
        if (ref == null) {
            return;
        }
        double f = (tv.pmean / (double)tv.cnt - len + 1) / (tv.pmean / (double)tv.cnt);
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

    private static final BedRowFormat DEFAULT_BED_ROW_FORMAT = new BedRowFormat(2, 6, 7, 9, 10, 12);
    private static final BedRowFormat CUSTOM_BED_ROW_FORMAT = new BedRowFormat(0, 1, 2, 3, 1, 2);

    public static class BedRowFormat {
        public final int chrColumn;
        public final int startColumn;
        public final int endColumn;
        public final int geneColumn;
        public final int thickStartColumn;
        public final int thickEndColumn;

        public BedRowFormat(int chrColumn, int startColumn, int endColumn, int geneColumn, int thickStartColumn, int thickEndColumn) {
            this.chrColumn = chrColumn;
            this.startColumn = startColumn;
            this.endColumn = endColumn;
            this.geneColumn = geneColumn;
            this.thickStartColumn = thickStartColumn;
            this.thickEndColumn = thickEndColumn;
        }

    }

    private static void ampVardict(List<List<Region>> segs, Map<String, Integer> chrs, String ampliconBasedCalling, String bam1, String sample, Configuration conf) throws IOException {
        Map<String, Integer> splice = new HashMap<>();
        int rlen = 0;
        for (List<Region> regions : segs) {
            List<Vars> vars = new ArrayList<>();
            Region rg = null;
            Map<Integer, List<Tuple2<Integer, Region>>> pos = new HashMap<>();
            int j = 0;
            for (Region region : regions) {
                rg = region;
                Tuple2<Integer, Vars> tpl = toVars(region, bam1, chrs, splice, ampliconBasedCalling, rlen, conf);
                rlen = tpl._1();
                vars.add(tpl._2());
                for (int p = region.istart; p <= region.iend; p++) {
                    List<Tuple2<Integer, Region>> list = pos.get(p);
                    if (list == null) {
                        list = new ArrayList<>();
                        pos.put(p, list);
                    }
                    list.add(Tuple2.newTuple(j, region));
                }
                j++;
            }
            for (Entry<Integer, List<Tuple2<Integer, Region>>> entry : pos.entrySet()) {
                final int p = entry.getKey();
                final List<Tuple2<Integer, Region>> v = entry.getValue();

                List<Tuple2<Var, String>> gvs = new ArrayList<>();
                List<Var> ref = new ArrayList<>();
                String nt = null;
                double maxaf = 0;
                String vartype = "SNV";
                boolean flag = false;
                Var vref;
                int nocov = 0;
                int maxcov = 0;
                Set<Integer> goodmap = new HashSet<>();
                List<Integer> vcovs = new ArrayList<>();
                for (Tuple2<Integer, Region> amps : v) {
                    final int amp = amps._1();
                    final String chr = amps._2().chr;
                    final int S = amps._2().start;
                    final int E = amps._2().end;

                    List<Var> l = vars.get(amp).var.get(p);
                    Var refAmpP = vars.get(amp).ref.get(p);
                    if (l != null && !l.isEmpty()) {
                        Var tv = l.get(0);
                        vcovs.add(tv.tcov);
                        vartype = varType(tv.refallele, tv.varallele);
                        if (isGoodVar(tv, refAmpP, vartype, splice, conf)) {
                            gvs.add(Tuple2.newTuple(tv, chr + ":" + S + "-" + E));
                            if (nt != null && !tv.n.equals(nt)) {
                                flag = true;
                            }
                            if (tv.freq > maxaf) {
                                maxaf = tv.freq;
                                nt = tv.n;
                                vref = tv;
                            }
                            goodmap.add(amp);
                            if (tv.tcov > maxcov) {
                                maxcov = tv.tcov;
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

                for (int t : vcovs) {
                    if (t < maxcov / 50) {
                        nocov++;
                    }
                }

                if (gvs.size() > 1) {
                    Collections.sort(gvs, new Comparator<Tuple2<Var, String>>() {
                        @Override
                        public int compare(Tuple2<Var, String> o1, Tuple2<Var, String> o2) {
                            return Double.compare(o2._1().freq, o1._1().freq);
                        }
                    });
                }
                if (ref.size() > 1) {
                    Collections.sort(ref, new Comparator<Var>() {
                        @Override
                        public int compare(Var o1, Var o2) {
                            return Integer.compare(o2.tcov, o1.tcov);
                        }
                    });
                }

                if (gvs.isEmpty()) { // Only referenece
                    if (conf.doPileup) {
                        if (!ref.isEmpty()) {
                            vref = ref.get(0);
                        } else {
                            System.err.println(join("\t", sample, rg.gene, rg.chr, p, p, "", "", 0, 0, 0, 0, 0, 0, "", 0,
                                    "0;0", 0, 0, 0, 0, 0, "", 0, 0, 0, 0, 0, 0, "", "", 0, 0,
                                    rg.chr + ":", +p + "-" + p, "", 0, 0, 0, 0));
                            continue;
                        }
                    } else {
                        continue;
                    }
                } else {
                    vref = gvs.get(0)._1();
                }
                if (flag) { // different good variants detected in different amplicons
                    String gdnt = gvs.get(0)._1().n;
                    int gcnt = 0;
                    for (Tuple2<Integer, Region> amps : v) {
                        int amp = amps._1();
                        Var var = getVarMaybe(vars.get(amp), p, varn, gdnt);
                        if (var != null && isGoodVar(var, getVarMaybe(vars.get(amp), p, VarsType.ref), vartype, splice, conf)) {
                            gcnt++;
                        }
                    }
                    if (gcnt == gvs.size()) {
                        flag = false;
                    }
                }

                List<Tuple2<Var, String>> badv = new ArrayList<>();
                for (Tuple2<Integer, Region> amps : v) {
                    int amp = amps._1();
                    Region reg = amps._2();
                    if (goodmap.contains(amp)) {
                        continue;
                    }
                    if (vref.sp >= reg.istart && vref.ep <= reg.iend) {

                        String regStr = reg.chr + ":" + reg.start + "-" + reg.end;

                        if (vars.get(amp).var.containsKey(p)) {
                            badv.add(Tuple2.newTuple(vars.get(amp).var.get(p).get(0), regStr));
                        } else if (vars.get(amp).ref.containsKey(p)) {
                            badv.add(Tuple2.newTuple(vars.get(amp).ref.get(p), regStr));
                        } else {
                            badv.add(Tuple2.newTuple((Var)null, regStr));
                        }
                    }
                }
                if (vartype.equals("Complex")) {
                    adjComplex(vref);
                }
                System.out.print(join("\t", sample, rg.gene, rg.chr, joinVar1(vref, "\t"), gvs.get(0)._2(), vartype, gvs.size(), gvs.size() + badv.size(), nocov, flag ? 1 : 0));
                if (conf.debug) {
                    System.out.print("\t" + vref.DEBUG);
                }
                if (conf.y) {
                    for (int gvi = 0; gvi < gvs.size(); gvi++) {
                        Tuple2<Var, String> tp = gvs.get(gvi);
                        System.out.print("\tGood" + gvi + " " + join(" ", joinVar2(tp._1(), " "), tp._2()));
                    }
                    for (int bvi = 0; bvi < badv.size(); bvi++) {
                        Tuple2<Var, String> tp = badv.get(bvi);
                        System.out.print("\tBad" + bvi + " " + join(" ", joinVar2(tp._1(), " "), tp._2()));
                    }
                }
                System.out.println();
            }
        }
    }

    private static String joinVar2(Var var, String delm) {
        StringBuilder sb = new StringBuilder();
        sb.append(var.tcov).append(delm);
        sb.append(var.cov).append(delm);
        sb.append(var.rfc).append(delm);
        sb.append(var.rrc).append(delm);
        sb.append(var.fwd).append(delm);
        sb.append(var.rev).append(delm);
        sb.append(var.genotype).append(delm);
        sb.append(var.freq).append(delm);
        sb.append(var.bias).append(delm);
        sb.append(var.pmean).append(delm);
        sb.append(var.pstd).append(delm);
        sb.append(var.qual).append(delm);
        sb.append(var.qstd).append(delm);
        sb.append(var.mapq).append(delm);
        sb.append(var.qratio).append(delm);
        sb.append(var.hifreq).append(delm);
        sb.append(var.extrafreq);
        return null;
    }

    private static String joinVar1(Var var, String delm) {
        StringBuilder sb = new StringBuilder();
        sb.append(var.sp).append(delm);
        sb.append(var.ep).append(delm);
        sb.append(var.refallele).append(delm);
        sb.append(var.varallele).append(delm);
        sb.append(var.tcov).append(delm);
        sb.append(var.cov).append(delm);
        sb.append(var.rfc).append(delm);
        sb.append(var.rrc).append(delm);
        sb.append(var.fwd).append(delm);
        sb.append(var.rev).append(delm);
        sb.append(var.genotype).append(delm);
        sb.append(var.freq).append(delm);
        sb.append(var.bias).append(delm);
        sb.append(var.pmean).append(delm);
        sb.append(var.pstd).append(delm);
        sb.append(var.qual).append(delm);
        sb.append(var.qstd).append(delm);
        sb.append(var.mapq).append(delm);
        sb.append(var.qratio).append(delm);
        sb.append(var.hifreq).append(delm);
        sb.append(var.extrafreq).append(delm);
        sb.append(var.shift3).append(delm);
        sb.append(var.msi).append(delm);
        sb.append(var.msint).append(delm);
        sb.append(var.nm).append(delm);
        sb.append(var.hicnt).append(delm);
        sb.append(var.hicov).append(delm);
        sb.append(var.leftseq).append(delm);
        sb.append(var.rightseq);
        return sb.toString();
    }


    private static void adjComplex(Var vref) {
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
            vref.leftseq = substr(vref.leftseq, 0, -n);
            vref.leftseq += substr(refnt, 0, n);
        }
        refnt = vref.refallele;
        varnt = vref.varallele;
        n = 1;
        while( refnt.length() - n > 0 && varnt.length() - n > 0 && substr(refnt, -n, 1).equals(substr(varnt, -n, 1))) {
            n++;
        }
        if (n > 1) {
            vref.ep -= n - 1;
            vref.refallele = substr(refnt, 0, 1 - n);
            vref.varallele = substr(varnt, 0, 1 - n);
            vref.rightseq = substr(refnt, 1 - n, n - 1) + substr(vref.rightseq, 0, 1 - n);
        }


    }


    static enum VarsType {varn,ref,var,}

    private static Var getVarMaybe(Vars vars, int key, VarsType type, Object... keys) {
        switch (type) {
            case var:
                if (vars.var.containsKey(key) && vars.var.get(key).size() > (Integer)keys[0]) {
                    return vars.var.get(key).get((Integer)keys[0]);
                }
            case varn:
                if (vars.varn.containsKey(key)) {
                    return vars.varn.get(key).get(keys[0]);
                }
            case ref:
                return vars.ref.get(key);
        }
        return null;
    }

    public static String join(String delim, Object... args) {
        if (args.length == 0) {
            return "";
        }
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < args.length; i++) {
            sb.append(args[i]);
            if (i + 1 != args.length) {
                sb.append(delim);
            }
        }
        return sb.toString();
    }

    private static boolean isGoodVar(Var vref, Var rref, String type, Map<String, Integer> splice, Configuration conf) {
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
            double d = vref.mapq + vref.refallele.length() + vref.varallele.length();
            if ((d - 2 < 5 && rref.mapq > 20)
                    || (1 + d) / (rref.mapq + 1) < 0.25d) {
                return false;
            }
        }

        int spl = splice.containsKey(vref.sp + "-" + vref.ep) ? splice.get(vref.sp + "-" + vref.ep) : 0;

        if (type.equals("Deletion") && spl != 0) {
            return false;
        }
        if(vref.qratio < conf.qratio) {
            return false;
        }
        if (vref.freq > 0.35d) {
            return true;
        }
        if(vref.mapq < conf.mapq) {
            return false;
        }
        if (vref.msi >= 13 && vref.freq <= 0.275d && vref.msint == 1) {
            return false;
        }
        if (vref.msi >= 8 && vref.freq <= 0.2d && vref.msint > 1) {
            return false;
        }
        if ( vref.bias.equals("2;1") && vref.freq < 0.25d ) {
            if (type == null || type.isEmpty() || type.equals("SNV") || (vref.refallele.length() <= 3 && vref.varallele.length() <= 3)) {
                return false;
            }
        }

        return true;
    }


    private static String varType(String ref, String var) {
        if (ref.length() == 1 && var.length() == 1) {
            return "SNV";
        } else if (ref.charAt(0) != var.charAt(0)) {
            return "Complex";
        } else if (ref.length() == 1 && var.length() > 1 && var.startsWith(ref)) {
            return "Insertion";
        } else if (ref.length() > 1 && var.length() == 1 && ref.startsWith(var)) {
            return "Deletion";
        }
        return "Complex";
    }

    public static List<String> globalFind(Pattern pattern, String string) {
        List<String> result = new LinkedList<>();
        Matcher matcher = pattern.matcher(string);
        while(matcher.find()) {
            result.add(matcher.group(1));
        }
        return result;

    }

    public static class ToVarsContext {
        Map<Integer, Character> ref = new HashMap<>();
        Map<Integer, Map<String, Variation>> hash = new HashMap<>();
        Map<Integer, Map<String, Variation>> iHash = new HashMap<>();
        Map<Integer, Integer> cov = new HashMap<>();
        Map<Integer, Sclip> sclip3 = new HashMap<>();
        Map<Integer, Sclip> sclip5 = new HashMap<>();
        Map<Integer, Map<String, Integer>> inc = new HashMap<>();
    }

    public static class Tuple2 <T1, T2> {

        private final T1 _f;
        private final T2 _s;

        public Tuple2(T1 f, T2 s) {
            _f = f;
            _s = s;
        }

        public T1 _1() {
            return _f;
        }

        public T2 _2() {
            return _s;
        }

        public static <T1, T2> Tuple2<T1, T2> newTuple(T1 f, T2 s) {
            return new Tuple2<T1, T2>(f, s);
        }
    }

    public static class Tuple3 <T1, T2, T3> {

        private final T1 _f;
        private final T2 _s;
        private final T3 _t;

        public Tuple3(T1 f, T2 s, T3 t) {
            _f = f;
            _s = s;
            _t = t;
        }

        public T1 _1() {
            return _f;
        }

        public T2 _2() {
            return _s;
        }

        public T3 _3() {
            return _t;
        }

        public static <T1, T2, T3> Tuple3<T1, T2, T3> newTuple(T1 f, T2 s, T3 t) {
            return new Tuple3<T1, T2, T3>(f, s, t);
        }
    }

    public static int sum(List<?> list) {
        int result = 0;
        for (Object object : list) {
            result += toInt(String.valueOf(object));
        }
        return result;
    }

    public static int toInt(String intStr) {
        return Integer.parseInt(intStr);
    }

    public static String substr(String string, int idx) {
        if(idx >= 0) {
            return string.substring(idx);
        } else {
            return string.substring(Math.max(0, string.length() + idx));
        }
    }

    public static String substr(String string, int begin, int len) {
        if (begin < 0) {
            begin = string.length() + begin;
        }
        if (len > 0) {
            return string.substring(begin, Math.min(begin + len, string.length()));
        } else if (len == 0) {
            return "";
        } else {
            int end = string.length() + len;
            if (end < begin) {
                return "";
            }
            return string.substring(begin, end);
        }
    }

    // TODO validate region format chr:start[-end][:gene]
    public static List<List<Region>> buildRegions(String region, final int numberNucleotideToExtend, final boolean zeroBased) {
        List<List<Region>> segs = new ArrayList<>();
        String[] split = region.split(":");
        String chr = split[0];
        String gene = split.length < 3 ? chr : split[2];
        String[] range = split[1].split("-");
        int start = toInt(range[0].replaceAll(",", ""));
        int end = range.length < 2 ? start : toInt(range[1].replaceAll(",", ""));
        start -= numberNucleotideToExtend;
        end += numberNucleotideToExtend;
        if (zeroBased)
            start++;
        if (start > end)
            start = end;
        segs.add(Collections.singletonList(new Region(chr, start, end, gene)));

        return segs;
    }
}
