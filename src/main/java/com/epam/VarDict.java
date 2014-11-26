package com.epam;

import static java.lang.String.format;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
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
        boolean zeroBased;
        String ampliconBasedCalling; //-a
        int columnForChromosome = -1;
        String sampleNameRegexp; // -n
        String sampleName; //-N
        String fasta;
        Bam bam;
        Double downsampling;
        boolean chromosomeNameIsNumber;
        Integer mappingQuality;//-Q
        boolean  nonPrimaryAlignment; //-F
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
        boolean debug = false; // -D
        double freq = 0.5; // -f and -p
        boolean  moveIndelsTo3 = false;


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

    public static void buildRegionsFromFile(Configuration conf) throws IOException {
        String a = conf.ampliconBasedCalling;
        List<String> segraw = new ArrayList<>();
        boolean zeroBased = conf.zeroBased;

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
                            if (a6 > a1 && a7 < a2) {
                                a = "10:0.95";
                                zeroBased = true;
                            }
                        } catch (NumberFormatException e) {
                            continue;
                        }
                    }
                    segraw.add(line);
                }

            }
        }
        List<List<Region>> segs = new LinkedList<>();

        if (a != null) {
            Map<String, List<Region>> tsegs = new HashMap<>();
            for (String string : segraw) {
                String[] split = string.split(conf.delimiter);
                String chr = split[0];
                int start = toInt(split[1]);
                int end = toInt(split[2]);
                String gene = split[3];
                int istart = toInt(split[6]);
                int iend = toInt(split[7]);
                if (zeroBased) {
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
            ampVardict(segs);
        } else {
            BedRowFormat format = DEFAULT_BED_ROW_FORMAT;
            for (String seg : segraw) {
                String[] splitA = seg.split(conf.delimiter);
                if (!conf.isColumnForChromosomeSet() && splitA.length == 4) {
                    try {
                        int a1 = toInt(splitA[1]);
                        int a2 = toInt(splitA[2]);
                        if (a1 <= a2) {
                            format = CUSTOM_BED_ROW_FORMAT;
                        }
                    } catch (NumberFormatException e) {
                    }
                }
                String chr = splitA[format.chrColumn];
                int cdss = toInt(splitA[format.startColumn]);
                int cdse = toInt(splitA[format.endColumn]);
                String gene = format.geneColumn < splitA.length ? splitA[format.geneColumn] : chr;

                String[] starts = splitA[format.thickStartColumn].split(","); // TODO why?
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
                    if (zeroBased)
                        s++;
                    cds.add(new Region(chr, s, e, gene));
                }
                segs.add(cds);
            }
        }
    }

    public static class Bam {
        private final String[] bamNames;
        private final String[] bams;

        public Bam(String value) {
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
    }

    public void finalCycle(List<List<Region>> segs, Map<String, Integer> chrs, Map<String, Integer> SPLICE, Configuration conf, String ampliconBasedCalling) throws IOException {
        for (List<Region> list : segs) {
            for (Region region : list) {
                if (conf.bam.hasBam2()) {
                    somdict(region, toVars(region, conf.bam.getBam1(), chrs, SPLICE, conf, ampliconBasedCalling), toVars(region, conf.bam.getBam2(), chrs, SPLICE, conf, ampliconBasedCalling));
                } else {
                    vardict(region, toVars(region, conf.bam.getBam1(), chrs, SPLICE, conf, ampliconBasedCalling));
                }
            }
        }
    }

    final static Pattern SAMPLE_PATTERN = Pattern.compile("([^\\/\\._]+).sorted[^\\/]*.bam");
    final static Pattern SAMPLE_PATTERN2 = Pattern.compile("([^\\/]+)[_\\.][^\\/]*bam");

    public static List<String> getSample(String bam1, String bam2, String sampleName, String regexp) {
        String sample = null;
        String samplem = "";

        if (sampleName != null) {
            sample = sampleName;
        } else if (regexp != null) {
            Pattern rn = Pattern.compile(regexp);
            Matcher m = rn.matcher(bam1);
            if (m.find()) {
                sample = m.group(1);
            }
        } else {
            Matcher matcher = SAMPLE_PATTERN.matcher(bam1);
            if (matcher.find()) {
                sample = matcher.group(1);
            } else {
                matcher = SAMPLE_PATTERN2.matcher(bam1);
                if (matcher.find()) {
                    sample = matcher.group(1);
                }
            }
        }

        if (bam2 != null) {
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
        }

        return Collections.unmodifiableList(Arrays.asList(sample, samplem));
    }

    private void somdict(Region region, Object vars, Object vars2) {
        // TODO Auto-generated method stub

    }

    private void vardict(Region region, Object vars) {
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

    public static class Sclip {
        Variation variation = new Variation();
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
    private static final Pattern D_M_D_ID = Pattern.compile("^(\\d+)M(\\d+)([ID])");
    private static final Pattern D_ID_D_M = Pattern.compile("(\\d+)([ID])(\\d+)M$");
    private static final Pattern D_M_D_DD_M_D_I_D_M_D_DD = Pattern.compile("^(.*?)(\\d+)M(\\d+)D(\\d)M(\\d+)I(\\d)M(\\d+)D");
    private static final Pattern D_M_D_DD_M_D_I_D_M_D_DD_prim = Pattern.compile("(\\d+)M(\\d+)D(\\d)M(\\d+)I(\\d)M(\\d+)D");
    private static final Pattern D_MIS = Pattern.compile("(\\d+)[MIS]");
    private static final Pattern D_MD = Pattern.compile("(\\d+)[MD]");
    private static final Pattern D_DD_D_M_D_DD_DI = Pattern.compile("(\\d+)D(\\d)M(\\d+)D(\\d+I)?");


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


    private static Object toVars(Region region, String bam, Map<String, Integer> chrs, Map<String, Integer> SPLICE, Configuration conf, String ampliconBasedCalling) throws IOException {
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
            if (conf.chromosomeNameIsNumber && chr .startsWith("chr")) {
                chr = region.getChr().substring("chr".length());
            }
            try (SamtoolsReader reader = new SamtoolsReader("faidx", "view", bami, chr + ":" + s_start + "-" + s_end)) {
                Map<String, Boolean> dup = new HashMap<>();
                int dupp = -1;
                String line;
                while ((line = reader.read()) != null) {
                    if (conf.isDownsampling() && RND.nextDouble() <= conf.downsampling) {
                        continue;
                    }
                    String[] a = line.split("\t");
                    Flag alignmentFlag = new Flag(toInt(a[1]));

                    assert a.length > 4;
                    final int Qmean = toInt(a[4]);
                    final String querySEQ = ampliconBasedCalling.length() > 9 ? a[9] : null;

                    if (conf.hasMappingQuality() && Qmean < conf.mappingQuality) { // ignore low mapping quality reads
                        continue;
                    }
                    if (conf.nonPrimaryAlignment && alignmentFlag.isNotPrimaryAlignment()) { // ignore "non-primary alignment"
                        continue;
                    }
                    if (alignmentFlag.isNotPrimaryAlignment() && "*".equals(querySEQ)) {
                        continue;
                    }

                    assert a.length > 6;
                    final String mrnm = a[6];
                    int position = toInt(a[3]);
                    if (conf.removeDuplicatedReads) {
                        if (position != dupp) {
                            dup.clear();
                        }
                        if (mrnm.equals("=") && a.length > 7) {
                            if (dup.containsKey(position + "-" + a[7])) {
                                continue;
                            }
                            dup.put(position + "-" + a[7], Boolean.TRUE);
                            dupp = position;
                        }
                    }

                    int nm = 0;
                    int idlen = 0;
                    String cigarStr = a[5];
                    for (String s : globalFind(INDEL, cigarStr)) {
                        idlen += toInt(s);
                    }
                    Matcher nmMatcher = NUMBER_MISMATCHES.matcher(line);
                    if (nmMatcher.find()) { // number of mismatches. Don't use NM since it includes gaps, which can be from indels
                        nm = toInt(nmMatcher.group(1)) - idlen;
                        if (nm > conf.mismatch) { // edit distance - indels is the # of mismatches
                            continue;
                        }
                    } else {
                        if (conf.y && !cigarStr.equals("*")) {
                            System.err.println("No XM tag for mismatches. " + line);
                            continue;
                        }
                    }
                    int n = 0; // keep track the read position, including softclipped
                    int p = 0; // keep track the position in the alignment, excluding softclipped
                    boolean dir = alignmentFlag.isReverseStrand();
                    List<String> segs = globalFind(MATCH_INSERTION, cigarStr); // Only match and insertion counts toward read length
                    List<String> segs2 = globalFind(SOFT_CLIPPED, cigarStr); // For total length, including soft-clipped bases
                    int rlen = sum(segs); //The read length for matched bases
                    int rlen2 = sum(segs2); //The total length, including soft-clipped bases
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
                        int rlen3 = sum(globalFind(ALIGNED_LENGTH, cigarStr)); // The total aligned length, excluding soft-clipped bases and insertions
                        int segstart = position;
                        int segend = segstart + rlen3 - 1;

                        if (cigarStr.matches("^(\\d+)S.*")) {
                            int ts1 = segstart > s_start ? segstart : s_start;
                            int te1 = segend < s_end ? segend : s_end;
                            if (Math.abs(ts1 - te1) / (double)(segend - segstart) > ovlp == false) {
                                continue;
                            }
                        } else if (cigarStr.matches("^.*(\\d+)S$")) {
                            int ts1 = segstart > s_start ? segstart : s_start;
                            int te1 = segend < s_end ? segend : s_end;
                            if (Math.abs(te1 - ts1) / (double)(segend - segstart) > ovlp == false) {
                                continue;
                            }

                        } else {
                          if (mrnm.equals("=") && a.length > 8) {
                              int isize = toInt(a[8]);
                              if (isize > 0) {
                                  segend = segstart + isize -1;
                              } else {
                                  segstart = toInt(a[7]);
                                  segend = segstart - isize - 1;
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
                    if (alignmentFlag.isUnmappedMate()) {
                        // to be implemented
                    } else {
                        if (mrnm.equals("=")) {
                            if (line.matches(".*\\tSA:Z:(\\S+).*")) {
                                if (alignmentFlag.isSupplementaryAlignment()) { // the supplementary alignment
                                    continue; // Ignore the supplmentary for now so that it won't skew the coverage
                                }
                            }
                        }

                    }

                    //Modify the CIGAR for potential mis-alignment for indels at the end of reads to softclipping and let VarDict's algorithm to figure out indels
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
                        cm = D_M_D_ID.matcher(cigarStr);
                        if (cm.find()) {
                            int tmid = toInt(cm.group(1));
                            if (tmid <= 8 ) {
                                String tslen = tmid + (cm.group(3).equals("I") ? toInt(cm.group(2)) : 0) + "S";
                                position = tmid + (cm.group(3).equals("D") ? toInt(cm.group(2)) : 0);
                                cigarStr = cm.replaceFirst(tslen);
                                flag = true;
                            }
                        }
                        cm = D_ID_D_M.matcher(cigarStr);
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
                                while (rdoff + rn < querySEQ.length()
                                        && Character.valueOf(querySEQ.charAt(rdoff + rn)).equals(ref.get(refoff + rn))) {
                                    rn++;
                                }
                                RDOFF += rn;
                                dlen -= rn;
                                tslen -= rn;
                                cigarStr = D_M_D_DD_M_D_I_D_M_D_DD_prim.matcher(cigarStr).replaceFirst(RDOFF + "M" + dlen + "D" + tslen + "I");
                                flag = true;
                            }
                        }
                        cm = D_DD_D_M_D_DD_DI.matcher(cigarStr);
                        if (cm.find()) {
                            int ilen = toInt(cm.group(2));
                            int dlen = ilen + toInt(cm.group(1)) + toInt(cm.group(3));
                            String istr = cm.group(4);
                            if (istr != null) {
                                ilen += toInt(istr.substring(0, istr.length() - 1));
                            }
                            cigarStr = cm.replaceFirst(dlen + "D" + ilen + "I");
                        }
                    }

                    List<String> cigar = new ArrayList<>();
                    Matcher cigarM = CIGAR_PAIR.matcher(cigarStr);
                    while (cigarM.find()) {
                        cigar.add(cigarM.group(1));
                        cigar.add(cigarM.group(2));
                    }

                    int start = position;
                    assert a.length > 10;
                    final String quality = a[10];

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
                                            && ref.get(start - 1).equals(querySEQ.charAt(m - 1))
                                            && quality.charAt(m - 1) - 33 > 10) {
                                        Variation variation = getVariation(hash, start - 1, ref.get(start - 1).toString());
                                        if (variation.cnt != 0) {
                                            variation.cnt = 0;
                                        }
                                        addCnt(variation, dir, m, quality.charAt(m - 1) - 33, Qmean, nm, conf.goodq);
                                        incCnt(cov, start - 1, 1);
                                        start--;
                                        m--;
                                    }
                                    if (m > 0) {
                                        int q = 0;
                                        int qn = 0;
                                        int lowqcnt = 0;
                                        for (int si = m - 1; si >= 0; si--) {
                                            if (querySEQ.charAt(si) == 'N') {
                                                break;
                                            }
                                            int tq = quality.charAt(si) - 33;
                                            if (tq < 7)
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
                                                Character ch = querySEQ.charAt(si);
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
                                                addCnt(variation, dir, si - (m - qn), quality.charAt(si) - 33, Qmean, nm, conf.goodq);
                                            }
                                            if (sclip.variation.cnt != 0)
                                                sclip.variation.cnt = 0;
                                            addCnt(sclip.variation, dir, m, q / qn, Qmean, nm, conf.goodq);
                                        }

                                    }
                                    m = toInt(cigar.get(ci));
                                } else if (ci == cigar.size() - 2) { // 3' soft clipped
                                    int qmean = quality.charAt(n) - 33;
                                    while (n < querySEQ.length()
                                            && ref.containsKey(start)
                                            && ref.get(start).equals(querySEQ.charAt(n))
                                            && qmean > 10) {

                                        Variation variation = getVariation(hash, start, ref.get(start).toString());
                                        if (variation.cnt != 0) {
                                            variation.cnt = 0;
                                        }
                                        addCnt(variation, dir, rlen2 - p, qmean, Qmean, nm, conf.goodq);
                                        incCnt(cov, start, 1);

                                        n++;
                                        start++;
                                        m--;
                                        p++;
                                    }
                                    if (querySEQ.length() - n > 0) {
                                        int q = 0;
                                        int qn = 0;
                                        int lowqcnt = 0;
                                        for (int si = 0; si < m; si++) {
                                            if ( querySEQ.charAt(n+si) == 'N' ) {
                                                break;
                                            }
                                            int tq = quality.charAt(n + si) - 33;
                                            if (tq < 7) {
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
                                                Character ch = querySEQ.charAt(n + si);
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
                                                addCnt(variation, dir, qn -  si, quality.charAt(n + si) - 33, Qmean, nm, conf.goodq);
                                            }
                                            if (sclip.variation.cnt != 0)
                                                sclip.variation.cnt = 0;
                                            addCnt(sclip.variation, dir, m, q / qn, Qmean, nm, conf.goodq);
                                        }

                                    }

                                }
                                n += m;
                                offset = 0;
                                start = toInt(a[3]);  // had to reset the start due to softclipping adjustment
                                continue;
                            case "H":
                                continue;
                            case "I": {
                                StringBuilder s = new StringBuilder(substr(querySEQ, n, m));
                                StringBuilder q = new StringBuilder(substr(quality, n, m));
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
                                    s.append("#").append(substr(querySEQ, n + m, ci2));
                                    q.append(substr(quality, n + m, ci2));
                                    s.append('^').append(cigar.get(ci + 5).equals("I") ? substr(querySEQ, n + m + ci2, ci4) : ci4);
                                    q.append(cigar.get(ci + 5).equals("I") ? substr(quality, n + m + ci2, ci4) : quality.charAt(n + m + ci2));
                                    multoffs += ci2 + (cigar.get(ci + 5).equals("D") ? ci4 : 0);
                                    multoffp += ci2 + (cigar.get(ci + 5).equals("I") ? ci4 : 0);
                                    ci += 4;
                                } else {
                                    if (cigar.size() > ci + 3 && cigar.get(ci + 3).contains("M")) {
                                        int ci2 = toInt(cigar.get(ci + 2));
                                        for (int vi = 0; vi < conf.vext && vi < ci2; vi++) {
                                            if (querySEQ.charAt(n + m + vi) == 'N') {
                                                break;
                                            }
                                            if (quality.charAt(n + m + vi) - 33 < conf.goodq) {
                                                break;
                                            }
                                            if (ref.containsKey(start + vi) && String.valueOf(querySEQ.charAt(n + m + vi)).equals(ref.get(start + vi))) {
                                                offset = vi + 1;
                                            }
                                        }
                                        if (offset != 0) {
                                            ss.append(substr(querySEQ, n+m, offset));
                                            q.append(substr(quality, n+m, offset));
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
                                    hv.Qmean += Qmean;
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
                                            && String.valueOf(querySEQ.charAt(n-1)).equals(ref.get(start - 1))) {

//                                        subCnt(getVariation(hash, start - 1, ref.get(start - 1 ).toString()), dir, tp, tmpq, Qmean, nm, conf);
                                        subCnt(getVariation(hash, start - 1, String.valueOf(querySEQ.charAt(n - 1))), dir, tp, quality.charAt(n - 1) - 33, Qmean, nm, conf);
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
                                        ttref.Qmean += Qmean;
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
                                    char q1 = quality.charAt(n - 1);
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

                                        s.append("#").append(substr(querySEQ, n, ci2));
                                        q.append(substr(quality, n, ci2));
                                        s.append('^').append(cigar.get(ci + 5).equals("I") ? substr(querySEQ, n + ci2, ci4) : ci4);
                                        q.append(cigar.get(ci + 5).equals("I") ? substr(quality, n + ci2, ci4) : "");
                                        multoffs += ci2 + (cigar.get(ci + 5).equals("D") ? ci4 : 0);
                                        multoffp += ci2 + (cigar.get(ci + 5).equals("I") ? ci4 : 0);
                                        ci += 4;
                                    } else if (cigar.size() > ci + 3 && cigar.get(ci + 3).equals("I")) {
                                       int ci2 = toInt(cigar.get(ci + 2));
                                       s.append("^").append(substr(querySEQ, n, ci2));
                                       q.append(substr(quality, n, ci2));
                                       multoffp += ci2;
                                       ci += 2;
                                    } else {
                                        if (cigar.size() > ci + 3 && cigar.get(ci + 3).contains("M")) {
                                            int ci2 = toInt(cigar.get(ci + 2));
                                            for (int vi = 0; vi < conf.vext && vi < ci2; vi++) {
                                                if (querySEQ.charAt(n + vi) == 'N') {
                                                    break;
                                                }
                                                if (quality.charAt(n + vi) - 33 < conf.goodq) {
                                                    break;
                                                }
                                                if (ref.containsKey(start + m + vi) && String.valueOf(querySEQ.charAt(n + vi)).equals(ref.get(start + m + vi))) {
                                                    offset = vi + 1;
                                                }
                                            }
                                            if (offset != 0) {
                                                ss.append(substr(querySEQ, n + m, offset));
                                                q.append(substr(quality, n + m, offset));
                                                for (int osi = 0; osi < offset; osi++) {
                                                    incCnt(cov, start + osi, 1);
                                                }
                                            }
                                        }
                                    }
                                    if ( offset > 0 ) {
                                        s.append("&").append(ss);
                                    }
                                    char q2 = quality.charAt(n + offset);
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
                                        hv.Qmean += Qmean;
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
                            String s = String.valueOf(querySEQ.charAt(n));
                            if (s.equals("N")) {
                                start++;
                                n++;
                                p++;
                                continue;
                            }
                            int q = quality.charAt(n) - 33;
                            int qbases = 1;
                            // for more than one nucleotide mismatch
                            StringBuilder ss = new StringBuilder();
                            // More than one mismatches will only perform when all nucleotides have quality > $GOODQ
                            // Update: Forgo the quality check.  Will recover later
                            while((start + 1) >= region.start
                                    && (start + 1) <= region.end && (i + 1) < m
                                    && ref.containsKey(start) && !ref.get(start).equals(querySEQ.charAt(n))) {

                                if (querySEQ.charAt(n + 1) == 'N' ) {
                                    break;
                                }
                                if ( ref.containsKey(start + 1) && !ref.get(start + 1).equals(querySEQ.charAt(n +  1))) {
                                    ss.append(querySEQ.charAt(n +  1));
                                    q += quality.charAt(n + 1) - 33;
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
                                    && (ss.length() > 0 || Character.valueOf(querySEQ.charAt(n)).equals(ref.get(start)))
                                    && quality.charAt(n) - 33 > conf.goodq) {
                                while (i + 1 < m) {
                                    s += querySEQ.charAt(n + 1);
                                    q += quality.charAt(n + 1) - 33;
                                    qbases++;
                                    i++;
                                    n++;
                                    p++;
                                    start++;
                                }
                                s = s.replace("&", "");
                                s = "-" + cigar.get(ci + 2) + "&" + s;

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
                                    q = q / qbases;
                                    if(hv.pstd == false && hv.pp != 0 && tp != hv.pp) {
                                        hv.pstd = true;
                                    }
                                    if(hv.qstd == false && hv.pq != 0 && tp != hv.pq) {
                                        hv.qstd = true;
                                    }
                                    hv.pmean += tp;
                                    hv.qmean += q;
                                    hv.Qmean += Qmean;
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
                                    if (s.contains("-")) {
                                        increment(dels5, start - qbases + 1, s);
                                        for (int qi = 0; qi < cigar.get(ci + 2).length(); qi++) {
                                            incCnt(cov, start + qi, 1);
                                        }
                                        start += cigar.get(ci + 2).length();
                                        ci += 2;
                                    }
                                }
                            }
                            if (operation.equals("I")) {
                                start++;
                            } else if (operation.equals("D")) {
                                n++;
                                p++;
                            }
                        }
                        offset = 0;

                    }

                }
            } //TODO abcall

        }
        adjMNP(hash, mnp, cov);

        if (conf.performLocalRealignment) {
            realigndel(hash, dels5, cov, sclip5, sclip3, ref, region.chr, chrs, bams, conf);
            realignins(hash, iHash, ins, cov, sclip5, sclip3, ref, region.chr, chrs, conf);
            realignlgdel(hash, cov, sclip5, sclip3, ref, region.chr, chrs, bams, conf);
            realignlgins(hash, iHash, cov, sclip5, sclip3, ref, region.chr, chrs, bams, conf);
            realignlgins30(hash, iHash, cov, sclip5, sclip3, ref, region.chr, chrs, bams, conf);
        }


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
                Var tvref = new Var();
                tvref.n = n;
                tvref.cov = cnt.cnt;
                tvref.fwd = fwd;
                tvref.rev = rev;
                tvref.bias = bias;
                tvref.freq = (fwd + rev) / (double)tcov;
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

                    Var tvref = new Var();
                    tvref.n = n;
                    tvref.cov = cnt.cnt;
                    tvref.fwd = fwd;
                    tvref.rev = rev;
                    tvref.bias = bias;
                    tvref.freq = (fwd + rev) / (double)tcov;
                    tvref.pmean = round(cnt.pmean / (double)cnt.cnt, 1);
                    tvref.pstd = cnt.pstd;
                    tvref.qual = vqual;
                    tvref.qstd = cnt.qstd;
                    tvref.mapq = mq;
                    tvref.qratio = round(hicnt / (locnt != 0 ? locnt : 0.5d), 3);
                    tvref.hifreq = hicov > 0 ? hicnt / (double)hicov : 0;
                    tvref.extrafreq = 0;
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
                    return Double.compare(o1.qual * o1.cov, o2.qual * o2.cov);
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
                        //TODO ???
                    }

                }
            }

        }

        return null;
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
        int bias;
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
            String[] bams,
            Configuration conf) throws IOException {

        List<Tuple2<Integer, Sclip>> tmp5 = new ArrayList<>();
        for (Entry<Integer, Sclip> ent5 : sclip5.entrySet()) {
            tmp5.add(Tuple2.newTuple(ent5.getKey(), ent5.getValue()));
        }
        Collections.sort(tmp5, COMP2); //  TODO sort {$a->[1] <=> $b->[1];} ????

        List<Tuple2<Integer, Sclip>> tmp3 = new ArrayList<>();
        for (Entry<Integer, Sclip> ent3 : sclip3.entrySet()) {
            tmp3.add(Tuple2.newTuple(ent3.getKey(), ent3.getValue()));
        }
        Collections.sort(tmp3, COMP2); //  TODO sort {$a->[1] <=> $b->[1];} ????

        for (Tuple2<Integer,Sclip> t5 : tmp5) {
            final int p5 = t5._1();
            final Sclip sc5v = t5._2();
            if (sc5v.used) {
                continue;
            }

            for (Tuple2<Integer,Sclip> t3 : tmp3) {
                final int p3 = t3._1();
                final Sclip sc3v = t3._2();
                if (sc3v.used) {
                    continue;
                }
                if (p3 < p5) {
                    continue;
                }
                if (p3 - p5 > 100) { // if they're too far away, don't even try
                    break;
                }
                final String seq5 = findconseq(sc5v);
                final String seq3 = findconseq(sc3v);
                if (seq5.length() <= 15 && seq3.length() <= 15) {
                    continue;
                }
                if (conf.y) {
                    System.err.printf("lgins30: %s %s %s %s \n", p3, p5, seq3, seq5);
                }
                Tuple3<Integer, Integer, Integer> tpl = find35match(seq5, seq3);
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
                StringBuilder tmp = new StringBuilder();
                if (conf.y) {
                    System.err.printf("lgins30: %s %s %s %s\n", p3, p5, ins.length(), ins);
                }
                if ( ins.length() <= p3 - p5 ) { // Tandem duplication
                    for(int p = p5; p < p5 + (p3 - p5 - ins.length())/2; p++) {
                        Character ch = ref.get(p);
                        if (ch != null) {
                            tmp.append(ch);
                        }
                    }
                } else {
                    for(int p = p5; p < p3; p++) {
                        Character ch = ref.get(p);
                        if (ch != null) {
                            tmp.append(ch);
                        }
                    }
                }
                ins = tmp + ins;
                sc3v.used = true;
                sc5v.used = true;
                final int bi = p5 - 1;
                final Variation iref = getVariation(iHash, bi, "+" + ins);
                iref.cnt = 0;
                iref.pstd = true;
                iref.qstd = true;
                Variation mref = getVariationMaybe(hash, bi, ref.get(bi));
                adjCnt(iref, sc3v.variation, mref);
                adjCnt(iref, sc5v.variation);
                if (ins.length() <= p3 - p5 && bams.length > 0 && p3 - p5 >= 5 && p3 - p5 < 75
                        && mref != null
                        && noPassingReads(chr, p5, p3, bams, conf)
                        && iref.cnt > 2 * mref.cnt) {

                    adjCnt(iref, mref, mref);
                }

                final String key = ins;
                Map<Integer, Map<String, Integer>> tins = new HashMap() {{
                    put(bi, new HashMap() {{
                            put("+" + key, iref.cnt);
                        }});
                }};
                realignins(hash, iHash, tins, cov, sclip5, sclip3, ref, chr, chrs, conf);
            }
        }
    }

    private static Tuple3<Integer, Integer, Integer> find35match(String seq5, String seq3) {
        final int longmm = 3;
        int max = 0;
        int b3 = 0;
        int b5 = 0;
        for (int i = 0; i < seq5.length() - 10; i++) {
            for (int j = 0; j < seq3.length() - 10; j++) {
                int nm = 0;
                int n = 0;
                while(n+j < seq3.length() && i+n < seq5.length()) {
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
                        && n - nm > 15
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
            adjCnt(iref, sc5v.variation);
            boolean rpflag = true;
            for (int i = 0; i < ins.length(); i++) {
                if (!isEquals(ref.get(bi + 1 + i), ins.charAt(i))) {
                    rpflag = false;
                    break;
                }
            }
            incCnt(cov, bi, sc5v.variation.cnt);
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
                    adjCnt(tvr, tv);
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
            if (rpflag && bams.length > 0 && ins.length() >= 5 && ins.length() < 75
                    && mref != null
                    && noPassingReads(chr, bi, bi + ins.length(), bams, conf)
                    && iref.cnt > 2 * mref.cnt) {
                adjCnt(iref, mref, mref);
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
            adjCnt(iref, sc3v.variation, lref);
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
            for (int ii = len - offset; ii < sc3v.seq.size(); ii++) {
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
                    adjCnt(vref, tv);
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
            if (rpflag && bams.length > 0 && ins.length() >= 5 && ins.length() < 75
                    && mref != null
                    && noPassingReads(chr, bi, bi + ins.length(), bams, conf)
                    && iref.cnt > 2 * mref.cnt) {

                adjCnt(iref, mref, mref);
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
                }
                m.add(seq.charAt(i + n));
                if (mm > maxmm) {
                    break;
                }
            }
            int mnt = m.size();
            if (mnt < 2) { // at least three different NT for overhang sequences, weeding out low complexity seq
                continue;
            }
            if ((mnt >= 3 && i + n >= seq.length() - 1 && i > 6 && mm / (double)i < 0.15)
                    || (mnt >= 2 && mm == 0 && i + n == seq.length() && n >= 20 && i >= 5)) {

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
                        return Tuple3.newTuple(p - 1 - extra.length(), insert.toString(), p - 1);
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
                        return Tuple3.newTuple(p + s, insert.toString(), p + s + extra.length());
                    } else if (i - mm > score) {
                        bi = p + s;
                        ins = insert.toString();
                        bi2 = p + s + extra.length();
                        score = i - mm;
                    }
                }
            }
        }

        return Tuple3.newTuple(bi, ins, bi2);
    }

    static final Comparator<Tuple2<Integer, Sclip>> COMP2 = new Comparator<Tuple2<Integer, Sclip>>() {
        @Override
        public int compare(Tuple2<Integer, Sclip> o1, Tuple2<Integer, Sclip> o2) {
            return Integer.compare(o1._2().variation.cnt, o2._2().variation.cnt);
        }
    };

    private static void realignlgdel(Map<Integer, Map<String, Variation>> hash,
            Map<Integer, Integer> cov,
            Map<Integer, Sclip> sclip5,
            Map<Integer, Sclip> sclip3,
            Map<Integer, Character> ref,
            String chr,
            Map<String, Integer> chrs,
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
            final int bp = findbp(seq, p - 5, ref, conf.indelsize, -1, chr, chrs);
            final int dellen = p - bp;
            if (bp == 0) {
                continue;
            }
            if (conf.y) {
                System.err.printf("Realignlgdel: $bp -$dellen \n", bp, dellen);
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
                incCnt(cov, tp, sc5v.variation.cnt);
            }
            adjCnt(tv, sc5v.variation);
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
                    adjCnt(getVariation(hash, bp, gt), sclip.variation, getVariation(hash, bp, ref.get(bp).toString()));
                } else {
                    adjCnt(getVariation(hash, bp, gt), sclip.variation);
                }
                for (int tp = bp; tp < bp + dellen; tp++) {
                    incCnt(cov, tp, sclip.variation.cnt);
                }
                for (int ip = bp + 1; ip < sc3p; ip++) {
                    Variation vv = getVariation(hash, ip, ref.get(dellen + ip).toString());
                    rmCnt(vv, sclip.variation);
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

            realigndel(hash, dels5, cov, sclip5, sclip3, ref, chr, chrs, bams, conf);
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
            int bp = findbp(seq, p + 5, ref, conf.indelsize, 1, chr, chrs);
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
                incCnt(cov, tp, sc3v.variation.cnt);
            }
            sc3v.variation.pmean += dellen * sc3v.variation.cnt;
            adjCnt(tv, sc3v.variation);
            sc3v.used = true;

            Map<Integer, Map<String, Integer>> dels5 = new HashMap<>();
            HashMap<String, Integer> map = new HashMap<>();
            map.put(gt, tv.cnt);
            dels5.put(bp, map);
            realigndel(hash, dels5, cov, sclip5, sclip3, ref, chr, chrs, bams, conf);
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
        if (vref.cnt < 0)
            vref.cnt = 0;
        if (vref.hicnt < 0)
            vref.hicnt = 0;
        if (vref.locnt < 0)
            vref.locnt = 0;
        if (vref.pmean < 0)
            vref.pmean = 0;
        if (vref.qmean < 0)
            vref.qmean = 0;
        if (vref.Qmean < 0)
            vref.Qmean = 0;
        if (vref.getDir(true) < 0)
            vref.addDir(true, -vref.getDir(true));
        if (vref.getDir(false) < 0)
            vref.addDir(false, -vref.getDir(false));
    }


    private static int findbp(String seq, int sp, Map<Integer, Character> ref, int dis, int dir, String chr, Map<String, Integer> chrs) {
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
                if (!Character.valueOf(seq.charAt(i)).equals(ref.get(sp + dir * n + dir * i))) {
                    mm++;
                }
                m.add(seq.charAt(i));
                if (mm > maxmm - n / 100) {
                    break;
                }
            }
            if (m.size() >= 3) {
                continue;
            }
            if (mm <= maxmm - n / 100 && i >= seq.length() - 2 && i >= 6 && mm / (double)i < 0.12) {
                int lbp = sp + dir * n - (dir < 0 ? dir : 0);
                if (mm == 0 && i == seq.length()) {
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
    private static final Pattern B_AMP_ATGC = Pattern.compile("^\\+([ATGC]+)");
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
            return Integer.compare(x1, x2);
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
            Matcher mtch = B_AMP_ATGC.matcher(vn);
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
            int wustart = p - 100 - vn.length() + 1;
            if (wustart <= 1) {
                wustart = 1;
            }
            String wupseq = joinRef(ref, wustart, p) + vn; // 5prime flanking seq
            wupseq = wupseq.replaceFirst("\\+", "").replaceFirst("&", "").replaceFirst("#", "");
            Integer tend = chrs.get(chr);
            int sanend = p + vn.length() + 100;
            if (tend != null && tend < sanend) {
                sanend = tend;
            }
            String sanpseq = vn + joinRef(ref, p + extra.length() + 1 + compm.length(), sanend); // 3prime flanking seq
            sanpseq = sanpseq.replaceFirst("\\+", "").replaceFirst("&", "").replaceFirst("#", "");
            Tuple3<List<Tuple3<String, Integer, Integer>>, List<Integer>, Integer> findMM3 = findMM3(ref, p + 1, sanpseq, insert.length() + compm.length(), sclip3); // mismatches, mismatch positions, 5 or 3 ends
            Tuple3<List<Tuple3<String, Integer, Integer>>, List<Integer>, Integer> findMM5 = findMM5(ref, p+extra.length()+compm.length(), wupseq, insert.length() + compm.length(), sclip5);

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

                adjCnt(vref, tv, lref);
                hash.get(mp).remove(mm);
                if (hash.get(mp).isEmpty()) {
                    hash.remove(mp);
                }
            }
            for (Integer sc5pp : sc5p) {
                Sclip tv = sclip5.get(sc5pp);
                if (tv != null && !tv.used) {
                    String seq = findconseq(tv);
                    if (!seq.isEmpty() && ismatch(seq, wupseq, -1)) {
                        if (conf.y) {
                            System.err.printf("ins5: %s %s $s %s %s %s\n", p, sc5pp, seq, wupseq, icnt, tv.variation.cnt);
                        }
                        if (sc5pp > p) {
                            incCnt(cov, p, tv.variation.cnt);
                        }
                        adjCnt(vref, tv.variation);
                        tv.used = true;
                    }
                }
            }
            for (Integer sc3pp : sc3p) {
                Sclip tv = sclip3.get(sc3pp);
                if (tv != null && !tv.used) {
                    String seq = findconseq(tv);
                    if (conf.y) {
                        System.err.printf("ins3: %s %s %s %s %s %s %s\n", p, sc3pp, seq, sanpseq, vn, icnt, tv.variation.cnt);
                    }
                    if (!seq.isEmpty() && ismatch(seq, substr(sanpseq, sc3pp - p - 1), 1)) {
                        if (conf.y) {
                            System.err.printf("ins3: %s %s %s %s %s %s Match\n", p, sc3pp, seq, vn, icnt, tv.variation.cnt);
                        }
                        if (sc3pp <= p) {
                            incCnt(cov, p, tv.variation.cnt);
                        }
                        Variation lref = null;
                        if (sc3pp <= p &&
                                hash.containsKey(p) &&
                                ref.containsKey(p) &&
                                hash.get(p).containsKey(ref.get(p).toString())) {

                            lref = hash.get(p).get(ref.get(p).toString());
                        }
                        adjCnt(vref, tv.variation, lref);
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
            Matcher mtch = ATGS_ATGS.matcher(vn);
            if (mtch.find()) {
                String tn = mtch.group(1);
                Variation tref = iHash.get(p).get(tn);
                if (tref != null) {
                    if (vref.cnt < tref.cnt) {
                        adjCnt(tref, vref, hash.get(p).get(ref.get(p).toString()));
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
            String[] bams,
            Configuration conf) throws IOException {

//        int longmm = 3; //Longest continued mismatches typical aligned at the end
        List<Object[]> tmp = fillTmp(dels5);

        for (Object[] objects : tmp) {
            Integer p = (Integer)objects[0];
            String vn = (String)objects[1];
            Integer dcnt = (Integer)objects[2];
            if (conf.y) {
                System.err.println(format("Realigndel for: %s %s %s", p, vn, dcnt));
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
            int sanend = p + 2 * dellen + 50;
            if (sanend > chrs.get(chr)) {
                sanend = chrs.get(chr);
            }
            String sanpseq = extrains + compm + extra + joinRef(ref, p + dellen + extra.length() + compm.length(), sanend); // 3' flanking seq
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
                    System.err.printf("Realigndel Adj: %s %s %s %s %s %s %s %s\n", mm, mp, me, nm3, nm5, p, tv.cnt, tv.qmean);
                }
                // Adjust ref cnt so that AF won't > 1
                if (mp > p && me == 5) {
                    incCnt(cov, p, tv.cnt);
                }
                Variation lref = (mp > p && me == 3) ? (hash.containsKey(p) && hash.get(p).containsKey(ref.get(p).toString()) ? hash.get(p).get(ref.get(p).toString()) : null) : null;
                adjCnt(vref, tv, lref);
                hash.get(mp).remove(mm);
                if (hash.get(mp).isEmpty() ) {
                    hash.remove(mp);
                }
            }

            for (Integer sc5pp : sc5p) {
                if (sclip5.containsKey(sc5pp) && !sclip5.get(sc5pp).used) {
                    Sclip tv = sclip5.get(sc5pp);
                    String seq = findconseq(tv);
                    if (conf.y) {
                        System.err.printf("Realigndel 5: %s %s %s %s %s %s %s\n",  sc5pp, seq, wupseq, tv.variation.cnt, dcnt, vn, p);
                    }
                    if (!seq.isEmpty() && ismatch(seq, wupseq, -1)) {
                        if (conf.y) {
                            System.err.printf("Realigndel 5: %s %s %s %s %s %s %s used\n",  sc5pp, seq, wupseq, tv.variation.cnt, dcnt, vn, p);
                        }
                        if (sc5pp > p) {
                            incCnt(cov, p, tv.variation.cnt);
                        }
                        adjCnt(vref, tv.variation);
                        sclip5.get(sc5pp).used = true;

                    }
                }
            }

            for (Integer sc3pp : sc3p) {
                if (sclip3.containsKey(sc3pp) && !sclip3.get(sc3pp).used) {
                    Sclip tv = sclip3.get(sc3pp);
                    String seq = findconseq(tv);
                    if(conf.y) {
                        System.err.printf("Realigndel 3: %s %s %s %s %s %s %s\n", sc3pp, seq, sanpseq, tv.variation.cnt, dcnt, vn, p);
                    }
                    if (!seq.isEmpty() && ismatch(seq, substr(sanpseq, sc3pp - p), 1)) {
                        if (conf.y) {
                            System.err.printf("Realigndel 3: %s %s %s %s %s %s %s used\n", sc3pp, seq, sanpseq, tv.variation.cnt, dcnt, vn, p);
                        }
                        if (sc3pp <=p ) {
                            incCnt(cov, p, tv.variation.cnt);
                        }
                        Variation lref = sc3pp <= p ? null : hash.get(p).get(ref.get(p).toString());
                        adjCnt(vref, tv.variation, lref);
                        sclip3.get(sc3pp).used = true;
                    }
                    if (bams != null
                            && sc3pp - p >= 5
                            && sc3pp - p < 75
                            && hash.containsKey(p) && hash.get(p).containsKey(ref.get(p).toString())
                            && noPassingReads(chr, p, sc3pp, bams, conf)
                            && vref.cnt > 2 * hash.get(p).get(ref.get(p).toString()).cnt) {
                        Variation hv = getVariation(hash, p, ref.get(p).toString());
                        adjCnt(vref, hv, hv);
                    }

                }
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
                Variation tref = hash.get(p).get(vn);
                if (tref != null) {
                    if (vref.cnt < tref.cnt) {
                        adjCnt(tref, vref);
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
            try(SamtoolsReader reader = new SamtoolsReader("view", bam, chr + ":" + s + "-" + e)) {
                String line;
                while ((line = reader.read()) != null) {
                    String[] a = line.split("\t");
                    int rs = toInt(a[3]);
                    if (a[5].matches("^(\\d+)M$")) {
                        int re = rs + toInt(substr(a[5], 0, -1));
                        if (re > e + 3 && rs < s -3) {
                            cnt++;
                        }
                    }

                }

            }
        }
        if (conf.y) {
            System.err.printf("Passing Read CNT: %s\n", cnt);
        }
        return cnt <= 0;
    }


    private static boolean ismatch(String seq1, String seq2, int dir) {
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
            Map<Integer, Integer> cov) {

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
                            adjCnt(vref, tref);
                            hash.get(p).remove(left);
                        }
                    }
                    if (hash.containsKey(p + i + 1) && hash.get(p + i + 1).containsKey(right)) {
                        Variation tref = getVariation(hash, p + i + 1, right);
                        if (tref.cnt < vref.cnt && tref.pmean / tref.cnt <= mnt.length() - i - 1) {
                            adjCnt(vref, tref);
                            incCnt(cov, p, tref.cnt);
                            hash.get(p + i + 1).remove(right);
                        }
                    }

                }

            }
        }

    }

    private static void adjCnt(Variation vref, Variation tv) {
        adjCnt(vref, tv, null);
    }

    private static void adjCnt(Variation vref, Variation tv, Variation ref) {
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
        if (ref.getDir(true) < 0)
            ref.addDir(true, -ref.getDir(true));
        if (ref.getDir(false) < 0)
            ref.addDir(false, -ref.getDir(false));
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

    private static void ampVardict(List<List<Region>> segs) {
        // TODO Auto-generated method stub

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
    public static List<Region> buildRegions(String region, final int numberNucleotideToExtend, final boolean zeroBased) {
        List<Region> segs = new ArrayList<>();
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
        segs.add(new Region(chr, start, end, gene));

        return segs;
    }
}
