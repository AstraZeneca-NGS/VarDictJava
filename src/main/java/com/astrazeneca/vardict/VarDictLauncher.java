package com.astrazeneca.vardict;

import com.astrazeneca.GlobalReadOnlyScope;
import com.astrazeneca.utils.Tuple.Tuple2;
import com.astrazeneca.utils.Tuple.Tuple3;
import com.astrazeneca.vardict.modes.AbstractMode;
import com.astrazeneca.vardict.modes.AmpliconMode;
import com.astrazeneca.vardict.modes.SimpleMode;
import com.astrazeneca.vardict.modes.SomaticMode;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static com.astrazeneca.GlobalReadOnlyScope.instance;
import static com.astrazeneca.utils.Tuple.tuple;
import static com.astrazeneca.utils.Utils.join;
import static com.astrazeneca.utils.Utils.toInt;

public class VarDictLauncher {

    //perl version: 15 str
    private final static Pattern SAMPLE_PATTERN = Pattern.compile("([^\\/\\._]+).sorted[^\\/]*.bam");

    //perl version: 21 str
    private final static Pattern SAMPLE_PATTERN2 = Pattern.compile("([^\\/]+)[_\\.][^\\/]*bam");

    private List<List<Region>> segments;

    public VarDictLauncher(Configuration conf) {
        initResources(conf);
    }

    public void start() throws IOException {
        printHeader();

        final Configuration conf = instance().conf;
        final AbstractMode mode;

        if (conf.regionOfInterest != null || instance().ampliconBasedCalling == null) {
            mode = conf.bam.hasBam2() ?
                    new SomaticMode(segments) :
                    new SimpleMode(segments);
        } else {
            mode = new AmpliconMode(segments);
        }

        if (instance().conf.threads == 1)
            mode.notParallel();
        else
            mode.parallel();
    }

    public void initResources(Configuration conf) {
        try {
            Tuple3<String, Boolean, List<String>> tpl = readBedFile(conf);
            String ampliconBasedCalling = tpl._1;
            Boolean zeroBased = tpl._2;
            List<String> segraw = tpl._3;

            initGlobalScope(conf, ampliconBasedCalling);
            RegionBuilder builder = new RegionBuilder(instance().chrLengths, conf);

            if (conf.regionOfInterest != null) {
                segments = builder.buildRegionFromConfiguration();
            } else {
                if (ampliconBasedCalling != null) {
                    segments = builder.buildAmpRegions(segraw, zeroBased != null ? zeroBased : false);
                } else {
                    segments = builder.buildRegions(segraw, zeroBased);
                }
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    private static void initGlobalScope(Configuration conf, String ampliconBasedCalling) throws IOException {
        Map<String, Integer> chrLengths = readChr(conf.bam.getBamX());
        Tuple2<String, String> samples;
        if ((conf.regionOfInterest != null) && (conf.bam.hasBam2())) {
            samples = getSampleNamesSomatic(conf);

        } else {
            samples = getSampleNames(conf);
        }

        //TODO Does it necessary? Porting BUG?? (Perl version uses config only)
        if (conf.regionOfInterest != null) {
            ampliconBasedCalling = conf.ampliconBasedCalling;
        }

        GlobalReadOnlyScope.init(conf, chrLengths, samples._1, samples._2, ampliconBasedCalling);
    }

    private void printHeader() {
        if (instance().conf.printHeader) {
            System.out.println(join("\t",
                    "Sample", "Gene", "Chr", "Start", "End", "Ref", "Alt", "Depth", "AltDepth", "RefFwdReads",
                    "RefRevReads", "AltFwdReads", "AltRevReads", "Genotype", "AF", "Bias", "PMean", "PStd",
                    "QMean", "QStd", "5pFlankSeq", "3pFlankSeq"));
        }
    }

    private static Tuple3<String, Boolean, List<String>> readBedFile(Configuration conf) throws IOException {
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

    //init GlobalScope
    private static Tuple2<String, String> getSampleNames(Configuration conf) {
        String sample = null;
        String samplem = "";

        if (conf.sampleName != null) {
            sample = conf.sampleName;
        } else {
            Pattern rn;
            if (conf.sampleNameRegexp == null || conf.sampleNameRegexp.isEmpty()) {
                rn = SAMPLE_PATTERN;
            } else {
                rn = Pattern.compile(conf.sampleNameRegexp);
            }

            Matcher matcher = rn.matcher(conf.bam.getBamRaw());
            if (matcher.find()) {
                sample = matcher.group(1);
            }
        }

        if (sample == null) {
            Matcher matcher = SAMPLE_PATTERN2.matcher(conf.bam.getBamRaw());
            if (matcher.find()) {
                sample = matcher.group(1);
            }
        }

        return tuple(sample, samplem);
    }

    //Init GlobalScope
    private static Tuple2<String, String> getSampleNamesSomatic(Configuration conf) {
        Tuple2<String, String> samples = getSampleNames(conf);
        String sample = samples._1;
        String samplem = samples._2;

        if (conf.sampleNameRegexp != null) {
            Pattern rn = Pattern.compile(conf.sampleNameRegexp);
            Matcher m = rn.matcher(conf.bam.getBam1());
            if (m.find()) {
                sample = m.group(1);
            }
            m = rn.matcher(conf.bam.getBam2());
            if (m.find()) {
                samplem = m.group(1);
            }
        } else {
            if (conf.sampleName != null) {
                String[] split = conf.sampleName.split("\\|");
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
}
