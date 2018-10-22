package com.astrazeneca.vardict;

import com.astrazeneca.vardict.collection.Tuple;
import com.astrazeneca.vardict.data.ReferenceResource;
import com.astrazeneca.vardict.data.Region;
import com.astrazeneca.vardict.data.scopedata.GlobalReadOnlyScope;
import com.astrazeneca.vardict.modes.AbstractMode;
import com.astrazeneca.vardict.modes.AmpliconMode;
import com.astrazeneca.vardict.modes.SimpleMode;
import com.astrazeneca.vardict.modes.SomaticMode;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static com.astrazeneca.vardict.data.Patterns.SAMPLE_PATTERN;
import static com.astrazeneca.vardict.data.Patterns.SAMPLE_PATTERN2;
import static com.astrazeneca.vardict.data.scopedata.GlobalReadOnlyScope.instance;
import static com.astrazeneca.vardict.Utils.join;
import static com.astrazeneca.vardict.Utils.toInt;
import static com.astrazeneca.vardict.collection.Tuple.tuple;
import static com.astrazeneca.vardict.data.Patterns.INTEGER_ONLY;

public class VarDictLauncher {
    private List<List<Region>> segments;

    ReferenceResource referenceResource;

    public VarDictLauncher(ReferenceResource referenceResource) {
        this.referenceResource = referenceResource;
    }

    public void start(Configuration config) {
        initResources(config);

        final Configuration conf = instance().conf;
        final AbstractMode mode;

        if (conf.regionOfInterest != null || instance().ampliconBasedCalling == null) {
            mode = conf.bam.hasBam2() ?
                    new SomaticMode(segments, referenceResource) :
                    new SimpleMode(segments, referenceResource);
        } else {
            mode = new AmpliconMode(segments, referenceResource);
        }

        if (instance().conf.threads == 1)
            mode.notParallel();
        else
            mode.parallel();
    }

    private void initResources(Configuration conf) {
        try {
            Map<String, Integer> chrLengths = readChr(conf.bam.getBamX());
            Tuple.Tuple2<String, String> samples;
            if ((conf.regionOfInterest != null) && (conf.bam.hasBam2())) {
                samples = getSampleNamesSomatic(conf);
            } else {
                samples = getSampleNames(conf);
            }

            RegionBuilder builder = new RegionBuilder(chrLengths, conf);
            String ampliconBasedCalling = null;

            if (conf.regionOfInterest != null) {
                segments = builder.buildRegionFromConfiguration();
            } else {
                Tuple.Tuple3<String, Boolean, List<String>> tpl = readBedFile(conf);
                ampliconBasedCalling = tpl._1;
                Boolean zeroBased = tpl._2;
                List<String> segraw = tpl._3;

                if (ampliconBasedCalling != null) {
                    segments = builder.buildAmpRegions(segraw, zeroBased != null ? zeroBased : false);
                } else {
                    segments = builder.buildRegions(segraw, zeroBased);
                }
            }

            GlobalReadOnlyScope.init(conf, chrLengths, samples._1, samples._2, ampliconBasedCalling);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    /**
     * Method reads BED file line by line and checks if an amplicon mode sets in Configuration or by BED file.
     *
     * @param conf Vardict Configuration (contains amplicon based calling parameter)
     * @return tuple of amplicon parameters (distance to the edges and overlap fraction), zero based parameter
     * and list of lines from BED file
     * @throws IOException
     */
    private static Tuple.Tuple3<String, Boolean, List<String>> readBedFile(Configuration conf) throws IOException {
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

    /**
     * Read map of chromosome lengths
     * @param bam BAM file name
     * @return Map of chromosome lengths. Key - chromosome name, value - length
     * @throws IOException if BAM/SAM file can't be opened
     */
    public static Map<String, Integer> readChr(String bam) throws IOException {
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

    private static Tuple.Tuple2<String, String> getSampleNames(Configuration conf) {
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

    private static Tuple.Tuple2<String, String> getSampleNamesSomatic(Configuration conf) {
        Tuple.Tuple2<String, String> samples = getSampleNames(conf);
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
