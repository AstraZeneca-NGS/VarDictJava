package com.epam;

import com.epam.VarDict.BedRowFormat;

public class Configuration {
    boolean printHeader; //-h
    String delimiter; // -d
    String bed; // -b
    int numberNucleotideToExtend; // -x
    Boolean zeroBased; // -z,  default true if set -R
    String ampliconBasedCalling; //-a
    int columnForChromosome = -1;
    BedRowFormat bedRowFormat;
    String sampleNameRegexp; // -n
    String sampleName; //-N
    String fasta; // -G
    BamNames bam;
    Double downsampling;
    boolean chromosomeNameIsNumber; // -C
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
    boolean  moveIndelsTo3 = false; //-3
    String samfilter = "0x500"; //-F
    String regionOfInterest; //-R chr:start[-end]
    int readPosFilter = 5; // -P The read position filter, default 5
    double qratio = 1.5; //-o
    double mapq = 0; // -O The minimun mean mapping quality to be considered, default 0
    boolean doPileup = false; // -p Do pileup regarless the frequency
    double lofreq = 0.05d; // -V The lowest frequency in normal sample allowed for a putative somatic mutations, default to 0.05
    int threads;

    public boolean isColumnForChromosomeSet() {
        return columnForChromosome >= 0;
    }

    public boolean isDownsampling() {
        return downsampling != null;
    }

    public boolean hasMappingQuality() {
        return mappingQuality != null;
    }

    public boolean isZeroBasedDefined() {
        return zeroBased != null;
    }

    public static class BamNames {
        private final String[] bamNames;
        private final String[] bams;
        private final String bamRaw;

        public BamNames(String value) {
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

}