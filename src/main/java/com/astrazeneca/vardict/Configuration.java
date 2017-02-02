package com.astrazeneca.vardict;

import htsjdk.samtools.ValidationStringency;

import com.astrazeneca.vardict.RegionBuilder.BedRowFormat;

public class Configuration {
    /**
     * Print a header row decribing columns
     */
    public boolean printHeader; //-h

    /**
     * The delimiter for split region_info
     */
    public String delimiter; // -d

    public String bed;

    /**
     * The number of nucleotide to extend for each segment
     */
    public int numberNucleotideToExtend; // -x

    /**
     * Indicate wehther is zero-based cooridates, as IGV does
     * When use -R option, it's set to false
     */
    public Boolean zeroBased; // -z,  default true if set -R

    /**
     * Indicate it's amplicon based calling.  Reads don't map to the amplicon will be skipped.  A read pair is considered belonging
     * the amplicon if the edges are less than int bp to the amplicon, and overlap fraction is at least float.  Default: 10:0.95
     */
    public String ampliconBasedCalling; //-a

    public int columnForChromosome = -1;

    public BedRowFormat bedRowFormat;

    /**
     * The regular expression to extract sample name from bam filenames.
     */
    public String sampleNameRegexp; // -descriptionString

    /**
     * The sample name to be used directly
     */
    public String sampleName; //-N

    /**
     * The the reference fasta
     */
    public String fasta; // -G

    /**
     * The indexed BAM file name(s)
     */
    public BamNames bam;

    /**
     * For downsampling fraction
     */
    public Double downsampling;

    public boolean chromosomeNameIsNumber; // -C

    /**
     * If set, reads with mapping quality less than INT will be filtered and ignored
     */
    public Integer mappingQuality;//-Q

    /**
     * Indicate to remove duplicated reads
     */
    public boolean removeDuplicatedReads; //-t

    /**
     * If set, reads with mismatches more than INT will be filtered and ignored
     */
    public int mismatch; //-m, default = 8

    public boolean y; //-y TODO ???

    /**
     * The phred score for a base to be considered a good call
     */
    public int goodq; // -q, default = 23

    public final int buffer = 200;

    /**
     * Extension of bp to look for mismatches after insersion or deletion
     */
    public int vext = 3; // -X, default 3

    /**
     * Trim bases after [INT] bases in the reads
     */
    public int trimBasesAfter = 0; // -T

    /**
     * Indicate whether to perform local realignment
     */
    public boolean performLocalRealignment; // -k, default false

    /**
     * The indel size
     */
    public int indelsize = 120; // -I, default 120

    /**
     * The cutoff to decide whether a position has read strand strandBiasFlag
     */
    public double bias = 0.05d;

    /**
     * The minimum reads for strandBiasFlag calculation
     */
    public int minb = 2; // -B, default 2.

    /**
     * The minimum # of variance reads
     */
    public int minr = 2; // -r, default 2 //TODO -p

    public boolean debug = false; // -D

    /**
     * The threshold for allele frequency
     */
    public double freq = 0.5; // -f and -p

    /**
     * Indicate to move indels to 3-prime if alternative alignment can be achieved
     */
    public boolean moveIndelsTo3 = false; //-3

    public String samfilter = "0x500"; //-F

    /**
     * chr:start[-end]
     */
    public String regionOfInterest; //-R

    /**
     * The read position filter
     */
    public int readPosFilter = 5; // -P default 5

    /**
     * The Qratio of (good_quality_reads)/(bad_quality_reads+0.5)
     */
    public double qratio = 1.5; //-o

    /**
     * The minimun mean mapping quality to be considered
     */
    public double mapq = 0; // -O  default 0

    /**
     * Do pileup regarless the frequency
     */
    public boolean doPileup = false; // -p

    /**
     * The lowest frequency in normal sample allowed for a putative somatic mutations
     */
    public double lofreq = 0.05d; // -V default to 0.05

    public final int lowqual = 10;

    public int minmatch = 0; // -M The minimum matches for a read to be considered

    public boolean outputSplicing = false; // -i Output splicing read counts

    public ValidationStringency validationStringency = ValidationStringency.LENIENT;

    /**
     * Threads count
     */
    public int threads;

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