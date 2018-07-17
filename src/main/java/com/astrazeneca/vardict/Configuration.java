package com.astrazeneca.vardict;

import htsjdk.samtools.ValidationStringency;

import com.astrazeneca.vardict.VarDict.BedRowFormat;

public class Configuration {
    /**
     * Print a header row describing columns
     */
    boolean printHeader; //-h
    /**
     * The delimiter for split region_info
     */
    String delimiter; // -d
    /**
     * Path to bed file with regions
     */
    String bed;
    /**
     * The number of nucleotide to extend for each segment
     */
    int numberNucleotideToExtend; // -x
    /**
     * Indicate whether is zero-based coordinates, as IGV does
     * When use -R option, it is set to false TODO: is it so? In perl it doesn't changed
     */
    Boolean zeroBased; // -z
    /**
     * Indicate it's amplicon based calling.  Reads don't map to the amplicon will be skipped.
     * A read pair is considered belonging the amplicon if the edges are less than int bp to the amplicon,
     * and overlap fraction is at least float. Default: 10:0.95
     */
    String ampliconBasedCalling; //-a

    int columnForChromosome = -1; //-c

    BedRowFormat bedRowFormat;
    /**
     * The regular expression to extract sample name from bam filenames.
     */
    String sampleNameRegexp; // -n
    /**
     * The sample name to be used directly
     */
    String sampleName; //-N
    /**
     * The reference fasta
     */
    String fasta; // -G
    /**
     * The indexed BAM file name(s)
     */
    BamNames bam; //-b
    /**
     * For downsampling fraction
     */
    Double downsampling; //-Z

    boolean chromosomeNameIsNumber; // -C
    /**
     * If set, reads with mapping quality less than INT will be filtered and ignored
     */
    Integer mappingQuality;//-Q
    /**
     * Indicate to remove duplicated reads
     */
    boolean removeDuplicatedReads; //-t
    /**
     * If set, reads with mismatches more than INT will be filtered and ignored
     */
    int mismatch; //-m, default = 8
    /**
     * Verbose mode. Will output variant calling process.
     */
    boolean y; //-y
    /**
     * The phred score for a base to be considered a good call
     */
    double goodq; // -q, default = 22.5
    /**
     * Extension of bp to look for mismatches after insersion or deletion
     */
    int vext = 3; // -X, default 3
    /**
     * Trim bases after [INT] bases in the reads
     */
    int trimBasesAfter = 0; // -T
    /**
     * Indicate whether to perform local realignment
     */
    boolean performLocalRealignment; // -k, default false
    /**
     * The indel size
     */
    int indelsize = 50; // -I, default 50
    /**
     * The cutoff to decide whether a position has read strand bias
     */
    double bias = 0.05d;
    /**
     * The minimum reads for bias calculation. Default: 2
     */
    int minb = 2; // -B
    /**
     * The minimum # of variance reads. Default: 2. If -p, it is set to 0.
     */
    int minr = 2; // -r

    /**
     * Debug mode. Will print some error messages and append full genotype at the end.
     */
    boolean debug = false; // -D
    /**
     * The threshold for allele frequency. If -p it is set to -1.
     */
    double freq = 0.01; // -f
    /**
     * Indicate to move indels to 3-prime if alternative alignment can be achieved.
     */
    boolean moveIndelsTo3 = false; //-3

    /**
     * The hexical to filter reads.
     */
    String samfilter = "0x500"; //-F
    /**
     * chr:start[-end]. If end is omitted, then a single position.
     */
    String regionOfInterest; //-R
    /**
     * The read position filter. Default: 5
     */
    int readPosFilter = 5; // -P
    /**
     * The Qratio of (good_quality_reads)/(bad_quality_reads+0.5)
     */
    double qratio = 1.5; // -o
    /**
     * The minimum mean mapping quality to be considered. Default: 0.
     */
    double mapq = 0; // -O
    /**
     * Do pileup regardless the frequency
     */
    boolean doPileup = false; // -p
    /**
     * The lowest allele frequency in normal sample allowed for a putative somatic mutations. Default: 0.05.
     */
    double lofreq = 0.05d; // -V

    /**
     * Any base with quality <=10 will be consider low quality in soft-clipped seq and extension will stop.
     */
    final int lowqual = 10;

    /**
     * The minimum matches for a read to be considered
     */
    int minmatch = 0; // -M
    /**
     *  Output splicing read counts
     */
    boolean outputSplicing = false; // -i

    /**
     * How strict to be when reading a SAM or BAM.
     */
    ValidationStringency validationStringency = ValidationStringency.LENIENT; // -VS
    
    /**
     * Include Ns in the total depth calculation.
     */
    boolean includeNInTotalDepth = false; // -K

    /**
     * Indicate unique mode, which when mate pairs overlap,
     * the overlapping part will be counted only once using forward read only.
     */
    boolean uniqueModeAlignmentEnabled = false; // -u

    /**
     * Indicate unique mode, which when mate pairs overlap,
     * the overlapping part will be counted only once using first read only.
     */
    boolean uniqueModeSecondInPairEnabled = false; // -UN

    /**
     * Threads count to use in multithreading mode
     */
    int threads; //-th

    /**
     * The larger seed size
     */
    int seed1 = 17;

    /**
     * The smaller seed size
     */
    int seed2 = 12;

    /**
     * The adaptor size
     */
    int adseed = 6;
    /**
     *
     * Indicate to turn off chimeric reads filtering.  Chimeric reads are artifacts from library construction,
     * where a read can be split into two segments, each will be aligned within 1-2 read length distance,
     * but in opposite direction.
     */
    public boolean chimeric = false; // --chimeric

    /**
     * Turn off structural variant calling when set to true
     */
    boolean disableSV = false; //-U

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
