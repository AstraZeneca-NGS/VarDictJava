package com.astrazeneca.vardict;

import com.astrazeneca.vardict.printers.PrinterType;
import htsjdk.samtools.ValidationStringency;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

public class Configuration {
    public static final String HG19 = "/ngs/reference_data/genomes/Hsapiens/hg19/seq/hg19.fa";
    public static final String HG38 = "/ngs/reference_data/genomes/Hsapiens/hg38/seq/hg38.fa";
    public static final String MM10 = "/ngs/reference_data/genomes/Mmusculus/mm10/seq/mm10.fa";

    /**
     * Print a header row describing columns
     */
    public boolean printHeader; //-h
    /**
     * The delimiter for split region_info
     */
    public String delimiter; // -d
    /**
     * Path to bed file with regions
     */
    public String bed;
    /**
     * The number of nucleotide to extend for each segment
     */
    public int numberNucleotideToExtend; // -x
    /**
     * Indicate whether is zero-based coordinates, as IGV does
     * When use -R option, it is set to false
     */
    public Boolean zeroBased; // -z
    /**
     * Indicate it's amplicon based calling.  Reads don't map to the amplicon will be skipped.
     * A read pair is considered belonging the amplicon if the edges are less than int bp to the amplicon,
     * and overlap fraction is at least float. Default: 10:0.95
     */
    public String ampliconBasedCalling; //-a

    public int columnForChromosome = -1; //-c

    public RegionBuilder.BedRowFormat bedRowFormat;
    /**
     * The regular expression to extract sample name from bam filenames.
     */
    public String sampleNameRegexp; // -n
    /**
     * The sample name to be used directly
     */
    public String sampleName; //-N
    /**
     * The reference fasta
     */
    public String fasta; // -G
    /**
     * The indexed BAM file name(s)
     */
    public BamNames bam; //-b
    /**
     * For downsampling fraction
     */
    public Double downsampling; //-Z

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
    /**
     * Verbose mode. Will output variant calling process.
     */
    public boolean y; //-y
    /**
     * The phred score for a base to be considered a good call
     */
    public double goodq; // -q, default = 22.5
    /**
     * Extension of bp to look for mismatches after insertion or deletion
     */
    public int vext = 2; // -X, default 2
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
    public int indelsize = 50; // -I, default 50
    /**
     * The cutoff to decide whether a position has read strand bias
     */
    public double bias = 0.05d;
    /**
     * The minimum reads for bias calculation. Default: 2 $minb
     */
    public int minBiasReads = 2; // -B
    /**
     * The minimum # of variance reads. Default: 2. If -p, it is set to 0.
     */
    public int minr = 2; // -r

    /**
     * Debug mode. Will print some error messages and append full genotype at the end.
     */
    public boolean debug = false; // -D
    /**
     * The threshold for allele frequency. If -p it is set to -1.
     */
    public double freq = 0.01; // -f
    /**
     * Indicate to move indels to 3-prime if alternative alignment can be achieved.
     */
    public boolean moveIndelsTo3 = false; //-3

    /**
     * The hexical to filter reads.
     */
    public String samfilter = "0x504"; //-F
    /**
     * chr:start[-end]. If end is omitted, then a single position.
     */
    public String regionOfInterest; //-R
    /**
     * The read position filter. Default: 5
     */
    public int readPosFilter = 5; // -P
    /**
     * The Qratio of (good_quality_reads)/(bad_quality_reads+0.5)
     */
    public double qratio = 1.5; // -o
    /**
     * The minimum mean mapping quality to be considered. Default: 0.
     */
    public double mapq = 0; // -O
    /**
     * Do pileup regardless the frequency
     */
    public boolean doPileup = false; // -p
    /**
     * The lowest allele frequency in normal sample allowed for a putative somatic mutations. Default: 0.05.
     */
    public double lofreq = 0.05d; // -V

    /**
     * Any base with quality &lt;=10 will be consider low quality in soft-clipped seq and extension will stop.
     */
    public static final int LOWQUAL = 10;

    /**
     * The minimum matches for a read to be considered
     */
    public int minmatch = 0; // -M
    /**
     *  Output splicing read counts
     */
    public boolean outputSplicing = false; // -i

    /**
     * How strict to be when reading a SAM or BAM.
     */
    public ValidationStringency validationStringency = ValidationStringency.LENIENT; // -VS
    
    /**
     * Include Ns in the total depth calculation.
     */
    public boolean includeNInTotalDepth = false; // -K

    /**
     * Indicate unique mode, which when mate pairs overlap,
     * the overlapping part will be counted only once using forward read only.
     */
    public boolean uniqueModeAlignmentEnabled = false; // -u

    /**
     * Indicate unique mode, which when mate pairs overlap,
     * the overlapping part will be counted only once using first read only.
     */
    public boolean uniqueModeSecondInPairEnabled = false; // -UN

    /**
     * Threads count to use in multithreading mode
     */
    public int threads; //-th

    /**
     * The larger seed size
     */
    public static final int SEED_1 = 17;

    /**
     * The smaller seed size
     */
    public static final int SEED_2 = 12;

    /**
     * The adaptor size
     */
    public static final int ADSEED = 6;
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
    public boolean disableSV = false; //-U

    /**
     * Turn on deleting of duplicate variants that can appear due to VarDict linear work on regions.
     */
    public boolean deleteDuplicateVariants = false;

    /**
     * Applying Fisher exact test on forward and reverse counts of variant.
     */
    public boolean fisher = false;

    /**
     * The minimum distance between two SV clusters in term of read length
     */
    public static final double MINSVCDIST = 1.5;
    /**
     * The minimum position in read for mapping quality in SV analysis
     */
    public static final int MINMAPBASE = 15;
    /**
     * The minimum distance between start position and end of SV structure in inter-chr translocation
     */
    public static final int MINSVPOS = 25;
    /**
     * Mean Insert size
     */
    public int INSSIZE = 300; //-w
    /**
     * Insert std
     */
    public int INSSTD = 100; //-W
    /**
     * Insert std amount
     */
    public int INSSTDAMT = 4; //-A
    /**
     * The minimum structural variant length to be presented using &lt;DEL&gt; &lt;DUP&gt; &lt;INV&gt; &lt;INS&gt;, etc.
     */
    public int SVMINLEN = 1000; //-L
    /**
     * Max Structure variant size to be called in realignment step
     */
    public static final int SVMAXLEN = 150000;
    /**
     * the flanking sequence length for SV
     */
    public static final int SVFLANK = 50;
    /**
     * The minimum mapping quality when structural variant is only supported by discordant pairs
     */
    public static final int DISCPAIRQUAL = 35;

    public static final int EXTENSION = 5000;

    public static final String DEFAULT_AMPLICON_PARAMETERS = "10:0.95";

    /**
     * Default reference extension $REFEXT
     */
    public int referenceExtension = 1200;

    /**
     * Default printer for variants - system.out
     */
    public PrinterType printerType = PrinterType.OUT;

    /**
     * Exception counter
     * */
    public AtomicInteger exceptionCounter = new AtomicInteger(0);

    /**
     * Maximum of exception to continue work
     */
    public static int MAX_EXCEPTION_COUNT = 10;

    /**
     * List of adaptor sequences
     */
    public List<String> adaptor = new ArrayList<>();

    /**
     * In CRISPR mode, the minimum amount in bp that a read needs to overlap with cutting site.
     */
    public int crisprFilteringBp = 0;
    /**
     * The genomic position that CRISPR/Cas9 suppose to cut, typically 3bp from the PAM NGG site and within the guide.
     */
    public int crisprCuttingSite = 0;

    /**
     * The variant frequency threshold to determine variant as good in case of monomer MSI
     */
    public double monomerMsiFrequency = 0.25d;  // -mfreq
    /**
     * The variant frequency threshold to determine variant as good in case of non-monomer MSI
     */
    public double nonMonomerMsiFrequency = 0.1d;  // -nmfreq

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
