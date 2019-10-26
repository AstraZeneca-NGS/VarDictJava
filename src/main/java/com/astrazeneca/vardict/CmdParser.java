package com.astrazeneca.vardict;

import com.astrazeneca.vardict.printers.PrinterType;
import htsjdk.samtools.ValidationStringency;
import org.apache.commons.cli.*;

import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

import static com.astrazeneca.vardict.RegionBuilder.DEFAULT_BED_ROW_FORMAT;

/**
 * Class to parse the parameters from the command line
 */
public class CmdParser {
    /**
     * Parses the array of command line parameters and fills configuration parameters.
     * @param args arguments from command line to be parsed
     * @return configuration with parameters from command line
     * @throws ParseException if parse can't be finished
     */
    public Configuration parseParams(String[] args) throws ParseException{
        Options options = buildOptions();

        CommandLineParser parser = new BasicParser();

        Configuration config = null;

        try{
            CommandLine cmd = parser.parse(options, args);
            if (cmd.getOptions().length == 0 || cmd.hasOption("H")) {
                help(options);
            }
            config = new CmdParser().parseCmd(cmd);
        } catch (MissingOptionException e) {
            List<?> missingOptions = e.getMissingOptions();
            System.err.print("Missing required option(s): ");
            for (Iterator<?> iterator = missingOptions.iterator(); iterator.hasNext(); ) {
                Object object = iterator.next();
                System.err.print(object);
                if (iterator.hasNext()) {
                    System.err.print(", ");
                }
            }
            System.err.println();
            help(options);
        }

        return config;
    }

    /**
     * For each parameter in CMD set the Configuration variable
     * @param cmd parsed CommandLine from apache CLI
     * @return configuration with parameters from command line
     * @throws ParseException
     */
    private Configuration parseCmd(CommandLine cmd) throws ParseException {
        Configuration config = new Configuration();

        String[] args = cmd.getArgs();

        if (args.length > 0) {
            config.bed = args[0];
        }
        config.printHeader = cmd.hasOption('h');
        config.chromosomeNameIsNumber = cmd.hasOption("C");
        config.debug = cmd.hasOption("D");
        config.removeDuplicatedReads = cmd.hasOption("t");
        config.moveIndelsTo3 = cmd.hasOption("3");
        config.samfilter = cmd.getOptionValue("F", "0x504");

        if (cmd.hasOption("z")) {
            config.zeroBased = 1 == getIntValue(cmd, "z", 1);
        }
        config.ampliconBasedCalling = cmd.getOptionValue("a");
        config.performLocalRealignment = 1 == getIntValue(cmd, "k", 1);

        config.fasta = setFastaFile(cmd);

        config.regionOfInterest = cmd.getOptionValue("R");
        config.delimiter = cmd.getOptionValue("d", "\t");
        config.sampleName = cmd.getOptionValue("N");

        if (cmd.hasOption('n')) {
            String regexp = cmd.getOptionValue('n');
            if (regexp.startsWith("/"))
                regexp = regexp.substring(1);
            if (regexp.endsWith("/"))
                regexp = regexp.substring(0, regexp.length() - 1);
            config.sampleNameRegexp = regexp;
        }
        config.bam = new Configuration.BamNames(cmd.getOptionValue("b"));

        int c_col = getColumnValue(cmd, "c", DEFAULT_BED_ROW_FORMAT.chrColumn);
        int S_col = getColumnValue(cmd, "S", DEFAULT_BED_ROW_FORMAT.startColumn);
        int E_col = getColumnValue(cmd, "E", DEFAULT_BED_ROW_FORMAT.endColumn);
        int s_col = getColumnValue(cmd, "s", DEFAULT_BED_ROW_FORMAT.thickStartColumn);
        int e_col = getColumnValue(cmd, "e", DEFAULT_BED_ROW_FORMAT.thickEndColumn);
        int g_col = getColumnValue(cmd, "g", DEFAULT_BED_ROW_FORMAT.geneColumn);

        if (cmd.hasOption("S") && !cmd.hasOption("s")) {
            s_col = S_col;
        }
        if (cmd.hasOption("E") && !cmd.hasOption("e")) {
            e_col = E_col;
        }

        config.bedRowFormat = new RegionBuilder.BedRowFormat(c_col, S_col, E_col, s_col, e_col, g_col);
        config.columnForChromosome = getColumnValue(cmd, "c", -1);

        config.numberNucleotideToExtend = getIntValue(cmd, "x", 0);
        config.freq = getDoubleValue(cmd, "f", 0.01d);
        config.minr = getIntValue(cmd, "r", 2);
        config.minBiasReads = getIntValue(cmd, "B", 2);
        if (cmd.hasOption("Q")) {
            config.mappingQuality = ((Number) cmd.getParsedOptionValue("Q")).intValue();
        }
        config.goodq = getDoubleValue(cmd, "q", 22.5);
        config.mismatch = getIntValue(cmd, "m", 8);
        config.trimBasesAfter = getIntValue(cmd, "T", 0);
        config.vext = getIntValue(cmd, "X", 2);
        config.readPosFilter = getIntValue(cmd, "P", 5);
        if (cmd.hasOption("Z")) {
            config.downsampling = getDoubleValue(cmd, "Z", 0);
        }
        config.qratio = getDoubleValue(cmd, "o", 1.5d);
        config.mapq = getDoubleValue(cmd, "O", 0);
        config.lofreq = getDoubleValue(cmd, "V", 0.05d);
        config.indelsize = getIntValue(cmd, "I", 50);


        if (cmd.hasOption("p")) {
            config.doPileup = true;
            config.freq = -1;
            config.minr = 0;
        }
        config.y = cmd.hasOption("y");
        config.outputSplicing = cmd.hasOption('i');

        if (cmd.hasOption('M')) {
            config.minmatch = ((Number) cmd.getParsedOptionValue("M")).intValue();
        }

        if (cmd.hasOption("VS")) {
            config.validationStringency = ValidationStringency.valueOf(cmd.getParsedOptionValue("VS").toString().toUpperCase());
        }

        config.includeNInTotalDepth = cmd.hasOption("K");
        config.chimeric = cmd.hasOption("chimeric");
        config.disableSV = cmd.hasOption("U");
        config.uniqueModeSecondInPairEnabled = cmd.hasOption("UN");
        config.uniqueModeAlignmentEnabled = cmd.hasOption("u");
        config.deleteDuplicateVariants = cmd.hasOption("deldupvar");

        config.INSSIZE = getIntValue(cmd, "w", 300);
        config.INSSTD = getIntValue(cmd, "W", 100);
        config.INSSTDAMT = getIntValue(cmd, "A", 4);
        config.SVMINLEN = getIntValue(cmd, "L", 1000);

        config.threads = Math.max(readThreadsCount(cmd), 1);
        config.referenceExtension = getIntValue(cmd, "Y", 1200);

        if (cmd.hasOption("adaptor")) {
            config.adaptor.addAll(Arrays.asList(cmd.getOptionValue("adaptor").split(",")));
        }

        config.fisher = cmd.hasOption("fisher");

        if (cmd.hasOption("DP")) {
            String defaultPrinter = cmd.getOptionValue("DP", PrinterType.OUT.name());
            switch(defaultPrinter) {
                case "OUT": config.printerType = PrinterType.OUT; break;
                case "ERR": config.printerType = PrinterType.ERR; break;
                default: config.printerType = PrinterType.OUT;
            }
        }

        config.crisprCuttingSite = getIntValue(cmd, "J", 0);
        config.crisprFilteringBp = getIntValue(cmd, "j", 0);

        config.monomerMsiFrequency = getDoubleValue(cmd, "mfreq", 0.25d);
        config.nonMonomerMsiFrequency = getDoubleValue(cmd, "nmfreq", 0.1d);
        return config;
    }

    /**
     * Reference can be set by using shorter command line option: hg19, hg38, mm10.
     * This paths are used only for AZ needs.
     * @param cmd parsed CommandLine from apache CLI
     * @return string contains path to fasta file
     */
    private String setFastaFile(CommandLine cmd) {
        String fasta = cmd.getOptionValue("G");
        if (fasta == null) {
            System.err.println("Reference file path wasn't set (option -G). Will be used the default fasta path.");
            fasta = Configuration.HG19;
        } else {
            switch (fasta) {
                case "hg19":
                    fasta = Configuration.HG19;
                    break;
                case "hg38":
                    fasta = Configuration.HG38;
                    break;
                case "mm10":
                    fasta = Configuration.MM10;
                    break;
            }
        }
        return fasta;
    }

    /**
     * Help information about options contains long and short option names and their descriptions
     * @return options from apache CLI
     */
    @SuppressWarnings("static-access")
    private Options buildOptions() {
        Options options = new Options();
        options.addOption("H", "?", false, "Print this help page");
        options.addOption("h", "header", false, "Print a header row describing columns");
        options.addOption("v", false, "VCF format output");
        options.addOption("i", "splice", false, "Output splicing read counts");
        options.addOption("p", false, "Do pileup regardless of the frequency");
        options.addOption("C", false, "Indicate the chromosome names are just numbers, such as 1, 2, not chr1, chr2 (deprecated)");
        options.addOption("D", "debug", false, "Debug mode.  Will print some error messages and append full genotype at the end.");
        options.addOption("t", "dedup", false, "Indicate to remove duplicated reads.  Only one pair with same start positions will be kept");
        options.addOption("3", false, "Indicate to move indels to 3-prime if alternative alignment can be achieved.");
        options.addOption("K", false, "Include Ns in the total depth calculation");
        options.addOption("u", false, "Indicate unique mode, which when mate pairs overlap, the overlapping part will be counted only once using forward read only.");
        options.addOption("UN", false, "Indicate unique mode, which when mate pairs overlap, the overlapping part will be counted only once using first read only.");
        options.addOption("chimeric", false, "Indicate to turn off chimeric reads filtering.");
        options.addOption("deldupvar", false, "Turn on deleting of duplicate variants. Variants in this mode are considered and outputted only if start position of variant is inside the region interest.");
        options.addOption("fisher", false, "Experimental feature: Changes R script (teststrandbias.R and testsomatic.) to Java implementation of Fisher exact test.");
        options.addOption("U", "nosv", false, "Turn off structural variant calling.");

        options.addOption(OptionBuilder.withArgName("bit")
                .hasArg(true)
                .withDescription("The hexical to filter reads using samtools. Default: 0x504 (filter 2nd alignments, unmapped reads and duplicates).  Use -F 0 to turn it off.")
                .withType(String.class)
                .isRequired(false)
                .create('F'));

        options.addOption(OptionBuilder.withArgName("0/1")
                .hasOptionalArgs(1)
                .withDescription("Indicate whether coordinates are zero-based, as IGV uses.  Default: 1 for BED file or amplicon BED file.\n"
                        + "Use 0 to turn it off. When using the -R option, it's set to 0")
                .withType(Number.class)
                .isRequired(false)
                .create('z'));

        options.addOption(OptionBuilder.withArgName("0/1")
                .hasOptionalArgs(1)
                .withDescription("Indicate whether to perform local realignment.  Default: 1.  Set to 0 to disable it.  For Ion or PacBio, 0 is recommended.")
                .withType(Number.class)
                .isRequired(false)
                .create('k'));

        options.addOption(OptionBuilder.withArgName("int:float")
                .hasArg(true)
                .withDescription("Indicate it's amplicon based calling.  Reads that don't map to the amplicon will be skipped.  A read pair is considered belonging "
                        + " to the amplicon if the edges are less than int bp to the amplicon, and overlap fraction is at least float.  Default: 10:0.95")
                .withType(Number.class)
                .isRequired(false)
                .withLongOpt("amplicon")
                .create('a'));

        options.addOption(OptionBuilder.withArgName("INT")
                .hasArg(true)
                .withDescription("The column for chromosome")
                .withType(Number.class)
                .isRequired(false)
                .create('c'));

        options.addOption(OptionBuilder.withArgName("Genome fasta")
                .hasArg(true)
                .withDescription("The reference fasta. Should be indexed (.fai).\n"
                        + "Defaults to: /ngs/reference_data/genomes/Hsapiens/hg19/seq/hg19.fa. " +
                        " Also short commands can be used to set path to: \n " +
                        "hg19 - /ngs/reference_data/genomes/Hsapiens/hg19/seq/hg19.fa\n" +
                        "hg38 - /ngs/reference_data/genomes/Hsapiens/hg38/seq/hg38.fa\n" +
                        "mm10 - /ngs/reference_data/genomes/Mmusculus/mm10/seq/mm10.fa" )
                .withType(String.class)
                .isRequired(false)
                .create('G'));

        options.addOption(OptionBuilder.withArgName("Region")
                .hasArg(true)
                .withDescription("The region of interest.  In the format of chr:start-end.  If end is omitted, then a single position.  No BED is needed.")
                .withType(String.class)
                .isRequired(false)
                .create('R'));

        options.addOption(OptionBuilder.withArgName("delemiter")
                .hasArg(true)
                .withDescription("The delimiter for split region_info, default to tab \"\\t\"")
                .withType(String.class)
                .isRequired(false)
                .create('d'));

        options.addOption(OptionBuilder.withArgName("regular_expression")
                .hasArg(true)
                .withDescription("The regular expression to extract sample name from BAM filenames.  Default to: /([^\\/\\._]+?)_[^\\/]*.bam/")
                .withType(String.class)
                .isRequired(false)
                .create('n'));

        options.addOption(OptionBuilder.withArgName("string")
                .hasArg(true)
                .withDescription("The sample name to be used directly.  Will overwrite -n option")
                .withType(String.class)
                .isRequired(false)
                .create('N'));

        options.addOption(OptionBuilder.withArgName("string")
                .hasArg(true)
                .withDescription("The indexed BAM file")
                .withType(String.class)
                .isRequired(true)
                .create('b'));

        options.addOption(OptionBuilder.withArgName("INT")
                .hasArg(true)
                .withDescription("The column for region start, e.g. gene start")
                .withType(Number.class)
                .isRequired(false)
                .create('S'));

        options.addOption(OptionBuilder.withArgName("INT")
                .hasArg(true)
                .withDescription("The column for region end, e.g. gene end")
                .withType(Number.class)
                .isRequired(false)
                .create('E'));

        options.addOption(OptionBuilder.withArgName("INT")
                .hasArg(true)
                .withDescription("The column for segment starts in the region, e.g. exon starts")
                .withType(Number.class)
                .isRequired(false)
                .create('s'));

        options.addOption(OptionBuilder.withArgName("INT")
                .hasArg(true)
                .withDescription("The column for segment ends in the region, e.g. exon ends")
                .withType(Number.class)
                .isRequired(false)
                .create('e'));

        options.addOption(OptionBuilder.withArgName("INT")
                .hasArg(true)
                .withDescription("The column for gene name, or segment annotation")
                .withType(Number.class)
                .isRequired(false)
                .create('g'));

        options.addOption(OptionBuilder.withArgName("INT")
                .hasArg(true)
                .withDescription("The number of nucleotide to extend for each segment, default: 0")
                .withType(Number.class)
                .isRequired(false)
                .create('x'));

        options.addOption(OptionBuilder.withArgName("double")
                .hasArg(true)
                .withDescription("The threshold for allele frequency, default: 0.01 or 1%")
                .withType(Number.class)
                .isRequired(false)
                .create('f'));

        options.addOption(OptionBuilder.withArgName("minimum reads")
                .hasArg(true)
                .withDescription("The minimum # of variant reads, default 2")
                .withType(Number.class)
                .isRequired(false)
                .create('r'));

        options.addOption(OptionBuilder.withArgName("INT")
                .hasArg(true)
                .withDescription("The minimum # of reads to determine strand bias, default 2")
                .withType(Number.class)
                .isRequired(false)
                .create('B'));

        options.addOption(OptionBuilder.withArgName("INT")
                .hasArg(true)
                .withDescription("If set, reads with mapping quality less than INT will be filtered and ignored")
                .withType(Number.class)
                .isRequired(false)
                .create('Q'));

        options.addOption(OptionBuilder.withArgName("INT")
                .hasArg(true)
                .withDescription("The phred score for a base to be considered a good call.  Default: 25 (for Illumina)\n"
                        + "For PGM, set it to ~15, as PGM tends to under estimate base quality.")
                .withType(Number.class)
                .isRequired(false)
                .create('q'));

        options.addOption(OptionBuilder.withArgName("INT")
                .hasArg(true)
                .withDescription("If set, reads with mismatches more than INT will be filtered and ignored.  Gaps are not counted as mismatches.\n"
                        + "Valid only for bowtie2/TopHat or BWA aln followed by sampe.  BWA mem is calculated as NM - Indels.  Default: 8,\n"
                        + "or reads with more than 8 mismatches will not be used.")
                .withType(Number.class)
                .isRequired(false)
                .create('m'));

        options.addOption(OptionBuilder.withArgName("INT")
                .hasArg(true)
                .withDescription("Trim bases after [INT] bases in the reads")
                .withType(Number.class)
                .isRequired(false)
                .withLongOpt("trim")
                .create('T'));

        options.addOption(OptionBuilder.withArgName("INT")
                .hasArg(true)
                .withDescription("Extension of bp to look for mismatches after insersion or deletion.  Default to 2 bp, or only calls when they're within 2 bp.")
                .withType(Number.class)
                .isRequired(false)
                .create('X'));

        options.addOption(OptionBuilder.withArgName("number")
                .hasArg(true)
                .withDescription("The read position filter.  If the mean variants position is less that specified, it's considered false positive.  Default: 5")
                .withType(Number.class)
                .isRequired(false)
                .create('P'));

        options.addOption(OptionBuilder.withArgName("double")
                .hasArg(true)
                .withDescription("For downsampling fraction.  e.g. 0.7 means roughly 70% downsampling.  Default: No downsampling.  Use with caution.  The\n"
                        + "downsampling will be random and non-reproducible.")
                .withType(Number.class)
                .isRequired(false)
                .withLongOpt("downsample")
                .create('Z'));

        options.addOption(OptionBuilder.withArgName("Qratio")
                .hasArg(true)
                .withDescription("The Qratio of (good_quality_reads)/(bad_quality_reads+0.5).  The quality is defined by -q option.  Default: 1.5")
                .withType(Number.class)
                .isRequired(false)
                .create('o'));

        options.addOption(OptionBuilder.withArgName("MapQ")
                .hasArg(true)
                .withDescription("The reads should have at least mean MapQ to be considered a valid variant.  Default: no filtering")
                .withType(Number.class)
                .isRequired(false)
                .create('O'));

        options.addOption(OptionBuilder.withArgName("freq")
                .hasArg(true)
                .withDescription("The lowest frequency in the normal sample allowed for a putative somatic mutation.  Defaults to 0.05")
                .withType(Number.class)
                .isRequired(false)
                .create('V'));

        options.addOption(OptionBuilder.withArgName("INT")
                .hasArg(true)
                .withDescription("The indel size.  Default: 50bp")
                .withType(Number.class)
                .isRequired(false)
                .create('I'));

        options.addOption(OptionBuilder
                .isRequired(false)
                .withLongOpt("verbose")
                .create('y'));

        options.addOption(OptionBuilder.withArgName("INT")
                .hasOptionalArg()
                .withDescription("Threads count.")
                .withType(Number.class)
                .isRequired(false)
                .create("th"));

        options.addOption(OptionBuilder.withArgName("INT")
                .hasArg(true)
                .withDescription("The minimum matches for a read to be considered. If, after soft-clipping, the matched bp is less than INT, then the "
                        + "read is discarded. It's meant for PCR based targeted sequencing where there's no insert and the matching is only the primers.\n"
                        + "Default: 0, or no filtering")
                .withType(Number.class)
                .isRequired(false)
                .create('M'));

        options.addOption(OptionBuilder.withArgName("STRICT | LENIENT | SILENT")
                .hasArg(true)
                .withDescription("How strict to be when reading a SAM or BAM.\n"
                        + "STRICT\t- throw an exception if something looks wrong.\n"
                        + "LENIENT\t- Emit warnings but keep going if possible.\n"
                        + "SILENT\t- Like LENIENT, only don't emit warning messages.\n"
                        + "Default: LENIENT")
                .withType(String.class)
                .isRequired(false)
                .create("VS"));

        options.addOption(OptionBuilder.withArgName("string")
                .hasArg(true)
                .withDescription("The printer type used for different outputs. Default: OUT (i.e. System.out).")
                .withType(String.class)
                .isRequired(false)
                .withLongOpt("default-printer")
                .create("DP"));

        options.addOption(OptionBuilder.withArgName("INT")
                .hasArg(true)
                .withDescription("The number of STD. A pair will be considered for DEL " +
                        "if INSERT > INSERT_SIZE + INSERT_STD_AMT * INSERT_STD.  Default: 4")
                .withType(Number.class)
                .isRequired(false)
                .create('A'));

        options.addOption(OptionBuilder.withArgName("INT")
                .hasArg(true)
                .withDescription("The insert size STD.  Used for SV calling.  Default: 100")
                .withType(Number.class)
                .isRequired(false)
                .withLongOpt("insert-std")
                .create('W'));

        options.addOption(OptionBuilder.withArgName("INT")
                .hasArg(true)
                .withDescription("The insert size.  Used for SV calling.  Default: 300")
                .withType(Number.class)
                .isRequired(false)
                .withLongOpt("insert-size")
                .create('w'));

        options.addOption(OptionBuilder.withArgName("INT")
                .hasArg(true)
                .withDescription("The minimum structural variant length to be presented using <DEL> <DUP> <INV> <INS>, etc. "
                        + "Default: 1000. Any indel, complex variants less than this will be spelled out with exact nucleotides.")
                .withType(Number.class)
                .isRequired(false)
                .create('L'));

        options.addOption(OptionBuilder.withArgName("INT")
                .hasArg(true)
                .withDescription("Extension of bp of reference to build lookup table. Default to 1200 bp." +
                        " Increase the number will slowdown the program. The main purpose is to call large indels with 1000 bp" +
                        " that can be missed by discordant mate pairs.")
                .withType(Number.class)
                .isRequired(false)
                .withLongOpt("ref-extension")
                .create('Y'));

        options.addOption(OptionBuilder.withArgName("string")
                .hasArg(true)
                .withDescription("Filter adaptor sequences so that they are not used in realignment. Multiple adaptors can be supplied " +
                        "by setting them with comma, like: --adaptor ACGTTGCTC,ACGGGGTCTC,ACGCGGCTAG .")
                .withType(String.class)
                .isRequired(false)
                .create("adaptor"));

        options.addOption(OptionBuilder.withArgName("CRISPR_cutting_site")
                .hasArg(true)
                .withDescription("The genomic position that CRISPR/Cas9 suppose to cut, typically 3bp from the PAM NGG site and within the guide.  For\n" +
                        "        CRISPR mode only.  It will adjust the variants (mostly In-Del) start and end sites to as close to this location as possible,\n" +
                        "        if there are alternatives. The option should only be used for CRISPR mode.")
                .withType(Number.class)
                .isRequired(false)
                .withLongOpt("crispr")
                .create('J'));

        options.addOption(OptionBuilder.withArgName("CRISPR_filtering_bp")
                .hasArg(true)
                .withDescription("In CRISPR mode, the minimum amount in bp that a read needs to overlap with cutting site.  If a read does not meet the criteria,\n" +
                        "        it will not be used for variant calling, since it is likely just a partially amplified PCR.  Default: not set, or no filtering")
                .withType(Number.class)
                .isRequired(false)
                .create('j'));

        options.addOption(OptionBuilder.withArgName("double")
                .hasArg(true)
                .withDescription("The variant frequency threshold to determine variant as good in case of monomer MSI. Default: 0.25")
                .withType(Number.class)
                .isRequired(false)
                .create("mfreq"));

        options.addOption(OptionBuilder.withArgName("double")
                .hasArg(true)
                .withDescription("The variant frequency threshold to determine variant as good in case of non-monomer MSI. Default: 0.1")
                .withType(Number.class)
                .isRequired(false)
                .create("nmfreq"));

        return options;
    }

    private int getIntValue(CommandLine cmd, String option, int defaultValue) throws ParseException {
            Object value = cmd.getParsedOptionValue(option);
            return  value == null ? defaultValue : ((Number)value).intValue();
    }

    private int getColumnValue(CommandLine cmd, String opt, int defaultValue) throws ParseException {
        Object value = cmd.getParsedOptionValue(opt);
        return  value == null ? defaultValue : ((Number)value).intValue() - 1;
    }

    private double getDoubleValue(CommandLine cmd, String opt, double defaultValue) throws ParseException {
        Object value = cmd.getParsedOptionValue(opt);
        return  value == null ? defaultValue : ((Number)value).doubleValue();
    }

    /**
     * Calculates possible count of threads to use. If -th set without value, it will be set to number of
     * available processors.
     * @param cmd parsed CommandLine from apache CLI
     * @return number of threads
     * @throws ParseException if option -th can't be read
     */
    private int readThreadsCount(CommandLine cmd) throws ParseException {
        int threads = 0;
        if (cmd.hasOption("th")) {
            Object value = cmd.getParsedOptionValue("th");
            if (value == null) {
                threads = Runtime.getRuntime().availableProcessors();
            } else {
                threads = ((Number)value).intValue();
            }
        }
        return threads;
    }

    private void help(Options options) {
        HelpFormatter formater = new HelpFormatter();
        formater.setOptionComparator(null);
        formater.printHelp(142, "vardict [-n name_reg] [-b bam] [-c chr] [-S start] [-E end] [-s seg_starts] [-e seg_ends] "
                        + "[-x #_nu] [-g gene] [-f freq] [-r #_reads] [-B #_reads] region_info",
                "VarDict is a variant calling program for SNV, MNV, indels (<50 bp), and complex variants.  It accepts any BAM format, either\n"+
                        "from DNA-seq or RNA-seq.  There are several distinct features over other variant callers.  First, it can perform local\n"+
                        "realignment over indels on the fly for more accurate allele frequencies of indels.  Second, it rescues softly clipped reads\n"+
                        "to identify indels not present in the alignments or support existing indels.  Third, when given the PCR amplicon information,\n"+
                        "it will perform amplicon-based variant calling and filter out variants that show amplicon bias, a common false positive in PCR\n"+
                        "based targeted deep sequencing.  Forth, it has very efficient memory management and memory usage is linear to the region of\n"+
                        "interest, not the depth.  Five, it can handle ultra-deep sequencing and the performance is only linear to the depth.  It has\n"+
                        "been tested on depth over 2M reads.  Finally, it has a build-in capability to perform paired sample analysis, intended for\n"+
                        "somatic mutation identification, comparing DNA-seq and RNA-seq, or resistant vs sensitive in cancer research.  By default,\n"+
                        "the region_info is an entry of refGene.txt from IGV, but can be any region or bed files.\nOptions:",
                options, "AUTHOR\n"
                        + ".       Written by Zhongwu Lai, AstraZeneca, Boston, USA\n\n"
                        + "REPORTING BUGS\n"
                        + ".       Report bugs to zhongwu@yahoo.com\n\n"
                        + "COPYRIGHT\n"
                        + ".       This is free software: you are free to change and redistribute it.  There is NO WARRANTY, to the extent permitted by law.");

        System.exit(0);
    }
}
