package com.epam;

import org.apache.commons.cli.*;

public class Main {

    // unless( getopts(
    // 'hHvtzypDCF3d:b:s:e:S:E:n:c:g:x:f:r:B:N:Q:m:T:q:Z:X:P:k:R:G:a:o:O:V:' ))
    // {
    // USAGE();
    // }
    // USAGE() if ( $opt_H );
    // USAGE() unless ( $opt_b );

    public static void main(String[] args) throws ParseException {
        Options options = buildOptions();
        CommandLineParser parser = new BasicParser();
        CommandLine cmd = parser.parse(options, args);
        if (cmd.getOptions().length == 0) {
            help(options);
        }
        new Main().run(cmd);
    }

    private void run(CommandLine cmd) {

    }

    private static Options buildOptions() {
        Options options = new Options();
        options.addOption("h", false, "Print a header row decribing columns");
        options.addOption("h", false, "Print a header row decribing columns");

        return options;
    }

    private static void help(Options options) {
        HelpFormatter formater = new HelpFormatter();
        formater.printHelp(128, "vardict", "VarDict is a variant calling program for SNV, MNV, indels (<120 bp), and complex variants.  It accepts any BAM format, either\n"+
    "from DNA-seq or RNA-seq.  There're several distinct features over other variant callers.  First, it can perform local\n"+
    "realignment over indels on the fly for more accurate allele frequencies of indels.  Second, it rescues softly clipped reads\n"+
    "to identify indels not present in the alignments or support existing indels.  Third, when given the PCR amplicon information,\n"+
    "it'll perform amplicon-based variant calling and filter out variants that show amplicon bias, a common false positive in PCR\n"+
    "based targeted deep sequencing.  Forth, it has very efficient memory management and memory usage is linear to the region of\n"+
    "interest, not the depth.  Five, it can handle ultra-deep sequencing and the performance is only linear to the depth.  It has\n"+
    "been tested on depth over 2M reads.  Finally, it has a build-in capability to perform paired sample analysis, intended for\n"+
    "somatic mutation identification, comparing DNA-seq and RNA-seq, or resistant vs sensitive in cancer research.  By default,\n"+
    "the region_info is an entry of refGene.txt from IGV, but can be any region or bed files.", options, "");
        System.exit(0);
    }

}
