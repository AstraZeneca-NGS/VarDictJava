#VarDictJava

##Introduction
VarDictJava is a variant discovery program written in Java and Perl. It is a partial Java port of [VarDict variant caller](https://github.com/AstraZeneca-NGS/VarDict). 

The original Perl VarDict is a sensitive variant caller for both single and paired sample variant calling from BAM files. VarDict implements several novel features such as amplicon bias aware variant calling from targeted
sequencing experiments, rescue of long indels by realigning bwa soft clipped reads and better scalability
than many other Java based variant callers. The Java port is around 10x faster than the original Perl implementation.

Please cite VarDict:

Lai Z, Markovets A, Ahdesmaki M, Chapman B, Hofmann O, McEwen R, Johnson J, Dougherty B, Barrett JC, and Dry JR.  VarDict: a novel and versatile variant caller for next-generation sequencing in cancer research. Nucleic Acids Res. 2016, pii: gkw227.

The link to is article can be accessed through: http://nar.oxfordjournals.org/cgi/content/full/gkw227?ijkey=Tk8eKQcYwNlQRNU&keytype=ref

Original coded by Zhongwu Lai 2014.

VarDictJava can run in single sample (see Single sample mode section), paired sample (see Paired variant calling section), or amplicon bias aware modes. As input, VarDictJava takes reference genomes in FASTA format, aligned reads in BAM format, and target regions in BED format.

##Requirements
1. JDK 1.8 or later
2. R language (uses /usr/bin/env R)
3. Perl (uses /usr/bin/env perl)
4. Internet connection to download dependencies using gradle.

##Getting started
###Getting source code
The VarDictJava source code is located at [https://github.com/AstraZeneca-NGS/VarDictJava](https://github.com/AstraZeneca-NGS/VarDictJava).

To load the project, execute the following command:

```
git clone --recursive https://github.com/AstraZeneca-NGS/VarDictJava.git
```

Note that the original VardDict project is placed in this repository as a submodule and its contents can be found in the sub-directory VarDict in VarDictJava working folder. So when you use `teststrandbias.R` and `var2vcf_valid.pl.` (see details and examples below), you have to add prefix VarDict: `VarDict/teststrandbias.R` and `VarDict/var2vcf_valid.pl.`

###Compiling
The project uses [Gradle](http://gradle.org/) and already includes a gradlew script.

To build the project, in the root folder of the project, run the following command:

```
./gradlew clean installApp 
```

To generate Javadoc, in the build/docs/javadoc folder, run the following command:

```
./gradlew clean javadoc
```

###Single sample mode

To run VarDictJava in single sample mode, use a BAM file specified without the `|` symbol and perform Steps 3 and 4 (see the Program workflow section) using `teststrandbias.R` and `var2vcf_valid.pl.`
The following is an example command to run in single sample mode:
  
```
AF_THR="0.01" # minimum allele frequency
<path_to_vardict_folder>/build/install/VarDict/bin/VarDict -G /path/to/hg19.fa -f $AF_THR -N sample_name -b /path/to/my.bam -z -c 1 -S 2 -E 3 -g 4 /path/to/my.bed | VarDict/teststrandbias.R | VarDict/var2vcf_valid.pl -N sample_name -E -f $AF_THR
```

VarDictJava can also be invoked without a BED file if the region is specified in the command line with `-R` option.
The following is an example command to run VarDictJava for a region (chromosome 7, position from 55270300 to 55270348, EGFR gene) with `-R` option:

```
<path_to_vardict_folder>/build/install/VarDict/bin/VarDict  -G /path/to/hg19.fa -f 0.001 -N sample_name -b /path/to/sample.bam  -z -R  chr7:55270300-55270348:EGFR | VarDict/teststrandbias.R | VarDict/var2vcf_valid.pl -N sample_name -E -f 0.001 >vars.vcf
```

In single sample mode, output columns contain a description and statistical info for variants in the single sample. See section Output Columns for list of columns in the output. 

###Paired variant calling

To run paired variant calling, use BAM files specified as `BAM1|BAM2` and perform Steps 3 and 4 (see the Program Workflow section) using `testsomatic.R` and `var2vcf_paired.pl`.

In this mode, the number of statistics columns in the output is doubled: one set of columns is for the first sample, the other - for second sample.

The following is an example command to run in paired mode:
```
AF_THR="0.01" # minimum allele frequency
<path_to_vardict_folder>/build/install/VarDict/bin/VarDict -G /path/to/hg19.fa -f $AF_THR -N tumor_sample_name -b "/path/to/tumor.bam|/path/to/normal.bam" -z -F -c 1 -S 2 -E 3 -g 4 /path/to/my.bed | VarDict/testsomatic.R | VarDict/var2vcf_paired.pl -N "tumor_sample_name|normal_sample_name" -f $AF_THR
```
###Running Tests
#### Integration testing
The list of integration test cases is stored in files in `testdata/intergationtestcases` directory.
To run all integration tests, the command is:
```$xslt
./gradlew test --tests com.astrazeneca.vardict.integrationtests.IntegrationTest 
```
##### User extension of testcases
Each file in `testdata/intergationtestcases` directory represents a test case with input data and expected output

**1.** Create a txt file in `testdata/intergationtestcases` folder. 

The file contains testcase input (of format described in [Test cases file format](Readme.md#CSVhead)) in the first line and expected output in the remaining file part.

**2.** Extend or create [thin-FASTA](Readme.md#thFastahead) in `testdata/fastas` folder.
 
**3.** Run tests.
##### <a name="CSVhead"></a>Test cases file format 
Each input file represents one test case input description. In the input file the first line consists of the following fields separated by `,` symbol:

*Required fields:*
- test case name
- reference name
- bam file name
- chromosome name
- start of region
- end of region

*Optional fields:*
- start of region with amplicon case
- end of region with amplicon case

*Parameters field:*
- the last filed can be any other command line parameters string

Example of first line of input file:
```
Amplicon,hg19.fa,Colo829-18_S3-sort.bam,chr1,933866,934466,933866,934466,-a 10:0.95 -D
Somatic,hg19.fa,Colo829-18_S3-sort.bam|Colo829-19_S4-sort.bam,chr1,755917,756517
Simple,hg19.fa,Colo829-18_S3-sort.bam,chr1,9922,10122,-p
```

##### <a name="thFastahead"></a>Thin-FASTA Format
Thin fasta is needed to store only needed for tests regions of real references to decrease disk usage. Each thin-FASTA file is `.csv` file, each line of which represent part of reference data with information of: 
- chromosome name
- start position of contig
- end position of contig
- and nucleotide sequence that corresponds to region

thin-FASTA example:
```
chr1,1,15,ATGCCCCCCCCCAAA
chr1,200,205,GCCGA
chr2,10,12,AC
```

>Note: VarDict expands given regions by 700bp to left and right (plus given value by `-x` option). 

##Program Workflow
The VarDictJava program follows the workflow:

1.	Get regions of interest from a BED file or the command line.
2.	For each segment:
	a.	Find all variants for this segment in mapped reads:
		i.	Optionally skip duplicated reads, low mapping-quality reads, and reads having a large number of mismatches.
		ii.	Skip a read if it does not overlap with the segment.
		iii.	Preprocess the CIGAR string for each read.
		iv.	For each position, create a variant. If a variant is already present, adjust its count using the adjCnt function.
	b. Realign some of the variants using special ad-hoc approaches.
	c. Calculate statistics for the variant, filter out some bad ones, if any.
	d. Assign a type to each variant.
	e. Output variants in an intermediate internal format (tabular). Columns of the table are described in the Output Columns section.     
          **Note**: To perform Steps 1 and 2, use Java VarDict.

3.	Perform a statistical test for strand bias using an R script.  
    **Note**: Use R script for this step.
4.	Transform the intermediate tabular format to VCF. Output the variants with filtering and statistical data.  
     **Note**: Use the Perl scripts `var2vcf_valid.pl` or `var2vcf_paired.pl` for this step.



##Program Options

- `-H`  
    Print help page
- `-h`   
    Print a header row decribing columns
- `-i`
    Output splicing read counts
- `-p`   
    Do pileup regarless the frequency
- `-C`    
    Indicate the chromosome names are just numbers, such as 1, 2, not chr1, chr2
- `-D`    
    Debug mode.  Will print some error messages and append full genotype at the end.
- `-t`   
    Indicate to remove duplicated reads.  Only one pair with identical start positions will be kept
- `-3`   
     Indicate to move indels to 3-prime if alternative alignment can be achieved.
- `-F bit`  
     The hexical to filter reads. Default: `0x500` (filter 2nd alignments and duplicates).  Use `-F 0` to turn it off.
- `-z 0/1`       
    Indicate whether the BED file contains zero-based cooridates, the same way as the Genome browser IGV does.  -z 1 indicates that coordinates in a BED file start from 0. -z 0 indicates that the coordinates start from 1. Default: `1` for a BED file or amplicon BED file.  Use `0` to turn it off. When using `-R` option, it is set to `0`
- `-a int:float`    
    Indicate it is amplicon based calling.  Reads that do not map to the amplicon will be skipped.  A read pair is considered to belong to the amplicon if the edges are less than int bp to the amplicon, and overlap fraction is at least float.  Default: `10:0.95`
- `-k 0/1`   
    Indicate whether to perform local realignment.  Default: `1` or yes.  Set to `0` to disable it.
- `-G Genome fasta`  
    The reference fasta.  Should be indexed (.fai).  Defaults to: `/ngs/reference_data/genomes/Hsapiens/hg19/seq/hg19.fa`
- `-R Region`  
    The region of interest.  In the format of chr:start-end.  If chr is not start-end but start (end is omitted), then it is a single position.  No BED is needed.
- `-d delimiter`  
    The delimiter for splitting `region_info`, defaults to tab `"\t"`
- `-n regular_expression`  
    The regular expression to extract sample names from bam filenames.  Defaults to: `/([^\/\._]+?)_[^\/]*.bam/`
- `-N string`   
    The sample name to be used directly.  Will overwrite `-n` option
- `-b string`   
    The indexed BAM file
- `-c INT`   
    The column for chromosome
- `-S INT`   
    The column for the region start, e.g. gene start
- `-E INT`  
    The column for the region end, e.g. gene end
- `-s INT`   
    The column for a segment starts in the region, e.g. exon starts
- `-e INT`  
    The column for a segment ends in the region, e.g. exon ends
- `-g INT`     
    The column for a gene name, or segment annotation
- `-x INT`   
    The number of nucleotides to extend for each segment, default: `0`
- `-f double`   
    The threshold for allele frequency, default: `0.05` or `5%`
- `-r minimum reads`   
    The minimum # of variance reads, default: `2`
- `-B INT`  
    The minimum # of reads to determine strand bias, default: `2`
- `-Q INT`  
    If set, reads with mapping quality less than INT will be filtered and ignored
- `-q INT`   
    The phred score for a base to be considered a good call.  Default: 25 (for Illumina). For PGM, set it to ~15, as PGM tends to underestimate base quality.
- `-m INT`   
    If set, reads with mismatches more than `INT` will be filtered and ignored.  Gaps are not counted as mismatches. Valid only for bowtie2/TopHat or BWA aln followed by sampe.  BWA mem is calculated as NM - Indels.  Default: 8, or reads with more than 8 mismatches will not be used.
- `-T INT`  
    Trim bases after `[INT]` bases in the reads
- `-X INT`   
    Extension of bp to look for mismatches after insersion or deletion.  Default to 3 bp, or only calls when they're within 3 bp.
- `-P number`  
    The read position filter.  If the mean variants position is less that specified, it is considered false positive.  Default: 5
- `-Z double`  
    For downsampling fraction,  e.g. `0.7` means roughly `70%` downsampling.  Default: No downsampling.  Use with caution.  The downsampling will be random and non-reproducible.
- `-o Qratio`  
    The `Qratio` of `(good_quality_reads)/(bad_quality_reads+0.5)`.  The quality is defined by `-q` option.  Default: `1.5`
- `-O MapQ`  
    The reads should have at least mean `MapQ` to be considered a valid variant.  Default: no filtering
- `-V freq`  
    The lowest frequency in a normal sample allowed for a putative somatic mutations.  Defaults to `0.05`
- `-I INT`  
    The indel size.  Default: 120bp
- `-M INT`
    The minimum matches for a read to be considered.  If, after soft-clipping, the matched bp is less than INT, then the 
    read is discarded.  It's meant for PCR based targeted sequencing where there's no insert and the matching is only the primers.
    Default: 0, or no filtering
- `-th [threads]`  
    If this parameter is missing, then the mode is one-thread. If you add the     -th parameter, the number of threads equals to the number of processor cores. The parameter -th threads sets the number of threads explicitly.
- `-VS STRICT | LENIENT | SILENT` 
    How strict to be when reading a SAM or BAM.
     `STRICT`   - throw an exception if something looks wrong.
     `LENIENT`  - Emit warnings but keep going if possible.
     `SILENT`   - Like `LENIENT`, only don't emit warning messages.
    Default: `LENIENT`

##Output columns

1. Sample - sample name
2. Gene - gene name from a BED file
3. Chr - chromosome name
4. Start - start position of the variation
5. End - end position of the variation
6. Ref - reference sequence
7. Alt - variant sequence
8. Depth - total coverage
9. AltDepth - variant coverage
10. RefFwdReads - reference forward strand coverage
11. RefRevReads - reference reverse strand coverage
12. AltFwdReads - variant forward strand coverage
13. AltRevReads - variant reverse strand coverage
14. Genotype - genotype description string
15. AF - allele frequency
16. Bias - strand bias flag
17. PMean - mean position in read
18. PStd - flag for read position standard deviation
19. QMean - mean base quality
20. QStd - flag for base quality standard deviation
23. QRATIO - ratio of high quality reads to low-quality reads
24. HIFREQ - variant frequency for high-quality reads
25. EXTRAFR - Adjusted AF for indels due to local realignment
26. SHIFT3 - No. of bases to be shifted to 3 prime for deletions due to alternative alignment
27. MSI - MicroSattelite. > 1 indicates MSI
28. MSINT - MicroSattelite unit length in bp
29. NM - average number of mismatches for reads containing the variant
30. HICNT - number of high-quality reads with the variant
31. HICOV - position coverage by high quality reads
21. 5pFlankSeq - neighboring reference sequence to 5' end 
22. 3pFlankSeq - neighboring reference sequence to 3' end
23. SEGMENT:CHR_START_END - position description
24. VARTYPE - variant type

License
-------

The code is freely available under the [MIT license](http://www.opensource.org/licenses/mit-license.html).

Contributors
------------

Java port of [VarDict](https://github.com/AstraZeneca-NGS/VarDict) implemented based on the original Perl version ([Zhongwu Lai](https://github.com/zhongwulai)) by:

- [Viktor Kirst](https://github.com/vkirst)

- [Zaal Lyanov](https://github.com/jabbarish)

