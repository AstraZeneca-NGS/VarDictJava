[![Bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](https://bioconda.github.io/recipes/vardict-java/README.html)
[![European Galaxy server](https://img.shields.io/badge/usegalaxy-.eu-brightgreen?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABgAAAASCAYAAABB7B6eAAAABGdBTUEAALGPC/xhBQAAACBjSFJNAAB6JgAAgIQAAPoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAACXBIWXMAAAsTAAALEwEAmpwYAAACC2lUWHRYTUw6Y29tLmFkb2JlLnhtcAAAAAAAPHg6eG1wbWV0YSB4bWxuczp4PSJhZG9iZTpuczptZXRhLyIgeDp4bXB0az0iWE1QIENvcmUgNS40LjAiPgogICA8cmRmOlJERiB4bWxuczpyZGY9Imh0dHA6Ly93d3cudzMub3JnLzE5OTkvMDIvMjItcmRmLXN5bnRheC1ucyMiPgogICAgICA8cmRmOkRlc2NyaXB0aW9uIHJkZjphYm91dD0iIgogICAgICAgICAgICB4bWxuczp0aWZmPSJodHRwOi8vbnMuYWRvYmUuY29tL3RpZmYvMS4wLyI+CiAgICAgICAgIDx0aWZmOlJlc29sdXRpb25Vbml0PjI8L3RpZmY6UmVzb2x1dGlvblVuaXQ+CiAgICAgICAgIDx0aWZmOkNvbXByZXNzaW9uPjE8L3RpZmY6Q29tcHJlc3Npb24+CiAgICAgICAgIDx0aWZmOk9yaWVudGF0aW9uPjE8L3RpZmY6T3JpZW50YXRpb24+CiAgICAgICAgIDx0aWZmOlBob3RvbWV0cmljSW50ZXJwcmV0YXRpb24+MjwvdGlmZjpQaG90b21ldHJpY0ludGVycHJldGF0aW9uPgogICAgICA8L3JkZjpEZXNjcmlwdGlvbj4KICAgPC9yZGY6UkRGPgo8L3g6eG1wbWV0YT4KD0UqkwAAAn9JREFUOBGlVEuLE0EQruqZiftwDz4QYT1IYM8eFkHFw/4HYX+GB3/B4l/YP+CP8OBNTwpCwFMQXAQPKtnsg5nJZpKdni6/6kzHvAYDFtRUT71f3UwAEbkLch9ogQxcBwRKMfAnM1/CBwgrbxkgPAYqlBOy1jfovlaPsEiWPROZmqmZKKzOYCJb/AbdYLso9/9B6GppBRqCrjSYYaquZq20EUKAzVpjo1FzWRDVrNay6C/HDxT92wXrAVCH3ASqq5VqEtv1WZ13Mdwf8LFyyKECNbgHHAObWhScf4Wnj9CbQpPzWYU3UFoX3qkhlG8AY2BTQt5/EA7qaEPQsgGLWied0A8VKrHAsCC1eJ6EFoUd1v6GoPOaRAtDPViUr/wPzkIFV9AaAZGtYB568VyJfijV+ZBzlVZJ3W7XHB2RESGe4opXIGzRTdjcAupOK09RA6kzr1NTrTj7V1ugM4VgPGWEw+e39CxO6JUw5XhhKihmaDacU2GiR0Ohcc4cZ+Kq3AjlEnEeRSazLs6/9b/kh4eTC+hngE3QQD7Yyclxsrf3cpxsPXn+cFdenF9aqlBXMXaDiEyfyfawBz2RqC/O9WF1ysacOpytlUSoqNrtfbS642+4D4CS9V3xb4u8P/ACI4O810efRu6KsC0QnjHJGaq4IOGUjWTo/YDZDB3xSIxcGyNlWcTucb4T3in/3IaueNrZyX0lGOrWndstOr+w21UlVFokILjJLFhPukbVY8OmwNQ3nZgNJNmKDccusSb4UIe+gtkI+9/bSLJDjqn763f5CQ5TLApmICkqwR0QnUPKZFIUnoozWcQuRbC0Km02knj0tPYx63furGs3x/iPnz83zJDVNtdP3QAAAABJRU5ErkJggg==)](https://usegalaxy.eu/root?tool_id=vardict_java)


# VarDictJava

## Introduction

VarDictJava is a variant discovery program written in Java and Perl. It is a Java port of [VarDict variant caller](https://github.com/AstraZeneca-NGS/VarDict). 

The original Perl VarDict is a sensitive variant caller for both single and paired sample variant calling from BAM files. VarDict implements several novel features such as amplicon bias aware variant calling from targeted
sequencing experiments, rescue of long indels by realigning bwa soft clipped reads and better scalability
than many other Java based variant callers. The Java port is around 10x faster than the original Perl implementation.

Please cite VarDict:

Lai Z, Markovets A, Ahdesmaki M, Chapman B, Hofmann O, McEwen R, Johnson J, Dougherty B, Barrett JC, and Dry JR.  VarDict: a novel and versatile variant caller for next-generation sequencing in cancer research. Nucleic Acids Res. 2016, pii: gkw227.

The link to is article can be accessed through: https://academic.oup.com/nar/article/44/11/e108/2468301?searchresult=1

Original coded by Zhongwu Lai 2014.

VarDictJava can run in single sample (see Single sample mode section), paired sample (see Paired variant calling section), or amplicon bias aware modes. As input, VarDictJava takes reference genomes in FASTA format, aligned reads in BAM format, and target regions in BED format.

## Requirements

1. JDK 1.8 or later
2. R language (uses /usr/bin/env R)
3. Perl (uses /usr/bin/env perl)
4. Internet connection to download dependencies using gradle.

To see the help page for the program, run
```
 <path_to_vardict_folder>/build/install/VarDict/bin/VarDict -H.
```
## Getting started

### Getting source code

The VarDictJava source code is located at [https://github.com/AstraZeneca-NGS/VarDictJava](https://github.com/AstraZeneca-NGS/VarDictJava).

To load the project, execute the following command:

```
git clone --recursive https://github.com/AstraZeneca-NGS/VarDictJava.git
```

Note that the original VarDict project is placed in this repository as a submodule and its contents can be found in the sub-directory VarDict in VarDictJava working folder. So when you use `teststrandbias.R` and `var2vcf_valid.pl.` (see details and examples below), you have to add prefix VarDict: `VarDict/teststrandbias.R` and `VarDict/var2vcf_valid.pl.`

### Compiling

The project uses [Gradle](http://gradle.org/) and already includes a gradlew script.

To build the project, in the root folder of the project, run the following command:

```
./gradlew clean installDist
```
Clean will remove all old files from build folder.

To generate Javadoc, in the build/docs/javadoc folder, run the following command. If you want to save content of `build` folder as it is (for example after building the project), run it without `clean` option:

```
./gradlew clean javadoc
```

To generate release version in the build/distributions folder as tar or zip archive, run the following command:
```
./gradlew distZip
``` 
or
```
./gradlew distTar
```

#### Distribution Package Structure
When the build command completes successfully, the `build/install/VarDict` folder contains the distribution package.

The distribution package has the following structure:
* `bin/` - contains the launch scripts
* `lib/` - has the jar file that contains the compiled project code and the jar files of the third-party libraries that the project uses.

You can move the distribution package (the content of the `build/install/VarDict` folder) to any convenient location.

Generated zip and tar releases will also contain scripts from VarDict Perl repository in `bin/` directory (`teststrandbias.R`, 
`testsomatic.R`, `var2vcf_valid.pl`, `var2vcf_paired.pl`).

You can add VarDictJava on PATH by adding this line to `.bashrc`:
```
export PATH=/path/to/VarDict/bin:$PATH
```
After that you can run VarDict by `Vardict` command instead of full path to `<path_to_vardict_folder>/build/install/VarDict/bin/VarDict`. 

#### Third-Party Libraries
Currently, the project uses the following third-party libraries:
* JRegex (http://jregex.sourceforge.net, BSD license) is a regular expressions library that is used instead of the 
standard Java library because its performance is much higher than that of the standard library.
* Commons CLI (http://commons.apache.org/proper/commons-cli, Apache License) – a library for parsing the command line.
* HTSJDK (http://samtools.github.io/htsjdk/) is an implementation of a unified Java library for accessing common file formats, such as SAM and VCF.
* Mockito and TestNG are the testing frameworks (not included in distribution, used only in tests).

### Single sample mode

To run VarDictJava in single sample mode, use a BAM file specified without the `|` symbol and perform Steps 3 and 4 
(see the Program workflow section) using `teststrandbias.R` and `var2vcf_valid.pl`.
The following is an example command to run in single sample mode with BED file.   
You have to set options `-c`, `-S`, `-E`, `-g` using number of columns in your BED file for chromosome, start, end
 and gene of region respectively:
  
```
AF_THR="0.01" # minimum allele frequency
<path_to_vardict_folder>/build/install/VarDict/bin/VarDict -G /path/to/hg19.fa -f $AF_THR -N sample_name -b /path/to/my.bam -c 1 -S 2 -E 3 -g 4 /path/to/my.bed | VarDict/teststrandbias.R | VarDict/var2vcf_valid.pl -N sample_name -E -f $AF_THR > vars.vcf
```

VarDictJava can also be invoked without a BED file if the region is specified in the command line with `-R` option.
The following is an example command to run VarDictJava for a region (chromosome 7, position from 55270300 to 55270348, EGFR gene) with `-R` option:

```
<path_to_vardict_folder>/build/install/VarDict/bin/VarDict  -G /path/to/hg19.fa -f 0.001 -N sample_name -b /path/to/sample.bam -R  chr7:55270300-55270348:EGFR | VarDict/teststrandbias.R | VarDict/var2vcf_valid.pl -N sample_name -E -f 0.001 > vars.vcf
```

In single sample mode, output columns contain a description and statistical info for variants in the single sample. 
See section Output Columns for list of columns in the output. 

### Paired variant calling

To run paired variant calling, use BAM files specified as `BAM1|BAM2` and perform Steps 3 and 4 
(see the Program Workflow section) using `testsomatic.R` and `var2vcf_paired.pl`.

In this mode, the number of statistics columns in the output is doubled: one set of columns is 
for the first sample, the other - for second sample.

The following is an example command to run in paired mode.  
You have to set options `-c`, `-S`, `-E`, `-g` using number of columns in your bed file for chromosome, start, 
 end and gene of region respectively:

```
AF_THR="0.01" # minimum allele frequency
<path_to_vardict_folder>/build/install/VarDict/bin/VarDict -G /path/to/hg19.fa -f $AF_THR -N tumor_sample_name -b "/path/to/tumor.bam|/path/to/normal.bam" -c 1 -S 2 -E 3 -g 4 /path/to/my.bed | VarDict/testsomatic.R | VarDict/var2vcf_paired.pl -N "tumor_sample_name|normal_sample_name" -f $AF_THR
```

### Amplicon based calling 
This mode is active if the BED file uses 8-column format and the -R option is not specified.

In this mode, only the first list of BAM files is used even if the files are specified as `BAM1|BAM2` - like for paired variant calling.

For each segment, the BED file specifies the list of positions as start and end positions (columns 7 and 8 of 
the BED file). The Amplicon based calling mode outputs a record for every position between start and end that has 
any variant other than the reference one (all positions with the `-p` option). For any of these positions, 
VarDict in amplicon based calling mode outputs the following:
* Same columns as in the single sample mode for the most frequent variant
* Good variants for this position with the prefixes `GOOD1`, `GOOD2`, etc.
* Bad variants for this position with the prefixes `BAD1`, `BAD2`, etc.

For this running mode, the `-a` option (default: `10:0.95`) specifies the criteria of discarding reads that are too 
far away from segments. A read is skipped if its start and end are more than 10 positions away from the segment 
ends and the overlap fraction between the read and the segment is less than 0.95.

### Running Tests
#### Integration testing

The list of integration test cases is stored in files in `testdata/intergationtestcases` directory.
To run all integration tests, the command is:

```$xslt
./gradlew test --tests com.astrazeneca.vardict.integrationtests.IntegrationTest 
```

The results of the tests can be viewed in the `build/reports/tests/index.html` file.

##### User extension of testcases

Each file in `testdata/intergationtestcases` directory represents a test case with input data and expected output

**1.** Create a txt file in `testdata/intergationtestcases` folder. 

The file contains testcase input (of format described in [Test cases file format](Readme.md#CSVhead)) in the first line and expected output in the remaining file part.

**2.** Extend or create [thin-FASTA](Readme.md#thFastahead) in `testdata/fastas` folder.
 
**3.** Run tests.

##### <a name="CSVhead"></a>Test cases file format 

Each input file represents one test case input description. In the input file the first line consists of the following fields separated by `,` symbol:

*Required fields:*
- test case name (Amplicon/Somatic/Simple mode)
- reference name
- bam file name
- chromosome name
- start of region
- end of region

*Optional fields:*
- start of region with amplicon case
- end of region with amplicon case

*Parameters field:*
- the last field can be any other command line parameters string

Example of first line of input file:

```
Amplicon,hg19.fa,Colo829-18_S3-sort.bam,chr1,933866,934466,933866,934466,-a 10:0.95 -D
Somatic,hg19.fa,Colo829-18_S3-sort.bam|Colo829-19_S4-sort.bam,chr1,755917,756517
Simple,hg19.fa,Colo829-18_S3-sort.bam,chr1,9922,10122,-p
```

##### <a name="testCoverageReporthead"></a>Test coverage Report

To build test coverage report run the following command:

```
./gradlew test jacocoTestReport 
```

Then HTML report could be found in `build/reports/jacoco/test/html/index.html`

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

>Note: VarDict expands given regions by 1200bp to left and right (plus given value by `-x` option). 

## Program Workflow

#### The main workflow
The VarDictJava program follows the workflow:

1. Get regions of interest from a BED file or the command line.
2. For each segment:
   1. Find all variants for this segment in mapped reads:
      1. Optionally skip duplicated reads, low mapping-quality reads, and reads having a large number of mismatches.
      2. Skip unmapped reads.
      3. Skip a read if it does not overlap with the segment.
      4. Preprocess the CIGAR string for each read (section CIGAR Preprocessing).
      5. For each position, create a variant. If a variant is already present, adjust its count using the adjCnt function.
      6. Combine SNVs into MNV or SNV with indel to complex if variants are located closer than `-X` option (default <=2 bases) and have good base qualities. 
   2. Find structural variants (optionally can be disabled by option `-U`).
   3. Realign some of the variants using realignment of insertions, deletions, large insertions, and large deletions using unaligned parts of reads 
   (soft-clipped ends). This step is optional and can be disabled using the `-k 0` switch.
   4. Apply variant filtering rules (hard filters) defined in [Variant filtering](Readme.md#VarFiltering). 
   5. Assign a type to each variant.
   6. Output variants in an intermediate internal format (tabular). Columns of the table are described in the Output Columns section.
	 
	 **Note**: To perform Steps 1 and 2, use Java VarDict.

3.	Perform a statistical test for strand bias using an R script.  
    **Note**: Use R scripts `teststrandbias.R` or `testsomatic.R` for this step.
4.	Transform the intermediate tabular format to VCF. Output the variants with filtering and statistical data.  
     **Note**: Use the Perl scripts `var2vcf_valid.pl` or `var2vcf_paired.pl` for this step. Be aware that `var2vcf_valid.pl` or `var2vcf_paired.pl` by default will output the variant with the highest AF on position: if few variants start at one position, only the highest will be added in VCF. To output all variants use `-A` option with these perl scripts. 
     
     
#### CIGAR Preprocessing (Initial Realignment)
Read alignment is specified in a BAM file as a CIGAR string. VarDict modifies this string (and alignment) in the following special cases:
* Soft clipping next to insertion/deletion is replaced with longer soft-clipping. 
The same takes place if insertion/deletion is separated from soft clipping by no more than 10 matched bases.
* Short matched sequence and insertion/deletion at the beginning/end are replaced by soft-clipping.
* Two close deletions and insertions are combined into one deletion and one insertion
* Three close deletions are combined in one
* Three close insertions/deletions are combined in one deletion or in insertion+deletion
* Two close deletions are combined into one
* Two close insertions/deletions are combined into one
* Mis-clipping at the start/end are changed to matched sequences

#### Variants
Simple variants (SNV, simple insertions, and deletions) are constructed in the following way:
* Single-nucleotide variation (SNV). VarDict inserts an SNV into the variants structure for every matched or 
mismatched base in the reads. If an SNV is already present in variants, VarDict adjusts its counts and statistics.
* Simple insertion variant. If read alignment shows an insertion at the position, VarDict inserts +BASES 
string into the variants structure. If the variant is already present, VarDict adjusts its count and statistics.
* Simple Deletion variant. If read alignment shows a deletion at the position, VarDict inserts -NUMBER 
into the variants structure. If the variant is already present, VarDict adjusts its count and statistics.
* Complex variants: VarDict also handles complex variants (for example, an insertion that is close to SNV or to deletion) 
using specialized ad-hoc methods.

Structural Variants are looked for after simple variants. VarDict supported DUP, INV and DEL structural variants.

#### Variant Description String
The description string encodes a variant for VarDict internal use. 

The following table describes Variant description string encoding:

String | Description 
 ----- | ---------- 
[ATGC] | for SNPs 
+[ATGC]+ | for insertions 
-[0-9]+ | for deletions
...#[ATGC]+ | for insertion/deletion variants followed by a short matched sequence
...^[ATGC]+ | something followed by an insertion
...^[0-9]+ | something followed by a deletion
...&amp;[ATGC]+ | for insertion/deletion variants followed by a matched sequence

#### <a name="VarFiltering"></a>Variant Filtering
A variant appears in the output if it satisfies the following criteria (in this order). 
If variant doesn't fit criteria on the step, it will be filtered out and the next steps won't be checked (except for the step 8, read the explanation below):
1. Frequency of the variant exceeds the threshold set by the `-f` option (default = 1%).
2. The minimum number of high-quality reads supporting variant is larger than the threshold set by the `-r` option (default = 2).
3. The mean position of the variant in reads is larger than the value set by the `-P` option (default = 5).
4. The mean base quality (phred score) for the variant is larger than the threshold set by the `-q` option (default = 22.5).
5. Variant frequency is more than 25% or reference allele does not have much better mapping quality than the variant.
6. Deletion variants are not located in the regions where the reference genome is missing.
7. The ratio of high-quality reads to low-quality reads is larger than the threshold specified by `-o` option (default=1.5).
8. Variant frequency exceeds 30%. If so, next steps won't be checked and variant considered as "good". Otherwise, other steps will be also checked.
9. Mean mapping quality exceeds the threshold set by the `-O` option (default: no filtering)
10. In the case of an MSI region, the variant size is less than 12 nucleotides for the non-monomer MSI or 15 for the monomer MSI. 
Variant frequency is more than 10% for the non-monomer MSI (or set by `--nmfreq` option) and 25% for the monomer MSI (or set by `--mfreq` option).
11. Variant has not "2;1" bias or variant frequency more than 20%. If both conditions aren't met, then variant mustn't be SNV and any of variants refallele or varallele lengths must be more than 3 nucleotides.

#### Bias flag explanation
Bias flag can take values [0-2];[0-2] (i.e. "0;2", "2;1" and separator can be another in paired and single VCF).
The first value refers to reads that support the reference allele, and the second to reads that support the variant allele.

0 - small total count of reads (less than 12 for the sum of forward and reverse reads)
1 - strand bias
2 - no strand bias

#### Variant classification in paired(somatic) analysis
In paired analysis, VarDict will classify each variant into the following types that are propagated 
into STATUS info tag after `var2vcf_paired.pl` script.

When both samples have coverage for the variant:

* Germline: detected in germline sample (pass all quality parameters)
* StrongSomatic: detected in tumor sample only
* LikelySomatic: the variant has at least one read support OR allele frequency < 5% (defined by –V option with default 0.05)
* StrongLOH: detected in germline sample only, opposite of StrongSomaitc
* LikelyLOH: detected in germline but either lost in tumor OR 20-80% in germline, but increased to 1-opt_V (95%).
* AFDiff: detected in tumor (pass quality parameters) and present in germline but didn’t pass quality parameters.

When only one sample has coverage for the variant:

* SampleSpecific: detected in tumor sample, but no coverage in germline sample (it’s more technical than biological, as it’s unlikely a tumor sample can gain a piece of sequence in reference that germline sample lacks).
* Deletion: detected in germline sample, but no coverage in tumor sample

These are only rough classification. You need to examine the p-value (after testsomatic.R script) to determine whether or not it's significant.
## Program Options
### VarDictJava options
- `-H|-?`  
    Print help page
- `-h|--header`   
    Print a header row describing columns
- `-i|--splice `
    Output splicing read counts
- `-p`   
    Do pileup regardless of the frequency
- `-C`    
    Indicate the chromosome names are just numbers, such as 1, 2, not chr1, chr2 (deprecated)
- `-D|--debug`    
    Debug mode.  Will print some error messages and append full genotype at the end.
- `-y|--verbose`   
    Verbose mode.  Will output variant calling process.
- `-t|--dedup`   
    Indicate to remove duplicated reads.  Only one pair with identical start positions will be kept
- `-3`   
    Indicate to move indels to 3-prime if alternative alignment can be achieved.
- `-K`
    Include Ns in the total depth calculation.
- `-F bit`  
    The hexical to filter reads. Default: `0x504` (filter unmapped reads, 2nd alignments and duplicates).  Use `-F 0` to turn it off.
- `-z 0/1`       
    Indicate whether the BED file contains zero-based coordinates, the same way as the Genome browser IGV does.  -z 1 indicates that coordinates in a BED file start from 0. -z 0 indicates that the coordinates start from 1. Default: `1` for a BED file or amplicon BED file (0-based).  Use `0` to turn it off. When using `-R` option, it is set to `0`
- `-a|--amplicon int:float`    
    Indicate it is amplicon based calling.  Reads that do not map to the amplicon will be skipped.  A read pair is considered to belong to the amplicon if the edges are less than int bp to the amplicon, and overlap fraction is at least float.  Default: `10:0.95`
- `-k 0/1`   
    Indicate whether to perform local realignment.  Default: `1` or yes.  Set to `0` to disable it.
- `-G Genome fasta`  
    The reference fasta.  Should be indexed (.fai).  Defaults to: `/ngs/reference_data/genomes/Hsapiens/hg19/seq/hg19.fa`  
    Also short commands can be used to set path to:  
    **hg19** - /ngs/reference_data/genomes/Hsapiens/hg19/seq/hg19.fa  
    **hg38** - /ngs/reference_data/genomes/Hsapiens/hg38/seq/hg38.fa  
    **mm10** - /ngs/reference_data/genomes/Mmusculus/mm10/seq/mm10.fa  
- `-R Region`  
    The region of interest.  In the format of chr:start-end.  If chr is not start-end but start (end is omitted), then it is a single position.  No BED is needed.
- `-d delimiter`  
    The delimiter for splitting `region_info`, defaults to tab `"\t"`
- `-n regular_expression`  
    The regular expression to extract sample names from bam filenames.  Defaults to: `/([^\/\._]+?)_[^\/]*.bam/`
- `-N string`   
    The sample name to be used directly.  Will overwrite `-n` option
- `-b string`   
    The indexed BAM file. Multiple BAM files can be specified with the “:” delimiter.
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
    The threshold for allele frequency, default: `0.01` or `1%`
- `-r minimum reads`   
    The minimum # of variance reads, default: `2`
- `-B INT`  
    The minimum # of reads to determine strand bias, default: `2`
- `-Q INT`  
    If set, reads with mapping quality less than INT will be filtered and ignored
- `-q double`   
    The phred score for a base to be considered a good call.  Default: 22.5 (for Illumina). For PGM, set it to ~15, as PGM tends to underestimate base quality.
- `-m INT`   
    If set, reads with mismatches more than `INT` will be filtered and ignored.  Gaps are not counted as mismatches. 
    Valid only for bowtie2/TopHat or BWA aln followed by sampe.  BWA mem is calculated as NM - Indels.  For STAR 
    you have to increase default if nM tag (for "paired" alignment) is presented in reads.Default: 8, or reads with more than 8 mismatches will not be used.
- `-T|--trim INT`  
    Trim bases after `[INT]` bases in the reads
- `-X INT`   
    Extension of bp to look for mismatches after insersion or deletion.  Default to 2 bp, or only calls when they are within 2 bp.
- `-P number`  
    The read position filter.  If the mean variants position is less that specified, it is considered false positive.  Default: 5
- `-Z|--downsample double`  
    For downsampling fraction,  e.g. `0.7` means roughly `70%` downsampling.  Default: No downsampling.  Use with caution.  The downsampling will be random and non-reproducible.
- `-o Qratio`  
    The `Qratio` of `(good_quality_reads)/(bad_quality_reads+0.5)`.  The quality is defined by `-q` option.  Default: `1.5`
- `-O MapQ`  
    The variant should has at least mean `MapQ` to be considered a valid variant.  Default: no filtering
- `-V freq`  
    The lowest frequency in a normal sample allowed for a putative somatic mutations. Used only in paired mode. Defaults to `0.05`
- `-I INT`  
    The indel size.  Default: 50bp. 
    Be cautious with -I option, especially in the amplicon mode, as amplicon sequencing is not a way 
    to find large indels. Increasing the search size might be slow and false positives may appear in low
    complexity regions. Increasing it to 200-300 bp is only recommend for hybrid capture sequencing.
- `-M INT`  
    The minimum matches for a read to be considered.  If, after soft-clipping, the matched bp is less than INT, then the 
    read is discarded.  It's meant for PCR based targeted sequencing where there's no insert and the matching is only the primers.
    Default: 0, or no filtering
- `-th [threads]`  
    If this parameter is missing, then the mode is one-thread. If you add the -th parameter, the number of threads 
    equals to the number of processor cores. The parameter -th threads sets the number of threads explicitly.
- `-VS STRICT | LENIENT | SILENT`   
    How strict to be when reading a SAM or BAM.
     `STRICT`   - throw an exception if something looks wrong.
     `LENIENT`  - Emit warnings but keep going if possible.
     `SILENT`   - Like `LENIENT`, only don't emit warning messages.
    Default: `LENIENT`
- `-u`  
    Indicate unique mode, which when mate pairs overlap, the overlapping part will be counted only once using **forward** read only.
    Default: unique mode disabled, all reads are counted.
- `-UN`  
    Indicate unique mode, which when mate pairs overlap, the overlapping part will be counted only once using **first** read only.
    Default: unique mode disabled, all reads are counted.
- `--chimeric`  
    Indicate to turn off chimeric reads filtering.  Chimeric reads are artifacts from library construction, 
    where a read can be split into two segments, each will be aligned within 1-2 read length distance,
    but in opposite direction. Default: filtering enabled
- `-U|--nosv`  
    Turn off structural variant calling
- `-L INT`   
   The minimum structural variant length to be presented using \<DEL\> \<DUP\> \<INV\> \<INS\>, etc. 
   Default: 1000. Any indel, complex variants less than this will be spelled out with exact nucleotides
- `-w|--insert-size INT` INSERT_SIZE  
   The insert size. Used for SV calling. Default: 300
- `-W|--insert-std INT` INSERT_STD  
   The insert size STD. Used for SV calling. Default: 100
- `-A INT` INSERT_STD_AMT  
   The number of STD. A pair will be considered for DEL if INSERT > INSERT_SIZE + INSERT_STD_AMT * INSERT_STD. Default: 4
- `-Y|--ref-extension INT`  
    Extension of bp of reference to build lookup table. Default to 1200 bp. Increasing the number will slow down the program.
    The main purpose is to call large indels with 1000 bp that can be missed by discordant mate pairs. 
- `--deldupvar`  
  Turn on deleting of duplicate variants in output that can appear due to VarDict linear work on regions. Variants in this mode are 
  considered and outputted only if start position of variant is inside the region interest.
- `-DP|--default-printer`   
    The printer type used for different outputs. Default: OUT (i.e. System.out).
- `--adaptor`  
    Filter adaptor sequences so that they are not used in realignment. Multiple adaptors can be supplied by setting them
     with comma, like:   
     --adaptor ACGTTGCTC,ACGGGGTCTC,ACGCGGCTAG 
- `-J|--crispr CRISPR_cutting_site`  
    The genomic position that CRISPR/Cas9 suppose to cut, typically 3bp from the PAM NGG site and within the guide.  For
   CRISPR mode only.  It will adjust the variants (mostly In-Del) start and end sites to as close to this location as possible,
    if there are alternatives. The option should only be used for CRISPR mode.
- `-j CRISPR_filtering_bp`  
    In CRISPR mode, the minimum amount in bp that a read needs to overlap with cutting site.  If a read does not meet the criteria,
    it will not be used for variant calling, since it is likely just a partially amplified PCR.  Default: not set, or no filtering  
- `--nmfreq`  
    The variant frequency threshold to determine variant as good in case of non-monomer MSI. Default: 0.1 
- `--mfreq`  
    The variant frequency threshold to determine variant as good in case of monomer MSI. Default: 0.25
- `--fisher`  
    EXPERIMENTAL FEATURE: to exclude R script from the VarDict pipeline we added this option to calculate pvalue and oddratio from Fisher Test. 
    It will decrease time processing on big samples because R script uses slow `textConnection` function.
   If you use this, do NOT run `teststrandbias.R` or `testsomatic.R` after Vardict, but use `var2vcf_valid.pl`
    or `var2vcf_paired.pl` after VarDictJava as usual.
    
### Important var2vcf_valid.pl options
The full list of options in VarDictPerl `var2vcf_valid.pl -h`
-  `-A`  
    Indicate to output all variants at the same position.  By default, only the variant with the highest allele frequency is converted to VCF.
-  `-S`  
    If set, variants that didn't pass filters will not be present in VCF file.
-  `-N` string  
    The sample name to be used directly.
-  `-f`  float  
    The minimum allele frequency.  Default to 0.02
    
### Important var2vcf_paired.pl options
The full list of options in VarDictPerl `var2vcf_paired.pl -h`
-  `-M`  
    If set, will increase stringency for candidate somatic: flag P0.01Likely and InDelLikely, and add filter P0.05
-  `-A`  
    Indicate to output all variants at the same position.  By default, only the variant with the highest allele frequency is converted to VCF.
-  `-S`  
    If set, variants that didn't pass filters will not be present in VCF file.
-  `-N` string  
    The sample name(s).  If only one name is given, the matched will be simply names as "name-match".  Two names are given separated by "|", such as "tumor|blood".
-  `-f`  float  
    The minimum allele frequency.  Default to 0.02
    
    
## Output columns
### Simple mode:
1. Sample - sample name
2. Gene - gene name from a BED file
3. Chr - chromosome name
4. Start - start position of the variation
5. End - end position of the variation
6. Ref - reference sequence
7. Alt - variant sequence
8. Depth (DP) - total coverage
9. AltDepth (VD) - variant coverage
10. RefFwdReads (REFBIAS) - reference forward strand coverage
11. RefRevReads (REFBIAS) - reference reverse strand coverage
12. AltFwdReads (VARBIAS) - variant forward strand coverage
13. AltRevReads (VARBIAS) - variant reverse strand coverage
14. Genotype - genotype description string
15. AF - allele frequency
16. Bias - strand bias flag
17. PMean - mean position in read
18. PStd - flag for read position standard deviation (1 if the variant is covered by at least 2 read segments with different positions, otherwise 0).
19. QMean - mean base quality
20. QStd - flag for base quality standard deviation
21. MAPQ - mapping quality
22. QRATIO (SN) - ratio of high quality reads to low-quality reads
23. HIFREQ (HIAF) - variant frequency for high-quality reads
24. EXTRAFR (ADJAF) - Adjusted AF for indels due to local realignment
25. SHIFT3 - No. of bases to be shifted to 3 prime for deletions due to alternative alignment
26. MSI - MicroSatellite. > 1 indicates MSI
27. MSILEN - MicroSatellite unit length in bp
28. NM - average number of mismatches for reads containing the variant
29. HICNT - number of high-quality reads with the variant
30. HICOV - position coverage by high quality reads
31. 5pFlankSeq (LSEQ) - neighboring reference sequence to 5' end 
32. 3pFlankSeq (RSEQ) - neighboring reference sequence to 3' end
33. SEGMENT:CHR_START_END - position description
34. VARTYPE - variant type
35. DUPRATE - duplication rate in fraction
36. SV splits-pairs-clusters: Splits (SPLITREAD) - No. of split reads supporting SV, Pairs (SPANPAIR) - No. of pairs supporting SV, 
Clusters - No. of clusters supporting SV 
37. CRISPR - only in crispr mode - how close to a CRISPR site is the variant

### Amplicon mode
In amplicon mode columns from #35 are changed to:  
(35) GoodVarCount (GDAMP) - number of good variants on amplicon  
(36) TotalVarCount (TLAMP) - number of good and bad variants on amplicon   
(37) Nocov (NCAMP) - number of variants on amplicon that has depth less than 1/50 of the max depth (they will be considered not working and thus not used).  
(38) Ampflag - if there are different good variants on different amplicons, it will be 1.

### Somatic mode
In somatic mode we have information from both samples:
1. Sample - sample name
2. Gene - gene name from a BED file
3. Chr - chromosome name
4. Start - start position of the variation
5. End - end position of the variation
6. Ref - reference sequence  
7. Alt - variant sequence  
   Fields from first sample:
8. Depth (DP) - total coverage
9. AltDepth (VD) - variant coverage
10. RefFwdReads (REFBIAS) - reference forward strand coverage
11. RefRevReads (REFBIAS) - reference reverse strand coverage
12. AltFwdReads (VARBIAS) - variant forward strand coverage
13. AltRevReads (VARBIAS) - variant reverse strand coverage
14. Genotype - genotype description string
15. AF - allele frequency
16. Bias - strand bias flag
17. PMean - mean position in read
18. PStd - flag for read position standard deviation
19. QMean - mean base quality
20. QStd - flag for base quality standard deviation
21. MAPQ - mapping quality
22. QRATIO (SN) - ratio of high quality reads to low-quality reads
23. HIFREQ (HIAF) - variant frequency for high-quality reads
24. EXTRAFR (ADJAF) - Adjusted AF for indels due to local realignment  
25. NM - average number of mismatches for reads containing the variant  
    Fields from second sample:
26. Depth - total coverage
27. AltDepth - variant coverage
28. RefFwdReads (REFBIAS) - reference forward strand coverage
29. RefRevReads (REFBIAS) - reference reverse strand coverage
30. AltFwdReads (VARBIAS) - variant forward strand coverage
31. AltRevReads (VARBIAS) - variant reverse strand coverage
32. Genotype - genotype description string
33. AF - allele frequency
34. Bias - strand bias flag
35. PMean - mean position in read
36. PStd - flag for read position standard deviation
37. QMean - mean base quality
38. QStd - flag for base quality standard deviation
39. MAPQ - mapping quality
40. QRATIO (SN) - ratio of high quality reads to low-quality reads
41. HIFREQ (HIAF) - variant frequency for high-quality reads
42. EXTRAFR (ADJAF) - Adjusted AF for indels due to local realignment  
43. NM - average number of mismatches for reads containing the variant  
    Common fields: 
44. SHIFT3 - No. of bases to be shifted to 3 prime for deletions due to alternative alignment
45. MSI - MicroSatellite. > 1 indicates MSI
46. MSILEN - MicroSatellite unit length in bp
47. 5pFlankSeq (LSEQ) - neighboring reference sequence to 5' end 
48. 3pFlankSeq (RSEQ) - neighboring reference sequence to 3' end
49. SEGMENT:CHR_START_END - position description
50. VarLabel - variant label due to type: StrongLOH, StrongSomatic...
51. VARTYPE - variant type
52. DUPRATE1 - duplication rate in fraction from first sample
53. SV_info1 - Splits - No. of split reads supporting SV, Pairs - No. of pairs supporting SV, 
Clusters - No. of clusters supporting SV from first sample
54. DUPRATE2 - duplication rate in fraction from second sample
55. SV_info2: Splits - No. of split reads supporting SV, Pairs - No. of pairs supporting SV, 
Clusters - No. of clusters supporting SV from second sample

### Input Files

#### BED File – Regions
VarDict uses 2 types of BED files for specifying regions of interest: 8-column and all others. 
The 8-column file format is used for targeted DNA deep sequencing analysis (amplicon based calling), amplicon analysis will 
try to start if BED with 8 columns was provided.
Otherwise you can start single and paired sample analysis by providing options `-c`, `-S`, `-E`, `-g` 
with number of columns for chromosome, start, end, gene of the region respectively.

All lines starting with #, browser, and track in a BED file are skipped. 
The column delimiter can be specified as the `-d` option (the default value is a tab “\t“).

The 8-column amplicon BED file format involves the following data:
* Chromosome name
* Region start position
* Region end position
* Gene name
* Score - not used by VarDict
* Strand - not used by VarDict
* Start position – VarDict starts outputting variants from this position
* End position – VarDict ends outputting variants from this position

For example 4-column BED file format involves the following data and VarDict must be start with `-c 1 -S 2 -E 3 -g 4` to
recognize it:
* Chromosome name
* Region start position
* Region end position
* Gene name

#### FASTA File - Reference Genome
The reference genome in FASTA format is read using HTSJDK library. 
For every invocation of the VarDict pipeline (usually 1 for a region in a BED file)
and for every BAM file, a part of the reference genome is extracted from the FASTA file. In some cases of Structural Variants finding
the reference can be reread in other regions. 

Region of FASTA extends and this extension can be regulated via the REFEXT variable (option `-Y INT`, default 1200 bp).

# Errors and warnings
Information about some of the errors and their causes is located in [wiki](https://github.com/AstraZeneca-NGS/VarDictJava/wiki)

# License
The code is freely available under the [MIT license](http://www.opensource.org/licenses/mit-license.html).

# Contributors
Java port of [VarDict](https://github.com/AstraZeneca-NGS/VarDict) implemented based on the original Perl version ([Zhongwu Lai](https://github.com/zhongwulai)) by:

- [Viktor Kirst](https://github.com/vkirst)

- [Zaal Lyanov](https://github.com/jabbarish)

