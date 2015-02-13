#VarDictJava

##Introduction
VarDictJava is a variant discovery program written in Java and Perl. It is a partial Java port of [VarDict variant caller](https://github.com/AstraZeneca-NGS/VarDict). 

Original Perl VarDict is a sensitive variant caller for both single and paired sample variant calling from BAM files. VarDict implements several novel features such as amplicon bias aware variant calling from targeted
sequencing experiments, rescue of long indels by realigning bwa soft clipped reads and better scalability
than Java based variant callers.

Original code by Zhongwu Lai 2014.

VarDictJava can run in [single sample](#singleSample), [paired sample](#pairedSample), or amplicon bias aware modes. As input, VarDictJava takes reference genomes in FASTA format, aligned reads in BAM format, and target regions in BED format.

##Requirements
1. JDK 1.7 or later
2. Samtools (must be in path)
3. R language (uses /usr/bin/env R)
4. Perl (uses /usr/bin/env perl)
3. Internet connection to download dependencies using gradle.

##Getting started
###Getting source code
The VarDictJava source code is located at [https://github.com/AstraZeneca-NGS/VarDictJava](https://github.com/AstraZeneca-NGS/VarDictJava).

To load the project, execute the following command:

```
git clone https://github.com/AstraZeneca-NGS/VarDictJava.git
```

###Compiling
The project uses [Gradle](http://gradle.org/) and already includes a gradlew script.

To build the project, in the root folder of the project, run the following command:

```
./gradlew clean installApp 
```

<a name="singleSample">
###Single sample mode
</a>

To run VarDictJava in single sample mode, use a BAM file specified without the `|` symbol and perform Steps 3 and 4 (see the [workflow](#programWorkflow)) using `teststrandbias.R` and `var2vcf_valid.pl.`
The following is an example command to run in single sample mode:
  
```
AF_THR="0.01" # minimum allele frequency
<path_to_vardict_folder>/bin/VarDict -G /path/to/hg19.fa -f $AF_THR -N sample_name -b /path/to/my.bam -z -c 1 -S 2 -E 3 -g 4 /path/to/my.bed | teststrandbias.R | var2vcf_valid.pl -N sample_name -E -f $AF_THR
```

VarDictJava can also be invoked without BED file if region is specified on command line with `-R` option.
The following is an example command to run VarDictJava for a region (chromosome 7, position from 55270300 to 55270348, EGFR gene) with `-R` option:

```
<path_to_vardict_folder>/bin/VarDict  -G /path/to/hg19.fa -f 0.001 -N sample_name -b /path/to/sample.bam  -z -R  chr7:55270300-55270348:EGFR | teststrandbias.R | var2vcf_valid.pl -N sample_name -E -f 0.001 >vars.vcf
```

In single sample mode, output columns contain a description and statistical info for variants in the single sample. See section [Output Columns](##Output Columns) for list of columns in the output. 

<a name="pairedSample">
###Paired variant calling
</a>

To run paired variant calling, use BAM files specified as `BAM1|BAM2` and perform Steps 3 and 4 (see the [workflow](##Program Workflow)) using `testsomatic.R` and `var2vcf_somatic.pl`.

In this mode, the number of statistics columns in the output is doubled: one set of columns is for the first sample, the other - for second sample.

The following is an example command to run in paired mode:
```
AF_THR="0.01" # minimum allele frequency
<path_to_vardict_folder>/bin/VarDict -G /path/to/hg19.fa -f $AF_THR -N tumor_sample_name -b "/path/to/tumor.bam|/path/to/normal.bam" -z -F -c 1 -S 2 -E 3 -g 4 /path/to/my.bed | testsomatic.R | var2vcf_somatic.pl -N "tumor_sample_name|normal_sample_name" -f $AF_THR
```

<a name="programWorkflow">
##Program Workflow
</a>
The VarDict program follows the workflow:

1.	Get regions of interest from a BED file or the command line.
2.	For each segment:
	1.	Find all variants for this segment in mapped reads:
		1.	Optionally skip duplicated reads, low mapping-quality reads, and reads having a large number of mismatches.
		2.	Skip a read if it does not overlap with the segment.
		3.	Preprocess the CIGAR string for each read.
		4.	For each position, create a variant. If a variant is already present, adjust its count using the adjCnt function.
	2. Realign some of the variants using special ad-hoc approaches.
	3. Calculate statistics for the variant, filter out some bad ones, if any.
	4. Assign a type to each variant.
	5. Output variants in an intermediate internal format (tabular). Columns of the table are described in [Output columns section](#outputColumns).     
          **Note**: To perform Steps 1 and 2, use the Java program VarDict-0.1.

3.	Perform a statistical test for strand bias using an R script.  
    **Note**: Use R script for this step.
4.	Transform the intermediate tabular format to VCF. Output the variants with filtering and statistical data.  
     **Note**: Use the Perl scripts `var2vcf_valid.pl` or `var2vcf_somatic.pl` for this step.



##Program Options

- `-H`  
    Print help page
- `-h`   
    Print a header row decribing columns
- `-p`   
    Do pileup regarless the frequency
- `-C`    
    Indicate the chromosome names are just numbers, such as 1, 2, not chr1, chr2
- `-D`    
    Debug mode.  Will print some error messages and append full genotype at the end.
- `-t`   
    Indicate to remove duplicated reads.  Only one pair with same start positions will be kept
- `-3`   
     Indicate to move indels to 3-prime if alternative alignment can be achieved.
- `-F bit`  
     The hexical to filter reads using samtools. Default: `0x500` (filter 2nd alignments and duplicates).  Use `-F 0` to turn it off.
- `-z 0/1`       
    Indicate whether is zero-based cooridates, as IGV does.  Default: `1` for BED file or amplicon BED file.  Use `0` to turn it off. When use `-R` option, it's set to `0`
- `-a int:float`    
    Indicate it's amplicon based calling.  Reads don't map to the amplicon will be skipped.  A read pair is considered belonging the amplicon if the edges are less than int bp to the amplicon, and overlap fraction is at least float.  Default: `10:0.95`
- `-k 0/1`   
    Indicate whether to perform local realignment.  Default: `1` or yes.  Set to `0` to disable it.
- `-G Genome fasta`  
    The reference fasta.  Should be indexed (.fai).  Default to: `/ngs/reference_data/genomes/Hsapiens/hg19/seq/hg19.fa`
- `-R Region`  
    The region of interest.  In the format of chr:start-end.  If end is omitted, then a single position.  No BED is needed.
- `-d delimiter`  
    The delimiter for split `region_info`, default to tab `"\t"`
- `-n regular_expression`  
    The regular expression to extract sample name from bam filenames.  Default to: `/([^\/\._]+?)_[^\/]*.bam/`
- `-N string`   
    The sample name to be used directly.  Will overwrite `-n` option
- `-b string`   
    The indexed BAM file
- `-c INT`   
    The column for chromosome
- `-S INT`   
    The column for region start, e.g. gene start
- `-E INT`  
    The column for region end, e.g. gene end
- `-s INT`   
    The column for segment starts in the region, e.g. exon starts
- `-e INT`  
    The column for segment ends in the region, e.g. exon ends
- `-g INT`     
    The column for gene name, or segment annotation
- `-x INT`   
    The number of nucleotide to extend for each segment, default: `0`
- `-f double`   
    The threshold for allele frequency, default: `0.05` or `5%`
- `-r minimum reads`   
    The minimum # of variance reads, default `2`
- `-B INT`  
    The minimum # of reads to determine strand bias, default `2`
- `-Q INT`  
    If set, reads with mapping quality less than INT will be filtered and ignored
- `-q INT`   
    The phred score for a base to be considered a good call.  Default: 25 (for Illumina). For PGM, set it to ~15, as PGM tends to under estimate base quality.
- `-m INT`   
    If set, reads with mismatches more than `INT` will be filtered and ignored.  Gaps are not counted as mismatches. Valid only for bowtie2/TopHat or BWA aln followed by sampe.  BWA mem is calculated as NM - Indels.  Default: 8, or reads with more than 8 mismatches will not be used.
- `-T INT`  
    Trim bases after `[INT]` bases in the reads
- `-X INT`   
    Extension of bp to look for mismatches after insersion or deletion.  Default to 3 bp, or only calls when they're within 3 bp.
- `-P number`  
    The read position filter.  If the mean variants position is less that specified, it's considered false positive.  Default: 5
- `-Z double`  
    For downsampling fraction.  e.g. `0.7` means roughly `70%` downsampling.  Default: No downsampling.  Use with caution.  The downsampling will be random and non-reproducible.
- `-o Qratio`  
    The `Qratio` of `(good_quality_reads)/(bad_quality_reads+0.5)`.  The quality is defined by `-q` option.  Default: `1.5`
- `-O MapQ`  
    The reads should have at least mean `MapQ` to be considered a valid variant.  Default: no filtering
- `-V freq`  
    The lowest frequency in normal sample allowed for a putative somatic mutations.  Default to `0.05`
- `-I INT`  
    The indel size.  Default: 120bp
- `-th threads`  
    Threads count. If omitted, number of threads is equal to number of processor cores.

<a name="outputColumns">
##Output columns
</a>

1. Sample - sample name
2. Gene - gene name from BED file
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
29. NM - average number of mismatches for reads containing variant
30. HICNT - number of high-quality reads with the variant
31. HICOV - position coverage by high quality reads
21. 5pFlankSeq - neighboring reference sequence to 5' end 
22. 3pFlankSeq - neighboring reference sequence to 3' end
23. SEGMENT:CHR_START_END - position description
24. VARTYPE - variant type