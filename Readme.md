#VarDict-java

This is a Java port of [VarDict variant caller](https://github.com/AstraZeneca-NGS/VarDict).

##Requirements
1. JDK 1.7 or later
2. Samtools (must be in `$PATH`)
2. Internet connection to download dependencies using gradle.

##Getting started
- To build:
```
cd path_to_vardict_folder
./gradlew clean installApp 
```

- Running in single sample mode:  
```
AF_THR="0.01" # minimum allele frequency
<path_to_vardict_folder>/bin/VarDict -G /path/to/hg19.fa -f 0.01 -N sample_name -b /path/to/my.bam -z -c 1 -S 2 -E 3 -g 4 /path/to/my.bed | teststrandbias.R | var2vcf_valid.pl -N sample_name -E -f 0.01
```

- Paired variant calling:
```
AF_THR="0.01" # minimum allele frequency
<path_to_vardict_folder>/bin/VarDict -G /path/to/hg19.fa -f $AF_THR -N tumor_sample_name -b "/path/to/tumor.bam|/path/to/normal.bam" -z -F -c 1 -S 2 -E 3 -g 4 /path/to/my.bed | testsomatic.R | var2vcf_somatic.pl -N "tumor_sample_name|normal_sample_name" -f $AF_THR
```


##Options

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

