#!/bin/bash -eu
set -o pipefail

#---
# Paths to VarDict java and perl. Change VARDICTJAVA_HOME and VARDICTPERL_HOME to your paths
#---
WORKSPACE="$HOME/IdeaProjects"
WORKSPACE="$HOME/workspace"
VARDICTJAVA_PROJECT="$WORKSPACE/VarDictJava"
VARDICTJAVA_HOME="$VARDICTJAVA_PROJECT/build/install/VarDict"
VARDICTJAVA="$VARDICTJAVA_HOME/bin/VarDict"

TESTS_DIR="$VARDICTJAVA_PROJECT/tests/integration"

VARDICTPERL_HOME="$VARDICTJAVA_PROJECT/VarDict"
VARDICTPERL="$TESTS_DIR/vardict.pl"
VARDICTPERL_R_PAIRED="$VARDICTPERL_HOME/testsomatic.R"
VARDICTPERL_VAR_PAIRED="$VARDICTPERL_HOME/var2vcf_paired.R"


# Parameters for Vardict
JAVA_THREADS=8
PARAMETERS="-c 1 -S 2 -E 3 -g 4 -f 0.001 -N abc"

# Multiallelic confirmed variants that aren't supported by Perl
CONFIRMED_DIFFERENCES_FILE="$TESTS_DIR/confirmed_differences.txt"

# File names and paths
DIR_INPUT="$TESTS_DIR/input"
DIR_OUTPUT="$TESTS_DIR/output"
VARDICT_OUT_JAVA="$DIR_OUTPUT/vardict.java.txt"
VARDICT_OUT_PERL="$DIR_OUTPUT/vardict.perl.txt"
VARDICT_OUT_SORT_JAVA="$DIR_OUTPUT/vardict.java.sort.txt"
VARDICT_OUT_SORT_PERL="$DIR_OUTPUT/vardict.perl.sort.txt"

FASTA_URL="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa"
FASTA=$(basename $FASTA_URL)
FASTA_BASE=$(basename $FASTA_URL .fa)
FASTA_PATH="$DIR_INPUT/$FASTA"

NORMAL_BAM_URL="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA12878/exome_alignment/NA12878.chrom20.ILLUMINA.bwa.CEU.exome.20121211.bam"
NORMAL_BAM=$(echo $NORMAL_BAM_URL | sed 's#.*/##')
NORMAL_BAM_PATH="$DIR_INPUT/$NORMAL_BAM"

TUMOR_BAM_URL="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA12889/exome_alignment/NA12889.chrom20.ILLUMINA.bwa.CEU.exome.20121211.bam"
TUMOR_BAM=$(echo $TUMOR_BAM_URL | sed 's#.*/##')
TUMOR_BAM_PATH="$DIR_INPUT/$TUMOR_BAM"

BED_URL="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/exome_pull_down_targets_phases1_and_2/20120518.consensus.annotation.bed"
BED=$(basename $BED_URL)
CHR="chr20"
BED_SPLIT=$BED.$CHR
BED_PATH="$DIR_INPUT/$BED"
BED_SPLIT_PATH="$DIR_INPUT/$BED_SPLIT"

#---
# Download data files
#---
echo Creating input directory
mkdir $DIR_INPUT || true
cd $DIR_INPUT

# Fasta downloading
echo Downloading fasta file
if [ ! -f "$FASTA_PATH" ]; then
	wget -nc $FASTA_URL.gz -O $FASTA_PATH.gz || true
	echo Unzipping fasta file
	gunzip -f $FASTA_PATH.gz
	echo Indexing FAST file
	samtools faidx $FASTA_PATH
fi

# BAM and BAI download
# Normal
echo Downloading normal BAM
wget -nc $NORMAL_BAM_URL -O $NORMAL_BAM_PATH || true
wget -nc $NORMAL_BAM_URL.bai -O $NORMAL_BAM_PATH.bai || true
# Tumor
echo Downloading tumor BAM
wget -nc $TUMOR_BAM_URL -O $TUMOR_BAM_PATH || true
wget -nc $TUMOR_BAM_URL.bai -O $TUMOR_BAM_PATH.bai || true

# BED download
echo Downloading BED
if [ ! -f "$BED_PATH" ]; then
	wget -nc $BED_URL -O $BED_PATH || true

	# Splitting BED file on chr20 for better performance in VarDict Perl
	cat $BED | grep $CHR > $BED_SPLIT
fi

#---
# Run both vardict versions and compare outputs
#---

echo Creating output directory
cd ..
mkdir $DIR_OUTPUT || true
cd $DIR_OUTPUT

# Run VarDict
echo Running VarDict java
time $VARDICTJAVA \
		-G $FASTA_PATH \
		$PARAMETERS \
		-b "$TUMOR_BAM_PATH|$NORMAL_BAM_PATH" \
		-th $JAVA_THREADS \
		$BED_SPLIT_PATH \
	| tee $VARDICT_OUT_JAVA \
	| sort \
	| grep -v -f $CONFIRMED_DIFFERENCES_FILE \
	> $VARDICT_OUT_SORT_JAVA

# Note: The '-F 0x504' flag can be deleted after Perl fix for filter unmapped reads by default
echo Running VarDict perl
time $VARDICTPERL \
		-G $FASTA_PATH \
		$PARAMETERS \
		-b "$TUMOR_BAM_PATH|$NORMAL_BAM_PATH" \
		-F 0x504 \
		$BED_SPLIT_PATH \
	| tee $VARDICT_OUT_PERL \
	| sort \
	| grep -v -f $CONFIRMED_DIFFERENCES_FILE \
	> $VARDICT_OUT_SORT_PERL

# Check if var files aren't empty
if [ ! -s "$VARDICT_OUT_PERL" ] || [ ! -s "$VARDICT_OUT_JAVA" ]; then
	echo "ERROR: Empty variant file/s"
	exit 1;
fi

# Compare differences
echo Compare differences 

DIFF_FILE="$DIR_OUTPUT/vardict.diff"
if diff $VARDICT_OUT_SORT_PERL $VARDICT_OUT_SORT_JAVA > $DIFF_FILE
then
	echo "OK: Raw VAR diff OK (no differences)";
else
	echo "ERROR: Raw VAR files have differences!"
	exit 1;
fi

#This part can be uncommented when .R and .pl scripts in vardict repositories will be updated.
#echo Running R script
#cat java_sorted.var | $VARDICTPERL_R_PAIRED > java_r.var
#cat perl_sorted.var | $VARDICTPERL_R_PAIRED > perl_r.var

#if [ ! -s "perl_r.var" ] || [ ! -s "java_r.var" ]; then 
#	echo "	Var files after R script are empty!" 
#	exit 1;
#fi
#echo Running Var2VCF script
#cat java_r.var | $VARDICTPERL_VAR_PAIRED > java.vcf
#cat perl_r.var | $VARDICTPERL_VAR_PAIRED > perl.vcf

#echo Running differences VCFs perl and java
#diff_vcf=$(diff perl.vcf java.vcf > diff_vcf.txt)
#ret2=$?
#if ["$ret2" = "0"]; then
#	echo "	VCF diff OK (no differences)";
#else 
#	echo "	VCF files have differences!"
#	exit 1;
#fi
