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

VARDICTPERL_HOME="$VARDICTJAVA_PROJECT/VarDict"
VARDICTPERL="$VARDICTPERL_HOME/vardict"
VARDICTPERL_R_PAIRED="$VARDICTPERL_HOME/testsomatic.R"
VARDICTPERL_VAR_PAIRED="$VARDICTPERL_HOME/var2vcf_paired.R"

TESTS_DIR="$VARDICTJAVA_PROJECT/tests/integration"

# Parameters for Vardict
JAVA_THREADS=8
PARAMETERS="-c 1 -S 2 -E 3 -g 4 -f 0.001 -N abc"

# Multiallelic confirmed variants that aren't supported by Perl
CONFIRMED_DIFFERENCES=\
"20\t24993259\t24993259\tG\tA|"\
"20\t31682933\t31682933\tG\tT|"\
"20\t31829197\t31829197\tC\tT|"\
"20\t44669022\t44669024\tTCC\tCTT|"\
"20\t36869637\t36869637\tC\tA|"\
"20\t50244187\t50244188\tGA\tAC"

# File names and paths
DIR_INPUT="$TESTS_DIR/input"
DIR_OUTPUT="$TESTS_DIR/output"

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
	-th $JAVA_THREADS \
	-b "$TUMOR_BAM_PATH|$NORMAL_BAM_PATH" \
	$BED_SPLIT_PATH  > java.var

#-F 0x504 flag can be deleted after Perl fix for filter unmapped reads by default
echo Running VarDict perl
time $VARDICTPERL \
	-G $FASTA_PATH \
	$PARAMETERS \
	-b "$TUMOR_BAM_PATH|$NORMAL_BAM_PATH" \
	-F 0x504 \
	$BED_SPLIT_PATH > perl.var

# Check if var files aren't empty
if [ ! -s "perl.var" ] || [ ! -s "java.var" ]; then
	echo "	Var files are empty!"
	exit 1;
fi

cat java.var | sort > java_sorted.var
cat perl.var | sort > perl_sorted.var

# Compare differences
echo Compare differences raw VARs perl and java versions
cat java_sorted.var | grep -Pv $CONFIRMED_DIFFERENCES > java_confirmed.var

diff_var=$(diff perl_sorted.var  java_confirmed.var > diff_var.txt)
ret1=$?
if [ "$ret1" = "0" ]; then
	echo "	Raw VAR diff OK (no differences)";
else
	echo "	Raw VAR files have differences!"
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
