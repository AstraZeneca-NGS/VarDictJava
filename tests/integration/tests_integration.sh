#!/bin/bash -eu
set -o pipefail

# Needs EBIvariation vcf-validator to be installed on path to validate VCF files
#---
# Config variables
#---

SCRIPT_DIR="$( cd "$( dirname "$0" )" && pwd )"
source $SCRIPT_DIR/config.sh

JAVA_THREADS=8
PARAMETERS="-c 1 -S 2 -E 3 -g 4 -f 0.001 -N abc"

# Output files
VARDICT_OUT_JAVA="$DIR_OUTPUT/vardict.java.txt"
VARDICT_OUT_JAVA_FISHER="$DIR_OUTPUT/vardict.java.fisher.txt"
VARDICT_OUT_PERL="$DIR_OUTPUT/vardict.perl.txt"
VARDICT_OUT_SORT_JAVA="$DIR_OUTPUT/vardict.java.sort.txt"
VARDICT_OUT_SORT_JAVA_FISHER="$DIR_OUTPUT/vardict.java.sort.fisher.txt"
VARDICT_OUT_R_JAVA="$DIR_OUTPUT/vardict.java.r.txt"
VARDICT_OUT_VCF_JAVA="$DIR_OUTPUT/vardict.java.vcf"
VARDICT_OUT_SORT_PERL="$DIR_OUTPUT/vardict.perl.sort.txt"
VARDICT_OUT_R_PERL="$DIR_OUTPUT/vardict.perl.r.txt"
VARDICT_OUT_VCF_PERL="$DIR_OUTPUT/vardict.perl.vcf"

# Download URLs and local file paths
FASTA_URL="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa"
FASTA=$(basename $FASTA_URL)
FASTA_BASE=$(basename $FASTA_URL .fa)
FASTA_PATH="$DIR_REFERENCE/$FASTA"

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
function download_files {
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
    cd ..
}


#---
# Compare expected and actual results.
# Usage: diff_results test_name expected actual
#---
function diff_results {
    TEST_NAME="$1"
    EXPECTED="$2"
    RESULTS="$3"

    DIFF_FILE="$DIR_OUTPUT/diff.$TEST_NAME.vcf"
    echo "Test $TEST_NAME: Compare results, creating diff file '$DIFF_FILE'"
    if diff <(sort $EXPECTED) <(sort $RESULTS) > $DIFF_FILE
    then
	    echo "Test $TEST_NAME: OK";
    else
	    echo "Test $TEST_NAME: ERROR";
	exit 1;
    fi
}

#Only for fisher test results from Java
function diff_results_with_epsilon {
    TEST_NAME="$1"
    EXPECTED="$2"
    RESULTS="$3"

    DIFF_FILE="$DIR_OUTPUT/diff.$TEST_NAME.var"
    DIFF_FILE_FISHER="$DIR_OUTPUT/diff.fisher.$TEST_NAME.var"
    EXP_COLUMNS="$DIR_OUTPUT/exp.columns.$TEST_NAME.var"
    RES_COLUMNS="$DIR_OUTPUT/res.columns.$TEST_NAME.var"
	COMBINED="$DIR_OUTPUT/combined.$TEST_NAME.var"

    echo "Test $TEST_NAME: Compare results without fisher test, creating diff file '$DIFF_FILE'"

	# Check first all columns except pvalue and oddratio columns, if there are differences, stop with error
	# Columns in paired mode to exclude:
	awk '{$26=""; $27=""; $46=""; $47=""; $60=""; $61=""; print $0}' $EXPECTED > $EXP_COLUMNS
	awk '{$26=""; $27=""; $46=""; $47=""; $60=""; $61=""; print $0}' $RESULTS > $RES_COLUMNS
    if diff <(sort $EXP_COLUMNS) <(sort $RES_COLUMNS) > $DIFF_FILE
    then
	    echo "Test without fisher columns $TEST_NAME: OK";
    else
	    echo "Test without fisher columns  $TEST_NAME: ERROR";
	    exit 1;
	fi

	# Check epsilon between oddratio and pvalue columns (must be leass then 1E-4)
	# Combine files to easier compare in awk
	paste $EXPECTED $RESULTS > $COMBINED

    # Check if differences between columns from both versions aren't more than epsilon
	awk 'function abs(x) {return ((x < 0.0)?-x : x)} \
	     eps=1E-4 \
	     {if((abs($26-$87)>eps) || (abs($27-$88)>eps) || (abs($46-$107)>eps) || \
	     (abs($47-$108)>eps) || (abs($60-$121)>eps) || (abs($61 - $122) > eps)) print;}' \
	     $COMBINED > $DIFF_FILE_FISHER

	if [ ! -s "$DIFF_FILE_FISHER" ]
    then
	    echo "Test with fisher columns $TEST_NAME: OK";
    else
	    echo "Test with fisher columns  $TEST_NAME: ERROR";
	    exit 1;
	fi
}


#---
# Run both vardict versions and compare outputs
#---
function run_perl_java {
    echo Creating output directory
    mkdir $DIR_OUTPUT || true

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
        > $VARDICT_OUT_SORT_JAVA

    echo Running VarDict java with fisher test
	time $VARDICTJAVA \
		    -G $FASTA_PATH \
		    $PARAMETERS \
		    --fisher \
		    -th $JAVA_THREADS \
		    -b "$TUMOR_BAM_PATH|$NORMAL_BAM_PATH" \
		    $BED_SPLIT_PATH \
		| tee $VARDICT_OUT_JAVA_FISHER \
        | sort \
        > $VARDICT_OUT_SORT_JAVA_FISHER

    echo Running VarDict perl
    time $VARDICTPERL \
            -G $FASTA_PATH \
            $PARAMETERS \
            -b "$TUMOR_BAM_PATH|$NORMAL_BAM_PATH" \
            $BED_SPLIT_PATH \
        | tee $VARDICT_OUT_PERL \
        | sort \
        > $VARDICT_OUT_SORT_PERL

    # Check if var files aren't empty
    if [ ! -s "$VARDICT_OUT_PERL" ] || [ ! -s "$VARDICT_OUT_JAVA" ]; then
        echo "ERROR: Empty variant file/s"
        exit 1;
    fi

    # Compare differences
    diff_results "VAR_JAVA_PERL" $VARDICT_OUT_SORT_PERL $VARDICT_OUT_SORT_JAVA

    echo Running R script
    cat $VARDICT_OUT_SORT_JAVA | $VARDICTPERL_R_PAIRED > $VARDICT_OUT_R_JAVA
    cat $VARDICT_OUT_SORT_PERL | $VARDICTPERL_R_PAIRED > $VARDICT_OUT_R_PERL

    if [ ! -s "$VARDICT_OUT_R_PERL" ] || [ ! -s "$VARDICT_OUT_R_JAVA" ]; then
        echo "	ERROR: Empty R variant file/s"
        exit 1;
    fi

    echo Running differences perl R and java fisher
	diff_results_with_epsilon "VAR_JAVA_FISHER_PERL_R"  $VARDICT_OUT_R_PERL $VARDICT_OUT_SORT_JAVA_FISHER

    echo Running Var2VCF script
    cat $VARDICT_OUT_R_JAVA | $VARDICTPERL_VAR_PAIRED > $VARDICT_OUT_VCF_JAVA
    cat $VARDICT_OUT_R_PERL | $VARDICTPERL_VAR_PAIRED > $VARDICT_OUT_VCF_PERL
    diff_results "VCF_JAVA_PERL" $VARDICT_OUT_VCF_PERL $VARDICT_OUT_VCF_JAVA

    echo Running Var2VCF script opt_M
    cat $VARDICT_OUT_R_JAVA | $VARDICTPERL_VAR_PAIRED -M > $VARDICT_OUT_VCF_JAVA
    cat $VARDICT_OUT_R_PERL | $VARDICTPERL_VAR_PAIRED -M > $VARDICT_OUT_VCF_PERL
    diff_results "VCF_JAVA_PERL_M" $VARDICT_OUT_VCF_PERL $VARDICT_OUT_VCF_JAVA
}

#---
# Run VCF integration testing for simple and paired and compare outputs with expected
#---
function vcf_testing {
    echo Running VCF integration testing

    RAWVARDICT_SIMPLE="$DIR_RAW_INPUT/raw.vardict.simple.var"
    RAWVARDICT_PAIRED="$DIR_RAW_INPUT/raw.vardict.paired.var"
    VARDICT_VCF_SIMPLE="$DIR_OUTPUT/vardict.simple.vcf"
    VARDICT_VCF_PAIRED="$DIR_OUTPUT/vardict.paired.vcf"
    VARDICT_VCF_SIMPLE_EXP="$DIR_EXPECTED/vardict.simple.vcf"
    VARDICT_VCF_PAIRED_EXP="$DIR_EXPECTED/vardict.paired.vcf"

    echo Running simple mode VCF integration testing
    cat $RAWVARDICT_SIMPLE | $VARDICTPERL_R_SIMPLE | $VARDICTPERL_VAR_SIMPLE > $VARDICT_VCF_SIMPLE
    diff_results "VCF_simple" $VARDICT_VCF_SIMPLE_EXP  $VARDICT_VCF_SIMPLE

    echo Running paired mode VCF integration testing
    cat $RAWVARDICT_PAIRED | $VARDICTPERL_R_PAIRED | $VARDICTPERL_VAR_PAIRED > $VARDICT_VCF_PAIRED
    diff_results "VCF_paired" $VARDICT_VCF_PAIRED_EXP  $VARDICT_VCF_PAIRED

    echo Running validation of VCF files
    $VCF_VALIDATOR -i $VARDICT_VCF_SIMPLE
    $VCF_VALIDATOR -i $VARDICT_VCF_PAIRED
}

download_files
run_perl_java
vcf_testing