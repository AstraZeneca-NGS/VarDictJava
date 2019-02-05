#!/bin/bash -eu
set -o pipefail

# Load config vars
SCRIPT_DIR="$( cd "$( dirname "$0" )" && pwd )"
source $SCRIPT_DIR/config.sh

REF_GENOME="$DIR_REFERENCE/hg19.fa"
COLO20="$TESTS_DIR/non_deterministic_behavior/Colo829.chr20.PerlNonDeterm.bam"
REGION="-R chr20:61659400-61659600"

N=5
for i in $(seq 1 $N); do 
	echo "Running VarDict Perl: $i / $N"
	time $VARDICTPERL \
			-G $REF_GENOME \
			-f 0.001 \
			-N abc \
			-b $COLO20 \
			$REGION \
		| sort \
		> $DIR_OUTPUT/vardictColo20.perl.$i.txt

 	echo "Running VarDict Java: $i / $N"
 	time $VARDICTJAVA \
 			-G $REF_GENOME \
 			-f 0.001 \
			-N abc \
 			-b $COLO20 \
 			$REGION \
 	| sort \
 	> $DIR_OUTPUT/vardictColo20.java.$i.txt

done

echo "Differences:"
for i in $(seq 1 $N); do 
	for j in $(seq 1 $N); do 
		echo "--------------------------------------------------------------------------------"
		echo
		echo
		echo "Difference $i-$j Perl:"
		diff $DIR_OUTPUT/vardictColo20.perl.$i.txt $DIR_OUTPUT/vardictColo20.perl.$j.txt || true
		echo
		echo
		echo "Difference $i-$j Java:"
		diff $DIR_OUTPUT/vardictColo20.java.$i.txt $DIR_OUTPUT/vardictColo20.java.$j.txt || true
	done
done
