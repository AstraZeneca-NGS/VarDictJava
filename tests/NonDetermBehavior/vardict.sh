#!/bin/bash
REF_GENOME="/ngs/fa/hg19.fa"
COLO20="Colo829.chr20.PerlNonDeterm.bam"
REGION="-R chr20:61659400-61659600"
#Flag to remove differences because Perl doesn't sort unmapped reads
UNMAPPED_FLAG="-F 0x504"
# Paths to Vardict

VARDICTJAVA_HOME="$HOME/IdeaProjects/VarDictJava"
VARDICTJAVA="$VARDICTJAVA_HOME/build/install/VarDict/bin/VarDict"
VARDICTPERL_HOME="$HOME/IdeaProjects/VarDict/vardict.pl"

echo Running VarDict perl for 5 times
for i in {1..5}; 
do 
echo Running VarDict perl $i time;
time $VARDICTPERL_HOME \
		-G $REF_GENOME \
		-f 0.001 -N abc \
		-b $COLO20 $UNMAPPED_FLAG \
		-c 1 -S 2 -E 3 -g 4 \
		$REGION | sort  > vardictColo20_$i.perl.txt;
done;

echo Running VarDict java
time $VARDICTJAVA \
		-G $REF_GENOME \
		-f 0.001 -N abc \
		-b $COLO20 $UNMAPPED_FLAG \
		-c 1 -S 2 -E 3 -g 4 \
		$REGION | sort  > vardictColo20.java.txt

echo Differences:
for i in {1..4}; 
do 
let j=$i+1;
echo Difference $i-$j time;
diff vardictColo20_$i.perl.txt vardictColo20_$j.perl.txt;
done;