#!/bin/bash -eu
set -o pipefail

#---
# Paths
#---
# Validation of VCF files with EBIvariation vcf-validator, must be on the PATH
VCF_VALIDATOR="vcf_validator"

WORKSPACE="$HOME/IdeaProjects"
WORKSPACE="$HOME/workspace"
VARDICTJAVA_PROJECT="$WORKSPACE/VarDictJava"
VARDICTJAVA_HOME="$VARDICTJAVA_PROJECT/build/install/VarDict"
VARDICTJAVA="$VARDICTJAVA_HOME/bin/VarDict"

TESTS_DIR="$VARDICTJAVA_PROJECT/tests/integration"

VARDICTPERL_HOME="$WORKSPACE/VarDict"
VARDICTPERL="$VARDICTPERL_HOME/vardict.pl"
VARDICTPERL_R_SIMPLE="$VARDICTPERL_HOME/teststrandbias.R"
VARDICTPERL_R_PAIRED="$VARDICTPERL_HOME/testsomatic.R"
VARDICTPERL_VAR_SIMPLE="$VARDICTPERL_HOME/var2vcf_valid.pl"
VARDICTPERL_VAR_PAIRED="$VARDICTPERL_HOME/var2vcf_paired.pl"

# File names and paths
DIR_INPUT="$TESTS_DIR/input"
DIR_RAW_INPUT="$TESTS_DIR/raw_input"
DIR_OUTPUT="$TESTS_DIR/output"
DIR_EXPECTED="$TESTS_DIR/expected"
DIR_REFERENCE="$TESTS_DIR/reference"

