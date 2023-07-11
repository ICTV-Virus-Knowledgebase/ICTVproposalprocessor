#!/usr/bin/env bash
#
# test run our docker image
#
TEST=proposalsTest3_binomial
TEST=proposalsTest_createNew
if [[ ! -z "$1" && "$1" != -* ]]; then TEST="$1"; shift; fi
echo TEST=$TEST
TEST_DIR=./testResultsDocker/$TEST
echo "# cleaning out $TEST_DIR/..."
echo "find testResultsDocker/$TEST -name '*new*' -o -name '*diff*' -exec rm {} +"
find testResultsDocker/$TEST -name '*new*' -o -name '*diff*' -exec rm {} +

sudo docker run -it \
	-v "$(pwd)/testData:/testData" \
	-v "$(pwd)/testResultsDocker:/testResults" \
	ictv_proposal_processor  \
	/merge_proposal_zips.R \
	--proposalsDir=/testData/$TEST \
	--outDir=/testResultsDocker/$TEST \
	--qcTsvRegression=QC.regression.new.tsv \
	-v \
	$*
