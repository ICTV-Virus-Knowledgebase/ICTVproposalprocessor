#!/usr/bin/env bash
#
# test run our docker image
#
TEST=proposalsTest3_binomial
TEST=msl38/proposalsTest_createNew
if [[ ! -z "$1" && "$1" != -* ]]; then TEST="$1"; shift; fi
echo TEST=$TEST
TEST_DIR=testData/$TEST
OUT_DIR=testResults/$TEST
echo "# cleaning out $OUT_DIR/..."
echo "find ./$OUT_DIR -name '*new*' -o -name '*diff*' -exec rm {} +"
find ./$OUT_DIR -name '*new*' -o -name '*diff*' -exec rm {} +

sudo docker run -it \
	-v "$(pwd)/testData:/testData" \
	-v "$(pwd)/testResults:/testResults" \
	ictv_proposal_processor  \
	/merge_proposal_zips.R \
	--proposalsDir=/$TEST_DIR \
	--outDir=/$OUT_DIR \
	--qcTsvRegression=QC.regression.new.tsv \
	-v \
	$*
