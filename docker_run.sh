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
echo ""

#
# target docker container naem
#
NAME_FILE=docker_container_name.txt
if [ ! -e "$NAME_FILE" ]; then 
	echo "ERROR: $NAME_FILE not found"
	exit 1
fi
echo READING $NAME_FILE
DOCKER_CONTAINER_NAME=$(egrep -v "^#" $NAME_FILE| sed 's/[ \t\r\n]*//g')
echo DOCKER_CONTAINER_NAME=$DOCKER_CONTAINER_NAME
echo ""

sudo docker run -it \
	-v "$(pwd)/testData:/testData" \
	-v "$(pwd)/testResults:/testResults" \
	${DOCKER_CONTAINER_NAME}  \
	/merge_proposal_zips.R \
	--proposalsDir=/$TEST_DIR \
	--outDir=/$OUT_DIR \
	--qcTsvRegression=QC.regression.new.tsv \
	-v \
	$*
