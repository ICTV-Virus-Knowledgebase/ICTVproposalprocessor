#!/usr/bin/env bash
#
# run and evaluate regression tests
# IN DOCKER
#
VER=""
if [ ! -z "$1" ]; then VER=":$1"; fi
TEST_DIR=testData
REPORT=regression_test.summary.txt
date > $REPORT
#
# scan for test directories
#
TESTS=$(cd $TEST_DIR; find . -maxdepth 1 -type d -name "proposal*" \! -name "*result*" -exec basename {} \;)
#TESTS=proposal2020

# 
# update docker image, just incase
#
echo "# Building docker image"
echo ./docker_build_image.sh
./docker_build_image.sh

#
# iterate
#
for TEST in $TESTS; do
    #
    # input/output for script
    #
    SRC_DIR=$TEST_DIR/$TEST
    DEST_DIR=${TEST_DIR}/results_docker/$TEST
    RESULTS=${DEST_DIR}/QC.regression.new.tsv
    BASELINE=${DEST_DIR}/QC.regression.tsv
    OUT=${DEST_DIR}/log.txt

    mkdir -p $DEST_DIR
    #
    # header
    #
    echo "#########################################"
    echo "###### $TEST "
    echo "#########################################"
    echo SRC_DIR=$SRC_DIR
    echo DEST_DIR=$DEST_DIR
    echo RESULTS=$RESULTS
    echo BASELINE=$BASELINE
    echo OUT=$OUT

    #
    # run script
    #
    echo "#" \
    sudo docker run -it \
	    -v "$(pwd)/testData:/testData":ro \
	    -v "$(pwd)/testData/results_docker:/testData/results_docker":rw \
            ictv_proposal_processor$VER  \
	    /merge_proposal_zips.R \
	    --proposalsDir=$SRC_DIR \
	    --outDir=$DEST_DIR \
	    --qcTsvRegression=$(basename $RESULTS) \
	    2>&1 >> $OUT
    sudo docker run -it \
	    -v "$(pwd)/testData:/testData":ro \
	    -v "$(pwd)/testData/results_docker:/testData/results_docker":rw \
            ictv_proposal_processor$VER  \
	    /merge_proposal_zips.R \
	    --proposalsDir=$SRC_DIR \
	    --outDir=$DEST_DIR \
	    --qcTsvRegression=$(basename $RESULTS) \
	    2>&1 > $OUT

    #
    # check output
    #    echo diff -q -w $RESULTS $BASELINE
    diff -q -w $RESULTS $BASELINE 2>&1 > /dev/null; RC=$?
    if [ $RC -eq "0" ]; then
	echo "PASS  $TEST" | tee -a $REPORT
    else
	echo "FAIL$RC $TEST" | tee -a $REPORT
    fi	
    
done
echo "#########################################"
echo "############### SUMMARY ################# "
echo "#########################################"
cat $REPORT
   
