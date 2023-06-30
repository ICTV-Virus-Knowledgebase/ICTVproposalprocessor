#!/usr/bin/env bash
#
# run and evaluate regression tests
#
#
TEST_DIR=testData

#
# scan for test directories
#
TESTS=$(cd $TEST_DIR; find . -type d -name "proposal*" \! -name "*result*" -depth 1 -exec basename {} +)
#TESTS=proposal2020

#
# iterate
#
for TEST in $TESTS; do
    #
    # input/output for script
    #
    SRC_DIR=$TEST_DIR/$TEST
    DEST_DIR=${TEST_DIR}/results/$TEST
    REGRESS=${DEST_DIR}/QC.regression.tsv
    OUT=${DEST_DIR}/log.txt

    mkdir -p $DEST_DIR
    #
    # header
    #
    echo "#########################################"
    echo "###### $TEST ######"
    echo "#########################################"
    echo SRC_DIR=$SRC_DIR
    echo DEST_DIR=$DEST_DIR
    echo REGRESS=$REGRESS
    echo OUT=$OUT

    #
    # run script
    #
    echo "#" \
        Rscript merge_proposal_zips.R \
	    --proposalsDir=$SRC_DIR \
	    --outDir=$DEST_DIR \
	    --qcTsvRegression=$(basename $REGRESS) \
	    2>&1 > $OUT
    Rscript merge_proposal_zips.R \
	    --proposalsDir=$SRC_DIR \
	    --outDir=$DEST_DIR \
	    --qcTsvRegression=$(basename $REGRESS) \
	    2>&1 >> $OUT

    #
    # check output
    #
    git status -s $REGRESS
done

   
