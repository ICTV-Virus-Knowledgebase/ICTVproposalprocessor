#!/usr/bin/env bash
#
# run and evaluate regression tests
#
#
TEST_DIR=testData
REPORT=regression_test.summary.txt
date > $REPORT
#
# scan for test directories
#
PAT="*"
if [ ! -z "$1" ]; then PAT="*$1*"; fi
TESTS=$(cd $TEST_DIR; find . -type d -name "$PAT" -name "proposal*" \! -name "*result*" -depth 1 -exec basename {} +)
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
    RESULTS=${DEST_DIR}/QC.regression.new.tsv
    RESULTSBASE=${DEST_DIR}/QC.regression.tsv
    RESULTSDIFF=${DEST_DIR}/QC.regression.diff
    LOG=${DEST_DIR}/log.new.txt
    LOGBASE=${DEST_DIR}/log.txt
    LOGDIFF=${DEST_DIR}/log.diff

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
    echo RESULTSBASE=$RESULTSBASE
    echo LOG=$LOG

    #
    # run script
    #
    echo "#" \
        Rscript merge_proposal_zips.R \
	    --proposalsDir=$SRC_DIR \
	    --outDir=$DEST_DIR \
	    --qcTsvRegression=$(basename $RESULTS) \
	    2>&1 | tee $LOG
    Rscript merge_proposal_zips.R \
	    --proposalsDir=$SRC_DIR \
	    --outDir=$DEST_DIR \
	    --qcTsvRegression=$(basename $RESULTS) \
	    2>&1 >> $LOG

    #
    # check output
    #
    echo "diff -yw -W 200 \<(cut -f 5- $RESULTS) \<(cut -f 5- $RESULTSBASE) \> $RESULTSDIFF" | tee $RESULTSDIFF
    diff -yw -W 200 <(cut -f 5- $RESULTS) <(cut -f 5- $RESULTSBASE) 2>&1 >> $RESULTSDIFF; RC=$?
    if [ $RC -eq "0" ]; then
	echo "PASS  $TEST" | tee -a $REPORT
    else
	echo "FAIL$RC $TEST" | tee -a $REPORT
    fi	
    #
    # check log
    #
    echo diff -yw -W 200 $LOG $LOGBASE \> $LOGDIFF | tee $LOGDIFF
    diff -yw -W 200 $LOG $LOGBASE 2>&1 >> $LOGDIFF; RC=$?
    if [ $RC -eq "0" ]; then
	echo "LOG_PASS  $TEST" | tee -a $REPORT
    else
	echo "LOG_FAIL$RC $TEST" | tee -a $REPORT
    fi
    echo "#-------------------------" | tee -a $REPORT
	
done
echo "#########################################"
echo "############### SUMMARY ################# "
echo "#########################################"
cat $REPORT
   
