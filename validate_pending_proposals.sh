#!/usr/bin/env bash
#
# QC proposals for new MSL (in ./pending)
#
# USAGE: ./validate_pending_proposals.sh [test_pattern] [container_version]
#
# On Linux, this runs the docker container
# On MacOS, this runs R directly. 
#

# which tests to run 
TEST_PAT="*revised"
#TEST_PAT="2023.017P*.xlsx"
echo TEST_PAT=$TEST_PAT

# pass-through args
MSL_NOTES="DRAFT: EC 56, Bari, Italy, August 2024; Email ratification February 2025 (MSL #40)"
SCRIPT_ARGS='--msl --mode=draft --newMslName=2024 '
if [ ! -z "$1" ]; then SCRIPT_ARGS="$SCRIPT_ARGS $*"; fi
echo SCRIPT_ARGS=$SCRIPT_ARGS

# which docker container to run
if [ "$(uname)" == "Linux" ]; then 
	CONTAINER=ictv_proposal_processor
	if [ ! -z "$1" ]; then CONTAINER="curtish/${CONTAINER}:$1"; shift; fi
	echo "CONTAINER=$CONTAINER"

	# 
	# update docker image, just incase
	#
	echo "# Building docker image"
	echo ./docker_build_image.sh
	./docker_build_image.sh
fi

#
# test cases location
# 
TEST_DIR=pending
echo TEST_DIR=$TEST_DIR
RESULTS_DIR=pending_results
echo RESULTS_DIR=$RESULTS_DIR

REPORT=pending.summary.txt
echo REPORT=$REPORT
(date; hostname) > $REPORT

#
# scan for test directories
#
#echo "#$ find $TEST_DIR -type d -name "$TEST_PAT" -exec basename {} \;"
#TESTS=$(find $TEST_DIR -type d -name "$TEST_PAT" -exec basename {} \;)
echo "#$ find $TEST_DIR  -name \"$TEST_PAT\" -exec basename {} \;"
TESTS=$(find $TEST_DIR -name "$TEST_PAT" -exec basename {} \;)
echo TESTS=$TEST

#
# iterate
#
for TEST in $TESTS; do
    #
    # input/output for script
    #
    SRC_DIR=$TEST_DIR/$TEST
    DEST_DIR=${RESULTS_DIR}/$TEST
    RESULTS=${DEST_DIR}/QC.regression.new.tsv
    RESULTSBASE=${DEST_DIR}/QC.regression.tsv
    RESULTSDIFF=${DEST_DIR}/QC.regression.diff
    RESULTSDWDIFF=${DEST_DIR}/QC.regression.dwdiff
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
    if [ -z "$CONTAINER" ]; then 
	   # update code/git version string
           ./version_git.sh
	    echo "#" \
	        Rscript merge_proposal_zips.R \
		    --proposalsDir=$SRC_DIR \
		    --outDir=$DEST_DIR \
		    --qcTsvRegression=$(basename $RESULTS) \
		    --qcTsvRegression=$(basename $RESULTS) \
		    --newMslNotes="$MSL_NOTES" \
		    $SCRIPT_ARGS \
		    2>&1 | tee $LOG
	    Rscript merge_proposal_zips.R \
		    --proposalsDir=$SRC_DIR \
		    --outDir=$DEST_DIR \
		    --qcTsvRegression=$(basename $RESULTS) \
		    --newMslNotes="$MSL_NOTES" \
		    $SCRIPT_ARGS \
		    2>&1 | tee -a $LOG
    else
	   echo "#" \
		sudo docker run -it \
		    -v "$(pwd)/$TEST_DIR:/testData":ro \
		    -v "$(pwd)/$RESULTS_DIR:/testResults":rw \
	            $CONTAINER  \
		    /merge_proposal_zips.R \
		    --proposalsDir=/testData/$TEST \
		    --outDir="/testResults/$TEST" \
		    --qcTsvRegression=$(basename $RESULTS) \
		    --newMslNotes="$MSL_NOTES" \
		    $SCRIPT_ARGS \
		    2>&1 | tee $LOG
	    sudo docker run -it \
		    -v "$(pwd)/$TEST_DIR:/testData":ro \
		    -v "$(pwd)/$RESULTS_DIR:/testResults":rw \
	            $CONTAINER  \
		    /merge_proposal_zips.R \
		    --proposalsDir=/testData/$TEST \
		    --outDir="/testResults/$TEST" \
		    --qcTsvRegression=$(basename $RESULTS) \
		    --newMslNotes="$MSL_NOTES" \
		    $SCRIPT_ARGS \
		    2>&1 >> $LOG
    fi	

    #
    # check output
    #
    echo "dwdiff --color  <(cut -f 5- $RESULTSBASE) <(cut -f 5- $RESULTS) #> $RESULTSDWDIFF" | tee $RESULTSDWDIFF
    dwdiff --color <(cut -f 5- $RESULTSBASE) <(cut -f 5- $RESULTS) 2>&1 >> $RESULTSDWDIFF; RC=$?
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
    # use "tail -n +3" to skip date/version/etc in first 2 lines
    #
    echo "diff -yw -W 200 \<(tail -n +3 $LOG) \<(tail -n +3 $LOGBASE) \> $LOGDIFF" | tee $LOGDIFF
    diff -yw -W 200 <(tail -n +3 $LOG) <(tail -n +3 $LOGBASE) 2>&1 >> $LOGDIFF; RC=$?
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
   
