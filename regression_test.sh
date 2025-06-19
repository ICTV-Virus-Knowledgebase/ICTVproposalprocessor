#!/usr/bin/env bash
#
# run and evaluate regression tests
#
# USAGE: ./regression_test.sh [test_pattern] [container_version]
#
# On Linux, this runs the docker container
# On MacOS, this runs R directly. 
#
# TO DO
# add -msl and diff msl.tsv vs ref MSL and add that to git.
#
# which tests to run 
TEST_PAT="*"
if [ ! -z "$1" ]; then TEST_PAT="*$1*"; shift; fi
echo TEST_PAT=$TEST_PAT

# which docker container to run
if [ "$(uname)" == "Linux" ]; then 
	CONTAINER=ictv_proposal_processor
	if [ ! -z "$1" ]; then CONTAINER="curtish/${CONTAINER}:$1"; shift; fi
	echo "CONTAINER=$CONTAINER"

	# 
	# update docker image, just incase
	#
	echo "# Building docker image"
	#echo "# SKIP ./docker_build_image.sh"
	./docker_build_image.sh
fi

#
# make sure version_git.txt is built
#
if [ ! -e "version_git.txt" ]; then
    echo "# BUILD version_git.txt"
    ./version_git.sh
fi

#
# test cases location
# 
MSL_DIR=current_msl
TEST_DIR=testData
echo TEST_DIR=$TEST_DIR
RESULTS_DIR=testResults
if [ ! -z "$CONTAINER" ]; then RESULTS_DIR=testResultsDocker; fi
echo RESULTS_DIR=$RESULTS_DIR

REPORT=QC.regression_test.summary.txt
echo REPORT=$REPORT
(date; hostname) > $REPORT

# linux cmdline colors
GREEN='\033[1;32m' # use bold version for legibility
#GREEN='\033[0;32m'
RED='\033[0;31m'
MAGENTA='\033[0;35m'
BLUE='\033[0;34m'
BLACK='\033[0;30m'
YELLOW='\033[0;33m'
CYAN='\033[0;36m'
WHITE='\033[0;37m'
RESET="\033[0m"


#
# scan for test directories
#
echo "# find $TEST_DIR -type d -name "$TEST_PAT" -name "proposal*" \! -name "*result*" | sed \"s|${TEST_DIR}.||\" "
TESTS=$(find $TEST_DIR -type d -name "$TEST_PAT" -name "proposal*" \! -name "*result*" | sed "s|${TEST_DIR}.||" )
echo TESTS=$TESTS

#
# iterate
#
for TEST in $TESTS; do
    #
    # input/output for script
    #
    TEST_MSL=.
    TEST_CASE=$TEST
    # support MSL-specific test cases
    if [[ $TEST == msl* ]]; then
	TEST_MSL=$(dirname $TEST)
	REF_MSL=current_msl/$TEST_MSL/taxonomy_node_export.utf8.txt
	TEST_MSL_NUM=$(echo $TEST_MSL | sed 's/msl//;s/v.*//;') # msl39v4 -> 39
	TEST_CASE=$(basename $TEST)
    fi
    
    SRC_DIR=$TEST_DIR/$TEST_MSL/$TEST_CASE
    DEST_DIR=${RESULTS_DIR}/$TEST_MSL/$TEST_CASE
    RESULTS=${DEST_DIR}/QC.regression.new.tsv
    RESULTSBASE=${DEST_DIR}/QC.regression.tsv
    RESULTSDIFF=${DEST_DIR}/QC.regression.diff
    DWDIFF_DELIM="" # delimiters for dwdiff (-P includes "-", which goes too far
    RESULTSDWDIFF=${DEST_DIR}/QC.regression.dwdiff
    RESULTSDWDIFFSHORT=${DEST_DIR}/QC.regression.sdwdiff
    LOGT=${DEST_DIR}/log.new.tmp
    LOG=${DEST_DIR}/log.new.txt
    LOGBASE=${DEST_DIR}/log.txt
    LOGDIFF=${DEST_DIR}/log.diff
    LOGDWDIFF=${DEST_DIR}/log.dwdiff
    MSL=${DEST_DIR}/msl.tsv
    MSLREFLOCAL=${DEST_DIR}/$TEST_MSL.tsv
    MSLDIFF=${DEST_DIR}/msl.vs.${TEST_MSL}.new.txt
    MSLDIFFBASE=${DEST_DIR}/msl.vs.${TEST_MSL}.txt
    MSLDIFFDIFF=${DEST_DIR}/msl.diff.txt
    
    mkdir -p $DEST_DIR
    #
    # header
    #
    echo -e "${BLUE}#########################################${RESET}"
    echo -e "${BLUE}###### $TEST ${RESET}"
    echo -e "${BLUE}#########################################${RESET}"
    echo SRC_DIR=$SRC_DIR
    echo DEST_DIR=$DEST_DIR
    echo RESULTS=$RESULTS
    echo RESULTSBASE=$RESULTSBASE
    echo MSLDIFFBASE=$MSLDIFFBASE
    echo LOG=$LOG

    #
    # run script
    #
    if [ -z "$CONTAINER" ]; then 
	    echo "#" \
	         Rscript merge_proposal_zips.R \
		    --refDir=$MSL_DIR/$TEST_MSL \
		    --proposalsDir=$SRC_DIR \
		    --outDir=$DEST_DIR \
		    --msl \
		    --qcTsvRegression=$(basename $RESULTS) \
		    '2>&1' | tee $LOG
	    Rscript merge_proposal_zips.R \
		    --refDir=$MSL_DIR/$TEST_MSL \
		    --proposalsDir=$SRC_DIR \
		    --outDir=$DEST_DIR \
		    --msl \
		    --qcTsvRegression=$(basename $RESULTS) \
		    1>> $LOG 2>&1
    else
	    echo "#" \
		sudo docker run -it \
		    -v "$(pwd)/${TEST_DIR}:/testData":ro \
		    -v "$(pwd)/${RESULTS_DIR}:/testResults":rw \
	            $CONTAINER  \
		    /merge_proposal_zips.R \
		    --refDir=current_msl/${TEST_MSL} \
		    --proposalsDir="testData/$TEST_MSL/$TEST_CASE" \
		    --outDir="testResults/$TEST_MSL/$TEST_CASE" \
		    --msl \
		    --qcTsvRegression=$(basename $RESULTS) \
		    2>&1 | tee $LOG
	    (sudo docker run -it \
		    -v "$(pwd)/${TEST_DIR}:/testData":ro \
		    -v "$(pwd)/${RESULTS_DIR}:/testResults":rw \
	            $CONTAINER  \
		    /merge_proposal_zips.R \
		    --refDir=current_msl/${TEST_MSL} \
		    --proposalsDir="testData/$TEST_MSL/$TEST_CASE" \
		    --outDir="testResults/$TEST_MSL/$TEST_CASE" \
		    --msl \
		    --qcTsvRegression=$(basename $RESULTS) \
		    ) 1>$LOGT 2>&1 
   	    # git rid of warnigns we only see inside docker
	    grep -v "to see the first 50" $LOGT >>$LOG
    fi	

    #
    # check output
    #
    if [[ ! -e $RESULTSBASE || ! -e $RESULTS ]]; then 
	echo -e "${ORANGE}*MISS  OUT  $TEST${RESET}" | tee -a $REPORT
    else
	# dwdiff (short): diff -w -u | dwdiff -u
        echo "diff -w -u <(cut -f 5- $RESULTSBASE) <(cut -f 5- $RESULTS) | dwdiff -u --delimiters='$DWDIFF_DELIM' --color #> $RESULTSDWDIFFSHORT" | tee $RESULTSDWDIFFSHORT
        diff -w -u <(cut -f 5- $RESULTSBASE) <(cut -f 5- $RESULTS) 2>&1  | dwdiff -u --delimiters="$DWDIFF_DELIM" --color >> $RESULTSDWDIFFSHORT
	# dwdiff 
        echo "dwdiff --delimiters='$DWDIFF_DELIM' --color  <(cut -f 5- $RESULTSBASE) <(cut -f 5- $RESULTS) #> $RESULTSDWDIFF" | tee $RESULTSDWDIFF
        dwdiff --delimiters="$DWDIFF_DELIM" --color <(cut -f 5- $RESULTSBASE) <(cut -f 5- $RESULTS) 2>&1 >> $RESULTSDWDIFF
	# diff -y
        echo "diff -yw --color -W 200 \<(cut -f 5- $RESULTSBASE) \<(cut -f 5- $RESULTS) \> $RESULTSDIFF" | tee $RESULTSDIFF
        diff -yw --color -W 200 <(cut -f 5- $RESULTSBASE) <(cut -f 5- $RESULTS) 2>&1 >> $RESULTSDIFF; RC=$?
        if [ $RC -eq "0" ]; then
            echo -e "${GREEN}ok     OUT  $TEST${RESET}" | tee -a $REPORT
        else
            echo -e "${RED}*FAIL  OUT  $TEST${RESET}" | tee -a $REPORT
        fi	
    fi

    #
    # taxonomy output - not working yet
    #
    # row order of msl.tsv vs msl39v4.tsv might be ONE of the issues
    #
    # get just the needed columns from REF_MSL, then strip MSL number off
    echo "cut -f 4-14,23-31 $REF_MSL | egrep \"^$TEST_MSL_NUM\\t\" > $MSLREFLOCAL" |tee $MSLREFLOCAL
    cut -f 4-14,23-31 $REF_MSL | egrep "^$TEST_MSL_NUM\\t" >> $MSLREFLOCAL
    echo "diff <(cut -f 2-19 $MSLREFLOCAL) <(cut -f 2- $MSL) 2>&1 >> $MSLDIFF" >> $MSLDIFF
    diff <(cut -f 2-19 $MSLREFLOCAL) <(cut -f 2- $MSL) 2>&1 >> $MSLDIFF
#    if [[ ! -e $MSLDIFFBASE || ! -e $MSLDIFF ]]; then 
#	echo "*MISS  MSL  $TEST" | tee -a $REPORT
#    else
#        echo "diff -u \<(cut -f 1- $MSLDIFFBASE) \<(cut -f 1- $MSLDIFF) | dwdiff -u --color \> $MSLDIFFDIFF" | tee $MSLDIFFDIFF
#        diff -u <(cut -f 1- $MSLDIFFBASE) <(cut -f 1- $MSLDIFF) | dwdiff -u --color 2>&1 >> $MSLDIFFDIFF; RC=$?
#        if [ $RC -eq "0" ]; then
#            echo "ok     MSL  $TEST" | tee -a $REPORT
#        else
#            echo "*FAIL  MSL  $TEST" | tee -a $REPORT
#        fi	
#    fi

    #
    # check log
    #
    # use "tail -n +3" to skip date/version/etc in first 2 lines
    #
    if [[ ! -e $LOGBASE || ! -e $LOG ]]; then 
	echo "*MISS  OUT  $TEST" | tee -a $REPORT
    else
	# official diff
        echo "diff -yw --color -W 200 \<(tail -n +3 $LOGBASE|sed -e 's/\[?25h//g'|awk '/^.+$/{print $0}') \<(tail -n +3 $LOG|sed -e 's/\[?25h//g'|awk '/^.+$/{print $0}') \> $LOGDIFF" | tee $LOGDIFF
        diff -yw --color -W 200 <(tail -n +3 $LOGBASE|sed -e 's/\[?25h//g'|awk '/^.+$/{print $0}') <(tail -n +3 $LOG|sed -e 's/\[?25h//g'|awk '/^.+$/{print $0}') 2>&1 >> $LOGDIFF; RC=$?
        if [ $RC -eq "0" ]; then
            echo -e "${GREEN}ok     LOG  $TEST${RESET}" | tee -a $REPORT
        else
            echo -e "${RED}*FAIL  LOG  $TEST${RESET}" | tee -a $REPORT
        fi
	# unofficial, prettier dwdiff
	echo "dwdiff --delimiters='$DWDIFF_DELIM' --color <(tail -n +3 $LOGBASE|sed -e 's/\[?25h//g'|awk '/^.+$/{print $0}') <(tail -n +3 $LOG|sed -e 's/\[?25h//g'|awk '/^.+$/{print $0}') 2>&1 #> $LOGDWDIFF" | tee $LOGDWDIFF
	dwdiff --delimiters="$DWDIFF_DELIM" --color <(tail -n +3 $LOGBASE|sed -e 's/\[?25h//g'|awk '/^.+$/{print $0}') <(tail -n +3 $LOG|sed -e 's/\[?25h//g'|awk '/^.+$/{print $0}') 2>&1 >> $LOGDWDIFF
    fi
    echo "#-------------------------" | tee -a $REPORT
	
done
echo -e "${BLUE}#########################################${RESET}"
echo -e "${BLUE}############### SUMMARY #################${RESET}"
echo -e "${BLUE}#########################################${RESET}"
cat $REPORT
   
