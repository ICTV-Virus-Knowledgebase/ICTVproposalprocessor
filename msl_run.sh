#!/usr/bin/env bash
#!/usr/bin/env bash
if [[ "$1" == -h* ]]; then 
	cat <<EOT
#----------------------------------------------------------------------
#
# Run over all proposals for an MSL, as preparing for a release
#
# USAGE: ./msl_run.sh [-help] [-verbose] [-tmi] [-docker] PROPOSAL_DIR [test_pattern] [container_version]
#
#----------------------------------------------------------------------
EOT
	exit 0
fi

#
#----------------------------------------------------------------------
# Argument processing
#----------------------------------------------------------------------
MSL_DIR=$(ls -1dtra current_msl/msl* | tail -1)
USE_DOCKER=0
PROPOSAL_DIR=
RESULTS_DIR=
TEST_PAT=
PASS_THROUGH_FLAGS=

#
# parse flags
#
while [[ "$1" == -* ]]; do

	# R or Docker
	if [[ "$1" == -d* ]]; then 
		USE_DOCKER=1
		shift 1
	fi
	echo "USE_DOCKER=$USE_DOCKER"

	# R pass through flags
	if [[ "$1" == -v* || "$1" == -t* ]]; then
		PASS_THROUGH_FLAGS="$PASS_THROUGH_FLAGS $1"
		shift
	fi
		
done

# input dir
PROPOSAL_DIR="$1"
shift
if [[ -z "$PROPOSAL_DIR" || ! -d "$PROPOSAL_DIR" ]]; then 
	echo "ERROR: MSL_DIR (input dir) is a required argument"
	echo "USAGE: ./msl_run.sh [-docker] MSL_DIR [test_pattern] [docker_container_version]"
	exit 1
fi
RESULTS_DIR=${PROPOSAL_DIR}_Results

# which tests to run 
if [ ! -z "$1" ]; then TEST_PAT="*$1*"; shift; fi
echo TEST_PAT=$TEST_PAT

# which docker container to run
if [ $USE_DOCKER -eq 1 ]; then 
	CONTAINER=ictv_proposal_processor
	if [ ! -z "$1" ]; then 
		CONTAINER="curtish/${CONTAINER}:$1"
		shift
	fi
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

#----------------------------------------------------------------------
# test cases location
#----------------------------------------------------------------------
# 
echo MSL_DIR=$MSL_DIR
echo PROPOSAL_DIR=$PROPOSAL_DIR
echo RESULTS_DIR=$RESULTS_DIR

REPORT=$0.qc_summary.txt
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
if [ -z "$TEST_PAT" ]; then
	TESTS=$PROPOSAL_DIR
else
	SPACES_IN_FILES=$(find $PROPOSAL_DIR \! -name '*result*' \! -name "*.zip" -name "* *")
	if [ ! -z "$SPACES_IN_FILES" ]; then
		echo "ERROR: filenames with spaces are not supported"
		echo "PLEASE CORRECT THE FOLLOWING:"
		find $PROPOSAL_DIR \! -name '*result*' \! -name "*.zip" -name "* *" 
		exit 1
	fi
	echo "# find $PROPOSAL_DIR \! -name '*result*' -name '$TEST_PAT'  "
	TESTS=$(find $PROPOSAL_DIR \! -name '*result*' -name "$TEST_PAT"  )
fi 
echo TESTS=$TESTS

cat <<EOT
#
# iterate
#
EOT
for TEST in $TESTS; do
    #
    # input/output for script
    #
    TEST_CASE=$(basename $TEST .xlsx)
    
    SRC_DIR=$TEST
    DEST_DIR=${RESULTS_DIR}/${TEST_CASE}

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
    if [ $USE_DOCKER -eq 0 ]; then 
	    echo "# RUN R Script"
	    echo "#" \
	         Rscript merge_proposal_zips.R \
		    --refDir=$MSL_DIR/$TEST_MSL \
		    --proposalsDir=$SRC_DIR \
		    --outDir=$DEST_DIR \
		    --msl \
		    --qcTsvRegression=$(basename $RESULTS) \
		    ${PASS_THROUGH_FLAGS} \
		    '2>&1' | tee $LOG
	    Rscript merge_proposal_zips.R \
		    --refDir=$MSL_DIR/$TEST_MSL \
		    --proposalsDir=$SRC_DIR \
		    --outDir=$DEST_DIR \
		    --msl \
		    --qcTsvRegression=$(basename $RESULTS) \
		    ${PASS_THROUGH_FLAGS} \
		    1>> $LOG 2>&1
    else
	    echo "# RUN via DOCKER image"
	    echo "#" \
		sudo docker run -it \
		    -v "$(pwd)/${PROPOSAL_DIR}:/testData":ro \
		    -v "$(pwd)/${RESULTS_DIR}:/testResults":rw \
	            $CONTAINER  \
		    /merge_proposal_zips.R \
		    --refDir=current_msl/${TEST_MSL} \
		    --proposalsDir="testData/$TEST_MSL/$TEST_CASE" \
		    --outDir="testResults/$TEST_MSL/$TEST_CASE" \
		    --msl \
		    --qcTsvRegression=$(basename $RESULTS) \
		    ${PASS_THROUGH_FLAGS} \
		    2>&1 | tee $LOG
	    (sudo docker run -it \
		    -v "$(pwd)/${PROPOSAL_DIR}:/testData":ro \
		    -v "$(pwd)/${RESULTS_DIR}:/testResults":rw \
	            $CONTAINER  \
		    /merge_proposal_zips.R \
		    --refDir=current_msl/${TEST_MSL} \
		    --proposalsDir="testData/$TEST_MSL/$TEST_CASE" \
		    --outDir="testResults/$TEST_MSL/$TEST_CASE" \
		    --msl \
		    --qcTsvRegression=$(basename $RESULTS) \
		    ${PASS_THROUGH_FLAGS} \
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
   
