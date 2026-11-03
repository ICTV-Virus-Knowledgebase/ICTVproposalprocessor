#!/usr/bin/env bash
#
# run and evaluate regression tests
#
# USAGE: ./regression_test_og.sh [--docker] [test_pattern] [container_version]
#
# On Linux, the script re-runs itself inside ./conda/r_env (conda run), so
# Rscript comes from that env, and writes to testResults.
# Use --docker to run the container and write to testResultsDocker.
# Other platforms use Rscript from PATH unless --docker is selected.
#
# TO DO
# add -msl and diff msl.tsv vs ref MSL and add that to git.
#
ORIG_ARGS=("$@") # for re-running inside the conda env
USE_DOCKER=0
if [[ "${1:-}" == "--docker" ]]; then
    USE_DOCKER=1
    shift
fi

# Which tests to run.
TEST_PAT="*"
if [[ -n "${1:-}" ]]; then TEST_PAT="*$1*"; shift; fi
echo TEST_PAT=$TEST_PAT

# Keep the legacy positional container-version form as an explicit Docker opt-in.
if [[ "$USE_DOCKER" -eq 0 && "$(uname)" == "Linux" && -n "${1:-}" ]]; then
    USE_DOCKER=1
fi

CONTAINER=""
R_COMMAND=()
if [[ "$USE_DOCKER" -eq 1 ]]; then
	CONTAINER=ictv_proposal_processor
	if [[ -n "${1:-}" ]]; then CONTAINER="curtish/${CONTAINER}:$1"; shift; fi
	echo "CONTAINER=$CONTAINER"
	echo "# Building docker image"
	./docker_build_image.sh || exit 2
else
	if [[ "$(uname)" == "Linux" ]]; then
		# Use the Rscript from ./conda/r_env. Wrapping each Rscript call in "conda run"
		# buffers R's output, writes its stderr before its stdout, and adds conda's own
		# ERROR line to the log whenever R exits non-zero. So instead, re-run this
		# whole script once inside the env, and call Rscript directly.
		R_ENV_BIN=./conda/r_env/bin
		if [[ ! -x "$R_ENV_BIN/Rscript" ]]; then
			echo "ERROR: $R_ENV_BIN/Rscript not found; the Linux R runner needs ./conda/r_env." >&2
			exit 2
		fi
		RSCRIPT_PATH=$(command -v Rscript 2>/dev/null)
		if [[ -z "$RSCRIPT_PATH" || "$(cd "$(dirname "$RSCRIPT_PATH")" && pwd -P)" != "$(cd "$R_ENV_BIN" && pwd -P)" ]]; then
			if [[ -n "${REGRESSION_IN_R_ENV:-}" ]]; then
				echo "ERROR: Rscript on PATH is not $R_ENV_BIN/Rscript, even inside conda run." >&2
				exit 2
			fi
			if ! command -v conda >/dev/null 2>&1; then
				echo "ERROR: conda is required for the Linux R runner." >&2
				exit 2
			fi
			REGRESSION_IN_R_ENV=1 exec conda run --no-capture-output -p ./conda/r_env bash "$0" "${ORIG_ARGS[@]}"
		fi
	fi
	R_COMMAND=(Rscript)
	if ! command -v "${R_COMMAND[0]}" >/dev/null 2>&1; then
		echo "ERROR: R runner '${R_COMMAND[0]}' was not found in PATH." >&2
		exit 2
	fi
fi

if [[ -n "${1:-}" ]]; then
	echo "USAGE: ./regression_test_og.sh [--docker] [test_pattern] [container_version]" >&2
	exit 2
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
if [[ -n "$CONTAINER" ]]; then RESULTS_DIR=testResultsDocker; fi
echo RESULTS_DIR=$RESULTS_DIR

REPORT=QC.regression_test.summary.txt
echo REPORT=$REPORT
(date; hostname) > $REPORT
FAILED=0

# linux cmdline colors
GREEN='\033[1;32m' # use bold version for legibility
#GREEN='\033[0;32m'
RED='\033[0;31m'
MAGENTA='\033[0;35m'
BLUE='\033[0;34m'
BLACK='\033[0;30m'
YELLOW='\033[0;33m'
ORANGE=$YELLOW # 8-color ANSI has no orange
CYAN='\033[0;36m'
WHITE='\033[0;37m'
RESET="\033[0m"

# exit status merge_proposal_zips.R should return for a test case (default 0)
expected_run_rc() {
	case "$1" in
		# no usable proposal files: NO_INPUT_FILES, quit(status=1)
		msl38/proposalsTest1Doc|msl38/proposalsTestEmpty) echo 1 ;;
		*) echo 0 ;;
	esac
}

# Text each comparison diffs. Each runs in a subshell with pipefail, so a failure
# at any step (e.g. an unreadable file) returns non-zero instead of producing
# empty or partial text that could compare as equal.
#   QC:  columns 5 onward, trailing carriage returns removed
#   log: line 3 onward (skips the command and VERSION lines), terminal cursor
#        codes and trailing carriage returns removed, then blank lines dropped;
#        carriage returns go first, so a CRLF blank line is dropped too
qc_compare_text() (
	set -o pipefail
	cut -f 5- "$1" | sed -e $'s/\r$//'
)
log_compare_text() (
	set -o pipefail
	tail -n +3 "$1" | sed -e $'s/\033\\[?25h//g' -e $'s/\r$//' | awk '/^.+$/{print $0}'
)

# compared text is written here, then diffed
CMP_DIR=$(mktemp -d "${TMPDIR:-/tmp}/regression_test_og.XXXXXX") || exit 2
trap 'rm -rf -- "$CMP_DIR"' EXIT

#
# scan for test directories
#
echo "# find $TEST_DIR -type d -name "$TEST_PAT" -name "proposal*" \! -name "*result*" | sed \"s|${TEST_DIR}.||\" "
TESTS=$(find $TEST_DIR -type d -name "$TEST_PAT" -name "proposal*" \! -name "*result*" | sed "s|${TEST_DIR}.||" )
echo TESTS=$TESTS
if [[ -z "$TESTS" ]]; then
	echo -e "${RED}*FAIL  no test cases in $TEST_DIR match '$TEST_PAT'${RESET}" | tee -a $REPORT
	exit 2
fi

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

    # clear this case's previous outputs, so a failed run can't leave stale ones behind;
    # if any can't be removed, fail the case rather than run it next to old results
    if ! rm -f -- "$RESULTS" "$LOG" "$LOGT" "$RESULTSDIFF" \
        "$RESULTSDWDIFF" "$RESULTSDWDIFFSHORT" "$LOGDIFF" "$LOGDWDIFF" "$MSLDIFF" \
        "$MSL" "$DEST_DIR/msl_load.sql" "$DEST_DIR/QC.summary.tsv" \
        "$DEST_DIR/QC.docx_summary.tsv" "$DEST_DIR/QC.pretty_summary.all.xlsx"; then
	echo -e "${RED}*FAIL  CLEAN  $TEST (could not remove previous outputs; case not run)${RESET}" | tee -a $REPORT
	FAILED=1
	echo "#-------------------------" | tee -a $REPORT
	continue
    fi

    #
    # run script
    #
    RUN_RC=0
    if [[ -z "$CONTAINER" ]]; then
	    CMD=("${R_COMMAND[@]}" merge_proposal_zips.R
		    "--refDir=$MSL_DIR/$TEST_MSL"
		    "--proposalsDir=$SRC_DIR"
		    "--outDir=$DEST_DIR"
		    --msl
		    "--qcTsvRegression=$(basename "$RESULTS")")
	    {
		    printf '#'
		    printf ' %q' "${CMD[@]}"
		    printf ' 2>&1\n'
	    } | tee "$LOG"
	    "${CMD[@]}" 1>> "$LOG" 2>&1
	    RUN_RC=$?
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
	    RUN_RC=$?
   	    # git rid of warnigns we only see inside docker
	    grep -v "to see the first 50" $LOGT >>$LOG
    fi	
    EXPECTED_RC=$(expected_run_rc "$TEST")
    if [[ "$RUN_RC" -ne "$EXPECTED_RC" ]]; then
	echo -e "${RED}*FAIL  RUN  $TEST (exit $RUN_RC, expected $EXPECTED_RC)${RESET}" | tee -a $REPORT
	FAILED=1
    fi

    #
    # check output
    #
    if [[ ! -e $RESULTSBASE || ! -e $RESULTS ]]; then 
	echo -e "${ORANGE}*MISS  OUT  $TEST${RESET}" | tee -a $REPORT
	FAILED=1
    elif ! qc_compare_text "$RESULTSBASE" > "$CMP_DIR/qc.base" || ! qc_compare_text "$RESULTS" > "$CMP_DIR/qc.new"; then
	echo -e "${RED}*FAIL  OUT  $TEST (could not read QC files to compare)${RESET}" | tee -a $REPORT
	FAILED=1
    else
	QC_CMP="<(qc_compare_text $RESULTSBASE) <(qc_compare_text $RESULTS)"
	# dwdiff (short): diff -u | dwdiff -u
        echo "diff -u $QC_CMP | dwdiff -u --delimiters='$DWDIFF_DELIM' --color #> $RESULTSDWDIFFSHORT" | tee $RESULTSDWDIFFSHORT
        diff -u --label baseline --label current "$CMP_DIR/qc.base" "$CMP_DIR/qc.new" 2>&1  | dwdiff -u --delimiters="$DWDIFF_DELIM" --color >> $RESULTSDWDIFFSHORT
	# dwdiff 
        echo "dwdiff --delimiters='$DWDIFF_DELIM' --color  $QC_CMP #> $RESULTSDWDIFF" | tee $RESULTSDWDIFF
        dwdiff --delimiters="$DWDIFF_DELIM" --color "$CMP_DIR/qc.base" "$CMP_DIR/qc.new" 2>&1 >> $RESULTSDWDIFF
	# diff -y
        echo "diff -y --color -W 200 $QC_CMP \> $RESULTSDIFF" | tee $RESULTSDIFF
        diff -y --color -W 200 "$CMP_DIR/qc.base" "$CMP_DIR/qc.new" 2>&1 >> $RESULTSDIFF; RC=$?
        if [ $RC -eq "0" ]; then
            echo -e "${GREEN}ok     OUT  $TEST${RESET}" | tee -a $REPORT
        else
            echo -e "${RED}*FAIL  OUT  $TEST${RESET}" | tee -a $REPORT
	    FAILED=1
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
	echo -e "${ORANGE}*MISS  LOG  $TEST${RESET}" | tee -a $REPORT
	FAILED=1
    elif ! log_compare_text "$LOGBASE" > "$CMP_DIR/log.base" || ! log_compare_text "$LOG" > "$CMP_DIR/log.new"; then
	echo -e "${RED}*FAIL  LOG  $TEST (could not read log files to compare)${RESET}" | tee -a $REPORT
	FAILED=1
    else
	LOG_CMP="<(log_compare_text $LOGBASE) <(log_compare_text $LOG)"
	# official diff
        echo "diff -y --color -W 200 $LOG_CMP \> $LOGDIFF" | tee $LOGDIFF
        diff -y --color -W 200 "$CMP_DIR/log.base" "$CMP_DIR/log.new" 2>&1 >> $LOGDIFF; RC=$?
        if [ $RC -eq "0" ]; then
            echo -e "${GREEN}ok     LOG  $TEST${RESET}" | tee -a $REPORT
        else
            echo -e "${RED}*FAIL  LOG  $TEST${RESET}" | tee -a $REPORT
	    FAILED=1
        fi
	# unofficial, prettier dwdiff
	echo "dwdiff --delimiters='$DWDIFF_DELIM' --color $LOG_CMP 2>&1 #> $LOGDWDIFF" | tee $LOGDWDIFF
	dwdiff --delimiters="$DWDIFF_DELIM" --color "$CMP_DIR/log.base" "$CMP_DIR/log.new" 2>&1 >> $LOGDWDIFF
    fi
    echo "#-------------------------" | tee -a $REPORT
	
done
echo -e "${BLUE}#########################################${RESET}"
echo -e "${BLUE}############### SUMMARY #################${RESET}"
echo -e "${BLUE}#########################################${RESET}"
cat $REPORT
exit "$FAILED"
