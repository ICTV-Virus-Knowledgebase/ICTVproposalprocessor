#!/usr/bin/env bash
#
# run and evaluate regression tests for merge_proposal_zips.py
#
# USAGE: ./regression_test_python.sh [test_pattern]
#
# Python counterpart of regression_test_og.sh: runs only merge_proposal_zips.py,
# writes to testResultsPython, and compares each case with the Python baselines
# in the same directory (QC.regression.tsv, log.txt). It never reads or writes
# testResults (the R baselines); comparing Python with R is a separate step.
#
# Interpreter: $PYTHON if set, else the python of the conda env
# ictv-proposal-validator-python (environment.python.yml), else python3 on PATH.
#
# Differences from regression_test_og.sh:
#   - each case runs on a temporary copy of its testData directory, because
#     merge_proposal_zips.py renames input files whose names contain spaces.
#     The copy's location is then removed from QC.regression.new.tsv and
#     log.new.txt, so they show testData/... paths, as an R run does.
#   - log comparison skips line 1 (the command) and the VERSION: line.
#   - no Docker mode, and no msl.tsv-vs-reference extract (unfinished in og).
#

# Which tests to run.
TEST_PAT="*"
if [[ -n "${1:-}" ]]; then TEST_PAT="*$1*"; shift; fi
echo TEST_PAT=$TEST_PAT

if [[ -n "${1:-}" ]]; then
	echo "USAGE: ./regression_test_python.sh [test_pattern]" >&2
	exit 2
fi

#
# find the python interpreter
#
PY_ENV_NAME=ictv-proposal-validator-python
if [[ -z "${PYTHON:-}" ]]; then
	CONDA_BIN="${CONDA_EXE:-$(command -v conda 2>/dev/null)}"
	for CANDIDATE in /opt/anaconda3/bin/conda "$HOME/miniconda3/condabin/conda"; do
		if [[ -z "$CONDA_BIN" && -x "$CANDIDATE" ]]; then CONDA_BIN="$CANDIDATE"; fi
	done
	if [[ -n "$CONDA_BIN" ]]; then
		PYTHON=$("$CONDA_BIN" run -n "$PY_ENV_NAME" python -c 'import sys; print(sys.executable)' 2>/dev/null | awk 'NF{last=$0} END{print last}')
	fi
	if [[ -z "$PYTHON" ]]; then
		echo "WARNING: conda env $PY_ENV_NAME not found; using python3 from PATH" >&2
		PYTHON=python3
	fi
fi
if ! command -v "$PYTHON" >/dev/null 2>&1; then
	echo "ERROR: Python interpreter '$PYTHON' was not found." >&2
	exit 2
fi
if ! "$PYTHON" -c 'import pandas, openpyxl' >/dev/null 2>&1; then
	echo "ERROR: '$PYTHON' cannot import pandas and openpyxl; set PYTHON or create the $PY_ENV_NAME env (environment.python.yml)." >&2
	exit 2
fi
PY_COMMAND=("$PYTHON" merge_proposal_zips.py)
# write stdout as it happens, so the log keeps stdout and stderr in order
export PYTHONUNBUFFERED=1
echo PYTHON=$PYTHON

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
RESULTS_DIR=testResultsPython
echo RESULTS_DIR=$RESULTS_DIR

REPORT=QC.regression_test_python.summary.txt
echo REPORT=$REPORT
(date; hostname; echo "PYTHON=$PYTHON ($("$PYTHON" --version 2>&1))") > $REPORT
if [[ -n "${MERGE_PROPOSAL_ZIPS_R_COMPAT:-}" ]]; then
	echo "MERGE_PROPOSAL_ZIPS_R_COMPAT=$MERGE_PROPOSAL_ZIPS_R_COMPAT (Python-only behavior switched off)" | tee -a $REPORT
fi
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

# exit status merge_proposal_zips.py should return for a test case (default 0)
expected_run_rc() {
	case "$1" in
		# no usable proposal files: NO_INPUT_FILES, exit status 1
		msl38/proposalsTest1Doc|msl38/proposalsTestEmpty) echo 1 ;;
		*) echo 0 ;;
	esac
}

# Text each comparison diffs. Each runs in a subshell with pipefail, so a failure
# at any step (e.g. an unreadable file) returns non-zero instead of producing
# empty or partial text that could compare as equal.
#   QC:  columns 5 onward, trailing carriage returns removed
#   log: line 2 onward (skips the command) without the VERSION: line, terminal
#        cursor codes and trailing carriage returns removed, then blank lines
#        dropped; carriage returns go first, so a CRLF blank line is dropped too
qc_compare_text() (
	set -o pipefail
	cut -f 5- "$1" | sed -e $'s/\r$//'
)
log_compare_text() (
	set -o pipefail
	tail -n +2 "$1" | sed -e $'s/\033\\[?25h//g' -e $'s/\r$//' | awk '/^.+$/ && !/^VERSION: /{print $0}'
)

# Remove the input copy's location ("$1", ends in "/") from a file, in place,
# so it records testData/... as if the case had run on testData directly.
strip_prefix() {
	local prefix=$1 file=$2 escaped
	[[ -e "$file" ]] || return 0
	escaped=$(printf '%s' "$prefix" | sed -e 's/[][\.*^$|/]/\\&/g')
	sed -e "s/$escaped//g" "$file" > "$file.strip.tmp" && mv -- "$file.strip.tmp" "$file"
}

# compared text is written to CMP_DIR, then diffed; input copies go in WORK_DIR.
# pwd -P: the path as the processor would resolve it (/var -> /private/var on macOS)
CMP_DIR=$(mktemp -d "${TMPDIR:-/tmp}/regression_test_python.XXXXXX") || exit 2
CMP_DIR=$(cd "$CMP_DIR" && pwd -P) || exit 2
trap 'rm -rf -- "$CMP_DIR"' EXIT
WORK_DIR=$CMP_DIR/work
mkdir -p "$WORK_DIR" || exit 2

#
# scan for test directories
#
echo "# find $TEST_DIR -type d -name "$TEST_PAT" -name "proposal*" \! -name "*result*" | sed \"s|${TEST_DIR}.||\" | sort"
TESTS=$(find $TEST_DIR -type d -name "$TEST_PAT" -name "proposal*" \! -name "*result*" | sed "s|${TEST_DIR}.||" | sort)
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
	TEST_CASE=$(basename $TEST)
    fi

    SRC_DIR=$TEST_DIR/$TEST_MSL/$TEST_CASE
    COPY_DIR=$WORK_DIR/$SRC_DIR
    DEST_DIR=${RESULTS_DIR}/$TEST_MSL/$TEST_CASE
    RESULTS=${DEST_DIR}/QC.regression.new.tsv
    RESULTSBASE=${DEST_DIR}/QC.regression.tsv
    RESULTSDIFF=${DEST_DIR}/QC.regression.diff
    DWDIFF_DELIM="" # delimiters for dwdiff (-P includes "-", which goes too far
    RESULTSDWDIFF=${DEST_DIR}/QC.regression.dwdiff
    RESULTSDWDIFFSHORT=${DEST_DIR}/QC.regression.sdwdiff
    LOG=${DEST_DIR}/log.new.txt
    LOGBASE=${DEST_DIR}/log.txt
    LOGDIFF=${DEST_DIR}/log.diff
    LOGDWDIFF=${DEST_DIR}/log.dwdiff
    MSL=${DEST_DIR}/msl.tsv

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
    echo LOG=$LOG

    # clear this case's previous outputs, so a failed run can't leave stale ones behind;
    # if any can't be removed, fail the case rather than run it next to old results
    if ! rm -f -- "$RESULTS" "$LOG" "$RESULTSDIFF" \
        "$RESULTSDWDIFF" "$RESULTSDWDIFFSHORT" "$LOGDIFF" "$LOGDWDIFF" \
        "$MSL" "$DEST_DIR/msl_load.sql" "$DEST_DIR/QC.summary.tsv" \
        "$DEST_DIR/QC.docx_summary.tsv" "$DEST_DIR/QC.pretty_summary.all.xlsx" \
        "$DEST_DIR/QC.pretty_summary.issues.xlsx"; then
	echo -e "${RED}*FAIL  CLEAN  $TEST (could not remove previous outputs; case not run)${RESET}" | tee -a $REPORT
	FAILED=1
	echo "#-------------------------" | tee -a $REPORT
	continue
    fi

    # independent copy of the inputs: the processor may rename files in it
    if ! { mkdir -p "$(dirname "$COPY_DIR")" && cp -RpL "$SRC_DIR" "$COPY_DIR"; }; then
	echo -e "${RED}*FAIL  COPY  $TEST (could not copy $SRC_DIR; case not run)${RESET}" | tee -a $REPORT
	FAILED=1
	echo "#-------------------------" | tee -a $REPORT
	continue
    fi

    #
    # run script
    #
    CMD=("${PY_COMMAND[@]}"
	    "--refDir=$MSL_DIR/$TEST_MSL"
	    "--proposalsDir=$COPY_DIR"
	    "--outDir=$DEST_DIR"
	    --msl
	    "--qcTsvRegression=$(basename "$RESULTS")")
    {
	    printf '#'
	    printf ' %q' "${CMD[@]}"
	    printf ' 2>&1\n'
    } > "$LOG"
    "${CMD[@]}" 1>> "$LOG" 2>&1
    RUN_RC=$?
    strip_prefix "$WORK_DIR/" "$LOG"
    strip_prefix "$WORK_DIR/" "$RESULTS"
    head -n 1 "$LOG"
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
	if command -v dwdiff >/dev/null 2>&1; then
	    # dwdiff (short): diff -u | dwdiff -u
	    echo "diff -u $QC_CMP | dwdiff -u --delimiters='$DWDIFF_DELIM' --color #> $RESULTSDWDIFFSHORT" | tee $RESULTSDWDIFFSHORT
	    diff -u --label baseline --label current "$CMP_DIR/qc.base" "$CMP_DIR/qc.new" 2>&1  | dwdiff -u --delimiters="$DWDIFF_DELIM" --color >> $RESULTSDWDIFFSHORT
	    # dwdiff
	    echo "dwdiff --delimiters='$DWDIFF_DELIM' --color  $QC_CMP #> $RESULTSDWDIFF" | tee $RESULTSDWDIFF
	    dwdiff --delimiters="$DWDIFF_DELIM" --color "$CMP_DIR/qc.base" "$CMP_DIR/qc.new" 2>&1 >> $RESULTSDWDIFF
	fi
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
    # check log
    #
    # log_compare_text skips the command line and the VERSION: line
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
	if command -v dwdiff >/dev/null 2>&1; then
	    echo "dwdiff --delimiters='$DWDIFF_DELIM' --color $LOG_CMP 2>&1 #> $LOGDWDIFF" | tee $LOGDWDIFF
	    dwdiff --delimiters="$DWDIFF_DELIM" --color "$CMP_DIR/log.base" "$CMP_DIR/log.new" 2>&1 >> $LOGDWDIFF
	fi
    fi
    echo "#-------------------------" | tee -a $REPORT

done
echo -e "${BLUE}#########################################${RESET}"
echo -e "${BLUE}############### SUMMARY #################${RESET}"
echo -e "${BLUE}#########################################${RESET}"
cat $REPORT
exit "$FAILED"
