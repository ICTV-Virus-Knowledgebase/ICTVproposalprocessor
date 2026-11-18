#!/usr/bin/env bash
#
# compare the R baselines with the Python baselines
#
# USAGE: ./regression_compare_r_python.sh [test_pattern]
#
#   R:      testResults/<msl>/<case>/QC.regression.tsv, log.txt
#   Python: testResultsPython/<msl>/<case>/QC.regression.tsv, log.txt
#
# Runs neither processor; it only reads the committed baselines. Comparison
# files go in testResultsPython/<msl>/<case>/ (nothing is written to testResults):
#   QC.R_vs_Python.diff     diff -y: R on the LEFT, Python on the RIGHT
#   QC.R_vs_Python.sdwdiff  diff -u | dwdiff -u: "-" is R, "+" is Python
#   QC.R_vs_Python.dwdiff   dwdiff: [-R-] {+Python+}
#   log.R_vs_Python.diff, log.R_vs_Python.dwdiff  the same for the logs
# The report goes to QC.regression_compare_r_python.summary.txt.
#
# Text compared, as in regression_test_og.sh / regression_test_python.sh:
#   QC:  columns 5 onward, trailing carriage returns removed
#   log: line 2 onward (skips the command) without the VERSION: line, terminal
#        cursor codes, trailing carriage returns and blank lines removed.
#        The results directory each implementation wrote to (testResults/ or
#        testResultsPython/) is replaced with RESULTS/ on both sides.
#
# Which tests to compare.
TEST_PAT="*"
if [[ -n "${1:-}" ]]; then TEST_PAT="*$1*"; shift; fi
echo TEST_PAT=$TEST_PAT

if [[ -n "${1:-}" ]]; then
	echo "USAGE: ./regression_compare_r_python.sh [test_pattern]" >&2
	exit 2
fi

#
# test cases location
#
TEST_DIR=testData
R_DIR=testResults
PY_DIR=testResultsPython
echo TEST_DIR=$TEST_DIR
echo "R_DIR=$R_DIR (R baselines)"
echo "PY_DIR=$PY_DIR (Python baselines)"

REPORT=QC.regression_compare_r_python.summary.txt
echo REPORT=$REPORT
(date; hostname; echo "R=$R_DIR  Python=$PY_DIR  (in diffs: R is LEFT / -, Python is RIGHT / +)") > $REPORT
FAILED=0
N_CASES=0; N_QC_SAME=0; N_QC_DIFF=0; N_LOG_SAME=0; N_LOG_DIFF=0; N_MISS=0

# linux cmdline colors
GREEN='\033[1;32m' # use bold version for legibility
RED='\033[0;31m'
BLUE='\033[0;34m'
YELLOW='\033[0;33m'
ORANGE=$YELLOW # 8-color ANSI has no orange
RESET="\033[0m"

# Text each comparison diffs. Each runs in a subshell with pipefail, so a failure
# at any step (e.g. an unreadable file) returns non-zero instead of producing
# empty or partial text that could compare as equal.
#   $1 = file, $2 = results directory that implementation wrote to
qc_compare_text() (
	set -o pipefail
	cut -f 5- "$1" | sed -e $'s/\r$//'
)
log_compare_text() (
	set -o pipefail
	tail -n +2 "$1" | sed -e $'s/\033\\[?25h//g' -e $'s/\r$//' -e "s#$2/#RESULTS/#g" \
		| awk '/^.+$/ && !/^VERSION: /{print $0}'
)
# data rows in a QC file (lines after the header)
qc_rows() {
	awk 'END{print (NR > 0 ? NR - 1 : 0)}' "$1"
}

# compared text is written here, then diffed
CMP_DIR=$(mktemp -d "${TMPDIR:-/tmp}/regression_compare_r_python.XXXXXX") || exit 2
trap 'rm -rf -- "$CMP_DIR"' EXIT

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
    N_CASES=$((N_CASES + 1))
    R_CASE=$R_DIR/$TEST
    PY_CASE=$PY_DIR/$TEST
    R_QC=$R_CASE/QC.regression.tsv
    PY_QC=$PY_CASE/QC.regression.tsv
    R_LOG=$R_CASE/log.txt
    PY_LOG=$PY_CASE/log.txt
    QCDIFF=$PY_CASE/QC.R_vs_Python.diff
    QCDWDIFF=$PY_CASE/QC.R_vs_Python.dwdiff
    QCDWDIFFSHORT=$PY_CASE/QC.R_vs_Python.sdwdiff
    LOGDIFF=$PY_CASE/log.R_vs_Python.diff
    LOGDWDIFF=$PY_CASE/log.R_vs_Python.dwdiff
    DWDIFF_DELIM="" # delimiters for dwdiff (-P includes "-", which goes too far

    #
    # header
    #
    echo -e "${BLUE}#########################################${RESET}"
    echo -e "${BLUE}###### $TEST ${RESET}"
    echo -e "${BLUE}#########################################${RESET}"
    echo "R      QC=$R_QC LOG=$R_LOG"
    echo "Python QC=$PY_QC LOG=$PY_LOG"

    mkdir -p "$PY_CASE"
    # clear this case's previous comparison, so a failed run can't leave a stale one behind
    if ! rm -f -- "$QCDIFF" "$QCDWDIFF" "$QCDWDIFFSHORT" "$LOGDIFF" "$LOGDWDIFF"; then
	echo -e "${RED}*FAIL  CLEAN  $TEST (could not remove previous comparison files)${RESET}" | tee -a $REPORT
	FAILED=1
	echo "#-------------------------" | tee -a $REPORT
	continue
    fi

    #
    # QC
    #
    if [[ ! -e $R_QC || ! -e $PY_QC ]]; then
	[[ -e $R_QC ]] || echo -e "${ORANGE}*MISS  QC   $TEST (no R baseline: $R_QC)${RESET}" | tee -a $REPORT
	[[ -e $PY_QC ]] || echo -e "${ORANGE}*MISS  QC   $TEST (no Python baseline: $PY_QC)${RESET}" | tee -a $REPORT
	N_MISS=$((N_MISS + 1))
	FAILED=1
    elif ! qc_compare_text "$R_QC" > "$CMP_DIR/qc.R" || ! qc_compare_text "$PY_QC" > "$CMP_DIR/qc.Python"; then
	echo -e "${RED}*FAIL  QC   $TEST (could not read QC files to compare)${RESET}" | tee -a $REPORT
	FAILED=1
    else
	QC_CMP="<(qc_compare_text $R_QC) <(qc_compare_text $PY_QC)"
	ROWS="R $(qc_rows "$R_QC") rows, Python $(qc_rows "$PY_QC") rows"
	if command -v dwdiff >/dev/null 2>&1; then
	    # dwdiff (short): diff -u | dwdiff -u
	    echo "# R=- Python=+ : diff -u $QC_CMP | dwdiff -u --delimiters='$DWDIFF_DELIM' --color" > $QCDWDIFFSHORT
	    diff -u --label "R $R_QC" --label "Python $PY_QC" "$CMP_DIR/qc.R" "$CMP_DIR/qc.Python" 2>&1 | dwdiff -u --delimiters="$DWDIFF_DELIM" --color >> $QCDWDIFFSHORT
	    # dwdiff
	    echo "# R=[-...-] Python={+...+} : dwdiff --delimiters='$DWDIFF_DELIM' --color $QC_CMP" > $QCDWDIFF
	    dwdiff --delimiters="$DWDIFF_DELIM" --color "$CMP_DIR/qc.R" "$CMP_DIR/qc.Python" 2>&1 >> $QCDWDIFF
	fi
	# diff -y
	echo "# R=LEFT Python=RIGHT : diff -y --color -W 200 $QC_CMP" > $QCDIFF
	diff -y --color -W 200 "$CMP_DIR/qc.R" "$CMP_DIR/qc.Python" 2>&1 >> $QCDIFF; RC=$?
	if [ $RC -eq "0" ]; then
	    echo -e "${GREEN}ok     QC   $TEST (R == Python; $ROWS)${RESET}" | tee -a $REPORT
	    N_QC_SAME=$((N_QC_SAME + 1))
	else
	    echo -e "${RED}*DIFF  QC   $TEST (R != Python; $ROWS)${RESET}" | tee -a $REPORT
	    N_QC_DIFF=$((N_QC_DIFF + 1))
	    FAILED=1
	fi
    fi

    #
    # log
    #
    if [[ ! -e $R_LOG || ! -e $PY_LOG ]]; then
	[[ -e $R_LOG ]] || echo -e "${ORANGE}*MISS  LOG  $TEST (no R baseline: $R_LOG)${RESET}" | tee -a $REPORT
	[[ -e $PY_LOG ]] || echo -e "${ORANGE}*MISS  LOG  $TEST (no Python baseline: $PY_LOG)${RESET}" | tee -a $REPORT
	N_MISS=$((N_MISS + 1))
	FAILED=1
    elif ! log_compare_text "$R_LOG" "$R_DIR" > "$CMP_DIR/log.R" || ! log_compare_text "$PY_LOG" "$PY_DIR" > "$CMP_DIR/log.Python"; then
	echo -e "${RED}*FAIL  LOG  $TEST (could not read log files to compare)${RESET}" | tee -a $REPORT
	FAILED=1
    else
	LOG_CMP="<(log_compare_text $R_LOG $R_DIR) <(log_compare_text $PY_LOG $PY_DIR)"
	echo "# R=LEFT Python=RIGHT : diff -y --color -W 200 $LOG_CMP" > $LOGDIFF
	diff -y --color -W 200 "$CMP_DIR/log.R" "$CMP_DIR/log.Python" 2>&1 >> $LOGDIFF; RC=$?
	if [ $RC -eq "0" ]; then
	    echo -e "${GREEN}ok     LOG  $TEST (R == Python)${RESET}" | tee -a $REPORT
	    N_LOG_SAME=$((N_LOG_SAME + 1))
	else
	    echo -e "${RED}*DIFF  LOG  $TEST (R != Python)${RESET}" | tee -a $REPORT
	    N_LOG_DIFF=$((N_LOG_DIFF + 1))
	    FAILED=1
	fi
	# unofficial, prettier dwdiff
	if command -v dwdiff >/dev/null 2>&1; then
	    echo "# R=[-...-] Python={+...+} : dwdiff --delimiters='$DWDIFF_DELIM' --color $LOG_CMP" > $LOGDWDIFF
	    dwdiff --delimiters="$DWDIFF_DELIM" --color "$CMP_DIR/log.R" "$CMP_DIR/log.Python" 2>&1 >> $LOGDWDIFF
	fi
    fi
    echo "#-------------------------" | tee -a $REPORT

done
echo "CASES=$N_CASES  QC: same=$N_QC_SAME differ=$N_QC_DIFF  LOG: same=$N_LOG_SAME differ=$N_LOG_DIFF  missing=$N_MISS" >> $REPORT
echo -e "${BLUE}#########################################${RESET}"
echo -e "${BLUE}############### SUMMARY #################${RESET}"
echo -e "${BLUE}#########################################${RESET}"
cat $REPORT
exit "$FAILED"
