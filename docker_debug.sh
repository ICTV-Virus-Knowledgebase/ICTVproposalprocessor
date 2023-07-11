#!/usr/bin/env bash
#
# test run our docker image
#
TEST=proposalsTest_createNew
if [ ! -z "$1" ]; then TEST="$1"; shift; fi
echo TEST=$TEST
SRC_DIR=testData/$TEST
TEST_DIR=testResultsDocker/$TEST
echo "# cleaning out $TEST_DIR/..."
echo "find testResultsDocker/$TEST -name '*new*' -o -name '*diff*' -exec rm {} +"
find testResultsDocker/$TEST -name '*new*' -o -name '*diff*' -exec rm {} +




cat <<EOT
# edit code (no VI or Nano in image)
awk '(NR==385){print "browser()"}{print \$0}' merge_proposal_zips.R > test2.R

# run
R
> source("test2.R")
EOT
sudo docker run -it \
	-v "$(pwd)/$SRC_DIR:/testData":ro \
	-v "$(pwd)/$TEST_DIR:/testResults" \
	ictv_proposal_processor  \
	bash


