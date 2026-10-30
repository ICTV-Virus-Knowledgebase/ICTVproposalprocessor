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

#
# target docker container naem
#
NAME_FILE=docker_container_name.txt
if [ ! -e "$NAME_FILE" ]; then 
	echo "ERROR: $NAME_FILE not found"
	exit 1
fi
echo READING $NAME_FILE
DOCKER_CONTAINER_NAME=$(egrep -v "^#" $NAME_FILE| sed 's/[ \t\r\n]*//g')
echo DOCKER_CONTAINER_NAME=$DOCKER_CONTAINER_NAME
echo ""



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
	${DOCKER_CONTAINER_NAME}  \
	bash


