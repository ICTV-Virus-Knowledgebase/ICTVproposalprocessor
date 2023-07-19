#!/usr/bin/env bash
#
# diff all log files
#
#
FILTER=""
if [ "$1" == "-q" ]; then
    FILTER='-I "VERSION:" -I "no host_source column"'
fi

echo "# diff -w  $FILTER"
for test_case in proposals*; do
    echo "#### $test_case ####"
    eval "diff -w $FILTER  $test_case/log.new.txt $test_case/log.txt"
done
