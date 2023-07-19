#!/usr/bin/env bash
#
# diff all regression files
#
#
echo '# diff -wy --suppress-common-lines '
for test_case in proposals*; do
    echo "#### $test_case ####"
    diff -wy --suppress-common-lines  <(cut -f 5- $test_case/QC.regression.new.tsv) <(cut -f 5- $test_case/QC.regression.tsv)
done
