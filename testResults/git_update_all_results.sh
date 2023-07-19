#!/usr/bin/env bash
#
# git add all current results & logs
#
#
for test_case in proposals*; do
    echo "#### $test_case ####"
    eval "pushd $test_case; ../git_update_results.sh; popd" 
done
