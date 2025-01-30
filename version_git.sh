#!/usr/bin/env bash
#
# get git hash and append to version.txt
#
echo "$(cat version.txt).$(date +"%Y%m%d").$(git log -1| head -1 | cut -c 8-14)" > version_git.txt
echo "VERSION: $(cat version.txt)  ($(cat version_git.txt))"
