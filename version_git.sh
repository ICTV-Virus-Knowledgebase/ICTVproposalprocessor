#!/usr/bin/env bash
#
# get git hash and append to version.txt
#
echo "$(cat version.txt).$(git log -1 --date=format:'%Y%m%d'| awk '/^Date:/{print $2}').$(git log -1| egrep "^commit " | cut -c 8-14)" > version_git.txt
echo "VERSION: $(cat version.txt)  ($(cat version_git.txt))"
