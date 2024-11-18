#!/usr/bin/env bash
#
# rename subdirs to NOT have spaces
#

echo "#"
echo "# dirs with spaces (level 1)"
echo "#"
find . -name "* *" -type d -maxdepth 1

echo "#"
echo '# fix with mv: change /[() ]+/ into underscore'
echo "#"
find . -name "* *" -type d -maxdepth 1 \
     | awk 'BEGIN{FS="\t"}{print $1;on=$1;nn=$1;x=gsub(/[() ]+/,"_",nn); cmd="mv \"" on "\" " nn; print "# " cmd; cmd | getline}'
   
