#!/usr/bin/env bash
#
# ../../fix_links.sh       # updates existing links
# OR
# ../../fix_links.sh -n    # create 2 expected links
# OR
# ./fix_links.sh -r        # recursive
#

# 
# RECURSE INTO ALL PROPOSAL* SUBDIRS
#
if [[ "$1" == -r* ]]; then 
	# clear that arg
	shift
	# recursive mode
	echo "RECURSIVE"
	DIRS=$(find . -maxdepth 2 -type d -name "proposals*")
	echo DIRS=$DIRS
	for DIR in $DIRS; do
		pushd $DIR
		../../$0 $*
		popd 
	done
	exit 0
fi

#
# MAIN WORK
#
#
# Update links for new msl## directory in path
FILES=$(find . -type l)
if [[ "$1" == -n* ]]; then
	FILES="QC.regression.tsv log.txt"
	echo "NEW: create $FILES"
fi

if [ -z "$FILES" ]; then 
	echo "ERROR: no symlinks to updated"
	exit 1
else
	echo FILES=$FILES
fi

# calc src dir
SRC_DIR=../../../testResults/$(basename $(dirname $PWD))/$(basename $PWD)
echo SRC_DIR=$SRC_DIR

for FILE in $FILES; do
	echo ln -fs $SRC_DIR/$FILE $FILE
	ln -fs $SRC_DIR/$FILE $FILE
	git add $FILE
done

echo "DONE"
