#!/bin/bash
#
# remove ".v#" from proposal names, and zip them up
#

SRC_DIR=proposals3
DEST_DIR=proposalsFinal

#
# scan for proposal codes
#
CODES=$(find $SRC_DIR -name "*docx" -o  -name "*xlsx" | grep -v Suppl | egrep -v "([#~])" | sed -e 's/.*\///;' | cut -d . -f 1-2 | sort | uniq )
echo "Found "$(echo $CODES | wc -w)" proposals in $SRC_DIR/"

#
# for each one
#
for CODE in $CODES; do

	# copy files to DEST, keeping top-level dir, removing ".v#[.fix]" 
	SRC_FILELIST=$(find $SRC_DIR -name "${CODE}.*docx" -o  -name "${CODE}.*xlsx" | egrep -v "([#~])")	
	for SRC_FILE in $SRC_FILELIST; do
		DEST_FILE=
	done
	# zip up files in DEST

done


