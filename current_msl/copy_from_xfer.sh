#!/usr/bin/env bash
#
# copy and rename files from xfer/ (export) directory
#
# copy_from_current DEST_DIR SRC_DIR
#
# The main export is done as part of the repository
#   https://github.com/ICTV-Virus-Knowledgebase/ICTVonlineDbLoad
#
# This copies the files from there, and renames them
#   TABLE.utf8.osx.txt -> TABLE.utf8.txt
#
SRC=../../ICTVonlineDbLoad/data
SRC2=../../xfer/prod/export_msl
DEST=$(ls -1dt ./msl* | head -1)

if [[ -z "$1" || -z "$2" ]]; then
    echo "SYNTAX: $0 DEST_DIR SRC_DIR"
    echo " "
    echo "Examples:"
    echo " "
    echo "$0 $DEST $SRC"
    echo "$0 $DEST $SRC2"
    exit 1
fi
DEST_DIR="$1"
SRC_DIR="$2"

echo DEST_DIR="$DEST_DIR"
echo SRC_DIR="$SRC_DIR"

if [[ ! -d "$DEST_DIR" ]]; then
       echo "ERROR: not a directory: $DEST_DIR"
       exit 1
fi
if [[ ! -d "$SRC_DIR" ]]; then
       echo "ERROR: not a directory: $SRC_DIR"
       exit 1
fi

for file in $DEST_DIR/*.utf8.txt; do
    # map file names
    SRC_FILE=$SRC_DIR/$(basename $file .txt).osx.txt

    # check if exits/diff
    if [ ! -e "$SRC_FILE" ]; then
	echo "MISS:  $file: source file does not exist: $SRC_FILE"
	continue
    fi
    
    # check if exits/diff
    diff -q "$SRC_FILE" "$file" 2>1 > /dev/null
    RC=$?
    if [ "$RC" -eq 0 ]; then
	echo "SAME:  $file = $SRC_FILE"
    else
	echo "UPDATE: $file <<< $SRC_FILE"
    fi

    # actually copy the file
    cp -a "$SRC_FILE" "$file"
done


	
    
	
