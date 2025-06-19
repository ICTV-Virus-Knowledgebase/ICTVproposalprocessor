#!/usr/bin/env bash
#
# clone ICTVdatabase repo, and pull out data files and copy here
#
REPO_URL="git@github.com:ICTV-Virus-Knowledgebase/ICTVdatabase.git"
TEMP=$(mktemp -d) 
REPO_DIR="$TEMP/ICTVdatabase"
REPO_DIR="/Users/curtish/Documents/ICTV/ICTVdatabase"

# clone from github
#git clone $REPO_URL $REPO_DIR

# copy here
for file in $(ls $REPO_DIR/data/*.utf8.txt); do
	dest_file=$(basename $file)
	echo "cp -a $file $dest_file"
	cp -a $file $dest_file
	git add $dest_file
done
