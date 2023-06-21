#!/usr/bin/env bash
#
# add git hash to version
./version_git.sh
VER_FILE=version_git.txt
if [ ! -e "$VER_FILE" ]; then 
	echo ERROR: $VER_FILE not found
	exit 1
fi
VER=$(cat $VER_FILE)

#
# Build image from Dockerfile
#

sudo docker build -t ictv_proposal_processor .


