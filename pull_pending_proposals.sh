#!/usr/bin/env bash
#
# ./pull_pending_proposals.sh [SERVER]
#
# Copy down the "pending" proposals from the 
# ICTV production Drupal server
#
# Use HTTP, since ssh/rsync doesn't work from test to prod

# Defaults
SERVER="ictv.global"
# RSYNC PATH
PENDING_DIR="/var/www/drupal/files/private/proposals/pending/All_Proposals_(Compressed_.zip_files)"
# HTTPS PATH
PENDING_DIR="/system/files/proposals/pending/All_Proposals_%28Compressed_.zip_files%29"
LOCAL_DIR="./pending"

ZIP_FILES="\
	Animal_DNA_viruses_and_retroviruses.zip \
	Animal_dsRNA_and_ssRNA-_viruses.zip \
	Archael_viruses.zip \
	Bacterial_viruses.zip \
	Fungal_and_Protist_viruses.zip \
	General_proposals.zip \
	Plant_viruses.zip \
"

# ARGS
if [[ ! -z "$1" ]]; then 
	SERVER="$1"
	shift
fi
echo "SERVER=     $SERVER"
echo "PENDING_DIR=$PENDING_DIR"
echo "LOCAL_DIR=  $LOCAL_DIR"
echo "ZIP_FILES=  $ZIP_FILES"

# COPY down ZIPs
for zip in $ZIP_FILES; do

	echo wget --mirror --no-parent --no-directories --directory-prefix "$LOCAL_DIR" "https://${SERVER}/${PENDING_DIR}/${zip}"
	wget --mirror --no-parent --no-directories --directory-prefix "$LOCAL_DIR" "https://${SERVER}/${PENDING_DIR}/${zip}"

	RC=$?
	if [ "$RC" -ne 0 ]; then 
		echo "ERROR: wget failed; RC=$RC"
		exit 1
	fi
done

# QC
ls -lstr $LOCAL_DIR/*.zip

# UNCOMPRESS
for zip in $ZIP_FILES; do
	echo "(cd $LOCAL_DIR; unzip $zip)"
	(cd $LOCAL_DIR; unzip $zip)

	RC=$?
	if [ "$RC" -ne 0 ]; then 
		echo "ERROR: unzip failed; RC=$RC"
		exit 1
	fi
done

# QC
ls -lstr $LOCAL_DIR/
