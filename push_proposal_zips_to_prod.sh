#!/usr/bin/env bash
#
# push_proposal_zips_to_prod.sh test|prod [pass-through-flags]
#
# copy the final per-proposal zip files
# to the apache /proposals/ directory on the
# production server
#
SRC_DIR="./pending/proposalsFinalZips/"
PROD_DEST="ictv.global:/var/www/drupal/files/ICTV/proposals/"
TEST_DEST="ubuntu@test.ictv.global:/var/www/drupal/files/ICTV/proposals/"

if [ "$1" == "test" ]; then DEST="$TEST_DEST"; shift 1
elif [ "$1" == "prod" ]; then DEST="$PROD_DEST"; shift 1
else
    echo "ERROR: unknown or missing args '$*'"
    echo " "
    echo "USAGE: $o test|prod [pass-through-flags]"
    exit 1
fi

echo "#  SRC=$SRC_DIR"
echo "# DEST=$DEST"
echo exec rsync $* -hav --progress --exclude='*~' --no-recursive $SRC_DIR/'20??.*.zip'  $SRC_DIR/'20??.*.doc*' "$DEST"
exec rsync $* -hav --progress --exclude='*~' --no-recursive $SRC_DIR/20??.*.zip $SRC_DIR/20??.*.doc* "$DEST"
