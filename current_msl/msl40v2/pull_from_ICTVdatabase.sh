#!/usr/bin/env bash
#
# clone ICTVdatabase repo, and pull out data files and copy here
#
set -euo pipefail

REPO_URL="git@github.com:ICTV-Virus-Knowledgebase/ICTVdatabase.git"
TEMP=$(mktemp -d) 
REPO_DIR="$TEMP/ICTVdatabase"

# Let the caller override it, otherwise derive it from the folder name.
# For current_msl/msl41v1 this becomes MSL41.v1
MSL_PREFIX="${MSL_PREFIX:-$(basename "$PWD" | sed -E 's/^msl([0-9]+)v([0-9]+)(_[0-9]{8})?$/MSL\1.v\2/')}"

# Find the newest tag that matches this MSL release.
TAG="$(
  git ls-remote --tags --refs "$REPO_URL" "refs/tags/${MSL_PREFIX}_*" \
    | awk '{print $2}' \
    | sed 's#refs/tags/##' \
    | sort -V \
    | tail -n 1
)"

if [ -z "$TAG" ]; then
  echo "ERROR: no tag found for prefix ${MSL_PREFIX}_" >&2
  exit 1
fi

git clone --branch "$TAG" --single-branch --depth 1 "$REPO_URL" "$REPO_DIR"

# copy here
for file in $(ls $REPO_DIR/data/*.utf8.txt); do
	dest_file=$(basename $file)
	echo "cp -a $file $dest_file"
	cp -a $file $dest_file
	git add $dest_file
done
