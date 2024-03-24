#!/usr/bin/env bash
#
# push final proposal files to box
#

DEST="box:/Virus Knowledgebase/taxonomy/ICTV_update/2024_updates/2024-03-17 MSL39/MSL39v5"

for SRC in Pending_Proposals proposalsFinal proposalsFinalZips; do
    echo rclone sync --progress "$SRC"  "$DEST/$SRC"
    rclone sync --progress "$SRC"  "$DEST/$SRC"
done
