#!/usr/bin/env bash
#
# tag and push our latest image to docker hub
#
# get version+git_hash
VER_FILE=version_git.txt
if [ ! -e "$VER_FILE" ]; then 
	echo "ERROR: $VER_FILE not found"
	exit 1
fi
VER_STR=$(cat $VER_FILE)

#
# tag
# 
echo sudo docker tag ictv_proposal_processor curtish/ictv_proposal_processor:$VER_STR
     sudo docker tag ictv_proposal_processor curtish/ictv_proposal_processor:$VER_STR
 
echo sudo docker tag ictv_proposal_processor curtish/ictv_proposal_processor:latest
     sudo docker tag ictv_proposal_processor curtish/ictv_proposal_processor:latest

#
# tag
#
echo sudo docker push curtish/ictv_proposal_processor:$VER_STR
#     sudo docker push curtish/ictv_proposal_processor:$VER_STR

echo sudo docker push curtish/ictv_proposal_processor:latest
#     sudo docker push curtish/ictv_proposal_processor:latest

#
# to pull instructions
#
echo "# to pull, run:"
echo "docker pull curtish/ictv_proposal_processor:$VER_STR"
echo "docker pull curtish/ictv_proposal_processor:latest"
