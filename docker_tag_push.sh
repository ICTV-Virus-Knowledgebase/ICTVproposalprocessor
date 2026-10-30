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
echo "RUNNING (local tagging):"
echo sudo docker tag ictv_proposal_processor ictvglobal/ictv_proposal_processor:$VER_STR
     sudo docker tag ictv_proposal_processor ictvglobal/ictv_proposal_processor:$VER_STR
 
echo sudo docker tag ictv_proposal_processor ictvglobal/ictv_proposal_processor:latest
     sudo docker tag ictv_proposal_processor ictvglobal/ictv_proposal_processor:latest

echo ""
echo "SUGGEST (local tagging):"
echo sudo docker tag ictv_proposal_processor ictvglobal/ictv_proposal_processor:prod

#
# tag
#
echo ""
echo "SUGGEST (push to docker hub):"
echo sudo docker push ictvglobal/ictv_proposal_processor:$VER_STR
#     sudo docker push ictvglobal/ictv_proposal_processor:$VER_STR
echo sudo docker push ictvglobal/ictv_proposal_processor:latest
#     sudo docker push ictvglobal/ictv_proposal_processor:latest
echo sudo docker push ictvglobal/ictv_proposal_processor:prod
#     sudo docker push ictvglobal/ictv_proposal_processor:prod

#
# to pull instructions
#
echo ""
echo "# to pull, run (on app.ictv.global, for example):"
echo "docker pull ictvglobal/ictv_proposal_processor:$VER_STR"
echo "docker pull ictvglobal/ictv_proposal_processor:latest"
echo "docker pull ictvglobal/ictv_proposal_processor:prod"
