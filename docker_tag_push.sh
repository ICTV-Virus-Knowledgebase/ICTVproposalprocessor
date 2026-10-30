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
echo READING $VER_FILE
VER_STR=$(cat $VER_FILE)
echo VER_STR=$VER_STR

# get docker user from docker login info
EXPECTED_USER="ictvglobal"
DOCKER_USER=$(sudo docker info 2> /dev/null | awk '/Username:/{print $2}')
if [ -z "$DOCKER_USER" ]; then
	echo "ERROR: could not recover Username: from 'sudo docker info'"
	exit 1
fi
if [ "$DOCKER_USER" != "$EXPECTED_USER" ]; then
	echo "ERROR: expected docker user '$EXPECTED_USER' is not logged in; instead '$DOCKER_USER' is currently logged in"
	exit 1
fi
echo DOCKER_USER=$DOCKER_USER
echo ""

#
# tag
# 
echo "RUNNING (local tagging):"
echo sudo docker tag ictv_proposal_processor ${DOCKER_USER}/ictv_proposal_processor:$VER_STR
     sudo docker tag ictv_proposal_processor ${DOCKER_USER}/ictv_proposal_processor:$VER_STR
 
echo sudo docker tag ictv_proposal_processor ${DOCKER_USER}/ictv_proposal_processor:latest
     sudo docker tag ictv_proposal_processor ${DOCKER_USER}/ictv_proposal_processor:latest

echo ""
echo "SUGGEST (local tagging):"
echo sudo docker tag ictv_proposal_processor ${DOCKER_USER}/ictv_proposal_processor:prod

#
# tag
#
echo ""
echo "SUGGEST (push to docker hub):"
echo sudo docker push ${DOCKER_USER}/ictv_proposal_processor:$VER_STR
#     sudo docker push ${DOCKER_USER}/ictv_proposal_processor:$VER_STR
echo sudo docker push ${DOCKER_USER}/ictv_proposal_processor:latest
#     sudo docker push ${DOCKER_USER}/ictv_proposal_processor:latest
echo sudo docker push ${DOCKER_USER}/ictv_proposal_processor:prod
#     sudo docker push ${DOCKER_USER}/ictv_proposal_processor:prod

#
# to pull instructions
#
echo ""
echo "# to pull, run (on app.ictv.global, for example):"
echo "docker pull ${DOCKER_USER}/ictv_proposal_processor:$VER_STR"
echo "docker pull ${DOCKER_USER}/ictv_proposal_processor:latest"
echo "docker pull ${DOCKER_USER}/ictv_proposal_processor:prod"
