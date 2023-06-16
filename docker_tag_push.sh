#!/usr/bin/env bash
#
# tag and push our latest image to docker hub
#
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!!!  need to add versions/real tagso   !!!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
sudo docker tag ictv_proposal_processor curtish/ictv_proposal_processor
sudo docker  push curtish/ictv_proposal_processor
