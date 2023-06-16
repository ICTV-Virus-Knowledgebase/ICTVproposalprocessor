#!/usr/bin/env bash
#
# remove all images with "processor" in the name
#
# with prejudice
#
sudo docker images -a | grep "processor" | awk '{print $3}' | xargs sudo docker rmi -f
