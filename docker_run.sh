#!/usr/bin/env bash
#
# test run our docker image
#
sudo docker run -it \
	-v "$(pwd)/proposalsTest:/proposalsTest":ro \
	-v "$(pwd)/results_new:/results" \
	ictv_proposal_processor  \
	/merge_proposal_zips.R -v
