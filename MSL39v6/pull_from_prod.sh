#!/usr/bin/env bash
#
# get proposals from ICTV's PROD Drupal server
#
# exclude anything with "2024", specifically the directories named
#   "Proposals to be considered for EC56 (2024)"
# -i ~/.ssh/ICTV-curtish.pem \
rsync \
      -hav \
      --progress \
      --exclude '*EC56*' \
      $* \
      ubuntu@ictv.global:/var/www/drupal/files/private/ICTV_files/Proposals/Pending_Proposals \
      .

