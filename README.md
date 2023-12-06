# MSL36_load

## WARNING: ./docker_build.sh broken

The last working image built was: 

	REPOSITORY                                       TAG             IMAGE ID       CREATED         SIZE
	curtish/ictv_proposal_processor                  v1.14.fca46f2   2ccf634e42d5   4 weeks ago     1.44GB

For some reason, our Dockerfile no longer builds a container that can load the R library DescTools. 

This *should* be solved by fixing Dockerfile to better specify versions, so we have a reproducible build. 

As a short-cut in the mean time, we load that image into a container, updaste the R script and version_git.txt, and save it as a new image:

```
## SRC
  * https://stackoverflow.com/questions/56551276/how-to-edit-docker-image
  * https://stackoverflow.com/questions/23935141/how-to-copy-docker-images-from-one-host-to-another-without-using-a-repository

ssh-ictv-test

# start container
ubuntu@test:/var/www/drupal/apps/proposal_validator/taxon_proposal_qc_and_load$ sudo docker run  -v $PWD:/SRC 2ccf634e42d5 sleep 500

# find container
ubuntu@test:/var/www/drupal/apps/proposal_validator/taxon_proposal_qc_and_load$ sudo docker ps
CONTAINER ID   IMAGE          COMMAND       CREATED              STATUS              PORTS     NAMES
fd985326b841   2ccf634e42d5   "sleep 500"   4 seconds ago        Up 3 seconds                  angry_lamarr

# connect into container
ubuntu@test:/var/www/drupal/apps/proposal_validator/taxon_proposal_qc_and_load$ sudo docker exec -it fd985326b841 bash

# update files in container
root@fd985326b841:/# cp /SRC/merge_proposal_zips.R /merge_proposal_zips.R 
root@fd985326b841:/# cp /SRC/version_git.txt /version_git.txt 

# write modified container to a new image
ubuntu@test:/var/www/drupal/apps/proposal_validator/taxon_proposal_qc_and_load$ sudo docker commit fd985326b841 ictv_proposal_processor/v1.15.119e2b3
sha256:801a872106bc81e9913e17a7cefe5a3bcd7ffc41c03a8f4ec33d8d48aaeb5772

# tag image
ubuntu@test:/var/www/drupal/apps/proposal_validator/taxon_proposal_qc_and_load$ ./docker_tag_push.sh 
sudo docker tag ictv_proposal_processor curtish/ictv_proposal_processor:v1.15.119e2b3
sudo docker tag ictv_proposal_processor curtish/ictv_proposal_processor:latest
# copy-paste-run these: 
# push image to hub.docker:
sudo docker push curtish/ictv_proposal_processor:v1.15.119e2b3
sudo docker push curtish/ictv_proposal_processor:latest
# to pull, run on target machine (APP):
docker pull curtish/ictv_proposal_processor:v1.15.119e2b3
docker pull curtish/ictv_proposal_processor:latest

# save image to disk/git
ubuntu@test:/var/www/drupal/apps/proposal_validator/taxon_proposal_qc_and_load$ mkdir -p images/curtish/ictv_proposal_processor
ubuntu@test:/var/www/drupal/apps/proposal_validator/taxon_proposal_qc_and_load$ sudo docker save curtish/ictv_proposal_processor:v1.14.fca46f2 | gzip -c > images/curtish/ictv_proposal_processor/v1.14.fca46f2.tar.gz

```
