#
# Docker image to run proposal.xlsx validation
#
# Includes reference data
#

# Ubunutu base
FROM ubuntu:22.04

# make sure installs dont hang on user input
ENV DEBIAN_FRONTEND=noninteractive


# install r-base and pre-requisitis
# installs R 3.6.3 (on ubuntu:20.04)
RUN set -e && apt-get -y update
RUN set -e && apt-get -y dist-upgrade 
RUN set -e && apt-get -y install --no-install-recommends --no-install-suggests apt-transport-https #=2.4.10
RUN set -e && apt-get -y install --no-install-recommends --no-install-suggests apt-utils #=2.4.10
RUN set -e && apt-get -y install --no-install-recommends --no-install-suggests ca-certificates
RUN set -e && apt-get -y install --no-install-recommends --no-install-suggests cmake #=3.22.1
RUN set -e && apt-get -y install --no-install-recommends --no-install-suggests curl #=7.81.0
RUN set -e && apt-get -y install --no-install-recommends --no-install-suggests g++ #=4:11.2.0
RUN set -e && apt-get -y install --no-install-recommends --no-install-suggests gcc #=4:11.2.0
RUN set -e && apt-get -y install --no-install-recommends --no-install-suggests gfortran #=4:11.2.0
RUN set -e && apt-get -y install --no-install-recommends --no-install-suggests git #=1:2.34.1
RUN set -e && apt-get -y install --no-install-recommends --no-install-suggests make #=4.3
RUN set -e && apt-get -y install --no-install-recommends --no-install-suggests libblas-dev #=3.10.0
RUN set -e && apt-get -y install --no-install-recommends --no-install-suggests libcurl4-gnutls-dev #=7.81.0
RUN set -e && apt-get -y install --no-install-recommends --no-install-suggests libfontconfig1-dev #=2.13.1
RUN set -e && apt-get -y install --no-install-recommends --no-install-suggests libfreetype6-dev #=2.11.1
RUN set -e && apt-get -y install --no-install-recommends --no-install-suggests libfribidi-dev #=1.0.8
RUN set -e && apt-get -y install --no-install-recommends --no-install-suggests libgit2-dev #=1.1.0
RUN set -e && apt-get -y install --no-install-recommends --no-install-suggests libharfbuzz-dev #=2.7.4
RUN set -e && apt-get -y install --no-install-recommends --no-install-suggests libiodbc2-dev #=3.52.9
RUN set -e && apt-get -y install --no-install-recommends --no-install-suggests libjpeg-dev #=8c
RUN set -e && apt-get -y install --no-install-recommends --no-install-suggests liblapack-dev #=3.10.0
RUN set -e && apt-get -y install --no-install-recommends --no-install-suggests libmariadb-dev #=1:10.6.12
RUN set -e && apt-get -y install --no-install-recommends --no-install-suggests libpng-dev #=1.6.37
RUN set -e && apt-get -y install --no-install-recommends --no-install-suggests libpq-dev #=14.9
RUN set -e && apt-get -y install --no-install-recommends --no-install-suggests libsqlite3-dev #=3.37.2
RUN set -e && apt-get -y install --no-install-recommends --no-install-suggests libssh-dev #=0.9.6
RUN set -e && apt-get -y install --no-install-recommends --no-install-suggests libssl-dev #=3.0.2
RUN set -e && apt-get -y install --no-install-recommends --no-install-suggests libtiff5-dev #=4.3.0
RUN set -e && apt-get -y install --no-install-recommends --no-install-suggests libxml2-dev #=2.9.13
RUN set -e && apt-get -y install --no-install-recommends --no-install-suggests locales #=2.35
RUN set -e && apt-get -y install --no-install-recommends --no-install-suggests pandoc #=2.9.2.1
RUN set -e && apt-get -y install --no-install-recommends --no-install-suggests pkg-config #=0.29.2
RUN set -e && apt-get -y install --no-install-recommends --no-install-suggests r-base #=4.1.2
RUN set -e && apt-get -y autoremove
RUN set -e && apt-get clean
RUN set -e && rm -rf /var/lib/apt/lists

# install several R packages
RUN /usr/bin/R -e "install.packages('optparse',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN /usr/bin/R -e "install.packages('data.table',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN /usr/bin/R -e "install.packages('yaml',dependencies=TRUE, repos='http://cran.rstudio.com/')"

#RUN R -e "install.packages('tidyverse',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN /usr/bin/R -e "install.packages('stringr',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN /usr/bin/R -e "install.packages('dplyr',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN /usr/bin/R -e "install.packages('tidyr',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN /usr/bin/R -e "install.packages('readr',dependencies=TRUE, repos='http://cran.rstudio.com/')"

# readxl 1.6.0 requires R >= 4.2; ubuntu:22.04 has R 4.1.2
RUN /usr/bin/R -e "install.packages(c('rematch','cellranger'), repos='http://cran.rstudio.com/')"
RUN /usr/bin/R -e "install.packages('https://cran.rstudio.com/src/contrib/Archive/readxl/readxl_1.5.0.1.tar.gz', repos=NULL, type='source')"
RUN /usr/bin/R -e "install.packages('writexl',dependencies=TRUE, repos='http://cran.rstudio.com/')"
# fails (not available) on R 3.6.3 / ubuntu:20.04
#RUN /usr/bin/R -e "install.packages('DescTools',dependencies=TRUE, repos='http://cran.rstudio.com/')"
# # read docx
RUN /usr/bin/R -e "install.packages('qdapTools',dependencies=TRUE, repos='http://cran.rstudio.com/')"

# install.packages() only warns on failure, so fail the build here
RUN /usr/bin/Rscript -e "pkgs <- c('data.table','dplyr','optparse','qdapTools','readr','readxl','stringi','stringr','tidyr','writexl','XML'); ok <- vapply(pkgs, requireNamespace, logical(1), quietly=TRUE); if (!all(ok)) stop('missing R packages: ', paste(pkgs[!ok], collapse=', '))"

# UTF-8 mode
RUN set -e \
      && locale-gen en_US.UTF-8 \
      && update-locale

#
# copy in our application
#
# do this as a git clone, instead!?!?
COPY merge_proposal_zips.R .
COPY abolished_accessions.tsv .
COPY version_git.txt .

#
# copy in reference data
#
RUN mkdir -p ./current_msl
COPY current_msl/ ./current_msl

# build cache of reference data
#RUN ./merge_proposal_zips.R --buildCache --ref ./current_msl

# what does ENTRYPOINT do exactly?
#ENTRYPOINT ["merge_proposal_zips.R"
CMD [ "./merge_proposal_zips.R" ]
