#
# Docker image to run proposal.xlsx validation
#
# Includes reference data
#

# Ubunutu base
FROM ubuntu:22.04

# make sure installs dont hang on user input
ENV DEBIAN_FRONTEND noninteractive


# install r-base and pre-requisitis
# installs R 3.6.3 (on ubuntu:20.04)
RUN set -e \
      && apt-get -y update \
      && apt-get -y dist-upgrade \
      && apt-get -y install --no-install-recommends --no-install-suggests \
        apt-transport-https apt-utils ca-certificates cmake curl g++ gcc \
        gfortran git make libblas-dev libcurl4-gnutls-dev libfontconfig1-dev \
        libfreetype6-dev libfribidi-dev libgit2-dev libharfbuzz-dev \
        libiodbc2-dev libjpeg-dev liblapack-dev libmariadb-dev libpng-dev \
        libpq-dev libsqlite3-dev libssh-dev libssl-dev libtiff5-dev \
        libxml2-dev locales pandoc pkg-config r-base \
      && apt-get -y autoremove \
      && apt-get clean \
      && rm -rf /var/lib/apt/lists/*

# install several R packages
RUN R -e "install.packages('optparse',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('data.table',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('yaml',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('tidyverse',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('readxl',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('writexl',dependencies=TRUE, repos='http://cran.rstudio.com/')"
# fails (not available) on R 3.6.3 / ubuntu:20.04
RUN R -e "install.packages('DescTools',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('qdapTools',dependencies=TRUE, repos='http://cran.rstudio.com/')"

# UTF-8 mode
RUN set -e \
      && locale-gen en_US.UTF-8 \
      && update-locale

#
# copy in our application
#
COPY merge_proposal_zips.R .

#
# copy in reference data
#
RUN mkdir -p ./current_msl
COPY current_msl/* ./current_msl

# build cache of reference data
#RUN ./merge_proposal_zips.R --buildCache --ref ./current_msl

# what does ENTRYPOINT do exactly?
#ENTRYPOINT ["merge_proposal_zips.R"
CMD [ "./merge_proposal_zips.R" ]
