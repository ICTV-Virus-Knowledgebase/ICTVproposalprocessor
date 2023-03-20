#
# Docker image to run proposal.xlsx validation
#
# Includes reference data
#

# Ubunutu base
FROM rocker/tidyverse:4

# make sure installs dont hang on user input
ENV DEBIAN_FRONTEND noninteractive


# install several R packages
#RUN R -e "install.packages('optparse',dependencies=TRUE, repos='http://cran.rstudio.com/')"
#RUN R -e "install.packages('data.table',dependencies=TRUE, repos='http://cran.rstudio.com/')"
#RUN R -e "install.packages('yaml',dependencies=TRUE, repos='http://cran.rstudio.com/')"
#RUN R -e "install.packages('readxl',dependencies=TRUE, repos='http://cran.rstudio.com/')"

RUN R -e "install.packages('writexl',dependencies=TRUE, repos='http://cran.rstudio.com/')"
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
