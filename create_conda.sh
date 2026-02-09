#!/usr/bin/env bash
#
# Create conda env for running R 
#
conda create -p ".conda/r_env" R

conda activate ./conda/r_env

# 
# install needed R packages
#
R -e "install.packages('optparse',dependencies=TRUE, repos='http://cran.rstudio.com/')"
R -e "install.packages('data.table',dependencies=TRUE, repos='http://cran.rstudio.com/')"
R -e "install.packages('yaml',dependencies=TRUE, repos='http://cran.rstudio.com/')"

#R -e "install.packages('tidyverse',dependencies=TRUE, repos='http://cran.rstudio.com/')"
R -e "install.packages('stringr',dependencies=TRUE, repos='http://cran.rstudio.com/')"
R -e "install.packages('dplyr',dependencies=TRUE, repos='http://cran.rstudio.com/')"
R -e "install.packages('readr',dependencies=TRUE, repos='http://cran.rstudio.com/')"

R -e "install.packages('readxl',dependencies=TRUE, repos='http://cran.rstudio.com/')"
R -e "install.packages('writexl',dependencies=TRUE, repos='http://cran.rstudio.com/')"

# fails (not available) on R 3.6.3 / ubuntu:20.04
#R -e "install.packages('DescTools',dependencies=TRUE, repos='http://cran.rstudio.com/')"

# read docx
R -e "install.packages('qdapTools',dependencies=TRUE, repos='http://cran.rstudio.com/')"
