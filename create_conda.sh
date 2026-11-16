#!/usr/bin/env bash
#
# Create the conda env for running the R scripts: ./conda/r_env
#
# USAGE: bash create_conda.sh   (from anywhere; the env goes in the repo root)
#
# All packages come from conda-forge, so nothing is compiled from CRAN source.
# R is pinned to 4.3.x: the testResults baselines were made with R 4.3.1 on
# Linux, but on macOS arm64 R 4.3.1 can't be installed with qdapTools (its
# dependency r-xml needs icu 75, and the r-base 4.3.1 builds need icu 72-73),
# so the solver picks 4.3.3 there.
#
# The Python env is separate: conda env create -f environment.python.yml
#
set -euo pipefail
cd "$(dirname "$0")"

conda create -p ./conda/r_env --override-channels -c conda-forge \
    "r-base=4.3" \
    r-optparse r-data.table r-stringr r-dplyr r-tidyr r-readr \
    r-readxl r-writexl \
    r-qdaptools \
    r-jsonlite r-rjson

# r-qdaptools: read docx
# r-jsonlite, r-rjson: only for the legacy load_json_msl.R

echo "# Created ./conda/r_env; check it with:"
echo "#   conda run -p ./conda/r_env Rscript --version"
