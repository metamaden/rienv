#!/usr/bin/bash

# Author Sean Maden

# Set up a conda environment for retained introns.

conda update conda
conda install python=3.7.0
conda install python=2.7.18
conda install -c bioconda bioconductor-biocinstaller

#---------
# main env 
#---------
conda activate r403
conda install -c conda-forge r-base
conda install -c conda-forge/label/gcc7 r-base

install.packages("BiocManager")
install.packages("IntEREst")

# fixes for openmp error on sirfinder
# conda install nomkl
# conda install -c conda-forge libgomp
conda install -c conda-forge openmp



conda install -c anaconda llvm-openmp # sirfinder dependency

conda install -c conda-forge r-devtools # to install from github
conda install -c conda-forge r-rfast
conda install -c bioconda bioconductor-interest

R
devtools::install_github("lbroseus/SIRFindeR", build_vignette = TRUE)