#!/usr/bin/bash

# Author: Sean Maden
#
# Set up conda environments for retained intron detection tools.
# 
# Notes:
# 
# * Setup is expedited by installing shared dependencies once and 
#   cloning the consensus environment for specific tools.
#
# * At any time, use `conda info --envs` to show all available 
#   environments.
# 

conda update conda
conda install python=3.7.0
conda install python=2.7.18
conda install r=3.5.1
conda install -c bioconda bioconductor-biocinstaller

#-------------------
# manage python envs
#-------------------

# python3 envs
conda create -n py370 python=3.7.0
# clone for preprocessing
conda create --name py370_preprocess --clone py370

# python 2 envs
conda create -n py2718 python=2.7.18
# clone for iread
conda create --name py2718_iread --clone py2718

#--------------
# manage r envs
#--------------

# make r v4### env
conda create -n r403; conda activate r403
conda install -c conda-forge r-base
conda install -c conda-forge/label/gcc7 r-base
conda deactivate
# clone r v4### env for interest and superintronic
conda create --name r403_interest --clone r403
conda create --name r403_superintronic --clone r403

# make r v3### env
conda create -n r351 r=3.5.1; conda attach r351
conda install -c conda-forge/label/gcc7 r-devtools=2.0.1
conda install -c conda-forge/label/gcc7 r-data.table=1.14.0
# clone r v3### env for sirfinder, kma
conda create --name r351_sirfinder --clone r351
conda create --name r351_kma --clone r351

#---------------------
# preprocess env setup
#---------------------
conda activate py370_preprocess
# dependencies
conda install -c bioconda samtools=1.3.1
conda install -c bioconda bowtie2=2.3.4.3
conda install -c bioconda star=2.7.6a

conda deactivate

#-------------------
# interest env setup
#-------------------
conda activate r403_interest
# install tool
conda install -c bioconda bioconductor-interest

conda deactivate

#------------------------
# superintronic env setup
#------------------------
conda activate r403_superintronic
# dependencies
conda install -c bioconda bioconductor-plyranges=1.10.0
conda install -c bioconda bioconductor-genomicfeatures=1.42.2
conda install -c conda-forge r-devtools=2.4.0
conda install -c conda-forge r-patchwork=1.1.1
conda install -c conda-forge r-ggplot2=3.3.3
# install superintronic
R
devtools::install_github("sa-lee/superintronic", build_vignette = FALSE)

conda deactivate

#----------------
# iread env setup
#----------------
conda activate py2718_iread
conda install perl=5.26.2
git clone https://github.com/genemine/iread/
# iread dependencies
conda install -c bioconda samtools=1.2
conda install -c bioconda bedops=2.4.20
conda install -c conda-forge argparse=1.4.0
conda install -c bioconda perl-parallel-forkmanager=2.02

#--------------------
# sirfinder env setup
#--------------------
conda activate r351_sirfinder
# dependencies
# conda install -c conda-forge/label/gcc7 r-devtools # to install from github
conda install -c conda-forge/label/gcc7 r-rfast=1.9.5
conda install -c conda-forge/label/gcc7 r-rlang=0.4.10

R 
devtools::install_github("lbroseus/SIRFindeR", build_vignette = TRUE)

# possible fixes for openmp error on sirfinder
# note: may only be applicable for R v4###
# conda install nomkl
# conda install -c conda-forge libgomp
# conda install -c conda-forge openmp
# conda install -c anaconda llvm-openmp # sirfinder dependency

#--------------
# kma env setup
#--------------
conda activate r351_kma
# dependencies
conda install -c conda-forge/label/gcc7 r-reshape2=1.4.3
conda install -c conda-forge/label/gcc7 r-dplyr=0.7.8
# get kma
R
devtools::install_github("https://github.com/adamtongji/kma")





