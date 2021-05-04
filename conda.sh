#!/usr/bin/bash

# Author Sean Maden
#
# Set up conda environments for retained intron tools.
#

conda update conda
conda install python=3.7.0
conda install python=2.7.18
conda install r=3.5.1
conda install -c bioconda bioconductor-biocinstaller

#--------------------
# manage python3 envs
#--------------------
conda create -n py370 python=3.7.0
# clone for preprocessing
conda create --name py370_preprocess --clone py370

conda info --envs # show all envs

#-------------------
# manage python2 envs
#-------------------
conda create -n py2718 python=2.7.18
# clone for iread
conda create --name py2718_iread --clone py2718

conda info --envs # show all envs

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
conda create -n r351 r=3.5.1
# clone r v3### env for sirfinder
conda create --name r351_sirfinder --clone r351

conda info --envs # show all envs

#---------------------
# preprocess env setup
#---------------------
conda activate py370_preprocess
# dependencies
conda install -c bioconda samtools
conda install -c bioconda bowtie2
conda install -c bioconda star

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
conda install -c bioconda bioconductor-plyranges
conda install -c bioconda bioconductor-genomicfeatures
conda install -c conda-forge r-devtools
conda install -c conda-forge r-patchwork
conda install -c conda-forge r-ggplot2
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
conda install -c conda-forge argparse
conda install -c bioconda perl-parallel-forkmanager

# test iread
# export PATH="$HOME/iread:$PATH"
chmod -R 755 ./iread/
cd ./iread
python iread.py data/mouse_test.bam meta/intron_mouse_3875.bed -o \
    tmp_output -t 62000000

#----------------------
# SIRFindeR environment
#----------------------
conda activate r351_sirfinder
# dependencies
conda install -c conda-forge r-devtools # to install from github
conda install -c conda-forge r-rfast

R
devtools::install_github("lbroseus/SIRFindeR", build_vignette = TRUE)

# fixes for openmp error on sirfinder
# conda install nomkl
# conda install -c conda-forge libgomp
# conda install -c conda-forge openmp
# conda install -c anaconda llvm-openmp # sirfinder dependency






