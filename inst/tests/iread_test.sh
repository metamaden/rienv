#!/usr/bin/bash

# Author: Sean Maden
# Test iREAD using readme example (https://github.com/genemine/iread/).
# Run after setting up iREAD env using conda.sh.

conda activate py2718_iread

#---------
# run test
#---------
export PATH="$HOME/iread:$PATH"
chmod -R 755 ./iread/
cd ./iread
python iread.py data/mouse_test.bam meta/intron_mouse_3875.bed -o \
    tmp_output -t 62000000