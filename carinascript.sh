#!/bin/bash

NUM_PROCESSORS=8
INPUT_DIR_ROOT="/Volumes/Data2/jtblair/20150303_100m_batchrun_2/batch_"
INPUT_FILE="25k_minbias.root"
OUTPUT_DIR="/Users/jtblair/histogram_test/"
OUTPUT_FILE_ROOT="histo_"
FIRST=36
LAST=1059
WRITE_FILE="/Users/jtblair/output6.txt"

nohup mpirun -np $NUM_PROCESSORS /Users/jtblair/utaustin/histoGenerator.o $INPUT_DIR_ROOT $INPUT_FILE $OUTPUT_DIR $OUTPUT_FILE_ROOT $NUM_PROCESSORS $FIRST $LAST > $WRITE_FILE &