#!/bin/bash

NUM_PROCESSORS=8
INPUT_DIR_ROOT="/Volumes/Data2/jtblair/20150303_100m_batchrun_2/batch_"
INPUT_FILE="25k_minbias.root"
OUTPUT_DIR="/Users/jtblair/invmass_test/"
OUTPUT_FILE_ROOT="histo_"
FIRST=0
LAST=3999
WRITE_FILE="/Users/jtblair/output12.txt"
RUN_FILE="/Users/jtblair/utaustin/invMassGenerator.o"

nohup mpirun -np $NUM_PROCESSORS $RUN_FILE $INPUT_DIR_ROOT $INPUT_FILE $OUTPUT_DIR $OUTPUT_FILE_ROOT $NUM_PROCESSORS $FIRST $LAST > $WRITE_FILE &
