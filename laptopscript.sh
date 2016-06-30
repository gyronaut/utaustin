#!/bin/bash

NUM_PROCESSORS=1
INPUT_DIR_ROOT="/Users/jtblair/Documents/minbias/20150303_100m_batchrun_2/batch_"
INPUT_FILE="25k_minbias.root"
OUTPUT_DIR="/Users/jtblair/K0/"
OUTPUT_FILE_ROOT="histo_"
FIRST=0
LAST=10
WRITE_FILE="/Users/jtblair/output_K0.txt"
RUN_FILE="/Users/jtblair/utaustin/K0histo.o"

$RUN_FILE $INPUT_DIR_ROOT $INPUT_FILE $OUTPUT_DIR $OUTPUT_FILE_ROOT $FIRST $LAST
