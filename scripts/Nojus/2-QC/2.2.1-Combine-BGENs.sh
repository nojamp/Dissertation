#!/bin/bash

BGENs_DIR=$1
OUTPUT_DIR=$2

module load gcc/10.5.0-vi5y
module load bgen/1.1.7-kdez

cat-bgen -g $BGENs_DIR/*.bgen -og $OUTPUT_DIR/combined.bgen