#!/bin/bash

INPUT=$1
INPUT_SAMPLE=$3
OUTPUT_DIR=$5

module add apps/plink2

plink2 \
  --bgen $INPUT_CHR23 ref-unknown \
  --sample $INPUT_SAMPLE \
  --check-sex \
  --out $OUTPUT_DIR/sex_check

awk '$3 != $4 {print $0}' $OUTPUT_DIR/sex_check.sexcheck > $OUTPUT_DIR/sex_mismatches.txt
