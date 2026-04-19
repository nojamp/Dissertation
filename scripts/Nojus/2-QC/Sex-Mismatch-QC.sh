#!/bin/bash

INPUT_BGEN=$1
INPUT_CHR23=$2
INPUT_SAMPLE=$3
OUTPUT_DIR=$4
OUTPUT_NAME=$5

module add apps/plink2

plink2 \
  --bgen $INPUT_CHR23 ref-first \
  --sample $INPUT_SAMPLE \
  --check-sex \
  --out $OUTPUT_DIR/sex_check

awk '$3 != $4 {print $0}' $OUTPUT_DIR/sex_check.sexcheck > $OUTPUT_DIR/sex_mismatches.txt

plink2 \
  --bfile $INPUT_BGEN \
  --sample $INPUT_SAMPLE
  --remove sex_mismatches.txt \
  --export bgen-1.3 id-delim=: \
  --out $OUTPUT_NAME