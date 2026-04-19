#!/bin/bash

INPUT_BGEN=$1
INPUT_SAMPLE=$2
REMOVE_SAMPLE=$3
OUTPUT=$4

module add apps/plink2

plink2 \
  --bfile $INPUT_BGEN \
  --sample $INPUT_SAMPLE
  --remove $REMOVE_SAMPLE \
  ----export bgen-1.3 id-delim=: \
  --out $OUTPUT