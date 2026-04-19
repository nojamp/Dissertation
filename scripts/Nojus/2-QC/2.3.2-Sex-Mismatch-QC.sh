#!/bin/bash

INPUT_BGEN=$1
INPUT_SAMPLE=$2
SEX_MISMATCH=$3
OUTPUT=$4

module add apps/plink2

plink2 \
  --bfile $INPUT_BGEN \
  --sample $INPUT_SAMPLE
  --remove $SEX_MISMATCH \
  ----export bgen-1.3 id-delim=: \
  --out $OUTPUT