#!/bin/bash

INPUT_BGEN=$1
INPUT_SAMPLE=$2
OUTPUT=$3

module add apps/plink2

# sample call rate > 97%
plink2 \
  --bgen $INPUT_BGEN ref-first \
  --sample $INPUT_SAMPLE \
  --mind 0.03 \
  --write-samples \
  --out $OUTPUT

# delete files that won't be used again (keeping only the sample file)
rm $INPUT_BGEN