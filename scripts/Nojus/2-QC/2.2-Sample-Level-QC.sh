#!/bin/bash

INPUT_BGEN=$1
INPUT_SAMPLE=$2
OUTPUT=$3

module add apps/plink2

plink2 \
  --bgen $INPUT_BGEN ref-first \
  --sample $INPUT_SAMPLE \
  --mind 0.03 \ 		            # sample call rate > 97%
##--king_cutoff 0.0884 \        # if wanting to filter out 2nd cousins or closer
  --export bgen-1.3 id-delim=:\
  --out $OUTPUT