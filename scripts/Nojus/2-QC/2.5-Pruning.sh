#!/bin/bash

INPUT_BGEN=$1
INPUT_SAMPLE=$2
OUTPUT=$3
OUTPUT_DIR=$4

module add apps/plink2

# Create LD-pruned set for PCA
plink2 \
  --bgen $INPUT_BGEN ref-unknown \
  --sample $INPUT_SAMPLE \
  --indep-pairwise 50 5 0.8 \
  --out $OUTPUT_DIR/ld_pruned
   
# Extract pruned variants
plink2 \
  --bgen $INPUT_BGEN ref-unknown \
  --sample $INPUT_SAMPLE \
  --extract $OUTPUT_DIR/ld_pruned.prune.in \
  --export bgen-1.3 id-delim=: \
  --out $OUTPUT