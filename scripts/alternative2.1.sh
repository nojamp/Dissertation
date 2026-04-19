#!/bin/bash

INPUT_BGEN=$1
INPUT_SAMPLE=$2
OUTPUT=$3

module add apps/plink2
CMD=(
  plink2
    --bgen $INPUT_BGEN ref-first
    --sample $INPUT_SAMPLE
    --maf 0.005                   # MAF >= 0.5%
    --geno 0.05                   # Variant call rate > 95%
    --hwe 1e-6                    # Hardy-Weinberg equilibrium p > 1e-6
    --rm-dup exclude-all          # Removes duplicates
    --export bgen-1.3 id-delim=:  # Save as new BGEN file
    --out $OUTPUT
)

"${CMD[@]}"