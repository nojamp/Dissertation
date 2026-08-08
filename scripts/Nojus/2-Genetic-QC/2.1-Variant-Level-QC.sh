#!/bin/bash

INPUT_BGEN=$1
INPUT_SAMPLE=$2
INPUT_ELIGIBLE_SAMPLE=$3
INPUT_CHR=$4
OUTPUT=$5

module add apps/plink2
CMD=(
  plink2
    --bgen $INPUT_BGEN ref-first snpid-chr
    --sample $INPUT_SAMPLE
    --lax-bgen-import
    --set-missing-var-ids @:#\$r,\$a
    --new-id-max-allele-len 100
    --keep $INPUT_ELIGIBLE_SAMPLE
    --maf 0.005                     # MAF >= 0.5%
    --geno 0.05                     # Variant call rate > 95%
    --hwe 1e-6                      # Hardy-Weinberg equilibrium p > 1e-6
    --make-pgen                     # Save as PLINK2 native files
    --out $OUTPUT
)

"${CMD[@]}"