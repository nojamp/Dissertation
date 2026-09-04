#!/bin/bash

INPUT_PREFIX=$1
IDs_THAT_PASSED_SLQC=$2
IDs_THAT_FAILED_KING=$3
SEX_CHECK=$4
OUTPUT=$5

grep "PROBLEM" $SEX_CHECK > $INPUT_PREFIX.txt

module add apps/plink2

CMD=(
    plink2
      --pfile $INPUT_PREFIX
      --keep $IDs_THAT_PASSED_SLQC     # keep samples that passed sample level QC
      --remove $INPUT_PREFIX.txt $IDs_THAT_FAILED_KING       # remove samples that failed sex concordance and King
      --make-pgen
      --out $OUTPUT
)

"${CMD[@]}"
