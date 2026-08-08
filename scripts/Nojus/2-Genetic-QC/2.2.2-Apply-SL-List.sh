#!/bin/bash

INPUT_PREFIX=$1
IDs_THAT_PASSED_SLQC=$2
SEX_CHECK=$3
OUTPUT=$4

grep "PROBLEM" $SEX_CHECK > $INPUT_PREFIX.txt

module add apps/plink2

CMD=(
    plink2
      --pfile $INPUT_PREFIX
      --keep $IDs_THAT_PASSED_SLQC  # keep samples that passed both sample level QCs
      --remove $INPUT_PREFIX.txt       # remove samples that failed sex concordance
      --make-pgen
      --out $OUTPUT
)

"${CMD[@]}"
