#!/bin/bash

INPUT_PREFIX=$1
IDs_THAT_PASSED_PhenoQC=$2
IDs_THAT_PASSED_SLQC=$3
SEX_CHECK=$4
OUTPUT=$5

# Combine the two passed IDs
COMBINED_KEEP="combined_keep.id"
awk 'NR==FNR {if(NR>1) a[$1,$2]; next} FNR==1 && /^#?FID/ {print; next} ($1,$2) in a' "$IDs_THAT_PASSED_PhenoQC" "$IDs_THAT_PASSED_SLQC" > "$COMBINED_KEEP"

grep "PROBLEM" $SEX_CHECK > $INPUT_PREFIX.txt

module add apps/plink2

CMD=(
    plink2
      --pfile $INPUT_PREFIX
      --keep $COMBINED_KEEP  # keep samples that passed both sample level QCs
      --remove $INPUT_PREFIX.txt       # remove samples that failed sex concordance
      --make-pgen
      --out $OUTPUT
)

"${CMD[@]}"
