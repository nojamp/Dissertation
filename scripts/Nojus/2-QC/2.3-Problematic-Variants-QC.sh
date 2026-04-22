#!/bin/bash

INPUT_PREFIX=$1
OUTPUT=$2

module add apps/plink2

# Remove multi-allelic variants
plink2 \
  --pfile $INPUT_PREFIX \
  --max-alleles 2 \
  --make-pgen \
  --out $OUTPUT.biallelic_only_temp

# Remove ambigous strand variants
awk 'NR>1 && (($4=="A" && $5=="T") || ($4=="T" && $5=="A") || \
    ($4=="G" && $5=="C") || ($4=="C" && $5=="G")) {print $3}' \
    $OUTPUT.biallelic_only_temp.pvar > $OUTPUT.ambiguous_snps.txt

plink2 \
  --pfile $OUTPUT.biallelic_only_temp \
  --exclude $OUTPUT.ambiguous_snps.txt \
  --make-pgen \
  --out $OUTPUT

# clean up temporary file
rm $OUTPUT.biallelic_only_temp.*  
