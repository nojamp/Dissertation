#!/bin/bash

INPUT_BGEN=$1
INPUT_SAMPLE=$2
OUTPUT_DIR=$3
OUTPUT=$4

module add apps/plink2

plink2 \
  --bgen $INPUT_BGEN ref-unknown \
  --sample $INPUT_SAMPLE \
  --max-alleles 2 \
  --export bgen-1.3 id-delim=: \
  --out $OUTPUT_DIR/biallelic_only_temp

# Export variant information including alleles
plink2 \
  --bgen $OUTPUT_DIR/biallelic_only_temp.bgen ref-unknown \
  --sample $OUTPUT_DIR/biallelic_only_temp.sample \
  --write-alleles \
  --out $OUTPUT_DIR/variant_info_temp

# This creates variant_info.alleles with format: CHR SNP POS A1 A2
awk '($4 == "A" && $5 == "T") || ($4 == "T" && $5 == "A") || \
     ($4 == "G" && $5 == "C") || ($4 == "C" && $5 == "G") \
     {print $2}' $OUTPUT_DIR/variant_info_temp.alleles > $OUTPUT_DIR/ambiguous_snps.txt

# Now exclude ambiguous variants
plink2 \
  --bgen $OUTPUT_DIR/biallelic_only_temp.bgen ref-unknown \
  --sample $OUTPUT_DIR/biallelic_only_temp.sample \
  --exclude $OUTPUT_DIR/ambiguous_snps.txt \
  --export bgen-1.3 id-delim=: \
  --out $OUTPUT

rm $OUTPUT_DIR/biallelic_only_temp*  
rm $OUTPUT_DIR/variant_info_temp.alleles
