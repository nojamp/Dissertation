#!/bin/bash

INPUT_BGEN=$1
INPUT_SAMPLE=$2
OUTPUT=$3

module add apps/plink2

plink2 \
  --bgen $INPUT_BGEN ref-unknown \
  --sample $INPUT_SAMPLE \
  --maf 0.005 \                     # MAF >= 0.5%
  --geno 0.05 \                     # Variant call rate > 95%
  --hwe 1e-6 \                      # Hardy-Weinberg equilibrium p > 1e-6
  --write-snplist allow-dups \      # records how many SNPs we have
  --make-bgen \
  --out $OUTPUT