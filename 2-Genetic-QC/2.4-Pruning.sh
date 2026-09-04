#!/bin/bash

INPUT_DIR=$1
OUTPUT_UNPRUNED=$2
OUTPUT_PRUNED=$3

module add apps/plink2

# merge all genetic data ready for pruning
# keep unpruned version for further analysis
ls $INPUT_DIR/ | sed 's/\.[^.]*$//' | uniq | grep "PV" | grep -v "ambiguous_snps" > $INPUT_DIR/file_names.txt

plink2 \
  --pmerge-list $INPUT_DIR/file_names.txt \
  --pmerge-list-dir $INPUT_DIR \
  --set-all-var-ids "@:#_\$r_\$a" \
  --new-id-max-allele-len 50 truncate \
  --rm-dup exclude-all \
  --make-pgen \
  --out $OUTPUT_UNPRUNED

# Create LD-pruned set for PCA
plink2 \
  --pfile $OUTPUT_UNPRUNED \
  --indep-pairwise 50 5 0.8 \
  --out $OUTPUT_PRUNED

# Extract pruned variants
plink2 \
  --pfile $OUTPUT_UNPRUNED \
  --extract $OUTPUT_PRUNED.prune.in \
  --make-pgen \
  --out $OUTPUT_PRUNED