#!/bin/bash

INPUT_DIR=$1
OUTPUT=$2

# create list of names of pgen files
ls $INPUT_DIR/ | sed 's/\.[^.]*$//' | uniq | grep "VL" > $INPUT_DIR/file_names.txt

module add apps/plink2

CMD=(
    plink2
      --pmerge-list $INPUT_DIR/file_names.txt         # merge pgens for sample level qc
      --pmerge-list-dir $INPUT_DIR
      --set-all-var-ids @:#\$r,\$a
      --new-id-max-allele-len 100 missing
      --mind 0.03                                     # filter samples with more than 3% missingness
      --king-cutoff 0.0884                            # filter for relatedness
      --write-samples                                 # write samples that pass the threshold
      --check-sex max-female-xf=0.2 min-male-xf=0.8   # check sex concordance using standard thresholds
      --delete-pmerge-result                          # remove merged pgens (we only need id outputs)
      --out $OUTPUT
)

"${CMD[@]}"
