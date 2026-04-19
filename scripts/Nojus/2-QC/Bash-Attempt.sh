#!/bin/bash

#SBATCH --job-name=QC
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --time=0:10:00
#SBATCH --mem=20GB
#SBATCH --array=21-23
#SBATCH --account=sscm037704

module add apps/plink2

RAW="../../data/raw/Nojus/genetic_data_raw"
QC="../../intermediate/2-Genetic-Data-QC"

# to handle file names that use 0X numbering system
PADDED_ID=$(printf "%02d" $SLURM_ARRAY_TASK_ID)

# Set up output variables
QC_VL="$QC/2.1-VL/VL_${PADDED_ID}"

QC_SL="$QC/2.2-SL/SL_${PADDED_ID}"
QC_SL_23="$QC/2.2-SL/SL_23"

QC_SEX_DIR="$QC/2.3-Sex-Mismatch"
QC_SEX_OUTPUT="$QC_SEX_DIR/Sex_Mismatch_${PADDED_ID}"

QC_PV_DIR="$QC/2.4-PV"
QC_PV_OUTPUT="$QC_PV_DIR/PV_${PADDED_ID}"

QC_PRUNING_DIR="$QC/2.5-Pruning"
QC_PRUNING_OUTPUT="$QC_PRUNING_DIR/Pruned_${PADDED_ID}"

# 2.1 Variant Level QC
. 2-Genotype-QC/2.1-Variant-Level-QC.sh $RAW/filtered_${PADDED_ID}.bgen $RAW/swapped.sample $QC_VL

# 2.2 Sample Level QC
. 2-Genotype-QC/2.2-Sample-Level-QC.sh $QC_VL.bgen $QC_VL.sample $QC_SL

# 2.3 Sex Mismatch QC
. 2-Genotype-QC/2.3-Sex-Mismatch-QC.sh $QC_SL.bgen $QC_SL_23.bgen $QC_SL.sample $QC_SEX_DIR $QC_SEX_OUTPUT

# 2.4 Problematic Variants QC
. 2-Genotype-QC/2.4-Problematic-Variants-QC.sh $QC_SEX_OUTPUT.bgen $QC_SEX_OUTPUT.sample $QC_PV_DIR $QC_PV_OUTPUT

# 2.5 Pruning
. 2-Genotype-QC/2.5-Pruning.sh $QC_PV_OUTPUT.bgen $QC_PV_OUTPUT.sample $QC_PRUNING_DIR $QC_PRUNING_OUTPUT
