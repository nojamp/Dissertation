# -*- mode: snakemake -*-

# Directories
RAW_DIR = "../../data/raw/Nojus"

QC_DIR = "../../data/intermediate/2-Genetic-Data-QC"
INTERMEDIATE_DIR = "../../data/intermediate"

RESULTS_DIR = "../../results"

OMICSPRED_DIR = "/user/work/cd24895/dissertation/working/scripts/Nojus/4-PRS/Nightingale"

PRS_TYPES = {
    "PGSC": {
        "traits": {
            "TC":  "PGS000677,PGS003138,PGS002783,PGS003819,PGS002424,PGP000561",
            "HDL": "PGS000686,PGS002781,PGS004156,PGS004043,PGS002450",
            "LDL": "PGS000688,PGS003978,PGS004981,PGS004974,PGS002409",
            "TG":  "PGS000699,PGS002784,PGS004342,PGS003802,PGS003148"
        },
        "out_suffix": "PGSC_Catalogue"
    },
    "OmicsPred": {
        "traits": {
            "TC":  f"{OMICSPRED_DIR}/OPGS003452.txt",
            "HDL": f"{OMICSPRED_DIR}/OPGS003520.txt",
            "LDL": f"{OMICSPRED_DIR}/OPGS003423.txt",
            "TG":  f"{OMICSPRED_DIR}/OPGS003486.txt"
        },
        "out_suffix": "Omics_Pred"
    }
}


# Extract metabolite names (assumes all sources have the same metabolites)
METABOLITES = list(PRS_TYPES["PGSC"]["traits"].keys())

# List of source keys (for expansion)
PRS_SOURCES = list(PRS_TYPES.keys())

# define chromosome patterns to parallelise process

CHROMOSOMES = list(range(1, 23+1))

def get_padded(wildcards):
  """Return padded number for input filename (01, 02, or 23)"""
  chrom = int(wildcards.chr)
  if chrom < 10:
      return f"0{chrom}"
  else:
      return str(chrom)

rule all:
  "The default rule"
  input:
        f"{RESULTS_DIR}/main_plot.png",
        f"{RESULTS_DIR}/prs_longtitudinal_analysis.csv",
        f"{RESULTS_DIR}/sensitivity_plot.png"

###################################################################################################
# (1) Recoding NAs
###################################################################################################

# 1.1 Recode NAs and create a sample sheet for eligible participants

rule pre_proccessing:
  "recode NAs, unify units, and select eligible individuals"
  input:
        pheno = f"{RAW_DIR}/B4891_Goudswaard_5Jan26.sav",
        omicsid = f"{RAW_DIR}/B4891_OmicsIDs_5Jan26.sav",
        sample = f"{RAW_DIR}/genetic_data_raw/swapped.sample"
  output:
        pheno = f"{INTERMEDIATE_DIR}/1-Pre-Proccessed/pheno_nas_recoded.csv",
        sampleid = f"{INTERMEDIATE_DIR}/1-Pre-Proccessed/eligible_sample.id",
        summary = f"{RESULTS_DIR}/summary_of_data.csv"
  shell:"""
  mkdir -p {INTERMEDIATE_DIR}/1-Pre-Proccessed
  Rscript 1-Pre-Proccessing/1-Pre-Proccessing.R {input.pheno} {input.omicsid} {input.sample} {output.pheno} {output.sampleid} {output.summary}
  """

###################################################################################################
# (2) Genetic QC + PCA
###################################################################################################

# 2.1 Variant-level-QC

rule variant_qc:
  "run variant-level QC"
  input:
        bgen = lambda w: f"{RAW_DIR}/genetic_data_raw/filtered_{get_padded(w)}.bgen",
        sample = f"{RAW_DIR}/genetic_data_raw/swapped.sample",
        eligible_sample = f"{INTERMEDIATE_DIR}/1-Pre-Proccessed/eligible_sample.id"
  output:
        pgen = f"{QC_DIR}/2.1-VL/VL_{{chr}}.pgen",
        psam = f"{QC_DIR}/2.1-VL/VL_{{chr}}.psam",
        pvar = f"{QC_DIR}/2.1-VL/VL_{{chr}}.pvar",
        log = f"{QC_DIR}/2.1-VL/VL_{{chr}}.log"
  params:
        out_prefix = f"{QC_DIR}/2.1-VL/VL_{{chr}}",
        chr_num = "{chr}"
  shell:"""
  mkdir -p {QC_DIR}/2.1-VL
  . 2-Genetic-QC/2.1-Variant-Level-QC.sh {input.bgen} {input.sample} {input.eligible_sample} {params.chr_num} {params.out_prefix}
  """

# 2.2 Sample-level QC

rule generate_slqc_list:
  "generate list of samples that passed sample-level qc and sex-check list"
  input:
        pgen = expand(f"{QC_DIR}/2.1-VL/VL_{{chr}}.pgen", chr=CHROMOSOMES),
        psam = expand(f"{QC_DIR}/2.1-VL/VL_{{chr}}.psam", chr=CHROMOSOMES),
        pvar = expand(f"{QC_DIR}/2.1-VL/VL_{{chr}}.pvar", chr=CHROMOSOMES)
  output:
        slqc_list = f"{QC_DIR}/2.2-SL/SL_QC.id",
        sm_list = f"{QC_DIR}/2.2-SL/SL_QC.sexcheck",
        king_list = f"{QC_DIR}/2.2-SL/SL_QC.king.cutoff.out.id"
  resources: mem_mb = 60000, threads = 20, runtime = 60
  params:
        in_dir = f"{QC_DIR}/2.1-VL",
        out_prefix = f"{QC_DIR}/2.2-SL/SL_QC"
  shell:"""
  mkdir -p {QC_DIR}/2.2-SL
  . 2-Genetic-QC/2.2.1-Generate-SLQC-List.sh {params.in_dir} {params.out_prefix}
  """

rule keep_slfiltered_samples:
  "keep samples that passed the filters"
  input:
        pgen = f"{QC_DIR}/2.1-VL/VL_{{chr}}.pgen",
        psam = f"{QC_DIR}/2.1-VL/VL_{{chr}}.psam",
        pvar = f"{QC_DIR}/2.1-VL/VL_{{chr}}.pvar",
        slqc_list = f"{QC_DIR}/2.2-SL/SL_QC.id",
        sm_list = f"{QC_DIR}/2.2-SL/SL_QC.sexcheck",
        king_list = f"{QC_DIR}/2.2-SL/SL_QC.king.cutoff.out.id"
  output:
        pgen = f"{QC_DIR}/2.2-SL/SL_{{chr}}.pgen",
        psam = f"{QC_DIR}/2.2-SL/SL_{{chr}}.psam",
        pvar = f"{QC_DIR}/2.2-SL/SL_{{chr}}.pvar",
        log = f"{QC_DIR}/2.2-SL/SL_{{chr}}.log"
  params:
        in_prefix = f"{QC_DIR}/2.1-VL/VL_{{chr}}",
        out_prefix = f"{QC_DIR}/2.2-SL/SL_{{chr}}"
  shell:"""
  . 2-Genetic-QC/2.2.2-Apply-SL-List.sh {params.in_prefix} {input.slqc_list} {input.king_list} {input.sm_list} {params.out_prefix}
  """

# 2.3 Remove problematic variants

rule problematic_variants_qc:
  "remove problematic variants"
  input:
        pgen = f"{QC_DIR}/2.2-SL/SL_{{chr}}.pgen",
        psam = f"{QC_DIR}/2.2-SL/SL_{{chr}}.psam",
        pvar = f"{QC_DIR}/2.2-SL/SL_{{chr}}.pvar"
  output:
        pgen = f"{QC_DIR}/2.3-PV/PV_{{chr}}.pgen",
        psam = f"{QC_DIR}/2.3-PV/PV_{{chr}}.psam",
        pvar = f"{QC_DIR}/2.3-PV/PV_{{chr}}.pvar",
        log = f"{QC_DIR}/2.3-PV/PV_{{chr}}.log"
  params:
        in_prefix = f"{QC_DIR}/2.2-SL/SL_{{chr}}",
        out_prefix = f"{QC_DIR}/2.3-PV/PV_{{chr}}"
  shell:"""
  mkdir -p {QC_DIR}/2.3-PV
  . 2-Genetic-QC/2.3-Problematic-Variants-QC.sh {params.in_prefix} {params.out_prefix}
  """

# 2.4 Pruning

rule pruning:
  "combine datasets and prune"
  input:
        pgen = expand(f"{QC_DIR}/2.3-PV/PV_{{chr}}.pgen", chr=CHROMOSOMES),
        psam = expand(f"{QC_DIR}/2.3-PV/PV_{{chr}}.psam", chr=CHROMOSOMES),
        pvar = expand(f"{QC_DIR}/2.3-PV/PV_{{chr}}.pvar", chr=CHROMOSOMES)
  output:
        pgen_unpruned = f"{QC_DIR}/2.4-Pruning/clean_unpruned.pgen",
        psam_unpruned = f"{QC_DIR}/2.4-Pruning/clean_unpruned.psam",
        pvar_unpruned = f"{QC_DIR}/2.4-Pruning/clean_unpruned.pvar",
        log_unpruned = f"{QC_DIR}/2.4-Pruning/clean_unpruned.log",
        pgen_pruned = f"{QC_DIR}/2.4-Pruning/clean_pruned.pgen",
        psam_pruned = f"{QC_DIR}/2.4-Pruning/clean_pruned.psam",
        pvar_pruned = f"{QC_DIR}/2.4-Pruning/clean_pruned.pvar",
        log_pruned = f"{QC_DIR}/2.4-Pruning/clean_pruned.log"
  resources: threads = 20, runtime = 25
  params:
        in_dir = f"{QC_DIR}/2.3-PV",
        unpruned_prefix = f"{QC_DIR}/2.4-Pruning/clean_unpruned",
        pruned_prefix = f"{QC_DIR}/2.4-Pruning/clean_pruned"
  shell:"""
  mkdir -p {QC_DIR}/2.4-Pruning
  . 2-Genetic-QC/2.4-Pruning.sh {params.in_dir} {params.unpruned_prefix} {params.pruned_prefix}
  """

# 2.5 PCA
rule pca:
  "PCA analysis of pruned variants"
  input:
        pgen_pruned = f"{QC_DIR}/2.4-Pruning/clean_pruned.pgen",
        psam_pruned = f"{QC_DIR}/2.4-Pruning/clean_pruned.psam",
        pvar_pruned = f"{QC_DIR}/2.4-Pruning/clean_pruned.pvar"
  output:
        eigenval = f"{QC_DIR}/2.5-PCA/genetic_pcs20.eigenval",
        eigenvec = f"{QC_DIR}/2.5-PCA/genetic_pcs20.eigenvec"
  resources: threads = 20, runtime = 90
  params:
        pruned_prefix = f"{QC_DIR}/2.4-Pruning/clean_pruned",
        pca_prefix = f"{QC_DIR}/2.5-PCA/genetic_pcs20"
  shell:"""
  mkdir -p {QC_DIR}/2.4-Pruning
  module add apps/plink2
  plink2 \
  --pfile {params.pruned_prefix} \
  --pca 20 \
  --out {params.pca_prefix}
  """

###################################################################################################
# (3) Metabolite Transformation
###################################################################################################

rule metab_transform:
  "Transform metabolites"
  input:
        pheno = f"{INTERMEDIATE_DIR}/1-Pre-Proccessed/pheno_nas_recoded.csv",
        eigenvec = f"{QC_DIR}/2.5-PCA/genetic_pcs20.eigenvec",
        eigenval = f"{QC_DIR}/2.5-PCA/genetic_pcs20.eigenval"
  output:
        pheno = f"{INTERMEDIATE_DIR}/3-Pre-Proccessed/pheno_transformed.csv",
        hist = f"{RESULTS_DIR}/transformed_metabolites.png"
  shell:"""
  mkdir -p {INTERMEDIATE_DIR}/3-Pre-Proccessed
  Rscript 3-Metabolite-Transformation /3-Metabolite-Transformation .R {input.pheno} {input.eigenvec} {input.eigenval} {output.pheno} {output.hist}
  """

###################################################################################################
# (4) PRS Calculation
###################################################################################################


def get_all_prs_inputs():
    """
    Return a flat list of all aggregated_scores.txt files.
    Order: for each source (PGSC, OmicsPred), for each metabolite (TC, HDL, LDL, TG).
    """
    files = []
    for source in PRS_SOURCES:
        suffix = PRS_TYPES[source]["out_suffix"]
        for m in METABOLITES:
            files.append(f"{INTERMEDIATE_DIR}/4-PRS/{m}/{suffix}/score/aggregated_scores.txt")
    return files

rule prs_calc_all:
    "Submit sbatch jobs for every (source, metabolite) combination"
    input:
        pgen = f"{QC_DIR}/2.4-Pruning/clean_unpruned.pgen",
        psam = f"{QC_DIR}/2.4-Pruning/clean_unpruned.psam",
        pvar = f"{QC_DIR}/2.4-Pruning/clean_unpruned.pvar"
    output:
        # This is a list of all aggregated score files (no wildcards)
        files = get_all_prs_inputs()
    resources:
        threads = 1,
        runtime = 300,
        mem_mb = 10
    params:
        # Build parallel arrays in the same order as `files`
        sources = " ".join([s for s in PRS_SOURCES for m in METABOLITES]),
        metabolites = " ".join([m for s in PRS_SOURCES for m in METABOLITES]),
        trait_lists = " ".join([
            PRS_TYPES[s]["traits"][m] for s in PRS_SOURCES for m in METABOLITES
        ]),
        out_suffixes = " ".join([
            PRS_TYPES[s]["out_suffix"] for s in PRS_SOURCES for m in METABOLITES
        ])
    shell:
        """
        # Convert params to bash arrays (order matches the input files)
        sources=({params.sources})
        metabolites=({params.metabolites})
        trait_lists=({params.trait_lists})
        out_suffixes=({params.out_suffixes})

        JID_FILE=$(mktemp)

        # Submit one sbatch per combination (source+metabolite)
        for i in ${{!sources[@]}}; do
            src="${{sources[$i]}}"
            meta="${{metabolites[$i]}}"
            traits="${{trait_lists[$i]}}"
            suffix="${{out_suffixes[$i]}}"
            out_dir="{QC_DIR}/4-PRS/$meta/$suffix/score"
            mkdir -p "$out_dir"
            sbatch --parsable \
                   --job-name="PRS-$meta-$src" \
                   src/run.sh \
                   --trait "$traits" \
                   --dir_out "$out_dir" >> "$JID_FILE" &
        done
        wait   # wait for all sbatch submissions to write to JID_FILE

        # Wait for each submitted job to finish (ignore failures)
        while read -r jid; do
            scontrol wait "$jid" || true
        done < "$JID_FILE"

        rm "$JID_FILE"

        # Touch all output files to satisfy Snakemake (even if some jobs failed)
        for i in ${{!sources[@]}}; do
            meta="${{metabolites[$i]}}"
            suffix="${{out_suffixes[$i]}}"
            touch "{QC_DIR}/4-PRS/$meta/$suffix/score/aggregated_scores.txt"
        done
        """

rule combine_all_prs:
    "Combine all PRS scores in one file"
    input:
        prs_files = get_all_prs_inputs()
    output:
        combined = f"{INTERMEDIATE_DIR}/4-PRS/all_aggregated_scores_combined.txt"
    params:
        metabolites = " ".join([m for s in PRS_SOURCES for m in METABOLITES]),
        prs_sources = " ".join([s for s in PRS_SOURCES for m in METABOLITES])
    shell:
        """
        metabolites=({params.metabolites})
        prs_sources=({params.prs_sources})
        input_files=({input.prs_files})

        out_file="{output.combined}"
        mkdir -p "$(dirname "$out_file")"

        # Find first non‑empty input file for the header
        first_file=""
        for f in "${{input_files[@]}}"; do
            if [ -s "$f" ]; then
                first_file="$f"
                break
            fi
        done

        if [ -z "$first_file" ]; then
            echo "ERROR: No non‑empty input files found." >&2
            exit 1
        fi

        header=$(head -n 1 "$first_file")
        echo -e "metabolite\\tprs_source\\t$header" > "$out_file"

        # Append data from each file with the two extra columns
        for i in "${{!input_files[@]}}"; do
            f="${{input_files[$i]}}"
            meta="${{metabolites[$i]}}"
            src="${{prs_sources[$i]}}"
            if [ -s "$f" ]; then
                tail -n +2 "$f" | awk -v meta="$meta" -v src="$src" '{{print meta "\\t" src "\\t" $0}}' >> "$out_file"
            else
                echo "Warning: $f is empty, skipping." >&2
            fi
        done
        """
###################################################################################################
# (5) Combine Datasets
###################################################################################################

rule combine_data:
  "combine pheno data and prs"
  input:
        pheno = f"{INTERMEDIATE_DIR}/3-Pre-Proccessed/pheno_transformed.csv",
        prs = f"{INTERMEDIATE_DIR}/4-PRS/all_aggregated_scores_combined.txt"
  output:
        agegroup = f"{INTERMEDIATE_DIR}/5-Combined-Data/combined_data_agegroup.csv",
        long = f"{INTERMEDIATE_DIR}/5-Combined-Data/combined_data_long.csv",
  shell:"""
  Rscript 5-Combine-Data/5-Combine-Data.R {input.pheno} {input.prs} {output.agegroup} {output.long}
  """

###################################################################################################
# (6) PRS Analysis
###################################################################################################

rule analysis_byagegroup:
  "analyse prs by age group"
  input: f"{INTERMEDIATE_DIR}/5-Combined-Data/combined_data_agegroup.csv"
  output: f"{RESULTS_DIR}/prs_performance_by_age_group.csv"
  resources: threads = 24, time = 120
  shell:"""
  Rscript 6-PRS-Analysis/6.1-PRS-Analysis-By-Age-Group.R {input} {output}
  """
rule sensitivity_analysis:
  "analyse longtitudinal subset data by age group"
  "analyse prs by age group"
  input: f"{INTERMEDIATE_DIR}/5-Combined-Data/combined_data_long.csv"
  output: f"{RESULTS_DIR}/prs_performance_by_age_group_sensitivity.csv"
  resources: threads = 24, time = 120
  shell:"""
  Rscript 6-PRS-Analysis/6.1-PRS-Analysis-By-Age-Group.R {input} {output}
  """

rule longtitudinal_analysis:
  "analyse each PRS longtitudinally"
  input: f"{INTERMEDIATE_DIR}/5-Combined-Data/combined_data_long.csv"
  output: f"{RESULTS_DIR}/prs_longtitudinal_analysis.csv"
  resources: threads = 12
  shell:"""
  Rscript 6-PRS-Analysis/6.2-PRS-Longtitudinal-Analysis.R {input} {output}
  """

###################################################################################################
# (7) Plots
###################################################################################################

rule plot_main:
  "plot results of prs by age group"
  input: f"{RESULTS_DIR}/prs_performance_by_age_group.csv"
  output: f"{RESULTS_DIR}/main_plot.png"
  shell:"""
  Rscript 7-Plots/7-PRS-By-Metabolite-By-Age-Group.R {input} {output}
  """
rule plot_sensitivity:
  "plot results of prs by age group"
  input: f"{RESULTS_DIR}/prs_performance_by_age_group_sensitivity.csv"
  output: f"{RESULTS_DIR}/sensitivity_plot.png"
  shell:"""
  Rscript 7-Plots/7-PRS-By-Metabolite-By-Age-Group.R {input} {output}
  """

# 0. Cleaning

rule clean1:
  "clean directory 1-Pre-Proccessed"
  shell:"""
  rm -rf {INTERMEDIATE_DIR}/1-Pre-Proccessed
  """
rule clean_VL:
  "clean directory 2.1-VL"
  shell:"""
  rm -rf {QC_DIR}/2.1-VL
  """
rule clean_SL:
  "clean directory 2.2-SL"
  shell:"""
  rm -rf {QC_DIR}/2.2-SL
  """
rule clean_PV:
  "clean directory 2.3-PV"
  shell:"""
  rm -rf {QC_DIR}/2.3-PV
  """
rule clean_pruning:
  "clean directory 2.4-Pruning"
  shell:"""
  rm -rf {QC_DIR}/2.4-Pruning
  """
rule clean_pca:
  "clean directory 2.5-PCA"
  shell:"""
  rm -rf {QC_DIR}/2.5-PCA
  """
rule clean_transformed:
  "clean directory 3-Pre-Proccessed"
  shell:"""
  rm -rf {INTERMEDIATE_DIR}/3-Pre-Proccessed
  """
rule clean_combined:
  "clean directory 5-Combined-Data"
  shell:"""
  rm -rf {INTERMEDIATE_DIR}/5-Combined-Data
  """
rule clean_results:
  "clean results directory"
  shell:"""
  rm -rf {RESULTS_DIR}
  """

rule clean_genetic:
  "clean all intermediate genetic QC directories (keep final results)"
  shell:"""
  snakemake -c1 clean1 clean_VL clean_SL clean_PV 
  """