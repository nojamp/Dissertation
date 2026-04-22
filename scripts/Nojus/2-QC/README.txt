Pipeline to run quality control (QC) on genetic data
To run type "snakemake --profile ." on a HPC with the correct credentials

## Bash files
2.1 - Imputs BGEN files, carries out sample-level QC, returns plink native files (pfiles)
2.2.1 - Merges pfiles, creates a list of samples that did not pass sample-level filtering,
and sex-mismatches. The merge is deleted to save space and individual files are used for parallelisation
2.2.2 - Uses the lists generated in 2.2.1 to remove samples from individual pfiles
2.3 - removes non-biallelic variants, and ambiguous variants
2.4 - merges pfiles, carries out ld-pruning. Returns both pfiles as both will be used downstream

## Snakefile
Contains a rule to run each of these bash files
Contains additional commands to clean data directories
Contains specific resource requests for certain rules (adjust as needed) 
