R script that log transforms, regresses, then performs inverse rank normalisation (IRN) on metabolites

It first combines the phenotype data, omics_ids, and PCs

Then for each study it log transforms metabolites, regresses them against age, sex, BMI, and top 10 PCs, and returns residuals 
It also then renames the study specific variables to simply age, and BMI
Also creates a new season variable
This is done for each study group as each one has differently coded variables

The different age groups are then combined and IRN is performed

Finally all the dataset splits are combined and then split by metabolite
There are extra variables study_group, and age_group denoting groupings
There are some duplicate IDs for some of age_group members as some had several screenings
