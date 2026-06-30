# Load libraries ----------------------------------------------------------

library(tidyverse)

# Load data ---------------------------------------------------------------

args = commandArgs(trailingOnly=TRUE)
INPUT_MET_TRANS_PLUS_PCS = args[1]
INPUT_METABOLITE_NAME = args[2]
INPUT_PRS_SCORES = args[3]

OUTPUT = args[4]

pheno_df <- read_csv(paste(INPUT_MET_TRANS_PLUS_PCS))
metabolite <- as.character(INPUT_METABOLITE_NAME)
prs_df <- read.table(paste(INPUT_PRS_SCORES), header = TRUE, sep = "\t")

# Standardise PRS scores and age ------------------------------------------

# normalise the scores
prs_std_df <- prs_df %>%
  group_by(PGS) %>%
  mutate(SUM_norm = (SUM - mean(SUM)) / sd(SUM)) %>%
  ungroup()

# transform age (standardise) for easier convergence
pheno_df_age_rescaled <- pheno_df %>% mutate(age_std = (age - mean(age))/sd(age))


# Combine and save --------------------------------------------------------

# put prs scores into compatible format ready to combine
prs_std_wide <- prs_std_df %>%
  pivot_wider(
    id_cols     = ID_sample,      # Column 1 = ID_sample
    names_from  = PGS,            # New column names come from unique values in PGS
    values_from = SUM_norm        # Fill those new columns with SUM_norm
  ) %>% mutate(
    gi_1000g_g0m_g1 = str_sub(ID_sample, 1, -2)
  )

# combine the data using the omics Ids
combined_data <- prs_std_wide %>% inner_join(pheno_df_age_rescaled, by = "gi_1000g_g0m_g1")

# save
combined_data %>% write_csv(paste(OUTPUT))
