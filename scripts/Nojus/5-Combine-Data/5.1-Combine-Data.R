###########################################################################
######## Combining PRS and Pheno Data #####################################
###########################################################################

# Load libraries ----------------------------------------------------------

library(tidyverse)

# Load data ---------------------------------------------------------------

args = commandArgs(trailingOnly=TRUE)
INPUT_PHENO = args[1]
INPUT_PRS = args[2]
OUTPUT = args[3]

# load data
dat_pheno <- read_csv(paste0(INPUT_PHENO))
dat_prs <- read.table(paste0(INPUT_PRS), header = TRUE, sep = "\t")

# only keep PRS ID and add metabolite code to the end of it
dat_prs <- dat_prs %>% mutate(
  PGS = str_remove(PGS, "_.*"),
  PGS = paste0(PGS, "_", metabolite)
)

# identify all unique prs cols
prs_cols <- unique(dat_prs$PGS)

# put prs scoreshead()# put prs scores into compatible form ready to combine
dat_prs_wide <- dat_prs %>%
  pivot_wider(
    id_cols = ID_sample,
    names_from = PGS,
    values_from = SUM
  ) %>% rename(IID = ID_sample)


# Combine datasets --------------------------------------------------------

# combine pheno data and prs
combined_data <- dat_pheno %>% inner_join(dat_prs_wide, by = "IID")

# standardise PRS
combined_data <- combined_data %>% mutate(
  across(all_of(prs_cols), ~ (.x - mean(.x, na.rm = TRUE)) / sd(.x, na.rm = TRUE))
)

# Plot histograms with density to check this has happened correctly
p1 <- dat_prs_wide %>%
  pivot_longer(prs_cols, names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value)) +
  geom_histogram(aes(y = after_stat(density)), bins = 15, fill = "steelblue", color = "white") +
  facet_wrap(~ variable, scales = "free") +
  theme_minimal() +
  labs(title = "Histograms of standardized variables", y = "Density")


# Save dataset ------------------------------------------------------------
combined_data %>% write_csv(paste0(OUTPUT))
