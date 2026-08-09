###########################################################################
######## Recoding NAs in Pheno Data #######################################
###########################################################################

# Summary -----------------------------------------------------------------

# 15678 individuals in original data set
# some NAs coded as NAs others as negative numbers
# 7260 missing data for all time points for all 4 metabolites of interest

# Load libraries ----------------------------------------------------------

library(tidyverse)
library(haven) # for reading sav files
library(naniar) # for dealing with NAs

set.seed(13)

# Load data ---------------------------------------------------------------

args = commandArgs(trailingOnly=TRUE)
INPUT_PHENO = args[1]
INPUT_OMICSID = args[2]
OUTPUT_PHENO = args[3]
OUTPUT_SAMPLEID = args[4]

# load pheno data
dat_pheno <- read_sav(paste0(INPUT_PHENO))
head(dat_pheno)

# Load in omicsid
omicsid <- read_sav(paste0(INPUT_OMICSID))
head(omicsid)

# Inspect data ------------------------------------------------------------

dim(dat_pheno)
# 15678 individuals with 176 variables

ls(dat_pheno)
# vars beggining with chol, Chol, or CHOL refer to cholesterol
# vars beggining HDL or hdl (NOT HDLC) refer to HDL
# vars beggining LDL or ldl (NOT LDLC) refer to LDL
# vars beggining trig, Trig, or TRIG refer to triglycerides

# function to summarise each metabolite individually from a data set
summarise_metabolites <- function(data, metabolite_name){
  data %>%  reframe(across(matches(paste0(metabolite_name,"_"), ignore.case = TRUE), summary))  %>%
    mutate(statistic = c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.", "NAs")) %>%
    select(statistic, everything()) %>%
    mutate(across(-statistic, ~round(., 2)))
}

# summarise  cholesterol
summarise_metabolites(dat_pheno, "chol")

# summarise HDL
summarise_metabolites(dat_pheno, "hdl")

# summarise LDL
summarise_metabolites(dat_pheno, "ldl")

# summarise triglycerides
summarise_metabolites(dat_pheno, "trig")

# Min for each is below 0 so must be recoded to be NA

# Recoding values <= zero as NAs ------------------------------------------

# all values for all variables should be positive
# recode all values that are <= 0 as NA using naniar package
dat_pheno_NA_recoded <- dat_pheno %>% replace_with_na_all(condition = ~.x <= 0)

# check that this has worked
summarise_metabolites(dat_pheno_NA_recoded, "chol")
summarise_metabolites(dat_pheno_NA_recoded, "hdl")
summarise_metabolites(dat_pheno_NA_recoded, "ldl")
summarise_metabolites(dat_pheno_NA_recoded, "trig")

# different number of NAs for different clinic visits across each metabolite (as expected as people drop out)
# similar number of NAs for same clinic visit regardless of metabolite (suggests same people had no data for all 4)

# Cleaning data with NAs in all 4 metabolites -----------------------------

# Find columns matching the four metabolites
metabolite_cols <- grep("chol_|hdl_|ldl_|trig_", names(dat_pheno), ignore.case = TRUE, value = TRUE)

# exclude columns for metabolite measures at 31 and 43 and cord --> not interested in these time points
metabolite_cols_ofinterest <- metabolite_cols[!grepl("31|43|cord", metabolite_cols)]

dim(dat_pheno_NA_recoded)[1] - (dat_pheno_NA_recoded %>% filter(if_any(all_of(metabolite_cols_ofinterest), ~!is.na(.))) %>% dim())[1]
# 7260 people completely missing data on all 4 metabolites of interest for all time points of interest

dat_pheno_complete_NAs_removed <- dat_pheno_NA_recoded %>%
  filter(if_any(all_of(metabolite_cols_ofinterest), ~!is.na(.)))

# Combine with omicsids ---------------------------------------------------

combined_data <- dat_pheno_complete_NAs_removed %>% inner_join(omicsid, by = c("cidB4891", "qlet")) %>%
                mutate( IID = paste0(gi_1000g_g0m_g1, qlet)) # for compatibility with genetic samples

# create file compatible with PLINK2
sample_id <- data.frame("#FID" = combined_data$IID, "IID" = combined_data$IID, check.names = FALSE)

# Save data ---------------------------------------------------------------

# save data
combined_data %>% write_csv(paste0(OUTPUT_PHENO))
sample_id %>% write_tsv(paste0(OUTPUT_SAMPLEID), quote = "none")
