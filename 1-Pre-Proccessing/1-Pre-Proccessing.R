###########################################################################
######## Recoding NAs, unifying units, selecting eligible individuals #####
###########################################################################

# Summary -----------------------------------------------------------------

# 15678 individuals in original data set
# some NAs coded as NAs others as negative numbers
# 7260 missing data for all time points for all 4 metabolites of interest

# Load libraries ----------------------------------------------------------

library(tidyverse)
library(haven) # for reading sav files
library(naniar) # for dealing with NAs
library(data.table) # for fread

set.seed(13)

# Load data ---------------------------------------------------------------

args = commandArgs(trailingOnly=TRUE)
INPUT_PHENO = args[1]
INPUT_OMICSID = args[2]
INPUT_SWAPPED_SAMPLE = args[3]
OUTPUT_PHENO = args[4]
OUTPUT_SAMPLEID = args[5]
OUTPUT_SUMMARY = args[6]

# load pheno data
dat_pheno <- read_sav(paste0(INPUT_PHENO))
head(dat_pheno)

# load omics ids
omicsid <- read_sav(paste0(INPUT_OMICSID))
head(omicsid)

# load swapped.sample
swapped_sample <- fread(paste0(INPUT_SWAPPED_SAMPLE), header = TRUE)
head(swapped_sample)

# Inspect data ------------------------------------------------------------

dim(dat_pheno)
# 15678 individuals with 176 variables

ls(dat_pheno)
# vars beggining with chol, Chol, or CHOL refer to cholesterol
# vars beggining HDL or hdl (NOT HDLC) refer to HDL
# vars beggining LDL or ldl (NOT LDLC) refer to LDL
# vars beggining trig, Trig, or TRIG refer to triglycerides

# select variables of interest
# f7, f9, fh, FJ, FK refer to variables for f7,f9,tf3,tf4,f24 respectively

id_cols <- c("cidB4891", "qlet")

membership_cols <- c("in_f07", "in_f09", "in_tf3", "in_tf4", "in_F24")

metabolite_cols <- dat_pheno %>% 
  select(matches("chol_|hdl_|ldl_|trig_", ignore.case = TRUE)) %>%          # match metabolite names
  select(matches("f7|f9|tf3|tf4|f24", ignore.case = TRUE)) %>% colnames()   # only for f7, f9, tf3, tf4, f24

age_cols <- c("f7003c", "f9003c", "fh0011a", "FJ003b", "FKAR0011")

sex_col <- "kz021"

bmi_cols <- c("f7ms026a", "f9ms026a", "fh3019", "FJMR022a", "FKMS1040")

month_cols <- c("f7001", "f9001", "fh0010a", "FJ002a", "FKAR0040")

cols_of_interest <- c(id_cols, membership_cols, metabolite_cols, age_cols, sex_col, bmi_cols, month_cols)

dat_pheno <- dat_pheno %>% select(all_of(cols_of_interest))

summary(dat_pheno)
# many NAs
# some variables have negative values when that is not possible

# Recoding values <=0 as NAs ----------------------------------------------

# using naniar package
dat_pheno_NA_recoded <- dat_pheno %>% replace_with_na_all(condition = ~.x <= 0)

# check all variables of interest are positive only
summary(dat_pheno_NA_recoded)
# check


# Unifying units ----------------------------------------------------------

dat_pheno_var_recoded <- dat_pheno_NA_recoded %>% 
  # f7003c, f9003c, fh0011a age in months but FJ003b, FKAR0011 age in years
  # convert f7003c, f9003c, fh0011a to age in years
  mutate(
    f7003c = f7003c / 12,
    f9003c = f9003c / 12,
    fh0011a = fh0011a /12
  ) %>%
  # convert month of attendance to season
  mutate(across(
    all_of(month_cols),
    ~ case_when(
      .x %in% c(12, 1, 2) ~ "Winter",
      .x %in% c(3, 4, 5)  ~ "Spring",
      .x %in% c(6, 7, 8)  ~ "Summer",
      .x %in% c(9, 10, 11) ~ "Autumn",
      TRUE ~ NA_character_
    ) %>%
      factor(levels = c("Spring", "Summer", "Autumn", "Winter"), ordered = TRUE)
  ))

# verify all variables behave as anticipated for all follow up groups
summary(dat_pheno_var_recoded)
# check


# Change data structure and define groups ---------------------------------

# create a lookup table
lookup_table <- tibble(
  focus = c("f7", "f9", "tf3", "tf4", "f24"),
  membership = membership_cols,
  age = age_cols,
  sex = c("kz021", "kz021", "kz021", "kz021", "kz021"),
  bmi = bmi_cols,
  season = month_cols,
  tc = grep("chol_", metabolite_cols, value = TRUE, ignore.case = TRUE),
  hdl = grep("hdl_", metabolite_cols, value = TRUE, ignore.case = TRUE),
  ldl = grep("ldl_", metabolite_cols, value = TRUE, ignore.case = TRUE),
  tg = grep("trig_", metabolite_cols, value = TRUE, ignore.case = TRUE)
)

# check it lines up
lookup_table

# create a row for each follow up
dat_pheno_long <- lookup_table %>%
  pmap_dfr(function(focus, membership, age, sex, bmi, season, tc, hdl, ldl, tg) {
    dat_pheno_var_recoded %>%
      select(cidB4891,
             qlet,
             all_of(membership), 
             all_of(age), 
             all_of(sex), 
             all_of(bmi), 
             all_of(season), 
             all_of(tc), 
             all_of(hdl), 
             all_of(ldl), 
             all_of(tg)) %>%
      # Keep only rows where the person attended the follow‑up
      filter(!is.na(.data[[membership]])) %>%
      # Add the time label
      mutate(focus = focus) %>%
      # Rename the variable columns to clean standard names
      rename(age    = all_of(age),
             sex    = all_of(sex),
             bmi    = all_of(bmi),
             season = all_of(season),
             tc     = all_of(tc),
             hdl    = all_of(hdl),
             ldl    = all_of(ldl),
             tg     = all_of(tg)) %>%
      # Drop the membership column (no longer needed)
      select(-all_of(membership))
  })

# check this has been done correctly
glimpse(dat_pheno_long)

# add group names
dat_pheno_long <- dat_pheno_long %>%
  mutate(
    age_group = case_when(
      focus %in% c("f7", "f9")    ~ "children",
      focus %in% c("tf3", "tf4")  ~ "teens",
      focus == "f24"              ~ "adults",
      TRUE                        ~ NA_character_
    ) %>%
      factor(levels = c("children", "teens", "adults"), ordered = TRUE)
  )


# Select for individuals with data in phenotype vars of interest ----------

# new metabolite col names
metabolite_cols <- c("tc", "hdl", "ldl", "tg")

# make sure each row has age, sex, bmi, season, and at least one of the metabolites
dat_pheno_eligible <- dat_pheno_long %>%
  drop_na(all_of(c("age", "sex", "bmi", "season"))) %>% 
  filter(if_any(all_of(metabolite_cols), ~ !is.na(.)))

dat_pheno_eligible %>% 
  select(cidB4891) %>% unique() %>% dim()
# 8257 with complete phenotype data 

# for longtitudinal analysis must have data on all 3 age groups
# create an indicator variable for this
dat_pheno_eligible <- dat_pheno_eligible %>% group_by(cidB4891) %>% mutate(
  all_three_times = as.integer(all(c("children","teens","adults") %in% age_group))
) %>% ungroup()

dat_pheno_eligible %>% 
  filter(all_three_times > 0 ) %>% 
  select(cidB4891) %>% unique() %>% dim()
# 1904 with eligible phenotype longtitudinal data


# Select for individuals with genetic data --------------------------------

# add column in omicsid to match swapped.sample
omicsid <- omicsid %>% mutate(ID_1 = paste0(gi_1000g_g0m_g1, qlet))

# combine the two so we can match to phenotype data
genetic_data <- omicsid %>% inner_join(swapped_sample, by = "ID_1") %>% select(cidB4891, qlet, gi_1000g_g0m_g1, ID_1)

# combine with pheno data
dat_eligible <- dat_pheno_eligible %>% inner_join(genetic_data, by = c("cidB4891", "qlet")) %>%
  rename(IID = ID_1) # for compatibility when joining pc data / prs later


dat_eligible %>% select(cidB4891) %>% unique() %>% dim()
# 6791 individuals with complete pheno data and genetic data

dat_eligible %>%
  group_by(age_group) %>%
  summarise(
    n_unique_cids = n_distinct(cidB4891),
    .groups = "drop"
  )
# 6132 in children
# 3770 in teens
# 2540 in adults

dat_eligible %>% 
  group_by(age_group, IID) %>%
  count(IID) %>% 
  filter(n>1) %>%
  ungroup() %>%
  group_by(age_group) %>%
  summarise(
    n_unique_cids = n_distinct(IID),
    .groups = "drop"
  )
# 3002 children with repeated measures
# 1872 teens with repeated measures
# 0 adults with repeated measures

dat_eligible %>% filter(all_three_times > 0) %>% select(cidB4891) %>% unique() %>% dim()
# 1718 individuals with pheno and genetic data in all 3 groups

# Justification for use of groups -----------------------------------------

dat_eligible %>%
  group_by(focus) %>%
  summarise(
    min = min(age, na.rm = TRUE),
    median = median(age, na.rm = TRUE),
    max = max(age, na.rm = TRUE),
    .groups = "drop"
  )
# groups f7 and f9 overlap in age
# groups tf3 and tf4 overlap in age

dat_eligible %>% group_by(cidB4891) %>% 
  filter(as.integer(all(c("f7", "f9","tf3", "tf4", "f24") %in% focus))>0) %>%
   ungroup() %>%
  select(cidB4891) %>% unique() %>% dim()
# only 682 individuals have eligible data for all 5 focus visits


# Summary of eligible data ------------------------------------------------

# function to create summary of numerical variables
numeric_summary <- function(data){
  res <- data %>%
    group_by(age_group) %>%
    summarise(
      across(
        c(age, bmi, tc, hdl, ldl, tg),
        list(
          mean = ~ mean(.x, na.rm = TRUE),
          sd = ~ sd(.x, na.rm = TRUE),
          n = ~ sum(!is.na(.x))
        ),
        .names = "{.col}_{.fn}"
      ),
      .groups = "drop"
    ) %>%
    pivot_longer(
      cols = -age_group,
      names_to = c("variable", "stat"),
      names_sep = "_",
      values_to = "value"
    ) %>%
    pivot_wider(
      id_cols = variable,
      names_from = c(age_group, stat), 
      values_from = value,
      names_sep = "_"
    ) %>%
    select(variable, starts_with("children"), starts_with("teens"), starts_with("adults"))
  
  return(res)
}

# function to create summary of sex
sex_summary <- function(data){
  res <- data %>%
    group_by(age_group, sex) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(age_group) %>%
    mutate(percent = n / sum(n) * 100) %>%
    pivot_wider(
      id_cols = sex,
      names_from = age_group,
      values_from = c(n, percent),
      names_glue = "{age_group}_{.value}"
    ) %>% 
    mutate(
      variable = paste0("sex:", sex)
    ) %>% select(-sex)
  
  return(res)
}

# function to create summary of season variable
season_summary <- function(data){
  res <- data %>%
    group_by(age_group, season) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(age_group) %>%
    mutate(percent = n / sum(n) * 100) %>%
    pivot_wider(
      id_cols = season,
      names_from = age_group,
      values_from = c(n, percent),
      names_glue = "{age_group}_{.value}"  # This puts age_group first, then the value name
    ) %>% 
    mutate(
      variable = paste0("season:", season)
    ) %>% select(-season)
  
  return(res)
}

# carry out a complete summary of data
complete_summary <- function(data){
  res <- bind_rows(
    numeric_summary(data),
    sex_summary(data),
    season_summary(data)
    )
  
  return(res)
}

# summarries of both datasets
summary_of_biggerset <- complete_summary(dat_eligible) %>%
                              mutate(dataset = "A")
summary_of_smallerset <- complete_summary(dat_eligible %>% filter(all_three_times>0)) %>%
                              mutate(dataset = "B")

# combine datasets
summary_of_all <- bind_rows(summary_of_biggerset, summary_of_smallerset)


# Save data ---------------------------------------------------------------

# create a sample file for use in the genetic qc pipeline
sample_id <- data.frame("#FID" = dat_eligible$IID, "IID" = dat_eligible$IID, check.names = FALSE)

# save data
dat_eligible %>% write_csv(paste0(OUTPUT_PHENO))
sample_id %>% write_tsv(paste0(OUTPUT_SAMPLEID), quote = "none")
summary_of_all %>% write_csv(paste0(OUTPUT_SUMMARY))
