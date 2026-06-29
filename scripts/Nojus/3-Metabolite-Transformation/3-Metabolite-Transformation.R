library(tidyverse)
library(RNOmni)

# log transform
# covariate adjustment
# inverse rank normalisation
# diagnostic plots


# Load data ---------------------------------------------------------------

args = commandArgs(trailingOnly=TRUE)
INPUT_PHENO_DATA = args[1]
INPUT_OMICS_IDS = args[2]
INPUT_PCS = args[3]
OUTPUT_DIR = args[4]

data <- read_csv(paste(INPUT_PHENO_DATA))
omics_ids <- read.csv(paste(INPUT_OMICS_IDS))

pcs <- read.table(paste(INPUT_PCS), header = FALSE)
colnames(pcs) <- c("FID", "IID", paste0("PC", 1:20))
pcs <- pcs %>%
  mutate(
    qlet = str_sub(FID, -1, -1),  # extract last character
    gi_1000g_g0m_g1 = str_sub(FID, 1, -2)                # remove last character from original
  )


ls(omics_ids)
head(pcs)
head(omics_ids$gi_1000g_g0m_g1)

length(pcs$FID)
length(omics_ids$gi_1000g_g0m_g1)

# dim(data)
dim(omics_ids)

# combine the data sets
omics_pcs <- inner_join(omics_ids, pcs, by = c("gi_1000g_g0m_g1", "qlet"))
data_post_qcs <- inner_join(data, omics_pcs, by = c("cidB4891", "qlet"))




dim(data_post_qcs)


# unselect all with _cord, 31, 43, LDLC, HDLC, SerumTG, ApoA1, ApoB

# Unselect (drop) columns matching any pattern
dat <- data_post_qcs %>%
  select(-matches("_cord|31|43|LDLC|HDLC|Serum|ApoA1|ApoB", ignore.case = TRUE))


ls(dat)


# Log transform then obtain residuals of metabolites ----------------------

# function to do this for each time point
pre_IRN_processing <- function(dat, study_var_codes, age_col, bmi_col, month_col) {
  
  dat %>%
    select(matches(paste0(study_var_codes, "|kz021|PC|cidb4891|gi_1000g_g0m_g1"), ignore.case = TRUE)) %>%     # select the study vars, sex, PCs, ids
    mutate(across(
      .cols = matches("CHOL|TRIG|LDL|HDL", ignore.case = TRUE),
      .fns = ~ residuals(
        lm(log(.x) ~ .data[[age_col]] + .data[[bmi_col]] + kz021 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,    # residuals of log transformed metabolites regressed on age, bmi, sex, PCs
           data = cur_data(),
           na.action = na.exclude)
      )
    )) %>%
    rename_with(
      .cols = matches("CHOL|TRIG|LDL|HDL", ignore.case = TRUE),
      .fn = ~ toupper(str_extract(., regex("CHOL|TRIG|LDL|HDL", ignore_case = TRUE)))    # rename the columns without reference to study for combining later
    ) %>%
    rename(age = !!sym(age_col), BMI = !!sym(bmi_col)) %>%    # code the season variable
    mutate(season = case_when(
      !!sym(month_col) %in% 3:5   ~ "Spring",
      !!sym(month_col) %in% 6:8   ~ "Summer",
      !!sym(month_col) %in% 9:11  ~ "Autumn",
      !!sym(month_col) %in% c(12,1,2) ~ "Winter"
    )) %>% mutate(age = if_else(age > 40, age / 12, age)) # make sure age is in years (all participants should be below 30 years of age and 7 year olds are about 84 months old so this is a safe threshold)
}

# use the function and add a column describing the group the data comes from
F7_res <- pre_IRN_processing(dat, "f7", "f7003c", "f7ms026a", "f7001") %>% mutate(study_group = "F7")
F9_res <- pre_IRN_processing(dat, "f9", "f9003c", "f9ms026a", "f9001") %>% mutate(study_group = "F9")
TF3_res <- pre_IRN_processing(dat, "TF3|fh", "fh7423", "fh3019", "fh0010a") %>% mutate(study_group = "TF3")
TF4_res <- pre_IRN_processing(dat, "TF4|FJ", "FJ003b", "FJMR022a", "FJ002a") %>% mutate(study_group = "TF4")
F24_res <- pre_IRN_processing(dat, "F24|FK", "FKAR0011", "FKMS1040", "FKAR0040") %>% mutate(study_group = "F24")


# Group by age ------------------------------------------------------------

children <- bind_rows(F7_res, F9_res) %>% mutate(age_group = "children")
teens <- bind_rows(TF3_res, TF4_res) %>% mutate(age_group = "teens")
adults <- F24_res %>% mutate(age_group = "adults")


# Rank Inverse Normalisation of metabolites -------------------------------

# function to perform RIN for each metabolite in a group
RIN_each_metab <- function(data) {
  # Define the four columns of interest
  cols <- c("CHOL", "TRIG", "LDL", "HDL")
  
  # Start with a copy of the original data
  result <- data
  
  for (col in cols) {
    # 1. Extract ID and the current column, remove NAs
    temp <- data %>%
      select(cidB4891, study_group, all_of(col)) %>%
      filter(!is.na(.data[[col]]))   # keep only non-NA rows for this column
    
    # 2. Apply Inverse Rank Norm transformation
    temp <- temp %>%
      mutate(transformed = RankNorm(.data[[col]])) %>%
      select(cidB4891, study_group, transformed)   # keep only ID, study group, and transformed values
    
    # 3. Left join the transformed values back to the main data
    #    Then replace the original column with the transformed one
    result <- result %>%
      left_join(temp, by = c("cidB4891", "study_group")) %>%
      mutate(!!sym(col) := coalesce(transformed, .data[[col]])) %>%
      select(-transformed)   # remove the temporary column
  }
  
  return(result)
}

# perform the RIN on the groups
children_transformed <- RIN_each_metab(children)
teens_transformed <- RIN_each_metab(teens)
adults_transformed <- RIN_each_metab(adults)

# check the transformations worked as intended
plot_transformed_met <- function(data, group_name){
  cols <- c("CHOL", "TRIG", "LDL", "HDL")
  for (col in cols) {
    data %>% pull(all_of(col)) %>% na.omit() %>% hist(main = paste("Histogram of transformed", col, group_name))
    print(data %>% select(all_of(col), cidB4891) %>% na.omit() %>% group_by(cidB4891) %>% sample_n(1) %>% ungroup() %>% dim())
  }
  print(data %>% select(all_of(cols)) %>% summary())
}

plot_transformed_met(children_transformed, "children")
plot_transformed_met(teens_transformed, "teens")
plot_transformed_met(adults_transformed, "adults")

# Combine the data then split by metabolite -------------------------------

combined_transformed <- bind_rows(children_transformed, teens_transformed, adults_transformed)

comb_chol <- combined_transformed %>% select(-TRIG, -LDL, -HDL) %>% filter(is.na(CHOL) == FALSE)
comb_trig <- combined_transformed %>% select(-CHOL, -LDL, -HDL) %>% filter(is.na(TRIG) == FALSE)
comb_ldl <- combined_transformed %>% select(-CHOL, -TRIG, -HDL) %>% filter(is.na(LDL) == FALSE)
comb_hdl <- combined_transformed %>% select(-CHOL, -TRIG, -LDL) %>% filter(is.na(HDL) == FALSE)


# Report data sizes for each metabolite and age group ---------------------

data_size <- function(data){
  print("total data")
  print(data %>% group_by(cidB4891) %>% sample_n(1) %>% ungroup() %>% dim())
  
  list_by_group <- list()
  
  for (i in c("children", "teens", "adults")) {
    list_with_met_data <- data %>% filter(age_group == i) %>% group_by(cidB4891) %>% sample_n(1) %>% ungroup() %>% pull(cidB4891)
    list_by_group[[i]] <- list_with_met_data
    print(i)
    print(length(list_with_met_data))
  }
  
  print("longtitudinal")
  print(length(intersect(list_by_group$children, list_by_group$teens) %>% intersect(list_by_group$adults)))

}

data_size(comb_chol)
data_size(comb_trig)
data_size(comb_ldl)
data_size(comb_hdl)


# Save the different groups -----------------------------------------------

comb_chol %>% write_csv(paste(OUTPUT_DIR,"/3-chol-data.csv"))
comb_trig %>% write_csv(paste(OUTPUT_DIR,"/3-trig-data.csv"))
comb_ldl %>% write_csv(paste(OUTPUT_DIR,"/3-ldl-data.csv"))
comb_hdl %>% write_csv(paste(OUTPUT_DIR,"/3-hdl-data.csv"))




