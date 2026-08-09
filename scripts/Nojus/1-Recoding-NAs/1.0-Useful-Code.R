# children born to the same mother at the same time (e.g. twins) have the same ID but different qlet
# there are 200 twins
dat_pheno_noNAs %>% count(cidB4891) %>% filter(n>1) %>% nrow()

# randomly select one of the IDs out of the duplicates
dat_pheno_noTwins <- dat_pheno_noNAs %>% group_by(cidB4891) %>% sample_n(1) %>% ungroup()

# Check duplicates successfully removed
dat_pheno_noTwins %>% count(cidB4891) %>% filter(n>1) %>% nrow()



#####################
# Recoding values outside 10sd as NAs -------------------------------------

# Function to identify values > 10 SD from mean
is_10sd_outlier <- function(x) {
  abs(x - mean(x, na.rm = TRUE)) > 10 * sd(x, na.rm = TRUE)
}

# recode values 10sd away from mean as NA
dat_pheno_NA_recoded <- dat_pheno_below_zero_NAs %>%
  mutate(across(all_of(metabolite_cols), ~ifelse(is_10sd_outlier(.), NA, .)))

# check if this has introduced new NAs
sum(is.na(dat_pheno_NA_recoded)) - sum(is.na(dat_pheno_below_zero_NAs))
# 15 new NAs


###############

dim(dat_pheno_NA_recoded)[1] - sum(apply(dat_pheno_NA_recoded[metabolite_cols_no3143], 1, anyNA))
# 139 people have data at all time points for all 4 metabolites (not sufficient)


#################################

dim(dat_pheno_complete_NAs_removed)
# 8315 individuals with at least one value in one of the metabolites of interest

# split up by metabolite
chol_cols <- grep("chol_", metabolite_cols_no3143, ignore.case = TRUE, value = TRUE)
hdl_cols <- grep("hdl_", metabolite_cols_no3143, ignore.case = TRUE, value = TRUE)
ldl_cols <- grep("ldl_", metabolite_cols_no3143, ignore.case = TRUE, value = TRUE)
trig_cols <- grep("trig_", metabolite_cols_no3143, ignore.case = TRUE, value = TRUE)

# function that returns how many individuals have at least one value in a metabolite
at_least_one_value <- function(named_cols){
  (dat_pheno_complete_NAs_removed %>% filter(if_any(all_of(named_cols), ~!is.na(.))) %>% dim())[1]
}

at_least_one_value(chol_cols)
at_least_one_value(hdl_cols)
at_least_one_value(ldl_cols)
at_least_one_value(trig_cols)
# 8313-8315 each


# Total sample available in each metabolite at each time point ------------

sample_size_individual_time <- function(named_cols){
  for (i in 1:length(named_cols)){
    one_time <- named_cols[i]
    sample_size <- sum(!is.na(dat_pheno_complete_NAs_removed[[one_time]]))
    print(paste(one_time, "has", sample_size, "number of samples"))
  }
}

sample_size_individual_time(chol_cols)
sample_size_individual_time(hdl_cols)
sample_size_individual_time(ldl_cols)
sample_size_individual_time(trig_cols)
# 868 for BBS otherwise at least 3000 sample size


# Longtitudinal sample size in groups -------------------------------------

# function to create 3 groups (child, teen, adult) for each metabolite
long_completeness_grouping <- function(named_cols, name, data){
  
  # seperate into groups
  vars_child <- named_cols[grepl("F7|F9|BBS", named_cols)]   # age 8-9
  vars_teen <- named_cols[grepl("TF3|TF4", named_cols)]      # age 15-17     
  vars_adult <- named_cols[grepl("F24", named_cols)]    # age 24-30
  
  # make new groups that are averages of age ranges
  dat_grouped <- data %>%
    mutate(
      !!paste0(name,"_group_child") := rowMeans(select(., all_of(vars_child)), na.rm = TRUE),
      !!paste0(name,"_group_teen") := rowMeans(select(., all_of(vars_teen)), na.rm = TRUE),
      !!paste0(name,"_group_adult") := rowMeans(select(., all_of(vars_adult)), na.rm = TRUE)
    )
  return( dat_grouped )
}

# add these groups to new data set
dat_grouped <- long_completeness_grouping(chol_cols, "chol", dat_pheno_complete_NAs_removed)
dat_grouped <- long_completeness_grouping(hdl_cols, "hdl", dat_grouped)
dat_grouped <- long_completeness_grouping(ldl_cols, "ldl", dat_grouped)
dat_grouped <- long_completeness_grouping(trig_cols, "trig", dat_grouped)

# check that the groups have been added correctly
groups <- names(dat_grouped)[grepl("group", names(dat_grouped))]
print(groups)

# seperate the group names out to check for longtitudinal sample sizes
chol_groups <- grep("chol", groups,ignore.case = TRUE, value = TRUE)
hdl_groups <- grep("hdl", groups,ignore.case = TRUE, value = TRUE)
ldl_groups <- grep("ldl", groups,ignore.case = TRUE, value = TRUE)
trig_groups <- grep("trig", groups,ignore.case = TRUE, value = TRUE)

# function to find complete cases
complete_individuals <- function(named_groups){
  sum(complete.cases(dat_grouped[named_groups]))
}

complete_individuals(chol_groups) # 1950
complete_individuals(hdl_groups) # 1563
complete_individuals(ldl_groups) # 1563
complete_individuals(trig_groups) # 1563

# over 1500 individuals with longtitudinal data for all age ranges when grouped


