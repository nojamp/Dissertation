###########################################################################
######## PRS Analysis for each age group ##################################
###########################################################################


# load libraries ----------------------------------------------------------

library(tidyverse)
library(performance)
library(parallel)

set.seed(13)

# load data ---------------------------------------------------------------

args = commandArgs(trailingOnly=TRUE)
INPUT = args[1]
OUTPUT = args[2]

dat <- read_csv(paste0(INPUT))


# Function to fit linear model --------------------------------------------

# fit linear model transformed metabolite ~ prs_score
# return a one row tibble of summary stats including coefficient of PRS and R2

fit_lin_model <- function(data, transformed_metabolite, prs_score){

  # build formula
  covariates <- c(prs_score)
  formula <- reformulate(covariates, response = transformed_metabolite)

  # fit model
  model <- lm(formula, data = data)
  prs_coef <- model$coefficients[[prs_score]]
  metrics <- as.data.frame(model_performance(model))

  # save into a one-row tibble
  result <- tibble(
    prs_coef    = prs_coef,
    r2          = metrics$R2,
    nobs        = nobs(model)
  )
}

# Function to carry out one bootstrap iteration ---------------------------

# returns the summary data of a linear model fitted to a bootstrapped sample
single_bootstrap <- function(data, transformed_metabolite, prs_score, unique_ids, n_unique){

  # sample ids
  sampled_ids <-  sample(unique_ids, size = n_unique, replace = TRUE)

  # For each sampled ID, grab exactly 1 random measurement (exclude repeats)
  boot_data <- data %>%
    filter(cidB4891 %in% sampled_ids) %>%         # Keep only sampled patients
    group_by(cidB4891) %>%                        # Group by patient
    slice_sample(n = 1) %>%                       # Randomly pick 1 row per patient
    ungroup()

  # fit linear model
  return(fit_lin_model(boot_data, transformed_metabolite, prs_score))
}


# Function to carry out n bootstraps --------------------------------------

# carry out n bootstraps and return a single row of bootstrapped values and confidence intervals
n_bootstrapped <- function(data, transformed_metabolite, prs_score, n_bootstraps){

  # get the unique IDs and number of unique ids
  unique_ids <- unique(data$cidB4891)
  n_unique <- length(unique_ids)

  # generate bootstrapped results
  results_list <- lapply(1:n_bootstraps, function(i) {
    single_bootstrap(data, transformed_metabolite, prs_score, unique_ids, n_unique)
  })

  # combine bootstrapped results into one data frame
  final_results <- bind_rows(results_list)

  # calculate the bootstrapped values and confidence values
  summary_row <- final_results %>%
    summarise(
      across(
        where(is.numeric),
        list(
          mean = ~mean(.x, na.rm = TRUE),
          lci  = ~quantile(.x, 0.025, na.rm = TRUE),
          uci  = ~quantile(.x, 0.975, na.rm = TRUE)
        ),
        .names = "{.col}_{.fn}"
      )
    ) %>%
    # Rename: remove the "_mean" suffix so the base metric stands alone
    rename_with(~ gsub("_mean$", "", .x), ends_with("_mean")) %>%
    mutate(prs_name = prs_score) %>% relocate(prs_name)

  return(summary_row)
}


# Carry out the analysis for each PRS score and age group -----------------

# Get all PRS columns
metab_cols <- c("tc", "tg", "hdl", "ldl")

all_prs_cols <- unique(unlist(lapply(metab_cols, function(m) {
  str_subset(names(dat), regex(paste0("_", m), ignore_case = TRUE))
})))

# Define age groups
age_groups <- c("children", "teens", "adults")

# Nested lapply: outer over PRS, inner over age groups
all_results <- bind_rows(mclapply(all_prs_cols, function(prs) {
  # Extract metabolite suffix from PRS column name
  metab <- str_extract(prs, "[^_]+$")
  # Find the corresponding transformed column (case‑insensitive)
  metab_trans <- grep(paste0(metab, "_transformed"), names(dat),
                      ignore.case = TRUE, value = TRUE)

  # Inner loop over age groups – return a tibble for each combination
  results_for_prs <- bind_rows(lapply(age_groups, function(i) {
    data_for_analysis <- dat %>% filter(age_group == i) %>% filter(!is.na(.data[[metab_trans]]))
    res <- n_bootstrapped(data_for_analysis, metab_trans, prs, 1000)
    res %>% mutate(prs_score = prs, age_group = i)
  }))
  results_for_prs
},mc.cores=12))


# add note about which metabolite each PRS is for
all_results <- all_results %>%
  mutate(
    metabolite = str_extract(prs_name, "[^_]+$") %>% toupper(),
    prs_name   = str_remove(prs_name, "_[^_]+$")
  )

# Save results ------------------------------------------------------------

all_results %>% write_csv(paste0(OUTPUT))