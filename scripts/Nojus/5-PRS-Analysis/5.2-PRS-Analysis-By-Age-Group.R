
# Libraries ---------------------------------------------------------------

library(tidyverse)
library(performance)

# Load data ---------------------------------------------------------------

args = commandArgs(trailingOnly=TRUE)
INPUT_DATA = args[1]
INPUT_METABOLITE_NAME = args[2]

OUTPUT = args[3]

data <- read_csv(paste(INPUT_DATA))
metabolite <- as.character(INPUT_METABOLITE_NAME)

# save unique time points and PRS scores
age_groups <- unique(data$age_group)
prs_scores <- data %>% select(matches("PGS")) %>% colnames()


# Function to fit linear model --------------------------------------------

# fit linear model transformed metabolite ~ prs_score + sex + P1:PC10
# return a one row tibble of summary stats including coefficient of PRS and R2
fit_lin_model <- function(data, metabolite, prs_score){
  
  # build formula
  covariates <- c(prs_score, "kz021", paste0("PC", 1:10))
  formula <- reformulate(covariates, response = metabolite)
  
  # fit model
  model <- lm(formula, data = data)
  prs_coef <- model$coefficients[[prs_score]]
  metrics <- as.data.frame(model_performance(model))
  
  # save into a one-row tibble
  result <- tibble(
    prs_coef    = prs_coef,
    r2          = metrics$R2,
    aic         = metrics$AIC,
    bic         = metrics$BIC,
    rmse        = metrics$RMSE,
    sigma       = metrics$Sigma
  )
}


# Function to carry out one bootstrap iteration ---------------------------

# returns the summary data of a linear model fitted to a bootstrapped sample
single_bootstrap <- function(data, metabolite, prs_score, unique_ids, n_unique){
  
  # sample ids
  sampled_ids <-  sample(unique_ids, size = n_unique, replace = TRUE)
  
  # For each sampled ID, grab exactly 1 random measurement (exclude repeats)
  boot_data <- data %>%
    filter(cidB4891 %in% sampled_ids) %>%         # Keep only sampled patients
    group_by(cidB4891) %>%                        # Group by patient
    slice_sample(n = 1) %>%                       # Randomly pick 1 row per patient
    ungroup()
  
  # fit linear model
  return(fit_lin_model(boot_data, metabolite, prs_score))
}


# Function to carry out n bootstraps --------------------------------------

# carry out n bootstraps and return a single row of bootstrapped values and confidence intervals 
n_bootstrapped <- function(data, metabolite, prs_score, n_bootstraps){
  
  # get the unique IDs and number of unique ids
  unique_ids <- unique(data$cidB4891)
  n_unique <- length(unique_ids)
  
  # generate bootstrapped results
  results_list <- lapply(1:n_bootstraps, function(i) {
    single_bootstrap(data, metabolite, prs_score, unique_ids, n_unique)
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

all_results <- tibble()

for (prs_score in prs_scores) {
  for (i in age_groups) {
    data_for_analysis <- combined_data  %>% filter(age_group == i)
    result <- n_bootstrapped(data_for_analysis, metabolite, prs_score, 1000) %>% 
      mutate(
        prs_score = prs_score,
        age_group = i
      )
    all_results <- bind_rows(all_results, result)
  }
}


# Save data ---------------------------------------------------------------

all_results %>% write_csv(paste(OUTPUT))