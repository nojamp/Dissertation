# Load libraries ----------------------------------------------------------

library(tidyverse)
library(splines)
library(glmmTMB)
library(performance)

# Load data ---------------------------------------------------------------

args = commandArgs(trailingOnly=TRUE)
INPUT_DATA = args[1]
METABOLITE_NAME = args[2]
OUTPUT_DATA = args[3]

data <- read_csv(paste(INPUT_DATA))

prs_scores <- data %>% select(matches("PGS")) %>% colnames()
metabolite <- as.character(METABOLITE_NAME)


# Function to fit Linear Mixed-Effects Models (LMEMs) ---------------------

fit_LMEM <- function(data, prs_score, metabolite, include_interaction=FALSE, REML_option=TRUE){
  
  interaction_symbol <- if (include_interaction) "*" else "+"
  
  formula_str <- paste(
    metabolite, "~",
    prs_score, interaction_symbol, "ns(age_std, df = 3)",
    "+ kz021",
    paste0("+PC", 1:10),
    "+ (1 + age_std | cidB4891)"
  )
  
  formula_obj <- as.formula(formula_str)
  
  return(glmmTMB(formula_obj, data = data, REML = REML_option))
}


# Function to fit LMEM main model and extract coefficients ----------------

summarise_LMEM_main <- function(data, prs_score, metabolite){
  main_model <- fit_LMEM(data = data, prs_score = prs_score, metabolite = metabolite, include_interaction = FALSE, REML_option = TRUE)
  
  coef_summary <- confint(main_model, parm=prs_score, estimate = TRUE)
  metrics <- as.data.frame(model_performance(main_model))

  # save into a one-row tibble
  result <- tibble(
    prs_name = prs_score,
    metabolite = metabolite,
    prs_coef_lci= coef_summary[1],
    prs_coef    = coef_summary[2],
    prs_coef_uci= coef_summary[3],
    R2_conditional= metrics$$R2_conditional,
    R2_marginal = metrics$R2_marginal,
    aic         = metrics$AIC,
    bic         = metrics$BIC,
    rmse        = metrics$RMSE,
    sigma       = metrics$Sigma
  )
  return(result)
}


# Function to fit LMEM with and without interaction and compare -----------

test_LMEM_interaction <- function(data, prs_score, metabolite){
  
  # fit LMEM with and without interaction
  # use REML_option = FALSE as looking at fixed effects
  no_int <- fit_LMEM(data = data, prs_score = prs_score, metabolite = metabolite, include_interaction = FALSE, REML_option = FALSE)
  with_int <- fit_LMEM(data = data, prs_score = prs_score, metabolite = metabolite, include_interaction = TRUE, REML_option = FALSE)
  
  # use anova test
  comparison <- anova(no_int, with_int)
  
  # save in single row tibble
  result <- tibble(
    prs_name = prs_score,
    interaction_test = comparison$`Pr(>Chisq)`,
  )
  
  return(result)
}


# Run analyses ------------------------------------------------------------

all_results <- tibble()
for (i in prs_scores) {
  print(i)
  result1 <- summarise_LMEM_main(data, i, metabolite)
  result2 <- test_LMEM_interaction(data, i, metabolite)
  
  result <- full_join(result1, result2, by = "prs_name")
  
  all_results <- bind_rows(all_results, result)
  
}


# Save results ------------------------------------------------------------

all_results %>% write_csv(paste(OUTPUT_DATA))
