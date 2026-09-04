###########################################################################
######## PRS Analysis Longtitudinals ######################################
###########################################################################

# Load libraries ----------------------------------------------------------

library(tidyverse)
library(splines)
library(glmmTMB)
library(performance)
library(parallel)

# Load data ---------------------------------------------------------------

args = commandArgs(trailingOnly=TRUE)
INPUT = args[1]
OUTPUT = args[2]

dat <- read_csv(paste0(INPUT))


# Function to fit LMEM models and return summaries ------------------------

fit_LMEMs <- function(data, prs_score, transformed_metabolite){

  # formula for LMEM no interaction
  formula1_str <- paste(
    transformed_metabolite, "~",
    prs_score,
    "+ns(age, df = 3)+ (1 + age | cidB4891)"
  )

  # formula for LMEM with interaction
  formula2_str <- paste(
    transformed_metabolite, "~",
    prs_score,
    "*ns(age, df = 3)+ (1 + age | cidB4891)"
  )

  # convert to formula objects
  formula1_obj <- as.formula(formula1_str)
  formula2_obj <- as.formula(formula2_str)

  # fit 3 models
  model1 <- glmmTMB(formula1_obj, data = data, REML = TRUE)
  model2 <- glmmTMB(formula1_obj, data = data, REML = FALSE)
  model3 <- glmmTMB(formula2_obj, data = data, REML = FALSE)

  # extract key stats from model
  coef_summary <- confint(model1, parm=prs_score, estimate = TRUE)
  metrics <- performance(model1)
  nobservations <- length(unique(data$cidB4891))

  # use F-test to compare model2 and model3
  comparison <- anova(model2, model3)
  p_value <- comparison$`Pr(>Chisq)`[2]

  # save in one row tibble
  result <- tibble(
    prs_name = prs_score,
    prs_coef_lci= coef_summary[1],
    prs_coef    = coef_summary[2],
    prs_coef_uci= coef_summary[3],
    R2_marginal = metrics$R2_marginal,
    R2_conditional = metrics$R2_conditional,
    number_observations = nobservations,
    interaction_pval = p_value
  )
  return(result)
}


# Run analysis ------------------------------------------------------------

# Get all PRS columns
metab_cols <- c("tc", "tg", "hdl", "ldl")

all_prs_cols <- unique(unlist(lapply(metab_cols, function(m) {
  str_subset(names(dat), regex(paste0("_", m), ignore_case = TRUE))
})))

# Define age groups
age_groups <- c("children", "teens", "adults")

all_results <- bind_rows(mclapply(all_prs_cols, function(prs) {
  metab <- str_extract(prs, "[^_]+$")
  # Find the corresponding transformed column (case‑insensitive)
  metab_trans <- grep(paste0(metab, "_transformed"), names(dat),
                      ignore.case = TRUE, value = TRUE)
  data_for_analysis <- dat %>%
    group_by(cidB4891) %>%
    # For each person, check that all three age groups appear among rows with non‑missing tc
    filter(all(c("children", "teens", "adults") %in% age_group[!is.na(.data[[metab_trans]])])) %>%
    ungroup()
  fit_LMEMs(data_for_analysis, prs, metab_trans)
},mc.cores = 12))


# add note about which metabolite each PRS is for
all_results <- all_results %>%
  mutate(
    metabolite = str_extract(prs_name, "[^_]+$") %>% toupper(),
    prs_name   = str_remove(prs_name, "_[^_]+$")
  )


# Save results ------------------------------------------------------------

all_results %>% write_csv(paste0(OUTPUT))