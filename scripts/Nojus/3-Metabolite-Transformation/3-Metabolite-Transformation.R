###########################################################################
######## Transform Metabolites ############################################
###########################################################################

# Load libraries ----------------------------------------------------------

library(tidyverse)
library(data.table)
library(RNOmni)

# Load data ---------------------------------------------------------------

args = commandArgs(trailingOnly=TRUE)
INPUT_PHENO = args[1]
INPUT_EIGENVEC = args[2]
INPUT_EIGENVAL = args[3]
OUTPUT_PHENO = args[4]
OUTPUT_HIST = args[5]

# load pheno data
dat_pheno <- read_csv(paste0(INPUT_PHENO))

# load pcs and their eigenvalues
pcs <- read.table(paste0(INPUT_EIGENVEC), header = FALSE)
eigenval <- read.table(paste0(INPUT_EIGENVAL), header = FALSE)
colnames(pcs) <- c("FID", "IID", paste0("PC", 1:20))


# Check PCs and combine with dat_pheno ------------------------------------

ggplot(pcs, aes(x = PC1, y = PC2)) +
  geom_point(alpha = 0.5) +
  labs(title = "Genetic Principal Components",
       subtitle = "Each point is one individual") +
  theme_minimal()

# calculate variance explained by each PC
variance_explained <- eigenval$V1 / sum(eigenval$V1) * 100

# Plot scree plot
plot(1:20, variance_explained[1:20], 
     type = "b",
     xlab = "Principal Component",
     ylab = "Variance Explained (%)",
     main = "Scree Plot of Genetic PCs")

# combine pheno and PC data
dat_pheno <- dat_pheno %>% inner_join(pcs, by = "IID")
colnames(dat_pheno)


# QC of metabolites -------------------------------------------------------

nas_before <- sum(is.na(dat_pheno))
print(nas_before)

# remove metabolites more than 10sd away from the mean
dat_pheno <- dat_pheno %>%
  group_by(age_group) %>%
  mutate(across(
    c(tc, hdl, ldl, tg),
    ~ {
      mean_val <- mean(.x, na.rm = TRUE)
      sd_val <- sd(.x, na.rm = TRUE)
      ifelse(abs(.x - mean_val) > 10 * sd_val, NA, .x)
    }
  )) %>%
  ungroup()

nas_after <- sum(is.na(dat_pheno))
print(nas_after)

nas_introduced <- nas_after - nas_before
print(nas_introduced)


# Transform Metabolites ---------------------------------------------------

# List of metabolite columns to transform
metabolites <- c("tc", "hdl", "ldl", "tg")

# function to transform metabolites
transform_metabolite <- function(data, metabolite_col, covars = c("age", "bmi", "sex", "season", 
                                                                  paste0("PC", 1:10))) {
  # Build the formula dynamically
  formula <- as.formula(paste0("I(log(", metabolite_col, ")) ~ ", 
                               paste(covars, collapse = " + ")))
  
  # Fit the model with na.action = na.exclude to get residuals with NAs
  model <- lm(formula, data = data, na.action = na.exclude)
  
  # Extract residuals (will have NAs for rows with missing in any predictor)
  res <- residuals(model)
  
  # Apply RankNorm only to non‑NA residuals
  res_non_na <- res[!is.na(res)]
  if (length(res_non_na) > 0) {
    res[!is.na(res)] <- RankNorm(res_non_na)
  }
  
  # Return the transformed residuals (same length as original data)
  return(res)
}

for (metabolite in metabolites) {
  new_col_name <- paste0(metabolite, "_transformed")
  dat_pheno[[new_col_name]] <- transform_metabolite(dat_pheno, metabolite)
}


# Plot the transformation -------------------------------------------------

# the col names of transformed metabolites
transformed_cols <- c("tc_transformed", "hdl_transformed", "ldl_transformed", "tg_transformed")

# create seperate datasets
data1 <- dat_pheno %>% mutate(set = "BiggerSet")
data2 <- dat_pheno %>% filter(all_three_times > 0) %>% mutate(set = "SmallerSet")

# plot metabolites to check the transformation was successful
# check the distributions allign
p <- bind_rows(data1, data2) %>%
  select(age_group, all_of(transformed_cols), set) %>%
  pivot_longer(
    cols = -c(age_group, set),
    names_to = "metabolite",
    values_to = "transformed_value"
  ) %>%
  mutate(metabolite = str_remove(metabolite, "_transformed")) %>%
  ggplot(aes(x = transformed_value, fill = set, color = set)) +
  geom_histogram(
    aes(y = after_stat(density)),   # <-- frequency (density) on y‑axis
    bins = 30,
    alpha = 0.4,
    position = "identity"
  ) +
  facet_grid(metabolite ~ age_group, scales = "free") +
  theme_bw() +
  labs(
    title = "Distribution of transformed metabolites by age group",
    x = "Transformed value (RankNorm of residuals)",
    y = "Density (frequency proportion)",
    fill = "Dataset",
    color = "Dataset"
  ) +
  theme(strip.text = element_text(face = "bold", size = 10))


# Save data ---------------------------------------------------------------

dat_pheno %>% write_csv(paste0(OUTPUT_PHENO))

ggsave(
  filename = paste0(OUTPUT_HIST),
  plot = p,
  width = 19,
  height = 15,
  units = "cm",
  dpi = 300
)