
# remove twins


# pc analysis

# season variable


# groups

# summary of vars

# summary of longtitudinal stuff


######

# 10sd remove

# metabolite transformation

# save


library(tidyverse)
library(data.table)
library(RNOmni)


# Load data ---------------------------------------------------------------

# load pheno data and omics ids
dat_pheno <- read_csv("pheno_nas_recoded.csv")
omics_id <- fread("clean_unpruned.psam", header = TRUE)

# load pcs
pcs <- read.table("genetic_pcs20.eigenvec", header = FALSE)
eigenval <- read.table("genetic_pcs20.eigenval", header = FALSE)
colnames(pcs) <- c("FID", "IID", paste0("PC", 1:20))

# only save the people with genetic data that passed the thresholds
dat_pheno <- dat_pheno %>% inner_join(omics_id, by = "IID")
dim(dat_pheno)
# 6780 people

# check no twins
dat_pheno %>% count(cidB4891) %>% filter(n>1) %>% nrow()


# Check PCs ---------------------------------------------------------------

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


# Pre-process data (unify units) ------------------------------------------

# vars for F7, F9, TF3, TF4, F24 are coded to start with f7, f9, fh, FJ and FK respectively
age_cols <- c("f7003c", "f9003c", "fh0011a", "FJ003b", "FKAR0011")
sex_col <- "kz021"
bmi_cols <- c("f7ms026a", "f9ms026a", "fh3019", "FJMR022a", "FKMS1040")
month_cols <- c("f7001", "f9001", "fh0010a", "FJ002a", "FKAR0040")

# metabolite cols
tc_cols <- dat_pheno %>% select(matches("chol_", ignore.case = TRUE)) %>% select(matches("f7|f9|tf3|tf4|f24", ignore.case = TRUE)) %>% colnames()
hdl_cols <- dat_pheno %>% select(matches("hdl_", ignore.case = TRUE)) %>% select(matches("f7|f9|tf3|tf4|f24", ignore.case = TRUE)) %>% colnames()
ldl_cols <- dat_pheno %>% select(matches("ldl_", ignore.case = TRUE)) %>% select(matches("f7|f9|tf3|tf4|f24", ignore.case = TRUE)) %>% colnames()
tg_cols <- dat_pheno %>% select(matches("trig_", ignore.case = TRUE)) %>% select(matches("f7|f9|tf3|tf4|f24", ignore.case = TRUE)) %>% colnames()
# vars for membership in each follow up group
membership_cols <- c("in_f07", "in_f09", "in_tf3", "in_tf4", "in_F24")

# convert age to all have same units
# f7003c, f9003c, fh0011a age in months
# FJ003b, FKAR0011 age in years
dat_pheno <- dat_pheno %>% 
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
dat_pheno %>% select(all_of(age_cols)) %>% summary()
membership_cols %>%
  map_dfr(~ dat_pheno %>%
            filter(!is.na(.data[[.x]])) %>%
            count(kz021) %>%
            mutate(prop = n / sum(n),
                   subgroup = .x),
          .id = NULL) %>%
  select(subgroup, kz021, n, prop) %>%
  print()
dat_pheno %>% select(all_of(bmi_cols)) %>% summary()
dat_pheno %>% select(all_of(month_cols)) %>% summary()
dat_pheno %>%
  select(all_of(month_cols)) %>%
  map(~ prop.table(table(.x)) * 100)

dat_pheno %>% select(all_of(tc_cols)) %>% summary()
dat_pheno %>% select(all_of(hdl_cols)) %>% summary()
dat_pheno %>% select(all_of(ldl_cols)) %>% summary()
dat_pheno %>% select(all_of(tg_cols)) %>% summary()

# Summarise constructed groups --------------------------------------------

lookup_table <- tibble(
  focus = c("f7", "f9", "tf3", "tf4", "f24"),
  membership = membership_cols,
  age = age_cols,
  sex = c("kz021", "kz021", "kz021", "kz021", "kz021"),
  bmi = bmi_cols,
  season = month_cols,
  tc = tc_cols,
  hdl = hdl_cols,
  ldl = ldl_cols,
  tg = tg_cols
)

# check it lines up
lookup_table

# create a row for each follow up

dat_long <- lookup_table %>%
  pmap_dfr(function(focus, membership, age, sex, bmi, season, tc, hdl, ldl, tg) {
    
    dat_pheno %>%
      # Select the ID, the membership column, and all variable columns for this time point
      select(cidB4891,
             qlet,
             gi_1000g_g0m_g1,
             IID,
             all_of(membership), 
             all_of(age), 
             all_of(sex), 
             all_of(bmi), 
             all_of(season), 
             all_of(tc), 
             all_of(hdl), 
             all_of(ldl), 
             all_of(tg)) %>%
      # Keep only rows where the person attended this follow‑up
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

# View the result
glimpse(dat_long)

# add group names
dat_long <- dat_long %>%
  mutate(
    age_group = case_when(
      focus %in% c("f7", "f9")    ~ "children",
      focus %in% c("tf3", "tf4")  ~ "teens",
      focus == "f24"              ~ "adults",
      TRUE                        ~ NA_character_
    ) %>%
      factor(levels = c("children", "teens", "adults"), ordered = TRUE)
  )


# summarise
numeric_summary <- dat_long %>%
  group_by(age_group) %>%
  summarise(
    across(
      c(age, bmi, tc, hdl, ldl, tg),
      list(
        n     = ~ sum(!is.na(.x)),
        mean  = ~ mean(.x, na.rm = TRUE),
        sd    = ~ sd(.x, na.rm = TRUE),
        median = ~ median(.x, na.rm = TRUE),
        min   = ~ min(.x, na.rm = TRUE),
        max   = ~ max(.x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

numeric_summary$tg_n

# For sex
sex_summary <- dat_long %>%
  group_by(age_group, sex) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(age_group) %>%
  mutate(percent = n / sum(n) * 100)

print(sex_summary)

# For season
season_summary <- dat_long %>%
  group_by(age_group, season) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(age_group) %>%
  mutate(percent = n / sum(n) * 100)

print(season_summary)

# how many repeats per age_group
dat_long %>% group_by(age_group) %>% count(cidB4891) %>% filter(n>1) %>% nrow()
dat_long %>%
  group_by(age_group, cidB4891) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(n > 1) %>%
  count(age_group, name = "participants_with_repeats")

# sample size for each
dat_long %>% filter(age_group == "children") %>% group_by(cidB4891) %>% sample_n(1) %>% ungroup() %>% dim()
dat_long %>% filter(age_group == "teens") %>% group_by(cidB4891) %>% sample_n(1) %>% ungroup() %>% dim()
dat_long %>% filter(age_group == "adults") %>% dim()

# QC of metabolites -------------------------------------------------------

# remove metabolites more than 10sd away from the mean
dat_long <- dat_long %>%
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


# Metabolite Transformation -----------------------------------------------

# only use data with non NAs in variables used for transformation
dat_long <- dat_long %>% filter(is.na(age) == FALSE, is.na(sex) == FALSE, is.na(bmi) == FALSE, is.na(season) == FALSE)

# add data on PCs
dat_long <- dat_long %>% inner_join(pcs, by = "IID")

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

# List of metabolite columns to transform
metabolites <- c("tc", "hdl", "ldl", "tg")

for (metabolite in metabolites) {
  new_col_name <- paste0(metabolite, "_transformed")
  dat_long[[new_col_name]] <- transform_metabolite(dat_long, metabolite)
}

# plot metabolites to check the transformation was successful
dat_long %>%
  select(age_group, all_of(transformed_cols)) %>%
  pivot_longer(
    cols = -age_group,
    names_to = "metabolite",
    values_to = "transformed_value"
  ) %>%
  mutate(metabolite = str_remove(metabolite, "_transformed")) %>%
  ggplot(aes(x = transformed_value)) +
  geom_histogram(
    aes(y = after_stat(density)),   # <-- frequency (density) on y‑axis
    bins = 30,
    fill = "steelblue",
    color = "black",
    alpha = 0.8
  ) +
  facet_grid(metabolite ~ age_group, scales = "free") +
  theme_bw() +
  labs(
    title = "Distribution of transformed metabolites by age group",
    x = "Transformed value (RankNorm of residuals)",
    y = "Density (frequency proportion)"
  ) +
  theme(strip.text = element_text(face = "bold", size = 10))


# Summaries of sample sizes -----------------------------------------------

# summary of sample sizes for each age group for each metabolite
dat_long %>%
  group_by(age_group) %>%
  summarise(across(ends_with("_transformed"), ~ sum(!is.na(.)))) %>%
  pivot_longer(-age_group, names_to = "metabolite", values_to = "n_non_na") %>%
  mutate(metabolite = str_remove(metabolite, "_transformed")) %>%
  pivot_wider(names_from = age_group, values_from = n_non_na) %>%
  arrange(metabolite)

# summary of sample sizes for each focus group
dat_long %>%
  # Select participant ID, focus, and transformed columns
  select(cidB4891, focus, ends_with("_transformed")) %>%
  # Pivot to long format (metabolite names and values)
  pivot_longer(
    cols = ends_with("_transformed"),
    names_to = "metabolite",
    values_to = "value"
  ) %>%
  # Clean metabolite names (remove "_transformed")
  mutate(metabolite = str_remove(metabolite, "_transformed")) %>%
  # Keep only non‑missing values
  filter(!is.na(value)) %>%
  # For each focus and metabolite, count distinct participants
  group_by(focus, metabolite) %>%
  summarise(n_participants = n_distinct(cidB4891), .groups = "drop") %>%
  # Pivot to wide format for easier comparison
  pivot_wider(
    names_from = focus,
    values_from = n_participants
  ) %>%
  arrange(metabolite)

# summary of longtitudinal sample sizes (metabolite values at all 3 age groups)
dat_long %>%
  # Keep only participant ID, age group, and transformed metabolites
  select(cidB4891, age_group, ends_with("_transformed")) %>%
  # Pivot to long format (one row per metabolite per person per age group)
  pivot_longer(
    cols = ends_with("_transformed"),
    names_to = "metabolite",
    values_to = "value"
  ) %>%
  # Clean the metabolite name (remove "_transformed" suffix)
  mutate(metabolite = str_remove(metabolite, "_transformed")) %>%
  # Remove rows where the value is missing (we only care about non-NA)
  filter(!is.na(value)) %>%
  # For each participant and metabolite, count how many distinct age groups have data
  group_by(cidB4891, metabolite) %>%
  summarise(
    n_age_groups_present = n_distinct(age_group),
    .groups = "drop"
  ) %>%
  # Keep only those who appear in all 3 age groups
  filter(n_age_groups_present == 3) %>%
  # Count how many participants meet this criterion per metabolite
  count(metabolite, name = "participants_with_all_3_ages")

# summary of longtitudinal sample sizes if using focus groups

# Define the total number of focus groups (5 in this case)
total_focus_groups <- dat_long$focus %>% n_distinct()

dat_long %>%
  # Keep participant ID, focus, and transformed metabolites
  select(cidB4891, focus, ends_with("_transformed")) %>%
  # Pivot to long format (one row per metabolite per person per focus)
  pivot_longer(
    cols = ends_with("_transformed"),
    names_to = "metabolite",
    values_to = "value"
  ) %>%
  # Clean the metabolite name
  mutate(metabolite = str_remove(metabolite, "_transformed")) %>%
  # Remove rows where the value is missing
  filter(!is.na(value)) %>%
  # For each participant and metabolite, count how many distinct focus groups have data
  group_by(cidB4891, metabolite) %>%
  summarise(
    n_focus_groups_present = n_distinct(focus),
    .groups = "drop"
  ) %>%
  # Keep only those who appear in ALL focus groups
  filter(n_focus_groups_present == total_focus_groups) %>%
  # Count how many participants meet this criterion per metabolite
  count(metabolite, name = "participants_with_all_5_focus_groups")


# summary of repeats in each group
dat_long %>%
  select(cidB4891, age_group, ends_with("_transformed")) %>%
  pivot_longer(
    cols = ends_with("_transformed"),
    names_to = "metabolite",
    values_to = "value"
  ) %>%
  mutate(metabolite = str_remove(metabolite, "_transformed")) %>%
  filter(!is.na(value)) %>%
  group_by(age_group, metabolite, cidB4891) %>%
  summarise(n_measurements = n(), .groups = "drop") %>%
  filter(n_measurements > 1) %>%
  count(age_group, metabolite, name = "participants_with_repeats") %>%
  pivot_wider(
    names_from = age_group,
    values_from = participants_with_repeats
  ) %>%
  arrange(metabolite)


# Save data ---------------------------------------------------------------

dat_long %>% write_csv()
