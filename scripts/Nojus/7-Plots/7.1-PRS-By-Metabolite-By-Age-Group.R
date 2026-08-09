###########################################################################
######## Plot of PRS performance for each age group #######################
###########################################################################

# load libraries ----------------------------------------------------------

library(tidyverse)
library(ggguides)


# Load data ---------------------------------------------------------------

args = commandArgs(trailingOnly=TRUE)
INPUT = args[1]
OUTPUT = args[2]

df <- read_csv(paste0(INPUT))

# Pre-process before plotting ---------------------------------------------

# put age at appropriate numeric position (median of age)
df$age_group_pos <- c(8.7,16.6,24.0)[match(df$age_group, c("children", "teens", "adults"))]

# Categorise (by hand) the PRS scores by method
PRS_CS_scores = c("PGS002781", "PGS003138", "PGS003148")
PRS_CSx_scores = c("PGS003978", "PGS004981")
ensemble_scores = c("PGS004156")
LDpred_scores = c("PGS004974")
LDpred2_scores = c("PGS004043","PGS003819", "PGS004333", "PGS004342", "PGS003802")
LDpred_family = c(LDpred_scores, LDpred2_scores)
PRS_family = c(PRS_CS_scores, PRS_CSx_scores)


df <- df %>% mutate(
  prs_name = str_extract(prs_name, "^PGS\\d+"),
  # create a new category variable for PRS method
  prs_method = case_when(
    prs_name %in% PRS_family ~ "PRS-CS Family",
    prs_name %in% LDpred_family ~ "LDpred Family",
    prs_name %in% ensemble_scores ~ "Ensemble"
    )
  )

# create a data frame to label PRSs on the plot
line_labels <- df %>% group_by(metabolite, prs_name) %>%
  filter(age_group_pos == max(age_group_pos)) %>%
  ungroup()

# create dataframe with max and min performance of OmicsPred
baseline <- data.frame(
  metabolite = c("TC", "TG", "HDL", "LDL"),
  max_r2 = c(0.087, 0.08, 0.104, 0.126),
  min_r2 = c(0.071, 0.061, 0.062, 0.062)
)

# Plot --------------------------------------------------------------------

# 4 plots side by side with identical axis of R^2 of each score for each age group
p <- ggplot(df, aes(x = age_group_pos, y = r2, color = prs_method, group = prs_name, shape = prs_method)) +

  # separate graph for each metabolite all in one row
  facet_wrap(~ fct_relevel(metabolite, "TC", "HDL", "LDL", "TG"),
           scales = "fixed", nrow = 1)  +
  # Shaded region with fill mapping
  geom_rect(
    data = baseline,
    aes(xmin = -Inf, xmax = Inf, ymin = min_r2, ymax = max_r2,
        fill = "Baseline range in\nOmicsPred"),   # <-- creates a new fill aesthetic
    alpha = 0.4,
    color = "gray30",
    linetype = "dashed",
    size = 0.3,
    inherit.aes = FALSE
  ) +

  # add points, error bars using 95% CIs, and join them with lines
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = r2_lci, ymax = r2_uci),
              position = position_dodge(width = 0.15),
              size = 0.8, alpha = 0.3) +
  geom_line(position = position_dodge(width = 0.15), linewidth = 0.9, alpha = 0.8) +

  # label x axis appropriately
  scale_x_continuous(
    breaks = c(8.5, 17.5, 24),
    labels = c("children", "teens", "adults"),
    limits = c(7,28) # create extra space for PRS name annotations
  ) +

  # manually assign color and shape
  # shape and color indicate the same family of methods used to develop the PRS
  # this is used to make the distinctions clearer
  scale_color_manual(
    name = "PRS Development Method:",
    values = c("PRS-CS Family" = "#009E73", "LDpred Family" = "#D55E00", "Ensemble" = "#56B4E9")
  ) +
  scale_shape_manual(
    name = "PRS Development Method:",   # <-- same title as color legend
    values = c("PRS-CS Family" = 16, "LDpred Family" = 17, "Ensemble" = 15)
  ) +

  # add labels for each PRS
#  geom_text(
#    data = line_labels,
#    aes(label = prs_name),
#    hjust = -0.15,
#    size = 1.5,
#    show.legend = FALSE,
#    color = "black"
#  ) +
  ylab(paste(expression(R^2), "of linear model")) + xlab(NULL) +
  scale_fill_manual(
    name = NULL,   # no title, or set to "Reference"
    values = c("Baseline range in\nOmicsPred" = "gray80"),
    guide = guide_legend(
      order = 2,                    # <-- place this legend after the color/shape one
      override.aes = list(
        alpha = 0.4,
        color = "gray30",
        linetype = "dashed",
        size = 0.8
      )
    )
  ) +
  guides(
    color = guide_legend(order = 1, override.aes = list(size = 4), position = "bottom"),
    shape = guide_legend(order = 1),   # shape legend will be merged with color,
    fill = guide_legend(order = 2, position = "bottom")
  ) + theme_gray()

# Save plot ---------------------------------------------------------------

ggsave(
  filename = paste0(OUTPUT),
  plot = p,
  width = 19,
  height = 15,
  units = "cm",
  dpi = 300
)