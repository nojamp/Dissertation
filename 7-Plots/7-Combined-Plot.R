###########################################################################
######## Plot of PRS performance for each age group #######################
###########################################################################

# load libraries ----------------------------------------------------------

library(tidyverse)
library(grid)
library(gridExtra)


# Load data ---------------------------------------------------------------

args = commandArgs(trailingOnly=TRUE)
INPUT = args[1]
OUTPUT = args[2]

df <- read_csv(paste0(INPUT))

# Pre-process before plotting ---------------------------------------------

# put age at appropriate numeric position (median of age)
df$age_group_pos <- c(8.7,16.5,24.0)[match(df$age_group, c("children", "teens", "adults"))]

# create dataframe with max and min performance of OmicsPred
baseline <- data.frame(
  metabolite = c("TC", "HDL","LDL", "TG"),
  max_r2 = c(0.087, 0.104, 0.126, 0.074),
  min_r2 = c(0.071, 0.062, 0.062, 0.061)
)

# rename
df <- df %>% rename(PRS = prs_name)

# set maximum y value for plot
y_max <- max(df$r2_uci) + 0.01

# set appropriately contrasting colours for each group of PRS
colour_df <- data.frame(
  colour = c(
    "#D62728", "#1F77B4", "#2CA02C",   # TC : red, blue, green
    "#FF7F0E", "#9467BD", "#17BECF",   # HDL: orange, purple, cyan
    "#F8766D", "#00BFC4", "#C51B7D",   # LDL: coral, teal, magenta
    "#8C564B", "#E377C2", "#000000"    # TG : brown, pink, black
  ),
  metabolite = rep(c("TC", "HDL", "LDL", "TG"), each = 3)
)

# Function to plot PRS for each metabolite --------------------------------

plot_by_age_group <- function(metabolite_name, data){
  
  # filter for metabolite
  plotting_data <- data %>% filter(metabolite == paste(metabolite_name))
  
  # order PRS by which one is top performing for children
  order_by_first <- plotting_data %>%
    filter(age_group == "children") %>%
    arrange(desc(r2)) %>%
    pull(PRS)
  
  plotting_data <- plotting_data %>% mutate(PRS = factor(PRS, levels = order_by_first))
  
  
  p <- ggplot(plotting_data, aes(x = age_group_pos, y = r2, colour = PRS, shape = PRS)) +
    
    # Shaded region with fill mapping
    geom_rect(
      data = baseline %>% filter(metabolite == paste(metabolite_name)),
      aes(xmin = -Inf, xmax = Inf, ymin = min_r2, ymax = max_r2),
      fill = "gray75",
      alpha = 0.4,
      color = "gray50",
      linetype = "dashed",
      size = 0.3,
      inherit.aes = FALSE
    )  +
    
    # add points, 95% CIs, and connect the points with lines
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = r2_lci, ymax = r2_uci),
                  position = position_dodge(width = 0.15),
                  linewidth = 0.8, alpha = 0.3) +
    geom_line(position = position_dodge(width = 0.15), linewidth = 0.9, alpha = 0.8) +
    
    # label x axis appropriately
    scale_x_continuous(
      breaks = c(8.7, 16.5, 24.0),
      labels = c("children", "teens", "adults"),
    ) +
    
    # use selected colours
    scale_color_manual(
      values = colour_df %>% filter(metabolite == paste(metabolite_name)) %>% select(colour) %>% pull()
    ) +
    
    # no y or x axis label since we will combine. alternatives provided commented out for single plots
    ylab(NULL) +    #ylab(paste(expression(R^2), "of linear model")) 
    xlab(NULL) +    #xlab("Age group")
    
    # set the limits to be the same for all plots (so they can be combined)
    scale_y_continuous(limits = c(0, y_max), expand = c(0,0) )+
    
    # put legend at the bottom
    guides(
      color = guide_legend(order = 1, override.aes = list(size = 4), position = "bottom", title = NULL, ncol = 1),
      shape = guide_legend(order = 1, title = NULL)   # shape legend will be merged with color
    )  + theme_bw()
  
  return(p)
}


# Plot for each metabolite ------------------------------------------------

# add title (metabolite name)
p_tc <- plot_by_age_group("TC", df) + ggtitle("TC") + theme(plot.title = element_text(hjust = 0.5))
p_hdl <- plot_by_age_group("HDL", df) + ggtitle("HDL-C") + theme(plot.title = element_text(hjust = 0.5))
p_ldl <- plot_by_age_group("LDL", df) + ggtitle("LDL-C") + theme(plot.title = element_text(hjust = 0.5))
p_tg <- plot_by_age_group("TG", df) + ggtitle("TG") + theme(plot.title = element_text(hjust = 0.5))


# Combined plot -----------------------------------------------------------

# combine all plots
# add y-axis label
combined_plot <- arrangeGrob(p_tc, p_hdl, p_ldl, p_tg, 
             nrow = 1,
             top = "PRS performance by metabolite by PRS",
             left = textGrob(expression(R^2 ~ "of linear model"), rot = 90, y = 0.6))

# add legend label
combined_plot <- grobTree(
  combined_plot,
  textGrob("PRS\nname", 
           x = 0.05, y = 0.10, 
           rot = 0, 
           gp = gpar(fontsize = 12))
)


# Save plot ---------------------------------------------------------------

ggsave(
  filename = paste0(OUTPUT),
  plot = combined_plot,
  width = 19,
  height = 15,
  units = "cm",
  dpi = 300
)
