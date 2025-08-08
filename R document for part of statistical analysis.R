# --- 0. Load Necessary Libraries ---
print("--- 0. Loading Libraries ---")
# Ensure these are installed: install.packages(c("readr", "dplyr", "ordinal", "ggplot2", "tidyr", "ggpattern", "rstatix"))
library(readr); library(dplyr); library(ordinal)
library(ggplot2); library(tidyr)
library(ggpattern); library(rstatix)

# --- 1. Load and Prepare Data ---
print("--- 1. Loading and Preparing Data ---")
Data <- readr::read_csv("data_for_aq_total_anova.csv")

# --- Data preparation code (unchanged) ---
iv1_col_name <- "Epilepsy_Condition_Num"
iv2_col_name <- "Depression_Condition_Num"
education_col_name <- "Q4" 
sds_item_names <- c("SDS_1", "SDS_2", "SDS_3", "SDS_4", "SDS_5", "SDS_6", "SDS_7")
aq_item_names <- c("AQ_Pity", "AQ_Danger", "AQ_Fear", "AQ_Blame", "AQ_Anger", "AQ_No_Help", "AQ_Coercion")
all_dvs <- c(aq_item_names, sds_item_names)
sds_name_map <- c("SDS_1" = "Neighbour", "SDS_2" = "Socialise", "SDS_3" = "Friends", "SDS_4" = "Work", "SDS_5" = "Marry", "SDS_6" = "Rent", "SDS_7" = "Hire")
sds_labels <- c("Definitely Willing", "Probably Willing", "Probably Unwilling", "Definitely Unwilling")
aq_name_map <- setNames(gsub("AQ_", "", aq_item_names), aq_item_names)
aq_levels <- as.character(1:9)
aq_fill_palette <- grey.colors(9, start = 0.95, end = 0.3); names(aq_fill_palette) <- aq_levels
aq_pattern_palette <- setNames(rep(c('none', 'stripe', 'crosshatch', 'pch'), length.out = 9), aq_levels)
sds_fill_palette <- grey.colors(4, start=0.95, end=0.3); names(sds_fill_palette) <- sds_labels
sds_pattern_palette <- c("Definitely Willing" = "none", "Probably Willing" = "stripe", "Probably Unwilling" = "crosshatch", "Definitely Unwilling" = "wave")
data_for_processing <- Data[, unique(c(all_dvs, iv1_col_name, iv2_col_name, education_col_name))]
data_for_processing$Epilepsy_Factor <- factor(data_for_processing[[iv1_col_name]], levels = c(0, 1), labels = c("Epilepsy Absent", "Epilepsy Present"))
data_for_processing$Depression_Factor <- factor(data_for_processing[[iv2_col_name]], levels = c(0, 1), labels = c("Depression Absent", "Depression Present"))
data_for_processing[[education_col_name]] <- as.numeric(data_for_processing[[education_col_name]])
data_for_processing$Condition <- factor(case_when(
  data_for_processing$Epilepsy_Factor == "Epilepsy Absent" & data_for_processing$Depression_Factor == "Depression Absent" ~ "Control",
  data_for_processing$Epilepsy_Factor == "Epilepsy Present" & data_for_processing$Depression_Factor == "Depression Absent" ~ "Epilepsy",
  data_for_processing$Epilepsy_Factor == "Epilepsy Absent" & data_for_processing$Depression_Factor == "Depression Present" ~ "Depression",
  data_for_processing$Epilepsy_Factor == "Epilepsy Present" & data_for_processing$Depression_Factor == "Depression Present" ~ "Combined"),
  levels = c("Control", "Epilepsy", "Depression", "Combined"))
print("--- Data Preparation Complete. ---")
aq_plot_list <- list()
sds_plot_list <- list()

# --- 2. Loop Through All DVs to Generate Individual Plots ---
options(contrasts = c("contr.sum", "contr.poly"))

for (dv_current in all_dvs) {
  print(paste("================ PROCESSING:", dv_current, "================"))
  
  # ========================================================================= #
  # --- PLOT CUSTOMIZATION SETTINGS ---
  bracket_start_y    <- 1.024   # How far above the bar the FIRST bracket starts.
  bracket_step_y     <- 0.045   # Vertical distance BETWEEN brackets. Smaller = tighter.
  plot_upper_limit_y <- 1.183    # Max height of the plot area. 1.4 = 140%.
  save_width         <- 5.5      # Saved image width in inches.
  save_height        <- 6      # Saved image height in inches. Taller for less squished look.
  # ========================================================================= #
  
  data_current_item <- na.omit(data_for_processing[, c(dv_current, "Condition", education_col_name)])
  
  is_aq <- dv_current %in% aq_item_names
  if (is_aq) {
    dv_levels_present <- sort(unique(data_current_item[[dv_current]])); dv_all_levels_for_legend <- aq_levels
    dv_fill_palette <- aq_fill_palette; dv_pattern_palette <- aq_pattern_palette; plot_title <- aq_name_map[dv_current]; legend_title <- "AQ Score"
    data_current_item[[dv_current]] <- factor(data_current_item[[dv_current]], levels = dv_levels_present, ordered = TRUE)
  } else {
    dv_levels_present <- sort(unique(data_current_item[[dv_current]])); dv_all_levels_for_legend <- sds_labels
    dv_labels_current <- sds_labels[dv_levels_present]; dv_fill_palette <- sds_fill_palette; dv_pattern_palette <- sds_pattern_palette
    plot_title <- sds_name_map[dv_current]; legend_title <- "SDS Response"
    data_current_item[[dv_current]] <- factor(data_current_item[[dv_current]], levels = dv_levels_present, labels = dv_labels_current, ordered = TRUE)
  }
  
  if (nrow(data_current_item) < 80) { print(paste("Skipping", dv_current, "due to low N.")); next }
  
  # --- Model Fitting and Data Pivoting ---
  formula_ordinal_str <- paste0("`", dv_current, "` ~ `", education_col_name, "` + Condition")
  model_clm <- clm(as.formula(formula_ordinal_str), data = data_current_item, link = "logit", Hess = TRUE)
  mean_education <- mean(data_current_item[[education_col_name]], na.rm = TRUE)
  pred_grid <- data.frame(Condition = levels(data_current_item$Condition), Q4 = round(mean_education, 0))
  predicted_probs_matrix <- predict(model_clm, newdata = pred_grid, type = "prob")$fit
  colnames(predicted_probs_matrix) <- levels(data_current_item[[dv_current]])
  predicted_probs_df <- cbind(pred_grid, predicted_probs_matrix)
  
  plot_data_long <- predicted_probs_df %>%
    pivot_longer(cols = -c(Condition, Q4), names_to = "Score_Category", values_to = "Probability") %>%
    mutate(
      Score_Category = factor(Score_Category, levels = (if(is_aq) as.character(dv_all_levels_for_legend) else dv_all_levels_for_legend)),
      Condition = factor(Condition, levels = c("Control", "Epilepsy", "Depression", "Combined")),
      # ADD: Create a column for the facet title, which will be the plot title
      facet_title = plot_title
    )
  
  # --- Significance Testing and Manual Bracket Creation ---
  print("... conducting post-hoc tests (NO adjustment) and building brackets manually ...")
  data_current_item$dv_numeric <- as.numeric(data_current_item[[dv_current]])
  stat_test <- dunn_test(data = data_current_item, formula = dv_numeric ~ Condition, p.adjust.method = "none")
  
  stat_test_sig <- stat_test %>%
    mutate(p.signif = case_when(p <= 0.0001 ~ "****", p <= 0.001  ~ "***", p <= 0.01   ~ "**", p <= 0.05   ~ "*", TRUE ~ "ns")) %>%
    filter(p.signif != "ns") %>%
    mutate(
      x = as.numeric(factor(group1, levels = levels(data_current_item$Condition))),
      xend = as.numeric(factor(group2, levels = levels(data_current_item$Condition))),
      label = p.signif
    )
  
  bracket_coords <- data.frame()
  label_coords <- data.frame()
  if(nrow(stat_test_sig) > 0) {
    y_pos_current <- bracket_start_y
    for (i in 1:nrow(stat_test_sig)) {
      row <- stat_test_sig[i, ]
      bracket_coords <- rbind(bracket_coords, data.frame(x = c(row$x, row$x, row$xend, row$xend), y = c(y_pos_current, y_pos_current + 0.01, y_pos_current + 0.01, y_pos_current), bracket_id = i))
      label_coords <- rbind(label_coords, data.frame(x = (row$x + row$xend) / 2, y = y_pos_current + 0.02, label = row$label))
      y_pos_current <- y_pos_current + bracket_step_y
    }
  }
  
  # --- Create the Plot ---
  p <- ggplot(plot_data_long, aes(x = Condition, y = Probability)) +
    geom_bar_pattern(
      stat = "identity", position = "fill", width = 0.825, color = "black",
      aes(fill = Score_Category, pattern = Score_Category),
      pattern_fill = 'black', pattern_angle = 45, pattern_density = 0.1, pattern_spacing = 0.02
    ) +
    # Add manual brackets if they exist
    {
      if(nrow(bracket_coords) > 0)
        list(
          geom_line(data = bracket_coords, aes(x = x, y = y, group = bracket_id), inherit.aes = FALSE),
          geom_text(data = label_coords, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 3.8, face = "bold")
        )
    } +
    # CHANGE: Use facet_wrap to create the title box
    facet_wrap(~facet_title) +
    scale_y_continuous(expand = c(0,0.03), breaks = seq(0, 1, 0.25), labels = scales::percent_format(accuracy = 1)) +
    coord_cartesian(ylim = c(0, plot_upper_limit_y), clip = "off") +
    scale_fill_manual(values = dv_fill_palette, name = legend_title, drop = FALSE) +
    scale_pattern_manual(values = dv_pattern_palette, name = legend_title, drop = FALSE) +
    labs(x = NULL, y = "Proportion of Responses") + # REMOVED: title argument
    theme_linedraw(base_size = 13) +
    theme(
      # ADD: Styling for the facet box (the "box thing")
      strip.background = element_rect(fill = "grey90", color = "black"),
      strip.text = element_text(face = "bold", size = 14, margin = margin(5,0,5,0), color = "black"),
      axis.title.y = element_text(face = "bold", size = 15, margin = margin(r = 10)),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 13, color = "black"),
      axis.text.y = element_text(color = "black", size = 13),
      legend.title = element_text(face= "bold", size = 13),
      legend.text = element_text(size=12),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 5.5, unit = "pt"), # REDUCED: top margin
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  if (is_aq) {
    aq_plot_list[[dv_current]] <- p
  } else {
    sds_plot_list[[dv_current]] <- p
  }
  print(p)
  ggsave(filename=paste0("Plot_", dv_current, "_Final_v3.png"), plot = p, 
         width = save_width, height = save_height, 
         dpi = 300, bg="white")
}

# --- 3. Assemble Individual Plots into Panelled Figures ---
# This section goes AFTER your main "for (dv_current in all_dvs)" loop.
# It assumes your loop has been modified to store the plot objects.

# First, let's slightly modify the loop to store the plots instead of just printing/saving
# (I'll provide the full corrected loop below for clarity)

# --- Full Corrected Script (Your code + patchwork assembly) ---

# --- 0. Load Necessary Libraries ---
# (Your library calls are fine, but ensure patchwork is loaded)
library(patchwork)

# --- 1. Load and Prepare Data ---
# (Your data prep code is fine)
# ... (paste your data prep code here) ...

# --- 2. Loop Through All DVs and STORE Plots in Lists ---
aq_plot_list <- list()
sds_plot_list <- list()

options(contrasts = c("contr.sum", "contr.poly"))

for (dv_current in all_dvs) {
  # ... (all your existing code inside the loop from 'print(paste...)' down to creating the plot 'p') ...
  # The only change is at the very end of the loop:
  # Instead of just print(p), we store the plot 'p' in the correct list
  
  if (is_aq) {
    aq_plot_list[[dv_current]] <- p
  } else {
    sds_plot_list[[dv_current]] <- p
  }
  
  # You can still print it if you want to see them one by one
  print(p) 
}

# --- 3. Assemble and Save the Panelled Figures ---
print("--- Assembling Panelled Figures ---")

# --- Figure 1: AQ Items ---
if (length(aq_plot_list) > 0) {
  # Combine plots. I'm using ncol=4 for a 2-row layout. Adjust as needed.
  # The `& theme(legend.position="none")` removes individual legends.
  panel_aq <- wrap_plots(aq_plot_list, ncol = 4) & theme(legend.position = "none")
  
  # Extract just one legend from the first plot
  legend_aq <- get_legend(aq_plot_list[[1]] + theme(legend.position = "right"))
  
  # Arrange the panel and the single legend
  figure1_aq_panel <- (panel_aq | legend_aq) +
    plot_layout(widths = c(4, 1)) # Give panel 4x the width of the legend
  
  print(figure1_aq_panel)
  ggsave(filename="Figure_AQ_Panel_Final.png", plot = figure1_aq_panel, 
         width = 11, height = 6, dpi = 300, bg="white")
}

# --- Figure 2: SDS Items ---
  if (length(sds_plot_list) > 0) {
  # Combine plots. ncol=3 for a 3x3 grid (with empty spots).
  panel_sds <- wrap_plots(sds_plot_list, ncol = 3) & theme(legend.position = "none")
  
  # Extract a legend
  legend_sds <- get_legend(sds_plot_list[[1]] + theme(legend.position = "right"))
  
  # Arrange the panel and the single legend
  figure2_sds_panel <- (panel_sds | legend_sds) +
  plot_layout(widths = c(3, 1)) 
  
  print(figure2_sds_panel)
  ggsave(filename="Figure_SDS_Panel_Final.png", plot = figure2_sds_panel, 
         width = 8, height = 7, dpi = 300, bg="white")
}