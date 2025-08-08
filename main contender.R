# encoding: UTF-8
################################################################################
# PROJECT TITLE: Investigating Layered Public Stigma (Epilepsy and Depression)
# SCRIPT PURPOSE: Final Statistical Analysis (CLM, LRT, Visualization) - AESTHETIC V4
# DATE: 2025-08-04
################################################################################

#===============================================================================
# PART 0: SETUP AND LIBRARY LOADING
#===============================================================================
start_time <- Sys.time()
print(paste("--- PART 0: SETUP --- Start Time:", start_time))

packages <- c("readr", "dplyr", "tidyr", "ggplot2", "ordinal", "stats",
              "ggpattern", "patchwork", "emmeans", "knitr", "broom", "stringr")

# Load libraries ensuring they are installed if necessary
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

lapply(packages, function(pkg) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
})

options(contrasts = c("contr.sum", "contr.poly"))

#===============================================================================
# PART 1: DATA PREPARATION
#===============================================================================
print("--- PART 1: DATA PREPARATION ---")

# --- 1.1 Load Data ---
load_data <- function() {
  # Determine file path
  if (file.exists("data_for_aq_total_anova.csv")) {
    file_path <- "data_for_aq_total_anova.csv"
  } else if (file.exists("data_final.csv")) {
    file_path <- "data_final.csv"
  } else {
    stop("ERROR: Data file not found.")
  }
  # Load data suppressing messages about column name repairs
  suppressMessages({
    Data <- readr::read_csv(file_path, show_col_types = FALSE)
  })
  print(paste("Loaded", file_path))
  return(Data)
}
Data_raw <- load_data()

# --- 1.2 Define Key Column Names and Aesthetics ---
IV1_NUM_NAME <- "Epilepsy_Condition_Num"
IV2_NUM_NAME <- "Depression_Condition_Num"
COV_NAME <- "Q4"

aq_item_names <- c("AQ_Pity", "AQ_Danger", "AQ_Fear", "AQ_Blame", "AQ_Anger", "AQ_No_Help", "AQ_Coercion")
sds_item_names <- paste0("SDS_", 1:7)
all_dvs <- c(aq_item_names, sds_item_names)

Data_raw[[COV_NAME]] <- as.numeric(Data_raw[[COV_NAME]])

# SDS Aesthetics
sds_name_map <- c("SDS_1" = "Neighbour", "SDS_2" = "Socialise", "SDS_3" = "Friends", "SDS_4" = "Work", "SDS_5" = "Marry", "SDS_6" = "Rent", "SDS_7" = "Hire")
sds_labels <- c("Definitely Willing", "Probably Willing", "Probably Unwilling", "Definitely Unwilling")
sds_fill_palette <- grey.colors(4, start=0.95, end=0.3); names(sds_fill_palette) <- sds_labels
sds_pattern_palette <- c("Definitely Willing" = "none", "Probably Willing" = "stripe", "Probably Unwilling" = "crosshatch", "Definitely Unwilling" = "circle")

# AQ Aesthetics
aq_name_map <- setNames(gsub("AQ_", "", aq_item_names), aq_item_names)
aq_name_map["AQ_No_Help"] <- "No Help"

aq_levels <- as.character(1:9)
aq_fill_palette <- grey.colors(9, start = 0.95, end = 0.3); names(aq_fill_palette) <- aq_levels
aq_pattern_palette <- setNames(rep(c('none', 'stripe', 'crosshatch', 'circle', 'weave'), length.out = 9), aq_levels)


# --- 1.3 Create Factor Variables and Filter ---
Data_analysis <- Data_raw %>%
  mutate(
    Epilepsy_Factor = factor(!!sym(IV1_NUM_NAME), levels = c(0, 1), labels = c("Epilepsy Absent", "Epilepsy Present")),
    Depression_Factor = factor(!!sym(IV2_NUM_NAME), levels = c(0, 1), labels = c("Depression Absent", "Depression Present")),
    Condition = case_when(
      Epilepsy_Factor == "Epilepsy Absent" & Depression_Factor == "Depression Absent" ~ "Control",
      Epilepsy_Factor == "Epilepsy Present" & Depression_Factor == "Depression Absent" ~ "Epilepsy",
      Epilepsy_Factor == "Epilepsy Absent" & Depression_Factor == "Depression Present" ~ "Depression",
      Epilepsy_Factor == "Epilepsy Present" & Depression_Factor == "Depression Present" ~ "Combined"
    )
  ) %>%
  mutate(Condition = factor(Condition, levels = c("Control", "Epilepsy", "Depression", "Combined"))) %>%
  # Apply AQ Blame Filter
  filter(!(Condition == "Control" & !is.na(AQ_Blame) & AQ_Blame > 2))

# --- 1.4 Final Plot Customization ---
# This section contains variables to control the final plot layouts and aesthetics.
# Change these values to alter the appearance of the saved figures.

# Layout Controls
NCOL <- 4
FINAL_PLOT_WIDTH <- 12
FINAL_PLOT_HEIGHT <- 12
SUBPLOT_ASPECT_RATIO <- 1.8
Y_AXIS_LOWER_PADDING <- 0.018 # e.g., 0.05 raises the bars off the x-axis by 5%

# Font and Theme Controls
BASE_FONT_SIZE <- 12
STRIP_TEXT_SIZE <- 11
STRIP_BG_COLOR <- "grey90"
X_AXIS_FONT_SIZE <- 10
LEGEND_TITLE_SIZE <- 11
LEGEND_TEXT_SIZE <- 10

# Bar & Pattern Controls
BAR_WIDTH <- 0.86
BAR_BORDER_COLOR <- "black"
PATTERN_FILL_COLOR <- "black"
PATTERN_DENSITY <- 0.05
PATTERN_SPACING <- 0.025

# Significance Bracket Controls
BRACKET_COLOR <- "black"
BRACKET_LINEWIDTH <- 0.5
BRACKET_LABEL_SIZE <- 4
BRACKET_LABEL_Y_OFFSET <- 0.027 # Controls height of asterisks from the bracket line


#===============================================================================
# PART 2: ANALYSIS LOOP & PRE-CALCULATION OF PLOT HEIGHTS (Two-Pass Approach)
#===============================================================================
print("--- PART 2: STARTING ANALYSIS LOOP (CLM) & PRE-CALCULATING PLOT HEIGHTS ---")

model_results_list <- list(); emm_results_list <- list()

# Define plot height parameters
bracket_start_y <- 1.02; bracket_step_y <- 0.06; base_upper_limit_y <- 1.03
max_plot_height <- base_upper_limit_y # Initialize max height

# --- Loop 1: Perform Analysis and Determine Maximum Plot Height Needed ---
for (dv_current in all_dvs) {
  
  # Data preparation
  data_current <- Data_analysis %>% select(all_of(c(dv_current, "Condition", "Epilepsy_Factor", "Depression_Factor", COV_NAME))) %>% na.omit()
  is_aq <- dv_current %in% aq_item_names
  levels_present_numeric <- sort(unique(data_current[[dv_current]]))
  
  if (is_aq) {
    data_current[[dv_current]] <- factor(data_current[[dv_current]], levels = levels_present_numeric, ordered = TRUE)
  } else {
    labels_present <- sds_labels[levels_present_numeric]
    data_current[[dv_current]] <- factor(data_current[[dv_current]], levels = levels_present_numeric, labels = labels_present, ordered = TRUE)
  }
  
  if (length(levels(data_current[[dv_current]])) < 2) next
  
  # Fit Models
  formula_full <- as.formula(paste0("`", dv_current, "` ~ `", COV_NAME, "` + Epilepsy_Factor * Depression_Factor"))
  model_full <- tryCatch(clm(formula_full, data = data_current, link = "logit", Hess = TRUE), error = function(e) NULL)
  formula_additive <- as.formula(paste0("`", dv_current, "` ~ `", COV_NAME, "` + Epilepsy_Factor + Depression_Factor"))
  model_additive <- tryCatch(clm(formula_additive, data = data_current, link = "logit"), error = function(e) NULL)
  if (is.null(model_full) || is.null(model_additive)) next
  
  # H2 Test (LRT) and PO Check
  lrt_result <- anova(model_additive, model_full)
  po_test_p <- tryCatch({
    nt <- suppressWarnings(ordinal::nominal_test(model_full))
    if (nrow(nt) == 0 || all(is.na(nt[, "Pr(>Chi)"]))) NA else min(nt[, "Pr(>Chi)"], na.rm = TRUE)
  }, error = function(e) NA)
  model_results_list[[dv_current]] <- data.frame(DV = dv_current, LRT_ChiSq = lrt_result$LR.stat[2], LRT_df = lrt_result$df[2], LRT_p = lrt_result$`Pr(>Chisq)`[2], PO_Assumption_p = po_test_p)
  
  # Post-Hoc
  formula_condition <- as.formula(paste0("`", dv_current, "` ~ `", COV_NAME, "` + Condition"))
  model_condition <- tryCatch(clm(formula_condition, data = data_current, link = "logit"), error = function(e) NULL)
  if (is.null(model_condition)) next
  mean_covariate <- mean(data_current[[COV_NAME]], na.rm = TRUE)
  emm <- emmeans(model_condition, specs = pairwise ~ Condition, adjust = "tukey", at = setNames(list(mean_covariate), COV_NAME))
  emm_summary <- summary(emm$contrasts)
  emm_results_list[[dv_current]] <- emm_summary
  
  # Determine required height
  contrasts_df <- as.data.frame(emm_summary) %>% filter(p.value <= 0.05)
  num_sig_contrasts <- nrow(contrasts_df)
  if (num_sig_contrasts > 0) {
    required_height <- bracket_start_y + (num_sig_contrasts * bracket_step_y)
    if (required_height > max_plot_height) {
      max_plot_height <- required_height
    }
  }
}

print(paste("Pre-calculation complete. Max plot height set to:", round(max_plot_height, 2)))

#===============================================================================
# PART 3: PLOTTING LOOP (Using Uniform Heights)
#===============================================================================
print("--- PART 3: GENERATING PLOTS WITH UNIFORM HEIGHTS ---")

plot_list_aq <- list(); plot_list_sds <- list()

# --- Loop 2: Generate Plots using the calculated Max Height ---
for (dv_current in all_dvs) {
  if (is.null(model_results_list[[dv_current]])) next
  print(paste("=== Plotting:", dv_current, "==="))
  
  # Data preparation
  data_current <- Data_analysis %>% select(all_of(c(dv_current, "Condition", "Epilepsy_Factor", "Depression_Factor", COV_NAME))) %>% na.omit()
  is_aq <- dv_current %in% aq_item_names
  levels_present_numeric <- sort(unique(data_current[[dv_current]]))
  if (is_aq) {
    all_possible_levels <- aq_levels; plot_title <- aq_name_map[dv_current]; legend_title <- "AQ Score"
    dv_fill_palette <- aq_fill_palette; dv_pattern_palette <- aq_pattern_palette
    data_current[[dv_current]] <- factor(data_current[[dv_current]], levels = levels_present_numeric, ordered = TRUE)
  } else {
    all_possible_levels <- sds_labels; plot_title <- sds_name_map[dv_current]; legend_title <- "SDS Response"
    dv_fill_palette <- sds_fill_palette; dv_pattern_palette <- sds_pattern_palette
    labels_present <- sds_labels[levels_present_numeric]
    data_current[[dv_current]] <- factor(data_current[[dv_current]], levels = levels_present_numeric, labels = labels_present, ordered = TRUE)
  }
  
  # Generate Predicted Probabilities
  formula_full <- as.formula(paste0("`", dv_current, "` ~ `", COV_NAME, "` + Epilepsy_Factor * Depression_Factor"))
  model_full <- clm(formula_full, data = data_current, link = "logit", Hess = TRUE)
  mean_covariate <- mean(data_current[[COV_NAME]], na.rm = TRUE)
  pred_grid <- expand.grid(Epilepsy_Factor = levels(data_current$Epilepsy_Factor), Depression_Factor = levels(data_current$Depression_Factor))
  pred_grid[[COV_NAME]] <- mean_covariate
  predicted_probs <- predict(model_full, newdata = pred_grid, type = "prob")$fit
  if (!is.data.frame(predicted_probs)) predicted_probs <- as.data.frame(predicted_probs)
  plot_data_long <- cbind(pred_grid, predicted_probs) %>%
    mutate(Condition = interaction(Epilepsy_Factor, Depression_Factor, sep=" / ")) %>%
    mutate(Condition = factor(Condition, levels = c("Epilepsy Absent / Depression Absent", "Epilepsy Present / Depression Absent", "Epilepsy Absent / Depression Present", "Epilepsy Present / Depression Present"), labels = c("Control", "Epilepsy", "Depression", "Combined"))) %>%
    pivot_longer(cols = -c(Epilepsy_Factor, Depression_Factor, !!sym(COV_NAME), Condition), names_to = "Score_Category", values_to = "Probability") %>%
    mutate(Score_Category = factor(Score_Category, levels = all_possible_levels), facet_title = plot_title)
  
  # Prepare Brackets
  emm_contrasts <- emm_results_list[[dv_current]]
  contrasts_df <- as.data.frame(emm_contrasts) %>%
    mutate(p.signif = case_when(p.value <= 0.0001 ~ "****", p.value <= 0.001 ~ "***", p.value <= 0.01 ~ "**", p.value <= 0.05 ~ "*", TRUE ~ "ns")) %>%
    filter(p.signif != "ns") %>%
    mutate(contrast = gsub(" ", "", contrast)) %>%
    separate(contrast, into = c("group1", "group2"), sep = "-")
  bracket_coords <- data.frame(); label_coords <- data.frame()
  if(nrow(contrasts_df) > 0) {
    y_pos_current <- bracket_start_y
    levels_map <- setNames(1:4, levels(plot_data_long$Condition))
    contrasts_df$x_start <- levels_map[contrasts_df$group1]; contrasts_df$x_end <- levels_map[contrasts_df$group2]
    if(!(any(is.na(contrasts_df$x_start)) || any(is.na(contrasts_df$x_end)))) {
      contrasts_df$dist <- abs(contrasts_df$x_start - contrasts_df$x_end)
      contrasts_df <- contrasts_df[order(contrasts_df$dist), ]
      for (i in 1:nrow(contrasts_df)) {
        row <- contrasts_df[i, ]; x_start <- row$x_start; x_end <- row$x_end
        bracket_coords <- rbind(bracket_coords, data.frame(x = c(x_start, x_start, x_end, x_end), y = c(y_pos_current, y_pos_current + 0.015, y_pos_current + 0.015, y_pos_current), bracket_id = i))
        label_coords <- rbind(label_coords, data.frame(x = (x_start + x_end) / 2, y = y_pos_current + BRACKET_LABEL_Y_OFFSET, label = row$p.signif))
        y_pos_current <- y_pos_current + bracket_step_y
      }
    }
  }
  
  # Generate the ggplot object using customization variables
  p <- ggplot(plot_data_long, aes(x = Condition, y = Probability)) +
    geom_bar_pattern(
      stat = "identity", position = "fill", width = BAR_WIDTH, color = BAR_BORDER_COLOR,
      aes(fill = Score_Category, pattern = Score_Category),
      pattern_fill = PATTERN_FILL_COLOR, pattern_angle = 45,
      pattern_density = PATTERN_DENSITY, pattern_spacing = PATTERN_SPACING
    ) +
    {
      if(nrow(bracket_coords) > 0)
        list(
          geom_line(data = bracket_coords, aes(x = x, y = y, group = bracket_id), inherit.aes = FALSE, color = BRACKET_COLOR, linewidth = BRACKET_LINEWIDTH),
          geom_text(data = label_coords, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = BRACKET_LABEL_SIZE, fontface = "bold")
        )
    } +
    facet_wrap(~facet_title) +
    scale_y_continuous(expand = c(Y_AXIS_LOWER_PADDING, 0.01), breaks = seq(0, 1, 0.25), labels = scales::percent_format(accuracy = 1)) +
    coord_cartesian(ylim = c(0, max_plot_height), clip = "off") +
    scale_fill_manual(values = dv_fill_palette, name = legend_title, drop = FALSE) +
    scale_pattern_manual(values = dv_pattern_palette, name = legend_title, drop = FALSE) +
    labs(x = NULL, y = NULL) +
    theme_linedraw(base_size = BASE_FONT_SIZE) +
    theme(
      aspect.ratio = SUBPLOT_ASPECT_RATIO,
      strip.background = element_rect(fill = STRIP_BG_COLOR, color = "black"),
      strip.text = element_text(face = "bold", size = STRIP_TEXT_SIZE, margin = margin(4,0,4,0), color = "black"),
      axis.text.y = element_blank(),
      axis.ticks.y = element_line(color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = X_AXIS_FONT_SIZE, color = "black"),
      legend.title = element_text(face= "bold", size = LEGEND_TITLE_SIZE),
      legend.text = element_text(size = LEGEND_TEXT_SIZE),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      plot.margin = margin(t = 15, r = 5, b = 5, l = 5, unit = "pt")
    )
  
  if (is_aq) { plot_list_aq[[dv_current]] <- p } else { plot_list_sds[[dv_current]] <- p }
}

#===============================================================================
# PART 4: FINAL OUTPUT COMPILATION
#===============================================================================
print("--- PART 4: COMPILING OUTPUTS ---")

# Formatting functions
format_p <- function(p) { if (is.na(p)) return("NA"); if (p < 0.001) return("<.001"); return(sprintf("%.3f", p)) }

# Helper function to extract legend
get_legend <- function(a.gplot){ tmp <- ggplot_gtable(ggplot_build(a.gplot)); leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box"); if (length(leg) == 0) return(NULL); legend <- tmp$grobs[[leg]]; return(legend) }

# Helper function to add Y-axis elements
add_y_axis_text <- function(plot) { plot + theme(axis.text.y = element_text(color = "black", size = 10)) }

# --- 4.1 & 4.2: Summary Tables ---
print("--- Table 1: Summary of Interaction Tests (LRT) and PO Assumption Checks ---")
if (length(model_results_list) > 0) {
  summary_table <- bind_rows(model_results_list)
  summary_table_formatted <- summary_table %>%
    mutate(LRT_ChiSq = sprintf("%.2f", LRT_ChiSq), LRT_p_formatted = sapply(LRT_p, format_p), PO_Assumption_p_formatted = sapply(PO_Assumption_p, format_p), Sig = ifelse(LRT_p < 0.05, "*", "")) %>%
    select(DV, LRT_ChiSq, LRT_df, LRT_p_formatted, Sig, PO_Assumption_p_formatted)
  print(knitr::kable(summary_table_formatted, col.names = c("DV", "LRT ChiSq (H2)", "df", "p-value (H2)", "Sig.", "p-value (PO Check)")))
}
print("\n--- Table 2: Test of H1 (Combined vs. Epilepsy Only) - Tukey Adjusted ---")
h1_results <- lapply(names(emm_results_list), function(dv) {
  res <- emm_results_list[[dv]]; h1_contrast <- res %>% filter(grepl("Epilepsy - Combined", contrast) | grepl("Combined - Epilepsy", contrast)); if (nrow(h1_contrast) == 0) return(NULL)
  estimate <- h1_contrast$estimate; z_ratio <- h1_contrast$z.ratio; if (grepl("Epilepsy - Combined", h1_contrast$contrast)) { estimate <- -estimate; z_ratio <- -z_ratio }
  data.frame(DV = dv, Estimate_LogOdds = estimate, SE = h1_contrast$SE, z.ratio = z_ratio, p.value = h1_contrast$p.value)
})
h1_table <- bind_rows(h1_results)
if (!is.null(h1_table) && nrow(h1_table) > 0) {
  h1_table_formatted <- h1_table %>% mutate(across(where(is.numeric), ~round(., 3)), p.value = sapply(p.value, format_p))
  print(knitr::kable(h1_table_formatted, caption = "H1: Stigma differences (Log-Odds Scale). Positive Estimate means higher stigma in Combined group."))
}

# --- 4.3 Assemble Panelled Figures (Aesthetic Adjustments Applied) ---
print("--- Assembling and Saving Panelled Figures (This may take time) ---")

# Function to assemble a final figure
assemble_figure <- function(plot_list, file_name, shared_y_title) {
  if (length(plot_list) == 0) return(NULL)
  
  # Extract legend from the first plot
  legend_grob <- get_legend(plot_list[[1]] + theme(legend.position = "right"))
  
  # Remove legend from all plots in the list
  plot_list <- lapply(plot_list, function(p) p + theme(legend.position = "none"))
  
  # Add Y-axis text to the leftmost plots
  for (i in seq_along(plot_list)) { 
    if ((i - 1) %% NCOL == 0) { 
      plot_list[[i]] <- add_y_axis_text(plot_list[[i]]) 
    } 
  }
  
  # Create a list of all grobs to be arranged, including the legend
  plots_for_layout <- c(plot_list)
  if (!is.null(legend_grob)) {
    plots_for_layout[[length(plots_for_layout) + 1]] <- legend_grob
  }
  
  # Use wrap_plots to arrange everything in a grid.
  main_panel <- wrap_plots(plots_for_layout, ncol = NCOL)
  
  # Create the shared y-axis label grob
  y_label_grob <- ggplot() + annotate(geom = "text", x = 1, y = 1, label = shared_y_title, angle = 90, fontface = "bold", size = 4.5) + theme_void()
  
  # Combine with the shared y-axis label
  final_figure <- (y_label_grob | main_panel) + 
    plot_layout(widths = c(0.5, 30))
  
  ggsave(filename = file_name, plot = final_figure, width = FINAL_PLOT_WIDTH, height = FINAL_PLOT_HEIGHT, dpi = 300, bg = "white")
}


# Create and Save Figures
print("Rendering Figure 1 (AQ Panel)...")
assemble_figure(plot_list_aq, file_name = "Figure_1_AQ_Panel_CLM_Final_V4.png", shared_y_title = "Proportion of Responses")

print("Rendering Figure 2 (SDS Panel)...")
assemble_figure(plot_list_sds, file_name = "Figure_2_SDS_Panel_CLM_Final_V4.png", shared_y_title = "Proportion of Responses")

#===============================================================================
# PART 5: VISUALIZE SIGNIFICANT INTERACTIONS
#===============================================================================
print("--- PART 5: GENERATING INTERACTION PLOTS ---")

# Identify DVs with significant interactions from the LRT results
significant_interactions_df <- bind_rows(model_results_list) %>% 
  filter(LRT_p < 0.05)

interaction_plot_list <- list()

if(nrow(significant_interactions_df) > 0) {
  # Create a single dataframe with all interaction data
  all_interaction_data <- do.call(rbind, lapply(significant_interactions_df$DV, function(dv_current) {
    print(paste("--- Calculating Interaction EMMs for:", dv_current, "---"))
    
    data_current <- Data_analysis %>% 
      select(all_of(c(dv_current, "Epilepsy_Factor", "Depression_Factor", COV_NAME))) %>% 
      na.omit()
    
    is_aq <- dv_current %in% aq_item_names
    levels_present_numeric <- sort(unique(data_current[[dv_current]]))
    if (is_aq) {
      data_current[[dv_current]] <- factor(data_current[[dv_current]], levels = levels_present_numeric, ordered = TRUE)
    } else {
      labels_present <- sds_labels[levels_present_numeric]
      data_current[[dv_current]] <- factor(data_current[[dv_current]], levels = levels_present_numeric, labels = labels_present, ordered = TRUE)
    }
    
    formula_full <- as.formula(paste0("`", dv_current, "` ~ `", COV_NAME, "` + Epilepsy_Factor * Depression_Factor"))
    model_full <- clm(formula_full, data = data_current, link = "logit", Hess = TRUE)
    
    mean_covariate <- mean(data_current[[COV_NAME]], na.rm = TRUE)
    interaction_emmeans <- emmeans(model_full, ~ Epilepsy_Factor * Depression_Factor, at = setNames(list(mean_covariate), COV_NAME))
    
    interaction_df <- as.data.frame(interaction_emmeans)
    
    plot_title <- ifelse(is_aq, aq_name_map[dv_current], sds_name_map[dv_current])
    interaction_df$plot_title <- plot_title
    
    return(interaction_df)
  }))
  
  # Ensure the facet titles are ordered correctly
  all_interaction_data$plot_title <- factor(all_interaction_data$plot_title, levels = unique(all_interaction_data$plot_title))
  
  # Create the combined plot
  p_interaction_panel <- ggplot(all_interaction_data, aes(x = Epilepsy_Factor, y = emmean, group = Depression_Factor, color = Depression_Factor)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3.5) +
    geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.1, linewidth = 0.8) +
    facet_wrap(~plot_title, ncol = NCOL) +
    labs(
      x = NULL,
      y = NULL,
      color = "Depression Condition"
    ) +
    theme_linedraw(base_size = BASE_FONT_SIZE) +
    theme(
      aspect.ratio = SUBPLOT_ASPECT_RATIO,
      strip.background = element_rect(fill = STRIP_BG_COLOR, color = "black"),
      strip.text = element_text(face = "bold", size = STRIP_TEXT_SIZE, margin = margin(4,0,4,0), color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = X_AXIS_FONT_SIZE, color = "black"),
      legend.title = element_text(face= "bold", size = LEGEND_TITLE_SIZE),
      legend.text = element_text(size = LEGEND_TEXT_SIZE),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      plot.margin = margin(t = 15, r = 5, b = 5, l = 5, unit = "pt"),
      legend.position = "bottom"
    )
  
  # Assemble and save the final figure
  print("--- Assembling and Saving Interaction Plot Panel ---")
  
  # Extract the legend and the panel
  legend_grob <- get_legend(p_interaction_panel)
  panel_grob <- p_interaction_panel + theme(legend.position = "none")
  
  # Create the shared y-axis label grob
  y_label_grob <- ggplot() + annotate(geom = "text", x = 1, y = 1, label = "Estimated Marginal Mean (Log-Odds Scale)", angle = 90, fontface = "bold", size = 4.5) + theme_void()
  
  # Combine panel and y-axis
  labeled_panel <- (y_label_grob | panel_grob) + plot_layout(widths = c(0.5, 30))
  
  # Add the legend at the bottom
  final_interaction_figure <- labeled_panel / legend_grob + plot_layout(heights = c(10, 1))
  
  ggsave(filename = "Figure_3_Interaction_Plots.png", plot = final_interaction_figure, 
         width = FINAL_PLOT_WIDTH, height = FINAL_PLOT_HEIGHT, dpi = 300, bg = "white")
  
} else {
  print("--- No significant interactions found to plot. ---")
}

#===============================================================================
# PART 6: SUMMARIZE INTERACTIONS WITH A FOREST PLOT
#===============================================================================
print("--- PART 6: GENERATING INTERACTION FOREST PLOT ---")

if(nrow(significant_interactions_df) > 0) {
  
  interaction_effects_df <- do.call(rbind, lapply(significant_interactions_df$DV, function(dv_current) {
    print(paste("--- Extracting Interaction Coefficient for:", dv_current, "---"))
    
    data_current <- Data_analysis %>% 
      select(all_of(c(dv_current, "Epilepsy_Factor", "Depression_Factor", COV_NAME))) %>% 
      na.omit()
    
    is_aq <- dv_current %in% aq_item_names
    levels_present_numeric <- sort(unique(data_current[[dv_current]]))
    if (is_aq) {
      data_current[[dv_current]] <- factor(data_current[[dv_current]], levels = levels_present_numeric, ordered = TRUE)
    } else {
      labels_present <- sds_labels[levels_present_numeric]
      data_current[[dv_current]] <- factor(data_current[[dv_current]], levels = levels_present_numeric, labels = labels_present, ordered = TRUE)
    }
    
    formula_full <- as.formula(paste0("`", dv_current, "` ~ `", COV_NAME, "` + Epilepsy_Factor * Depression_Factor"))
    model_full <- clm(formula_full, data = data_current, link = "logit", Hess = TRUE)
    
    # Use broom to get a tidy summary with confidence intervals
    tidy_model <- broom::tidy(model_full, conf.int = TRUE)
    
    # Find the interaction term
    interaction_term <- tidy_model %>% 
      filter(grepl("Epilepsy_Factor.*:Depression_Factor.*", term))
    
    # Get the clean name for the DV
    dv_clean_name <- ifelse(is_aq, aq_name_map[dv_current], sds_name_map[dv_current])
    interaction_term$DV <- dv_clean_name
    
    return(interaction_term)
  }))
  
  # Create the forest plot
  forest_plot <- ggplot(interaction_effects_df, aes(x = estimate, y = reorder(DV, estimate))) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
    geom_point(size = 3) +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
    labs(
      title = "Summary of Interaction Effects",
      x = "Interaction Effect (Log-Odds Ratio)",
      y = "Dependent Variable"
    ) +
    theme_linedraw(base_size = BASE_FONT_SIZE) +
    theme(
      panel.grid.major.y = element_line(color = "grey90", linetype = "dotted")
    )
  
  print("--- Saving Interaction Forest Plot ---")
  ggsave("Figure_4_Interaction_Forest_Plot.png", plot = forest_plot, width = 8, height = 6, dpi = 300, bg = "white")
  
} else {
  print("--- No significant interactions to summarize in a forest plot. ---")
}


end_time <- Sys.time()
print(paste("--- ANALYSIS SCRIPT COMPLETE --- Total Runtime:", round(difftime(end_time, start_time, units = "mins"), 2), "minutes"))


# Create an empty list to store the EMMs for each DV
all_emmeans_list <- list()

# Loop through all 14 dependent variables
for (dv_current in all_dvs) {
  # Skip if the analysis failed for this DV in the main script
  if (is.null(model_results_list[[dv_current]])) next
  
  print(paste("Calculating EMMs for:", dv_current))
  
  # Re-run the data prep for the current DV
  data_current <- Data_analysis %>%
    select(all_of(c(dv_current, "Condition", COV_NAME))) %>%
    na.omit()
  
  is_aq <- dv_current %in% aq_item_names
  levels_present_numeric <- sort(unique(data_current[[dv_current]]))
  if (is_aq) {
    data_current[[dv_current]] <- factor(data_current[[dv_current]], levels = levels_present_numeric, ordered = TRUE)
  } else {
    labels_present <- sds_labels[levels_present_numeric]
    data_current[[dv_current]] <- factor(data_current[[dv_current]], levels = levels_present_numeric, labels = labels_present, ordered = TRUE)
  }
  
  # Re-fit the 'Condition' model
  formula_condition <- as.formula(paste0("`", dv_current, "` ~ `", COV_NAME, "` + Condition"))
  model_condition <- clm(formula_condition, data = data_current, link = "logit")
  
  # Calculate the estimated marginal means
  mean_covariate <- mean(data_current[[COV_NAME]], na.rm = TRUE)
  emm <- emmeans(model_condition, specs = ~ Condition, at = setNames(list(mean_covariate), COV_NAME))
  
  # Store the summary in the list
  emm_summary_df <- as.data.frame(summary(emm))
  emm_summary_df$DV <- dv_current
  all_emmeans_list[[dv_current]] <- emm_summary_df
}

# Combine all the EMMs into a single data frame
all_emmeans_df <- bind_rows(all_emmeans_list)

# Save the combined data frame to a CSV file
write.csv(all_emmeans_df, "all_condition_means.csv", row.names = FALSE)

print("Successfully saved all_condition_means.csv")


# This code tests if the 'Combined' group has significantly more stigma than the 'Epilepsy' group

# Create an empty list to store the results
h1_test_list <- list()

# Loop through all 14 dependent variables
for (dv_current in all_dvs) {
  # Skip if the analysis failed for this DV in the main script
  if (is.null(model_results_list[[dv_current]])) next
  
  print(paste("Testing H1 for:", dv_current))
  
  # Re-run the data prep for the current DV
  data_current <- Data_analysis %>%
    select(all_of(c(dv_current, "Condition", COV_NAME))) %>%
    na.omit()
  
  is_aq <- dv_current %in% aq_item_names
  levels_present_numeric <- sort(unique(data_current[[dv_current]]))
  if (is_aq) {
    data_current[[dv_current]] <- factor(data_current[[dv_current]], levels = levels_present_numeric, ordered = TRUE)
  } else {
    labels_present <- sds_labels[levels_present_numeric]
    data_current[[dv_current]] <- factor(data_current[[dv_current]], levels = levels_present_numeric, labels = labels_present, ordered = TRUE)
  }
  
  # Re-fit the 'Condition' model
  formula_condition <- as.formula(paste0("`", dv_current, "` ~ `", COV_NAME, "` + Condition"))
  model_condition <- clm(formula_condition, data = data_current, link = "logit")
  
  # Calculate the estimated marginal means
  mean_covariate <- mean(data_current[[COV_NAME]], na.rm = TRUE)
  emm <- emmeans(model_condition, specs = pairwise ~ Condition, at = setNames(list(mean_covariate), COV_NAME))
  
  # Extract the specific contrast: Epilepsy vs. Combined
  h1_contrast <- as.data.frame(summary(emm$contrasts)) %>%
    filter(contrast == "Epilepsy - Combined")
  
  # Store the result
  if(nrow(h1_contrast) > 0){
    h1_test_list[[dv_current]] <- data.frame(
      DV = dv_current,
      Estimate = h1_contrast$estimate,
      SE = h1_contrast$SE,
      z.ratio = h1_contrast$z.ratio,
      p.value = h1_contrast$p.value
    )
  }
}

# Combine all results into a single data frame and format it
h1_results_table <- bind_rows(h1_test_list)
h1_results_table$p.value_formatted <- sapply(h1_results_table$p.value, format_p)
h1_results_table$Significant <- ifelse(h1_results_table$p.value < 0.05, "Yes", "No")


# Print the final table
print(knitr::kable(h1_results_table, caption = "Hypothesis 1 Test: Is 'Combined' stigma greater than 'Epilepsy' stigma?"))








install.packages("VGAM")
library(ordinal)   # for clm, scale_test, nominal_test
library(VGAM)      # for partial / generalized ordered logit

## 1. Inspect category counts
dat<-Data_analysis

table(dat$AQ_Fear, dat$Condition)
table(dat$AQ_Anger, dat$Condition)

## 2. Re-fit each item separately with full Hessian
m_fear  <- clm(AQ_Fear  ~ Condition * Facet, data = dat,
               link = "logit", Hess = TRUE, nAGQ = 7)

m_anger <- clm(AQ_Anger ~ Condition * Facet, data = dat,
               link = "logit", Hess = TRUE, nAGQ = 7)

## 3. Run PO score tests
scale_test(m_fear)     # tests parallel (slope) assumption
nominal_test(m_fear)   # tests equal intercepts assumption

scale_test(m_anger)
nominal_test(m_anger)

## 4. If the tests still return NA…
##    A) Collapse sparse categories (example: 5→4 levels)
dat$AQ_Fear4  <- collapseCategories(dat$AQ_Fear, 4)   # helper you write
dat$AQ_Anger4 <- collapseCategories(dat$AQ_Anger, 4)

##    …then refit and re-test.

##    B) Or switch to a partial proportional-odds model
m_fear_pp <- clm(AQ_Fear ~ Condition * Facet,
                 nominal = ~ Condition,      # lets Condition vary by threshold
                 data = dat, Hess = TRUE)

##    C) Or use a generalized ordered logit (no PO constraint at all)
m_fear_gologit <- vglm(AQ_Fear ~ Condition * Facet,
                       family = cumulative(parallel = FALSE),
                       data = dat)

## 5. Compare fit statistics
anova(m_fear, m_fear_pp)           # LR test for PO vs PPO
AIC(m_fear, m_fear_pp, m_fear_gologit)


library(ordinal)   # clm / nominal_test
library(VGAM)      # vglm for gologit
library(dplyr)
library(tidyr)

# ---------- helper to run one item ----------
fit_and_test_PO <- function(df, dv, cov = "Q4",
                            epi = "Epilepsy_Factor",
                            dep = "Depression_Factor",
                            alpha = 0.05,
                            collapse = TRUE,    # collapse sparse tail categories?
                            min_cell = 3)       # minimum count per (DV × Condition)
{
  dat <- df %>% select(all_of(c(dv, cov, epi, dep))) %>% na.omit()
  
  # Optional: merge highest categories until every cell ≥ min_cell
  if (collapse) {
    repeat {
      crosstab <- table(dat[[dv]], interaction(dat[[epi]], dat[[dep]]))
      too_small <- which(rowSums(crosstab < min_cell) > 0)
      if (length(too_small) == 0 || length(levels(dat[[dv]])) <= 3) break
      hi <- max(too_small)                          # highest sparse category
      levels(dat[[dv]])[hi-1] <- paste(levels(dat[[dv]])[hi-1],
                                       levels(dat[[dv]])[hi], sep = "/")
      dat[[dv]] <- factor(as.character(dat[[dv]])) # drop unused level
    }
  }
  
  # Make sure DV is ordered factor
  dat[[dv]] <- factor(dat[[dv]], ordered = TRUE)
  
  # Fit full interaction model
  fmla <- as.formula(sprintf("%s ~ %s + %s * %s",
                             dv, cov, epi, dep))
  mod  <- tryCatch(clm(fmla, data = dat, link = "logit", Hess = TRUE),
                   error = function(e) NULL)
  if (is.null(mod)) return(list(dv = dv, status = "clm failed"))
  
  # Score test of proportional odds
  po_p <- tryCatch({
    nt <- suppressWarnings(ordinal::nominal_test(mod))
    if (nrow(nt) == 0 || all(is.na(nt[,"Pr(>Chi)"]))) NA else
      min(nt[,"Pr(>Chi)"], na.rm = TRUE)
  }, error = function(e) NA)
  
  out <- list(dv = dv, PO_p = po_p)
  
  if (is.na(po_p) || po_p < alpha) {
    # --> fit partial or full non-parallel model
    mod_pp <- tryCatch(
      clm(fmla, nominal = ~ 1, data = dat, link = "logit", Hess = TRUE),
      error = function(e) NULL
    )
    if (!is.null(mod_pp)) {
      aic <- AIC(mod); aic_pp <- AIC(mod_pp)
      out$PPO_fit <- paste0("PPO OK; ΔAIC = ", round(aic - aic_pp, 1))
      out$final_model <- mod_pp
    } else {
      # fallback to fully generalised model (no parallel constraint)
      mod_go <- vglm(fmla,
                     family = cumulative(link = "logit",
                                         parallel = FALSE, reverse = FALSE),
                     data = dat)
      out$PPO_fit <- "PPO failed; used gologit"
      out$final_model <- mod_go
    }
  } else {
    out$PPO_fit <- "PO satisfied"
    out$final_model <- mod
  }
  out
}

# ---------- run it for AQ_Fear and AQ_Anger ----------
fixes <- lapply(c("AQ_Fear", "AQ_Anger"), fit_and_test_PO, df = Data_analysis)

fixes