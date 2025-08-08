# ================================================================
# Epilepsy × Depression Stigma Study — Full Analysis Pipeline (R)
# ================================================================
# Outputs: /results/figures/*  and  /results/tables/*
# Requires: data_final.csv with the columns seen in your glimpse()
# ----------------------------------------------------------------

# ---- 0) Packages ------------------------------------------------
pkgs <- c(
  "tidyverse","emmeans","ordinal","sandwich","car","effectsize",
  "broom","patchwork","scales","psych","DiagrammeR","DiagrammeRsvg","rsvg"
)
to_install <- pkgs[!sapply(pkgs, require, character.only = TRUE)]
if (length(to_install)) {
  install.packages(to_install, dependencies = TRUE)
  invisible(lapply(to_install, require, character.only = TRUE))
}
set.seed(123)

# ---- 1) Paths & Data -------------------------------------------
data_file  <- "data_final.csv"
results_dir <- "results"
fig_dir     <- file.path(results_dir, "figures")
tab_dir     <- file.path(results_dir, "tables")

dir.create(results_dir, showWarnings = FALSE)
dir.create(fig_dir,     showWarnings = FALSE)
dir.create(tab_dir,     showWarnings = FALSE)

dat <- readr::read_csv(data_file, col_types = cols())

# ---- 2) Factors & covariate centering --------------------------
stopifnot(all(c("Epilepsy_Condition_Num","Depression_Condition_Num",
                "Education_Numeric_Cov","AQ_Total","SDST_Total") %in% names(dat)))

dat <- dat %>%
  mutate(
    Epilepsy   = factor(if_else(Epilepsy_Condition_Num == 1, "Present", "Absent"),
                        levels = c("Absent","Present")),
    Depression = factor(if_else(Depression_Condition_Num == 1, "Present", "Absent"),
                        levels = c("Absent","Present")),
    Condition  = factor(case_when(
      Epilepsy == "Absent"  & Depression == "Absent"  ~ "Control",
      Epilepsy == "Present" & Depression == "Absent"  ~ "Epilepsy",
      Epilepsy == "Absent"  & Depression == "Present" ~ "Depression",
      TRUE                                              ~ "Combined"
    ), levels = c("Control","Epilepsy","Depression","Combined")),
    Education_c = scale(Education_Numeric_Cov, center = TRUE, scale = FALSE)[,1]
  )

options(contrasts = c("contr.sum","contr.poly"))  # for Type III

# Item names
aq_items  <- c("AQ_Pity","AQ_Danger","AQ_Fear","AQ_Blame","AQ_Anger","AQ_No_Help","AQ_Coercion")
sds_items <- paste0("SDS_", 1:7)
stopifnot(all(aq_items %in% names(dat)), all(sds_items %in% names(dat)))

# ---- 3) Scale reliability (alpha) -------------------------------
aq_alpha  <- psych::alpha(dat[aq_items],  check.keys = FALSE)
sds_alpha <- psych::alpha(dat[sds_items], check.keys = FALSE)

rel_tbl <- tibble::tibble(
  Scale = c("AQ","SDS"),
  Alpha = c(aq_alpha$total$raw_alpha, sds_alpha$total$raw_alpha),
  CI_L  = c(aq_alpha$total$raw_alpha_ci[1], sds_alpha$total$raw_alpha_ci[1]),
  CI_U  = c(aq_alpha$total$raw_alpha_ci[2], sds_alpha$total$raw_alpha_ci[2]),
  n_items = c(length(aq_items), length(sds_items))
)
readr::write_csv(rel_tbl, file.path(tab_dir,"Scale_Reliability.csv"))

# ---- 4) Convert items to ordered for CLMs -----------------------
dat[aq_items]  <- lapply(dat[aq_items],  function(x) factor(x, levels = 1:9, ordered = TRUE))
dat[sds_items] <- lapply(dat[sds_items], function(x) factor(x, levels = 1:4, ordered = TRUE))

# ---- 5) Totals: Robust ANCOVAs (HC3) + Effect sizes -------------
fit_AQ   <- lm(AQ_Total   ~ Epilepsy * Depression + Education_c, data = dat)
fit_SDST <- lm(SDST_Total ~ Epilepsy * Depression + Education_c, data = dat)

# HC3 inference
aq_anova_hc3   <- car::Anova(fit_AQ,   type = 3, white.adjust = "hc3")
sdst_anova_hc3 <- car::Anova(fit_SDST, type = 3, white.adjust = "hc3")
readr::write_csv(broom::tidy(aq_anova_hc3),   file.path(tab_dir,"ANCOVA_AQ_Total_TypeIII_HC3.csv"))
readr::write_csv(broom::tidy(sdst_anova_hc3), file.path(tab_dir,"ANCOVA_SDST_Total_TypeIII_HC3.csv"))

# Standard Type III (for SS to compute partial eta^2)
aq_anova_ss   <- car::Anova(fit_AQ,   type = 3)
sdst_anova_ss <- car::Anova(fit_SDST, type = 3)
aq_eta   <- effectsize::eta_squared(aq_anova_ss,   partial = TRUE)
sdst_eta <- effectsize::eta_squared(sdst_anova_ss, partial = TRUE)
readr::write_csv(as_tibble(aq_eta),   file.path(tab_dir,"EffectSizes_AQ_Total_eta2.csv"))
readr::write_csv(as_tibble(sdst_eta), file.path(tab_dir,"EffectSizes_SDST_Total_eta2.csv"))

# EMMs (Education_c = 0)
emm_AQ   <- emmeans(fit_AQ,   ~ Epilepsy * Depression, at = list(Education_c = 0))
emm_SDST <- emmeans(fit_SDST, ~ Epilepsy * Depression, at = list(Education_c = 0))
readr::write_csv(as.data.frame(emm_AQ),   file.path(tab_dir,"EMMs_AQ_Total.csv"))
readr::write_csv(as.data.frame(emm_SDST), file.path(tab_dir,"EMMs_SDST_Total.csv"))

pairs_AQ   <- as.data.frame(contrast(emm_AQ,   "pairwise", adjust = "tukey"))
pairs_SDST <- as.data.frame(contrast(emm_SDST, "pairwise", adjust = "tukey"))
readr::write_csv(pairs_AQ,   file.path(tab_dir,"Pairs_AQ_Total_Tukey.csv"))
readr::write_csv(pairs_SDST, file.path(tab_dir,"Pairs_SDST_Total_Tukey.csv"))

# ---- 6) CLMs per item with robust formula building --------------
# Helper: identify p-value column in nominal_test output
get_pcol <- function(x) grep("^Pr\\(", colnames(x), value = TRUE)[1]

# Safely fit a CLM for one item:
fit_item_clm <- function(item) {
  # Subset to variables needed and drop unused levels on the response
  df <- dat %>% dplyr::select(all_of(c(item, "Epilepsy","Depression","Education_c")))
  df[[item]] <- droplevels(df[[item]])
  # If the response has <2 levels after dropping, skip
  if (nlevels(df[[item]]) < 2) {
    return(list(model = NULL, skipped = TRUE, reason = "Only one response level after droplevels()",
                violated = NA, violators = character(0)))
  }
  # Build formula robustly
  rhs <- c("Epilepsy*Depression", "Education_c")
  form <- stats::reformulate(termlabels = rhs, response = item)
  # Fit proportional odds model
  m <- tryCatch(ordinal::clm(form, data = df, link = "logit"),
                error = function(e) e)
  if (inherits(m, "error")) {
    return(list(model = NULL, skipped = TRUE, reason = paste("clm fit error:", m$message),
                violated = NA, violators = character(0)))
  }
  # PO test
  po <- tryCatch(ordinal::nominal_test(m), error = function(e) NULL)
  violated <- FALSE; violators <- character(0)
  if (!is.null(po)) {
    pcol <- get_pcol(po)
    if (!is.na(pcol)) {
      violators <- rownames(po)[po[[pcol]] < 0.05]
      violators <- setdiff(violators, "(Intercept)")
      if (length(violators)) {
        violated <- TRUE
        m2 <- tryCatch(
          ordinal::clm(form, data = df, link = "logit",
                       nominal = stats::reformulate(violators)),
          error = function(e) e
        )
        if (!inherits(m2, "error")) m <- m2  # fallback succeeded
      }
    }
  }
  list(model = m, skipped = FALSE, reason = "", violated = violated, violators = violators)
}

all_items <- c(aq_items, sds_items)
item_models <- setNames(vector("list", length = length(all_items)), all_items)
po_flags <- tibble::tibble(Item = character(), Skipped = logical(),
                           Reason = character(), PO_violated = logical(),
                           Violators = character())

for (it in all_items) {
  res <- fit_item_clm(it)
  item_models[[it]] <- res$model
  po_flags <- bind_rows(po_flags,
                        tibble::tibble(Item = it,
                                       Skipped = isTRUE(res$skipped),
                                       Reason = res$reason,
                                       PO_violated = isTRUE(res$violated),
                                       Violators = paste(res$violators, collapse = "; ")))
}
readr::write_csv(po_flags, file.path(tab_dir,"PO_Tests_ItemLevel.csv"))

# ---- 7) Item-level contrasts (E+D vs E; E+D vs D) ---------------
get_two_contrasts <- function(m) {
  # emmeans on latent (logit) scale
  em <- emmeans(m, ~ Epilepsy * Depression, mode = "latent")
  emdf <- as.data.frame(em)
  iC <- with(emdf, which(Epilepsy=="Present" & Depression=="Present"))
  iE <- with(emdf, which(Epilepsy=="Present" & Depression=="Absent"))
  iD <- with(emdf, which(Epilepsy=="Absent"  & Depression=="Present"))
  v_C_vs_E <- rep(0, nrow(emdf)); v_C_vs_E[iC] <- 1; v_C_vs_E[iE] <- -1
  v_C_vs_D <- rep(0, nrow(emdf)); v_C_vs_D[iC] <- 1; v_C_vs_D[iD] <- -1
  cc <- contrast(em, method = list(Combined_vs_Epilepsy = v_C_vs_E,
                                   Combined_vs_Depression = v_C_vs_D))
  out <- summary(cc, infer = c(TRUE, TRUE))
  out$OR      <- exp(out$estimate)
  out$OR_low  <- exp(out$lower.CL)
  out$OR_high <- exp(out$upper.CL)
  out
}

item_contrasts <- purrr::map_df(names(item_models), function(it) {
  m <- item_models[[it]]
  if (is.null(m)) {
    tibble::tibble(contrast = c("Combined_vs_Epilepsy","Combined_vs_Depression"),
                   estimate = NA_real_, SE = NA_real_, df = NA_real_,
                   z.ratio = NA_real_, p.value = NA_real_,
                   lower.CL = NA_real_, upper.CL = NA_real_,
                   OR = NA_real_, OR_low = NA_real_, OR_high = NA_real_,
                   Item = it, Family = ifelse(it %in% aq_items,"AQ","SDS"))
  } else {
    s <- get_two_contrasts(m)
    s$Item <- it
    s$Family <- ifelse(it %in% aq_items, "AQ", "SDS")
    as_tibble(s)
  }
})

# FDR within family × contrast (ignore NAs)
item_contrasts <- item_contrasts %>%
  group_by(Family, contrast) %>%
  mutate(p_FDR = ifelse(is.na(p.value), NA_real_, p.adjust(p.value, method = "BH"))) %>%
  ungroup()

readr::write_csv(item_contrasts, file.path(tab_dir,"ItemLevel_Contrasts_ORs_FDR.csv"))

# ---- 8) Table 1 — Sample characteristics -----------------------
tab1_counts <- dat %>% count(Condition, name = "n")
readr::write_csv(tab1_counts, file.path(tab_dir,"Table1_CellSizes.csv"))

if ("Education_Factor" %in% names(dat)) {
  tab1_edu <- dat %>% count(Condition, Education_Factor) %>%
    group_by(Condition) %>% mutate(pct = round(100*n/sum(n),1)) %>% ungroup()
  readr::write_csv(tab1_edu, file.path(tab_dir,"Table1_EducationByCell.csv"))
} else {
  tab1_edu <- dat %>% 
    group_by(Condition) %>% 
    summarise(n = dplyr::n(),
              mean_Edu = mean(Education_Numeric_Cov, na.rm = TRUE),
              sd_Edu   = sd(Education_Numeric_Cov,   na.rm = TRUE),
              .groups = "drop")
  readr::write_csv(tab1_edu, file.path(tab_dir,"Table1_EducationNumericByCell.csv"))
}

# ---- 9) Figures -------------------------------------------------

# (a) Totals EMM plots
mk_cond_df <- function(emm_df) {
  emm_df %>%
    mutate(Condition = case_when(
      Epilepsy=="Absent" & Depression=="Absent" ~ "Control",
      Epilepsy=="Present" & Depression=="Absent" ~ "Epilepsy",
      Epilepsy=="Absent" & Depression=="Present" ~ "Depression",
      TRUE ~ "Combined"
    )) %>%
    mutate(Condition = factor(Condition, levels = c("Control","Epilepsy","Depression","Combined")))
}

emm_AQ_df   <- as.data.frame(emm_AQ)
emm_SDST_df <- as.data.frame(emm_SDST)

emm_AQ_plot   <- mk_cond_df(emm_AQ_df)
emm_SDST_plot <- mk_cond_df(emm_SDST_df)

p_AQ <- ggplot(emm_AQ_plot, aes(Condition, emmean)) +
  geom_col(width = 0.65, colour = "black", fill = "grey80") +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  labs(y = "Adjusted mean AQ_Total", x = NULL, title = "AQ_Total by condition") +
  theme_bw()

p_SDST <- ggplot(emm_SDST_plot, aes(Condition, emmean)) +
  geom_col(width = 0.65, colour = "black", fill = "grey80") +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  labs(y = "Adjusted mean SDST_Total", x = NULL, title = "SDST_Total by condition") +
  theme_bw()

p_totals <- p_AQ + p_SDST + patchwork::plot_annotation(
  title = "Estimated Marginal Means (95% CI) — Totals adjusted for Education"
)
ggsave(file.path(fig_dir,"Fig1_Totals_EMM.png"), p_totals, width = 10, height = 5, dpi = 300)
ggsave(file.path(fig_dir,"Fig1_Totals_EMM.svg"), p_totals, width = 10, height = 5)

# (b) Forest plot of ORs (skip NA rows)
forest_df <- item_contrasts %>%
  filter(!is.na(OR)) %>%
  select(Item, Family, contrast, OR, OR_low, OR_high, p_FDR) %>%
  mutate(Item = factor(Item, levels = c(aq_items, sds_items)),
         contrast = dplyr::recode(contrast,
                                  Combined_vs_Epilepsy   = "Combined vs Epilepsy",
                                  Combined_vs_Depression = "Combined vs Depression"))

p_forest <- ggplot(forest_df, aes(x = OR, y = Item, shape = contrast)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_point(position = position_dodge(width = 0.6), size = 2) +
  geom_errorbar(aes(xmin = OR_low, xmax = OR_high),
                position = position_dodge(width = 0.6), width = 0.25) +
  scale_x_log10() +
  facet_grid(Family ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Odds ratio (log scale)", y = NULL,
       title = "Comorbid vs single-condition contrasts across items") +
  theme_bw() +
  theme(legend.position = "bottom", legend.title = element_blank())

ggsave(file.path(fig_dir,"Fig2_Forest_ORs.png"), p_forest, width = 7.5, height = 7.5, dpi = 300)
ggsave(file.path(fig_dir,"Fig2_Forest_ORs.svg"), p_forest, width = 7.5, height = 7.5)

# (c) Stacked Likert predicted probabilities per item × condition
edu0  <- 0  # Education_c at mean
newdat <- expand.grid(Epilepsy = c("Absent","Present"),
                      Depression = c("Absent","Present"),
                      Education_c = edu0)

pred_probs_one <- function(m, item) {
  if (is.null(m)) return(NULL)
  pr <- predict(m, newdata = newdat, type = "prob")$fit
  out <- cbind(newdat, as.data.frame(pr))
  out$Item <- item
  out
}
pred_list <- purrr::map(names(item_models), ~ pred_probs_one(item_models[[.x]], .x))
pred_all  <- dplyr::bind_rows(pred_list)

if (!is.null(pred_all)) {
  cat_cols <- setdiff(names(pred_all), c("Epilepsy","Depression","Education_c","Item"))
  prob_long <- pred_all %>%
    tidyr::pivot_longer(cols = all_of(cat_cols), names_to = "Category", values_to = "Prob") %>%
    mutate(Category = readr::parse_number(Category) %>% factor(ordered = TRUE),
           Condition = case_when(
             Epilepsy=="Absent" & Depression=="Absent" ~ "Control",
             Epilepsy=="Present" & Depression=="Absent" ~ "Epilepsy",
             Epilepsy=="Absent" & Depression=="Present" ~ "Depression",
             TRUE ~ "Combined"
           ),
           Scale = if_else(Item %in% aq_items, "AQ", "SDS"))
  
  prob_sds <- prob_long %>% filter(Scale=="SDS")
  prob_aq  <- prob_long %>% filter(Scale=="AQ")
  
  p_sds <- ggplot(prob_sds, aes(x = Condition, y = Prob, fill = Category)) +
    geom_bar(stat = "identity", position = "fill", width = 0.8, colour = "black") +
    facet_wrap(~ Item, ncol = 4) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_grey(start = 0.95, end = 0.25) +
    labs(y = "Proportion of responses", x = NULL,
         title = "Predicted response distributions (SDS)") +
    theme_bw()
  
  p_aq <- ggplot(prob_aq, aes(x = Condition, y = Prob, fill = Category)) +
    geom_bar(stat = "identity", position = "fill", width = 0.8, colour = "black") +
    facet_wrap(~ Item, ncol = 4) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_grey(start = 0.95, end = 0.25) +
    labs(y = "Proportion of responses", x = NULL,
         title = "Predicted response distributions (AQ)") +
    theme_bw()
  
  ggsave(file.path(fig_dir,"Fig3_Stacked_SDS.png"), p_sds, width = 10, height = 6, dpi = 300)
  ggsave(file.path(fig_dir,"Fig4_Stacked_AQ.png"),  p_aq,  width = 12, height = 6, dpi = 300)
  ggsave(file.path(fig_dir,"Fig3_Stacked_SDS.svg"), p_sds, width = 10, height = 6)
  ggsave(file.path(fig_dir,"Fig4_Stacked_AQ.svg"),  p_aq,  width = 12, height = 6)
} else {
  message("No predicted-probability plots generated because all items were skipped (unlikely).")
}

# ---- 10) Analysis flowchart (Mermaid -> PDF) --------------------
fc <- DiagrammeR::mermaid("
graph TD;
  A[Load data & define factors] --> B[Reliability (alpha)];
  B --> C[Totals: robust ANCOVA (HC3, Type III)];
  C --> D[EMMs & Tukey pairwise];
  B --> E[Items: CLM (ordinal logit)];
  E --> F[PO test (nominal_test)];
  F --> G{PO violated?};
  G -- No --> H[Keep proportional odds model];
  G -- Yes --> I[Partial PO (nominal terms)];
  H --> J[Contrasts: Combined vs Epilepsy & vs Depression];
  I --> J;
  D --> K[Fig.1 Totals EMM plot];
  J --> L[Fig.2 Forest of ORs];
  J --> M[Fig.3–4 Stacked predicted probabilities];
  K --> N[Export tables & figures];
  L --> N;
  M --> N;
  N --> O[Done];
", width = "100%")

svg_code <- DiagrammeRsvg::export_svg(fc)
writeLines(svg_code, file.path(fig_dir,"Analysis_Flowchart.svg"))
rsvg::rsvg_pdf(file.path(fig_dir,"Analysis_Flowchart.svg"),
               file.path(fig_dir,"Analysis_Flowchart.pdf"))

# ---------------------- End of Script ----------------------------
