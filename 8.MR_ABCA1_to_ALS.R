library(TwoSampleMR)
library(MRPRESSO)
library(tidyverse)
library(data.table)

genes <- readRDS("./ML_features_select/common_genes.rds")
load("Data_preprocessing/disc_cohort.rda")
setwd("MR_analysis")

# Map gene SYMBOLs to ENSEMBL IDs
genes <- disc_cohort_mapping[disc_cohort_mapping$SYMBOL %in% genes,]$ENSEMBL
genes

# OpenGWAS API token (if set)
Sys.getenv('OPENGWAS_JWT')

# Initialize result storage
Allres <- data.frame()
j <- 1
out_id <- 'ebi-a-GCST90027163'

# Main loop over genes
for (i in genes) {
  print(paste("Processing gene", i, "(", j, "/", length(genes), ")"))
  
  tryCatch({
    # 1) Extract exposure instruments from OpenGWAS (cis-eQTL per gene)
    exp <- extract_instruments(
      outcomes = paste0('eqtl-a-', i),
      p1 = 5e-8,
      clump = TRUE,
      p2 = 5e-8,
      r2 = 0.001,
      kb = 10000
    )
    
    if (is.null(exp) || nrow(exp) == 0) {
      message("No valid instruments for gene ", i)
      next
    }
    
    # 2) Extract outcome data for the selected outcome
    out <- extract_outcome_data(snps = exp$SNP, outcomes = out_id)
    if (is.null(out) || nrow(out) == 0) {
      message("No outcome data for gene ", i)
      next
    }
    
    # 3) Harmonise exposure and outcome
    dat <- harmonise_data(exposure_dat = exp, outcome_dat = out)
    dat <- subset(dat, mr_keep == TRUE)
    if (nrow(dat) == 0) {
      message("No harmonised SNPs for gene ", i)
      next
    }
    
    # 4) Primary MR analysis (GagnonMR wrapper)
    res <- GagnonMR::primary_MR_analysis(dat = dat)
    
    # Identify p-value column name dynamically
    pval_col <- if ("pval" %in% names(res)) "pval" else "p"
    if (!pval_col %in% names(res)) {
      stop("No p-value column found in results")
    }
    res_pval <- res[[pval_col]]
    
    # Keep only significant results
    if (is.na(res_pval) || res_pval > 0.05) {
      message("Non-significant result (p = ", res_pval, ") for gene ", i)
      next
    }
    
    # 5) Save per-gene objects
    save(exp, out, dat, res, file = paste0('exp.out.dat.res', i, "_mrResult.rdata"))
    
    # 6) Aggregate into global results
    Allres <- bind_rows(Allres, res)
    print(paste("Success:", i))
    
  }, error = function(e) {
    message("Error in gene ", i, ": ", e$message)
  })
  
  j <- j + 1
}


# ───────────────────────────────────────────────────────────────────────────────
# ABCA1 → ALS (ebi-a-GCST90027163, van Rheenen W, European)
# ───────────────────────────────────────────────────────────────────────────────

library(TwoSampleMR)
library(MRPRESSO)
library(dplyr)
library(readr)
library(ggplot2)
library(stringr)

# 1) Prepare exposure and outcome data
exposure_dat <- extract_instruments(outcomes = 'eqtl-a-ENSG00000165029', p1 = 5e-8)
outcome_dat  <- extract_outcome_data(snps = exposure_dat$SNP, outcomes = 'ebi-a-GCST90027163')

# Harmonise alleles
harmonised_dat <- harmonise_data(exposure_dat, outcome_dat, action = 1)

# ── Instrument strength (exposure side)
# 1) Per-SNP F statistics
harmonised_dat <- harmonised_dat %>%
  mutate(F_stat = (beta.exposure^2) / (se.exposure^2))

# 2) I^2_GX for MR-Egger NOME check
calc_I2GX <- function(b, se) {
  ok <- is.finite(b) & is.finite(se) & (se > 0)
  b  <- b[ok]; se <- se[ok]
  k  <- length(b)
  if (k < 3) return(NA_real_)
  w  <- 1 / (se^2)
  mu <- sum(w * b) / sum(w)
  Q  <- sum(w * (b - mu)^2)
  if (!is.finite(Q) || Q <= 0) return(0)
  max(0, (Q - (k - 1)) / Q)
}
I2_GX_val <- calc_I2GX(harmonised_dat$beta.exposure, harmonised_dat$se.exposure)

# 3) Summary metrics
F_summary <- harmonised_dat %>%
  summarize(
    exposure_id   = first(id.exposure),
    outcome_id    = first(id.outcome),
    nsnp          = n(),
    mean_F        = mean(F_stat, na.rm = TRUE),
    median_F      = median(F_stat, na.rm = TRUE),
    prop_F_lt_10  = mean(F_stat < 10, na.rm = TRUE),
    min_F         = min(F_stat, na.rm = TRUE),
    max_F         = max(F_stat, na.rm = TRUE),
    I2_GX         = I2_GX_val
  )
print(F_summary)

# 4) Export instrument-strength tables
write_csv(
  harmonised_dat %>% select(SNP, beta.exposure, se.exposure, F_stat),
  "Instrument_strength_per_SNP_ABCA1_eQTL.csv"
)
write_csv(
  F_summary,
  "Instrument_strength_summary_ABCA1_eQTL.csv"
)

# 5) Visualise distribution of per-SNP F statistics
p <- ggplot(harmonised_dat, aes(x = F_stat)) +
  geom_density(fill = "#EEDC82", alpha = 0.6, color = "#E5C494", size = 1) +
  geom_vline(xintercept = 10, linetype = "dashed", color = "#EE6363", size = 1) +
  annotate("text", x = 12, y = 0.02, label = "Threshold F = 10",
           color = "#EE6363", hjust = 0, vjust = -1, size = 4) +
  geom_rug(aes(x = F_stat), sides = "b", alpha = 0.9, color = "#4682B4", size = 0.9) +
  labs(
    title    = "Distribution of per-SNP F statistics (exposure side)",
    subtitle = paste0(
      "nsnp = ", nrow(harmonised_dat),
      "; mean F = ", round(F_summary$mean_F, 1),
      "; median F = ", round(F_summary$median_F, 1),
      "; I²GX = ", round(F_summary$I2_GX, 3)
    ),
    x = "F statistic", y = "Density"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title      = element_text(face = "bold", size = 16),
    plot.subtitle   = element_text(size = 12, color = "grey40")
  )
p
ggsave("Fstat_density_ABCA1_eQTL.pdf", p, width = 7, height = 5)

# 6) MR using multiple estimators
mr_results <- mr(
  harmonised_dat,
  method_list = c(
    "mr_ivw",
    "mr_raps",
    "mr_ivw_mre",
    "mr_egger_regression",
    "mr_weighted_median",
    "mr_weighted_mode"
  )
)
print(mr_results)

# 7) Pleiotropy (Egger intercept)
pleio <- mr_pleiotropy_test(harmonised_dat)
p_pleio_egger <- pleio$pval[1]
p_pleio_egger

# 8) Heterogeneity (Q tests)
het <- mr_heterogeneity(harmonised_dat)
p_Q_ivw   <- het %>% dplyr::filter(grepl("Inverse variance weighted$", method)) %>% dplyr::pull(Q_pval)
p_Q_egger <- het %>% dplyr::filter(method == "MR Egger") %>% dplyr::pull(Q_pval)
p_Q_ivw
p_Q_egger

# 9) MR-PRESSO
mr_presso_result <- mr_presso(
  BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
  SdOutcome = "se.outcome",    SdExposure = "se.exposure",
  OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
  data = harmonised_dat, NbDistribution = 1000, SignifThreshold = 0.05
)
print(mr_presso_result$`Main MR results`)
print(mr_presso_result$`Global Test`$Pvalue)
if(length(mr_presso_result$`Outlier SNPs`) > 0){
  print(mr_presso_result$`Outlier SNPs`)
}
if(length(mr_presso_result$`Outlier SNPs`) > 0){
  outliers <- mr_presso_result$`Outlier SNPs`$`MR-PRESSO Outlier Test`$Outlier
  harmonised_dat_no_outliers <- harmonised_dat[ !harmonised_dat$SNP %in% outliers, ]
  mr_no_outlier_res <- mr(harmonised_dat_no_outliers, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))
  print(mr_no_outlier_res)
}

# 10) Steiger directionality
steiger_res <- directionality_test(harmonised_dat)
print(steiger_res)

# 11) MR scatter plot export
mr_scatter <- mr_scatter_plot(mr_results, harmonised_dat)

# 12) Polished scatter with IVW line
library(ggplot2)
options(repr.plot.width=5, repr.plot.height=5)

df <- mr_scatter[[1]]$data
res <- mr_results[1, ]
beta <- round(res$b, 2)
se <- round(res$se, 2)
or <- round(exp(res$b), 2)
pval <- signif(res$pval, 2)

p1 <- ggplot(df, aes(x = beta.exposure, y = beta.outcome)) +
  geom_errorbar(aes(ymin = beta.outcome - se.outcome,
                    ymax = beta.outcome + se.outcome),
                width = 0, color = "#E5C494") +
  geom_errorbarh(aes(xmin = beta.exposure - se.exposure,
                     xmax = beta.exposure + se.exposure),
                 height = 0, color = "#E5C494") +
  geom_point(color = "#EEDC82", size = 2.2) +
  geom_abline(intercept = 0, slope = beta, linetype = "dashed",
              color = "#EE6363", size = 1) +
  annotate("text",
           fontface = 'bold',
           x = min(df$beta.exposure) - 0.02 ,
           y = max(df$beta.outcome) + 0.015,
           hjust = 0, vjust = 1,
           size = 4.2, family = "sans",
           label = paste0(
             "Inverse variance weighted\n",
             expression(beta), ": ", beta, "\n",
             "SE: ", se, "\n",
             "OR: ", or, "\n",
             "P value: ", pval
           )) +
  labs(
    x = "SNP effect on ABCA1 expression",
    y = "SNP effect on ALS risk"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
p1
ggsave(p1, file="van Rheenen W ABCA1_MR_scatter_Plot.pdf", width = 5, height = 5)

# 13) Egger intercept table
egger_test <- mr_pleiotropy_test(harmonised_dat)
print(egger_test)

# 14) Leave-One-Out (LOO)
loo_result <- mr_leaveoneout(harmonised_dat)
print(loo_result)
mr_leaveoneout_plot(loo_result)

str(loo_result)
library(openxlsx)
write.xlsx(loo_result, file='van Rheenen W ABCA1_MR多效性残差与异常值检验.xlsx')

# 15) Custom LOO forest
library(ggplot2)
library(dplyr)

df <- loo_result %>%
  mutate(
    lower = b - 1.96 * se,
    upper = b + 1.96 * se,
    sig = case_when(
      SNP == "All" ~ "Main",
      lower > 0 | upper < 0 ~ "Significant",
      TRUE ~ "Nonsignificant"
    ),
    SNP = ifelse(SNP == "All", "All–IVW", SNP),
    SNP = factor(SNP, levels = rev(SNP))
  )

color_scheme <- c("Main" = "#EE6363",
                  "Significant" = "#EEB4B4",
                  "Nonsignificant" = "#999999")

p <- ggplot(df, aes(x = b, y = SNP)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.4) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.3, color = "gray60") +
  geom_point(aes(color = sig), shape = 22, fill = "white", size = 3.5, stroke = 1.2) +
  scale_color_manual(values = color_scheme) +
  labs(
    x = expression(beta~"(95% CI)"),
    y = NULL,
    title = "MR Leave-One-Out Sensitivity Analysis: ABCA1 → ALS"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 10),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "none"
  ) +
  xlim(min(df$lower) - 0.01, max(df$upper) + 0.01)
p
ggsave("van Rheenen W ABCA1_MR_leaveoneout_forest_final.pdf", p, width = 2.5, height = 5)

# Cochran’s Q (TwoSampleMR helper)
mr_heterogeneity(harmonised_dat)

# Single-SNP MR estimates
singlesnp_res <- mr_singlesnp(harmonised_dat)
singlesnp_res

# Funnel plot data
funnel_plot <- mr_funnel_plot(singlesnp_res)
funnel_plot[[1]]

library(ggplot2)

# Take first 11 rows for display
singlesnp_plot <- singlesnp_res[1:11, ]

# Reference IVW beta for line
ivw_beta <- mr_results$b[mr_results$method == "Inverse variance weighted"]

# Inverse SE for y-axis
singlesnp_plot$inv_se <- 1 / singlesnp_plot$se

p1 <- ggplot(singlesnp_plot, aes(x = b, y = inv_se)) +
  geom_point(shape = 16, size = 3, color = "#EE6363", alpha = 1) +
  geom_vline(xintercept = ivw_beta, linetype = "dashed", color = "grey30", linewidth = 1) +
  labs(
    x = expression(beta~"(IVW)"),
    y = expression("1/SE (IVW)")
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.title.x = element_text(face = "bold", size = 13),
    axis.title.y = element_text(face = "bold", size = 13),
    axis.text = element_text(size = 11)
  )
p1
ggsave(p1, file = 'van Rheenen W ABCA1漏斗图.pdf', height = 5, width = 5)

# OR table for all MR methods
or_results <- data.frame(
  Method = mr_results$method,
  Beta = mr_results$b,
  OR = exp(mr_results$b),
  Lower_CI = exp(mr_results$b - 1.96 * mr_results$se),
  Upper_CI = exp(mr_results$b + 1.96 * mr_results$se),
  P_value = mr_results$pval
)
print(or_results)

library(openxlsx)
write.xlsx(or_results, file = 'van Rheenen W ABCA1_MR_结果表格.xlsx')

# Directionality
directionality_test(harmonised_dat)


# ───────────────────────────────────────────────────────────────────────────────
# ABCA1 → ALS (ebi-a-GCST005647, Nicolas A, European)
# ───────────────────────────────────────────────────────────────────────────────

exposure_dat <- extract_instruments(outcomes = 'eqtl-a-ENSG00000165029', p1 = 5e-8)
outcome_dat  <- extract_outcome_data(snps = exposure_dat$SNP, outcomes = 'ebi-a-GCST005647')

harmonised_dat <- harmonise_data(exposure_dat, outcome_dat, action = 1)

mr_results <- mr(
  harmonised_dat,
  method_list = c(
    "mr_ivw",
    "mr_raps",
    "mr_ivw_mre",
    "mr_egger_regression",
    "mr_weighted_median",
    "mr_weighted_mode"
  )
)
print(mr_results)

pleio <- mr_pleiotropy_test(harmonised_dat)
p_pleio_egger <- pleio$pval[1]
p_pleio_egger

mr_presso_result <- mr_presso(
  BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
  SdOutcome = "se.outcome",    SdExposure = "se.exposure",
  OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
  data = harmonised_dat, NbDistribution = 1000, SignifThreshold = 0.05
)
print(mr_presso_result$`Main MR results`)
print(mr_presso_result$`Global Test`$Pvalue)
if(length(mr_presso_result$`Outlier SNPs`) > 0){
  print(mr_presso_result$`Outlier SNPs`)
}
if(length(mr_presso_result$`Outlier SNPs`) > 0){
  outliers <- mr_presso_result$`Outlier SNPs`$`MR-PRESSO Outlier Test`$Outlier
  harmonised_dat_no_outliers <- harmonised_dat[ !harmonised_dat$SNP %in% outliers, ]
  mr_no_outlier_res <- mr(harmonised_dat_no_outliers, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))
  print(mr_no_outlier_res)
}

het <- mr_heterogeneity(harmonised_dat)
p_Q_ivw   <- het %>% dplyr::filter(grepl("Inverse variance weighted$", method)) %>% dplyr::pull(Q_pval)
p_Q_egger <- het %>% dplyr::filter(method == "MR Egger") %>% dplyr::pull(Q_pval)
p_Q_ivw
p_Q_egger

steiger_res <- directionality_test(harmonised_dat)
print(steiger_res)

mr_scatter <- mr_scatter_plot(mr_results, harmonised_dat)

library(ggplot2)
options(repr.plot.width=5, repr.plot.height=5)

df <- mr_scatter[[1]]$data
res <- mr_results[1, ]
beta <- round(res$b, 2)
se <- round(res$se, 2)
or <- round(exp(res$b), 2)
pval <- signif(res$pval, 2)

p1 <- ggplot(df, aes(x = beta.exposure, y = beta.outcome)) +
  geom_errorbar(aes(ymin = beta.outcome - se.outcome,
                    ymax = beta.outcome + se.outcome),
                width = 0, color = "#E5C494") +
  geom_errorbarh(aes(xmin = beta.exposure - se.exposure,
                     xmax = beta.exposure + se.exposure),
                 height = 0, color = "#E5C494") +
  geom_point(color = "#EEDC82", size = 2.2) +
  geom_abline(intercept = 0, slope = beta, linetype = "dashed",
              color = "#EE6363", size = 1) +
  annotate("text",
           fontface = 'bold',
           x = min(df$beta.exposure) - 0.02 ,
           y = max(df$beta.outcome) + 0.015,
           hjust = 0, vjust = 1,
           size = 4.2, family = "sans",
           label = paste0(
             "Inverse variance weighted\n",
             expression(beta), ": ", beta, "\n",
             "SE: ", se, "\n",
             "OR: ", or, "\n",
             "P value: ", pval
           )) +
  labs(
    x = "SNP effect on ABCA1 expression",
    y = "SNP effect on ALS risk"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
p1
ggsave(p1, file="Nicolas A ABCA1_MR_scatter_Plot.pdf", width = 5, height = 5)

# Sensitivity: Egger intercept
egger_test <- mr_pleiotropy_test(harmonised_dat)
print(egger_test)

# LOO
loo_result <- mr_leaveoneout(harmonised_dat)
print(loo_result)

library(dplyr)
library(ggplot2)
library(openxlsx)

# Keep the “All” row
loo_all <- loo_result %>% filter(SNP == "All") %>% distinct()
b_all   <- loo_all$b[1]
p_all   <- loo_all$p[1]

# Deduplicate non-All rows by (SNP, b, se)
loo_clean <- loo_result %>%
  filter(SNP != "All") %>%
  distinct(SNP, b, se, .keep_all = TRUE) %>%
  arrange(SNP)

# Compute CI and simple influence markers
loo_clean <- loo_clean %>%
  mutate(
    ci_lo      = b - 1.96*se,
    ci_hi      = b + 1.96*se,
    delta_b    = b - b_all,
    rel_change = ifelse(abs(b_all) > 0, delta_b/abs(b_all), NA_real_),
    flips_sig  = xor(p < 0.05, p_all < 0.05)
  )

# Add back the “All” row
loo_clean_with_all <- bind_rows(
  loo_clean,
  loo_all %>% mutate(ci_lo = b - 1.96*se, ci_hi = b + 1.96*se,
                     delta_b = NA_real_, rel_change = NA_real_, flips_sig = NA)
)

df <- loo_clean_with_all %>%
  mutate(
    lower = b - 1.96 * se,
    upper = b + 1.96 * se,
    sig = case_when(
      SNP == "All" ~ "Main",
      lower > 0 | upper < 0 ~ "Significant",
      TRUE ~ "Nonsignificant"
    ),
    SNP = ifelse(SNP == "All", "All–IVW", SNP),
    SNP = factor(SNP, levels = rev(SNP))
  )

color_scheme <- c("Main" = "#EE6363",
                  "Significant" = "#EEB4B4",
                  "Nonsignificant" = "#999999")

p <- ggplot(df, aes(x = b, y = SNP)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.4) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.3, color = "gray60") +
  geom_point(aes(color = sig), shape = 22, fill = "white", size = 3.5, stroke = 1.2) +
  scale_color_manual(values = color_scheme) +
  labs(
    x = expression(beta~"(95% CI)"),
    y = NULL,
    title = "MR Leave-One-Out Sensitivity Analysis: ABCA1 → ALS"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 10),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "none"
  ) +
  xlim(min(df$lower) - 0.01, max(df$upper) + 0.01)
p
ggsave("Nicolas A ABCA1_MR_leaveoneout_forest_final.pdf", p, width = 2.5, height = 5)

str(loo_result)
library(openxlsx)
write.xlsx(loo_clean_with_all, file='Nicolas A ABCA1_MR多效性残差与异常值检验.xlsx')

mr_presso_result[[1]]
mr_presso_result[[2]]

# Cochran’s Q
mr_heterogeneity(harmonised_dat)

# Single-SNP MR
singlesnp_res <- mr_singlesnp(harmonised_dat)

# Funnel data
funnel_plot <- mr_funnel_plot(singlesnp_res)
funnel_plot[[1]]

library(ggplot2)
singlesnp_plot <- singlesnp_res[1:46, ]

ivw_beta <- mr_results$b[mr_results$method == "Inverse variance weighted"]

singlesnp_plot$inv_se <- 1 / singlesnp_plot$se

p1 <- ggplot(singlesnp_plot, aes(x = b, y = inv_se)) +
  geom_point(shape = 16, size = 3, color = "#EE6363", alpha = 1) +
  geom_vline(xintercept = ivw_beta, linetype = "dashed", color = "grey30", linewidth = 1) +
  labs(
    x = expression(beta~"(IVW)"),
    y = expression("1/SE (IVW)")
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.title.x = element_text(face = "bold", size = 13),
    axis.title.y = element_text(face = "bold", size = 13),
    axis.text = element_text(size = 11)
  )
p1
ggsave(p1, file='Nicolas A ABCA1漏斗图.pdf', height=5, width=5)

or_results <- data.frame(
  Method = mr_results$method,
  Beta = mr_results$b,
  OR = exp(mr_results$b),
  Lower_CI = exp(mr_results$b - 1.96 * mr_results$se),
  Upper_CI = exp(mr_results$b + 1.96 * mr_results$se),
  P_value = mr_results$pval
)
print(or_results)
library(openxlsx)
write.xlsx(or_results, file='Nicolas A ABCA1_MR_结果表格.xlsx')

directionality_test(harmonised_dat)


# ───────────────────────────────────────────────────────────────────────────────
# ABCA1 → ALS (ebi-a-GCST90013429, Iacoangeli A, Mixed)
# ───────────────────────────────────────────────────────────────────────────────

exposure_dat <- extract_instruments(outcomes = 'eqtl-a-ENSG00000165029', p1 = 5e-8)
outcome_dat  <- extract_outcome_data(snps = exposure_dat$SNP, outcomes = 'ebi-a-GCST90013429')

harmonised_dat <- harmonise_data(exposure_dat, outcome_dat, action = 1)

mr_results <- mr(
  harmonised_dat,
  method_list = c(
    "mr_ivw",
    "mr_raps",
    "mr_ivw_mre",
    "mr_egger_regression",
    "mr_weighted_median",
    "mr_weighted_mode"
  )
)
print(mr_results)

pleio <- mr_pleiotropy_test(harmonised_dat)
p_pleio_egger <- pleio$pval[1]
p_pleio_egger 

mr_presso_result <- mr_presso(
  BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
  SdOutcome = "se.outcome",    SdExposure = "se.exposure",
  OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
  data = harmonised_dat, NbDistribution = 1000, SignifThreshold = 0.05
)
print(mr_presso_result$`Main MR results`)
print(mr_presso_result$`Global Test`$Pvalue)
if(length(mr_presso_result$`Outlier SNPs`) > 0){
  print(mr_presso_result$`Outlier SNPs`)
}
if(length(mr_presso_result$`Outlier SNPs`) > 0){
  outliers <- mr_presso_result$`Outlier SNPs`$`MR-PRESSO Outlier Test`$Outlier
  harmonised_dat_no_outliers <- harmonised_dat[ !harmonised_dat$SNP %in% outliers, ]
  mr_no_outlier_res <- mr(harmonised_dat_no_outliers, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))
  print(mr_no_outlier_res)
}

het <- mr_heterogeneity(harmonised_dat)
p_Q_ivw   <- het %>% dplyr::filter(grepl("Inverse variance weighted$", method)) %>% dplyr::pull(Q_pval)
p_Q_egger <- het %>% dplyr::filter(method == "MR Egger") %>% dplyr::pull(Q_pval)
p_Q_ivw
p_Q_egger

steiger_res <- directionality_test(harmonised_dat)
print(steiger_res)


# ───────────────────────────────────────────────────────────────────────────────
# ABCA1 → ALS (ebi-a-GCST004901, Benyamin B, East Asian)
# ───────────────────────────────────────────────────────────────────────────────

exposure_dat <- extract_instruments(outcomes = 'eqtl-a-ENSG00000165029', p1 = 5e-8)
outcome_dat  <- extract_outcome_data(snps = exposure_dat$SNP, outcomes = 'ebi-a-GCST004901')

harmonised_dat <- harmonise_data(exposure_dat, outcome_dat, action = 1)

mr_results <- mr(
  harmonised_dat,
  method_list = c(
    "mr_ivw",
    "mr_raps",
    "mr_ivw_mre",
    "mr_egger_regression",
    "mr_weighted_median",
    "mr_weighted_mode"
  )
)
print(mr_results)

pleio <- mr_pleiotropy_test(harmonised_dat)
p_pleio_egger <- pleio$pval[1]
p_pleio_egger 

mr_presso_result <- mr_presso(
  BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
  SdOutcome = "se.outcome",    SdExposure = "se.exposure",
  OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
  data = harmonised_dat, NbDistribution = 1000, SignifThreshold = 0.05
)
print(mr_presso_result$`Main MR results`)
print(mr_presso_result$`Global Test`$Pvalue)
if(length(mr_presso_result$`Outlier SNPs`) > 0){
  print(mr_presso_result$`Outlier SNPs`)
}
if(length(mr_presso_result$`Outlier SNPs`) > 0){
  outliers <- mr_presso_result$`Outlier SNPs`$`MR-PRESSO Outlier Test`$Outlier
  harmonised_dat_no_outliers <- harmonised_dat[ !harmonised_dat$SNP %in% outliers, ]
  mr_no_outlier_res <- mr(harmonised_dat_no_outliers, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))
  print(mr_no_outlier_res)
}

het <- mr_heterogeneity(harmonised_dat)
p_Q_ivw   <- het %>% dplyr::filter(grepl("Inverse variance weighted$", method)) %>% dplyr::pull(Q_pval)
p_Q_egger <- het %>% dplyr::filter(method == "MR Egger") %>% dplyr::pull(Q_pval)
p_Q_ivw
p_Q_egger

steiger_res <- directionality_test(harmonised_dat)
print(steiger_res)
