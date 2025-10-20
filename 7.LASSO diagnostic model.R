library(tidyverse)
library(glmnet)
source("msvmRFE.R")
library(VennDiagram)
library(sigFeature)
library(e1071)
library(caret)
library(randomForest)

setwd("Diag_models/LASSO")

genes <- readRDS("../../ML_features_select/common_genes.rds")

# -------------------------------
# Train data
# -------------------------------
load("../../Data_preprocessing/train_cohort.rda")
ls()

diagnosis <- train_cohort_clin$diagnosis
datExpr   <- as.data.frame(cbind(diagnosis, t(train_cohort_expr)))
datExpr   <- datExpr[, c("diagnosis", genes)]
str(datExpr)

# Rename the first column to 'group'
colnames(datExpr)[1] <- "group"
datExpr$group <- as.factor(datExpr$group)

# Convert all expression columns (except the first) to numeric
datExpr[, -1] <- lapply(datExpr[, -1], function(x) as.numeric(as.character(x)))

# Check structure after conversion
str(datExpr)

# -------------------------------
# Feature selection with LASSO
# -------------------------------
# glmnet input format:
#   rows = samples, columns = genes
#   first col = class labels; remaining cols = expression matrix
x <- as.matrix(datExpr[, -1]) %>% scale()
y <- ifelse(datExpr$group == "CON", 0, 1)  # encode groups as 0/1

options(repr.plot.width = 5, repr.plot.height = 5)
fit <- glmnet(x, y, family = "binomial", alpha = 1, lambda = NULL)
plot(fit, xvar = "dev", label = FALSE)

pdf("A_lasso.pdf", width = 5, height = 5)
plot(fit, xvar = "dev", label = FALSE)
dev.off()

# Cross-validation to choose optimal lambda
lasso_model <- cv.glmnet(
  x, y, family = "binomial", alpha = 1,
  type.measure = "auc",
  nfold = 10
)

# Best lambda (also can use lambda.1se for fewer features)
best_lambda <- lasso_model$lambda.1se

pdf(file = "LASSO.Coefficient_Paths_Plot.pdf", height = 5, width = 5)
plot(lasso_model)
dev.off()

plot(lasso_model)
best_lambda

# Fit final LASSO with best lambda
lasso_fit <- glmnet(x, y, family = "binomial", alpha = 1, lambda = best_lambda)

# Inspect coefficients
coef(lasso_fit)
saveRDS(lasso_fit, file = "lasso_model.rds")

# In-sample prediction and ROC/AUC
predictions <- predict(lasso_fit, x, type = "response")

library(pROC)
roc_curve <- roc(y, predictions)
plot(roc_curve)
auc(roc_curve)

str(lasso_fit)

# Save train ROC object
train_ROC <- roc_curve

# -------------------------------
# Test data
# -------------------------------
load("../../Data_preprocessing/test_cohort.rda")
diagnosis <- test_cohort_clin$diagnosis
datExpr   <- as.data.frame(cbind(diagnosis, t(test_cohort_exp)))
datExpr   <- datExpr[, c("diagnosis", genes)]
str(datExpr)

# Rename the first column to 'group'
colnames(datExpr)[1] <- "group"
datExpr$group <- as.factor(datExpr$group)

# Convert all expression columns (except the first) to numeric
datExpr[, -1] <- lapply(datExpr[, -1], function(x) as.numeric(as.character(x)))

# Check structure after conversion
str(datExpr)

# Use training scaling parameters for test set
x_test <- scale(datExpr[-1], center = attr(x, "scaled:center"), scale = attr(x, "scaled:scale"))
y_test <- ifelse(datExpr$group == "CON", 0, 1)  # encode groups as 0/1

# Predict and evaluate on test set
prediction_test <- predict(lasso_fit, x_test, type = "response")

library(pROC)
roc_curve <- roc(y_test, prediction_test)
plot(roc_curve)
auc(roc_curve)

test_ROC <- roc_curve

# -------------------------------
# External validation data
# -------------------------------
load("../../Data_preprocessing/valid_cohort.rda")
diagnosis <- valid_cohort_clin$diagnosis
datExpr   <- as.data.frame(cbind(diagnosis, t(valid_cohort_expr)))
datExpr   <- datExpr[, c("diagnosis", genes)]
str(datExpr)

# Rename the first column to 'group'
colnames(datExpr)[1] <- "group"
datExpr$group <- as.factor(datExpr$group)

# Convert all expression columns (except the first) to numeric
datExpr[, -1] <- lapply(datExpr[, -1], function(x) as.numeric(as.character(x)))

# Keep ALS and CON only
datExpr_filtered <- datExpr[datExpr$group %in% c("ALS", "CON"), ]

# Inspect filtered structure
str(datExpr_filtered)

# Use training scaling parameters for validation set
x_val <- scale(datExpr_filtered[-1], center = attr(x, "scaled:center"), scale = attr(x, "scaled:scale"))
y_val <- ifelse(datExpr_filtered$group == "CON", 0, 1)  # encode groups as 0/1

# Predict and evaluate on validation set
prediction_val <- predict(lasso_fit, x_val, type = "response")

library(pROC)
roc_curve <- roc(y_val, prediction_val)
plot(roc_curve)
auc(roc_curve)

validation_ROC <- roc_curve

# -------------------------------
# Plot combined ROC curves (Train/Test/Validation)
# -------------------------------
library(pROC)
library(ggplot2)
library(dplyr)

# Extract ROC curves into data frames
train_data <- data.frame(
  specificity = train_ROC$specificities,
  sensitivity = train_ROC$sensitivities,
  dataset     = "Train"
)

test_data <- data.frame(
  specificity = test_ROC$specificities,
  sensitivity = test_ROC$sensitivities,
  dataset     = "Test"
)

validation_data <- data.frame(
  specificity = validation_ROC$specificities,
  sensitivity = validation_ROC$sensitivities,
  dataset     = "Validation"
)

# Combine all
roc_data <- bind_rows(train_data, test_data, validation_data)

# Compute AUCs
train_auc      <- auc(train_ROC)
test_auc       <- auc(test_ROC)
validation_auc <- auc(validation_ROC)

# Legend labels with AUCs
labels <- c(
  paste("Train (ROC =", round(train_auc, 3), ")"),
  paste("Test (ROC =", round(test_auc, 3), ")"),
  paste("Validation (ROC =", round(validation_auc, 3), ")")
)

# Factor order for legend and color mapping
roc_data$dataset <- factor(roc_data$dataset, levels = c("Train", "Test", "Validation"))

# Plot
p <- ggplot(roc_data, aes(x = specificity, y = sensitivity, color = dataset)) +
  geom_line(size = 1.2) +
  geom_segment(aes(x = 0, y = 1, xend = 1, yend = 0),
               linetype = "dashed", color = "gray") +
  labs(
    title    = "Receiver Operating Characteristic (ROC)",
    subtitle = "ALS vs CON",
    x = "Specificity",
    y = "Sensitivity",
    color = "LASSO Model:"
  ) +
  scale_x_reverse(limits = c(1, 0), breaks = seq(0, 1, 0.25), expand = c(0.02, 0.02)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0.02, 0.02)) +
  scale_color_manual(
    values = c("Train" = "#CDB38B", "Test" = "#8DB6CD", "Validation" = "#CD5555"),
    labels = labels
  ) +
  theme_minimal() +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold", size = 10),
    plot.subtitle   = element_text(hjust = 0.5, face = "bold", size = 8, color = "black"),
    axis.title      = element_text(face = "bold", size = 8),
    axis.text       = element_text(color = "black", size = 6),
    axis.ticks      = element_line(size = 1),
    axis.ticks.length = unit(-0.2, "cm"),
    panel.border    = element_rect(color = "black", fill = NA, size = 1),
    panel.grid      = element_blank(),
    legend.title    = element_text(face = "bold", size = 8, color = "black"),
    legend.text     = element_text(face = "bold", size = 8, color = "black"),
    legend.position = c(0.80, 0.15)
  )

print(p)
ggsave(p, file = "LASSO_ROC.Plot.pdf", height = 5, width = 5)

# -------------------------------
# Summarize ROC metrics
# -------------------------------
library(pROC)
library(dplyr)

# Helper to extract metrics from a pROC::roc object
metrics_from_roc <- function(roc_obj, dataset) {
  auc_val <- as.numeric(auc(roc_obj))
  best    <- coords(
    roc_obj, x = "best", best.method = "youden",
    ret = c("threshold", "sensitivity", "specificity", "accuracy"),
    transpose = FALSE
  )
  
  n_cases    <- length(roc_obj$cases)
  n_controls <- length(roc_obj$controls)
  
  tibble(
    dataset           = dataset,
    N                 = n_cases + n_controls,
    `Case proportion` = n_cases / (n_cases + n_controls),
    AUC               = auc_val,
    Accuracy          = as.numeric(best["accuracy"]),
    Sensitivity       = as.numeric(best["sensitivity"]),
    Specificity       = as.numeric(best["specificity"]),
    Threshold         = as.numeric(best["threshold"]),
    Direction         = roc_obj$direction
  )
}

summarize_three_rocs <- function(train_ROC, test_ROC, validation_ROC) {
  bind_rows(
    metrics_from_roc(train_ROC, "Train"),
    metrics_from_roc(test_ROC,  "Test"),
    metrics_from_roc(validation_ROC, "Validation")
  ) %>%
    mutate(across(c(AUC, Accuracy, Sensitivity, Specificity, Threshold, `Case proportion`),
                  ~ round(.x, 3)))
}

metrics_tbl <- summarize_three_rocs(train_ROC, test_ROC, validation_ROC)
print(metrics_tbl)
readr::write_csv(metrics_tbl, "LASSO_metrics_summary.csv")

# -------------------------------
# Plot non-zero LASSO coefficients (feature importance)
# -------------------------------
coef(lasso_fit)

beta  <- as.matrix(lasso_fit$beta)
genes <- rownames(beta)
coefs <- as.numeric(beta)

nonzero_idx <- which(coefs != 0)

coefs_df <- data.frame(
  gene = genes[nonzero_idx],
  coef = coefs[nonzero_idx]
) %>%
  mutate(direction = ifelse(coef > 0, "Positive", "Negative")) %>%
  arrange(coef) %>%
  mutate(gene = factor(gene, levels = gene))  # preserve order on y-axis

options(repr.plot.width = 5, repr.plot.height = 4)
p1 <- ggplot(coefs_df, aes(x = coef, y = gene, fill = direction)) +
  geom_col(width = 0.6) +
  # Label placement: positives on the left, negatives on the right (relative to bars)
  geom_text(
    aes(
      label = gene,
      x = ifelse(coef > 0, 0 - 0.05, 0 + 0.05)
    ),
    hjust = ifelse(coefs_df$coef > 0, 1, 0),
    color = ifelse(coefs_df$direction == "Positive", "#CD5555", "#00688B"),
    size = 6,
    fontface = ifelse(coefs_df$gene == "ABCA1", "bold", "plain")
  ) +
  scale_fill_manual(values = c("Positive" = "#CD5555", "Negative" = "#00688B")) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.5, color = "black") +
  scale_x_continuous(expand = expansion(mult = c(0.2, 0.2))) +
  labs(
    title = "Feature Importance in LASSO Diagnostic Model",
    x     = "Model Coefficient (β)",
    y     = NULL
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid       = element_blank(),
    panel.border     = element_rect(fill = NA, color = "black"),
    axis.text.y      = element_blank(),
    axis.ticks.y     = element_blank(),
    axis.title.x     = element_text(size = 14, face = "bold"),
    plot.title       = element_text(hjust = 0.5, face = "bold", size = 15),
    legend.position  = "none"
  )

p1
ggsave(p1, file = "LASSO模型特征重要性.pdf", height = 4, width = 5)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(readr)
  library(pROC)
  library(mgcv)
  library(grid)
})

# ============== Utilities ==============
# Clamp probabilities to (eps, 1-eps) to avoid qlogis overflow at 0/1
.safe_prob <- function(p, eps = 1e-6) pmin(pmax(as.numeric(p), eps), 1 - eps)
# Prefer cairo_pdf for better embedding; fall back to pdf
.pdf_device <- function() if (capabilities("cairo")) cairo_pdf else pdf

# Compute Eavg/Emax (for numeric summary only; not for plotting)
.eavg_emax <- function(y, p, bins = 10){
  tibble(y = as.integer(y), p = .safe_prob(p)) |>
    mutate(bin = ntile(p, bins)) |>
    group_by(bin) |>
    summarise(p_mean = mean(p), y_mean = mean(y), .groups = "drop") |>
    summarise(Eavg = mean(abs(y_mean - p_mean)),
              Emax = max(abs(y_mean - p_mean))) |>
    as.list() |> unlist()
}

# ============== Binning + Wilson CI ==============
# Wilson 95% CI for bin-level observed event rate
.wilson_ci <- function(k, n, conf = 0.95){
  z <- qnorm(1 - (1 - conf)/2)
  ph <- k/n
  den <- 1 + z^2/n
  cen <- (ph + z^2/(2*n)) / den
  half <- z * sqrt((ph*(1-ph) + z^2/(4*n))/n) / den
  c(lower = pmax(0, cen - half), upper = pmin(1, cen + half))
}

# Bin data by equal frequency; compute mean p, mean y, and Wilson CI
.bin_summary <- function(y, p, bins = 10){
  tibble(y = as.integer(y), p = .safe_prob(p)) |>
    mutate(bin = ntile(p, bins)) |>
    group_by(bin) |>
    summarise(n = n(), p_mean = mean(p), y_sum = sum(y), y_mean = mean(y), .groups = "drop") |>
    rowwise() |>
    mutate(ci = list(.wilson_ci(y_sum, n)),
           ci_l = ci[[1]], ci_u = ci[[2]]) |>
    ungroup() |> select(-ci)
}

# ============== Post-hoc calibrators (Platt / Isotonic) ==============
# Platt: logit(Pr) = a + b * logit(p)
platt_fit <- function(y, p){
  y <- as.integer(y); p <- .safe_prob(p)
  fit <- glm(y ~ qlogis(p), family = binomial())
  list(method = "platt", a = unname(coef(fit)[1]), b = unname(coef(fit)[2]))
}
platt_apply <- function(p, model){
  p <- .safe_prob(p)
  plogis(model$a + model$b * qlogis(p))
}

# Isotonic: monotone nonparametric; expose via linear/constant interpolation
isotonic_fit <- function(y, p){
  y <- as.integer(y); p <- .safe_prob(p)
  o   <- order(p)
  iso <- isoreg(p[o], y[o])
  xy  <- tibble(x = iso$x, y = iso$yf) |>
    group_by(x) |> summarise(y = mean(y), .groups = "drop") |>
    arrange(x)
  list(method = "isotonic", x = xy$x, yf = xy$y)
}
isotonic_apply <- function(p_new, model, iso_interp = c("linear","constant")){
  iso_interp <- match.arg(iso_interp)
  p_new <- .safe_prob(p_new)
  as.numeric(approx(x = model$x, y = model$yf, xout = p_new,
                    method = if (iso_interp == "linear") "linear" else "constant",
                    rule = 2)$y)
}

# ============== Logistic calibration line (blue) ==============
# Fit logit(Pr) = a + b * logit(p), then map onto probability scale over a grid
logistic_cal_line <- function(y, p, grid_n = 200){
  y <- as.integer(y); p <- .safe_prob(p)
  fit <- glm(y ~ qlogis(p), family = binomial())
  a <- unname(coef(fit)[1]); b <- unname(coef(fit)[2])
  xg <- seq(1e-6, 1 - 1e-6, length.out = grid_n)
  yg <- plogis(a + b * qlogis(xg))
  list(df = data.frame(x = xg, y = yg, series = "Logistic calibration"),
       a = a, b = b)
}

# ============== LOESS smoothed curve + 95% band (fallback to binomial GAM) ==============
.smooth_band <- function(y, p, span = 0.75, grid_n = 200){
  y <- as.integer(y); p <- .safe_prob(p)
  ok <- is.finite(y) & is.finite(p); y <- y[ok]; p <- p[ok]
  
  # Small-sample fallback to constant mean
  if (length(y) < 10 || length(unique(p)) < 2){
    prev <- mean(y); grid <- seq(1e-6, 1-1e-6, length.out = grid_n)
    se0  <- sqrt(pmax(prev*(1-prev)/max(1, length(y)), 1e-12))
    return(tibble(x = grid, y = rep(prev, grid_n),
                  y_l = pmax(0, prev - 1.96*se0),
                  y_u = pmin(1, prev + 1.96*se0)))
  }
  
  xgrid <- seq(max(min(p),1e-6), min(max(p),1-1e-6), length.out = grid_n)
  
  # Try LOESS with standard errors
  loess_try <- try({
    fit <- loess(y ~ p, span = span, family = "gaussian",
                 control = loess.control(surface = "direct"))
    pr  <- predict(fit, newdata = data.frame(p = xgrid), se = TRUE)
    list(y = pmin(pmax(pr$fit, 0), 1), se = pmax(pr$se.fit, 1e-12), ok = TRUE)
  }, silent = TRUE)
  
  if (!inherits(loess_try, "try-error") && isTRUE(loess_try$ok) &&
      all(is.finite(loess_try$y)) && all(is.finite(loess_try$se))){
    yhat <- loess_try$y; se <- loess_try$se
    return(tibble(x = xgrid, y = yhat,
                  y_l = pmax(0, yhat - 1.96*se),
                  y_u = pmin(1, yhat + 1.96*se)))
  }
  
  # Fallback: binomial GAM with logit link
  k <- min(10, max(4, floor(length(unique(p))/3) + 2))
  g  <- mgcv::gam(y ~ s(p, k = k), family = binomial(), method = "REML")
  pr <- predict(g, newdata = data.frame(p = xgrid), type = "link", se.fit = TRUE)
  eta <- pr$fit; se <- pmax(pr$se.fit, 1e-12)
  tibble(x = xgrid,
         y   = plogis(eta),
         y_l = plogis(eta - 1.96*se),
         y_u = plogis(eta + 1.96*se))
}

# ============== Nonparametric curve y(x) for bootstrap aggregation ==============
.nonparam_curve <- function(y, p, xgrid, span = 0.75){
  y <- as.integer(y); p <- .safe_prob(p)
  ok <- is.finite(y) & is.finite(p); y <- y[ok]; p <- p[ok]
  if (length(y) < 10 || length(unique(p)) < 2){
    return(rep(mean(y), length(xgrid)))
  }
  loess_try <- try({
    fit <- loess(y ~ p, span = span, family = "gaussian",
                 control = loess.control(surface = "direct"))
    pr  <- predict(fit, newdata = data.frame(p = xgrid), se = FALSE)
    pmin(pmax(as.numeric(pr), 0), 1)
  }, silent = TRUE)
  if (!inherits(loess_try, "try-error") && all(is.finite(loess_try))){
    return(loess_try)
  }
  k <- min(10, max(4, floor(length(unique(p))/3) + 2))
  g  <- mgcv::gam(y ~ s(p, k = k), family = binomial(), method = "REML")
  pr <- predict(g, newdata = data.frame(p = xgrid), type = "response")
  pmin(pmax(as.numeric(pr), 0), 1)
}

# ============== Frequency spikes (histogram-like rug) ==============
.freq_spikes <- function(p, bins = 40, max_height = 0.10){
  p <- .safe_prob(p)
  brks <- seq(0, 1, length.out = bins + 1)
  cnt  <- as.numeric(table(cut(p, breaks = brks, include.lowest = TRUE)))
  xmid <- (head(brks, -1) + tail(brks, -1))/2
  if (max(cnt) == 0) return(tibble(x = xmid, h = 0))
  tibble(x = xmid, h = cnt / max(cnt) * max_height)
}

# ============== Single calibration plot (raw or calibrated) ==============
# Overlays: Ideal / Logistic / LOESS / binned ± CI; frequency spikes at bottom
calib_plot_one <- function(y, p, setname, tag,
                           out_dir = "calib_plots", file_stub = NULL,
                           bins = 10, span = 0.75, dpi = 300, add_auc = TRUE,
                           metrics_pos = c("tr","tl","br","bl","none"),
                           metrics_alpha = 0.85, metrics_size = 3.6,
                           freq_bins = 40, freq_height = 0.10,
                           ci_alpha = 0.40) {
  metrics_pos <- match.arg(metrics_pos)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  if (is.null(file_stub)) file_stub <- gsub("\\s+","_", paste0(setname, "_", tag))
  
  y <- as.integer(y); p <- .safe_prob(p)
  
  # Metrics
  auc_val <- if (add_auc && length(unique(y)) == 2L) as.numeric(pROC::roc(y, p, quiet = TRUE)$auc) else NA_real_
  brier   <- mean((p - y)^2)
  ab_fit  <- glm(y ~ qlogis(p), family = binomial())
  a_hat   <- unname(coef(ab_fit)[1]); b_hat <- unname(coef(ab_fit)[2])
  
  df_bins <- .bin_summary(y, p, bins)
  df_lo   <- .smooth_band(y, p, span = span) |>
    mutate(series = "Nonparametric calibration (LOESS)",
           band   = "LOESS 95% CI")
  lc      <- logistic_cal_line(y, p)
  df_ideal <- tibble(x = c(0,1), y = c(0,1), series = "Ideal")
  df_freq  <- .freq_spikes(p, bins = freq_bins, max_height = freq_height)
  
  # Colors/linetypes
  col_map <- c("Ideal" = "grey40",
               "Logistic calibration" = "#2C7FB8",
               "Nonparametric calibration (LOESS)" = "#1F9E89",
               "Binned mean ± 95% CI" = "#E68613")
  lty_map  <- c("Ideal"="22","Logistic calibration"="solid","Nonparametric calibration (LOESS)"="solid")
  fill_map <- c("LOESS 95% CI" = "#C8E6DE")
  breaks_lines_points <- c("Ideal","Logistic calibration","Nonparametric calibration (LOESS)","Binned mean ± 95% CI")
  
  # Metrics panel position
  pos_tbl <- list(
    tl = list(x = 0.03, y = 0.97, h = 0, v = 1),
    tr = list(x = 0.97, y = 0.97, h = 1, v = 1),
    bl = list(x = 0.03, y = 0.03, h = 0, v = 0),
    br = list(x = 0.97, y = 0.03, h = 1, v = 0)
  )
  metrics_txt <- paste0(
    if (!is.na(auc_val)) paste0("AUC   ", sprintf('%.3f', auc_val), "\n") else "",
    "Brier   ",     sprintf('%.3f', brier),   "\n",
    "Intercept   ", sprintf('%.3f', a_hat),   "\n",
    "Slope   ",     sprintf('%.3f', b_hat),   "\n",
    {
      eavg <- mean(abs(df_bins$y_mean - df_bins$p_mean))
      emax <- max(abs(df_bins$y_mean - df_bins$p_mean))
      paste0("Eavg/Emax   ", sprintf('%.3f', eavg), "/", sprintf('%.3f', emax))
    }
  )
  
  plt <- ggplot() +
    # Frequency spikes
    geom_segment(data = df_freq,
                 aes(x = x, xend = x, y = 0, yend = h),
                 inherit.aes = FALSE, linewidth = 0.7, alpha = 0.55, colour = "grey45", show.legend = FALSE) +
    # LOESS band
    geom_ribbon(data = df_lo, aes(x = x, ymin = y_l, ymax = y_u, fill = band),
                alpha = ci_alpha, colour = NA, show.legend = TRUE) +
    # Lines
    geom_line(data = df_ideal,
              aes(x = x, y = y, colour = series, linetype = series), linewidth = 1.1, show.legend = TRUE) +
    geom_line(data = lc$df,
              aes(x = x, y = y, colour = series, linetype = series), linewidth = 1.35, show.legend = TRUE) +
    geom_line(data = df_lo,
              aes(x = x, y = y, colour = series, linetype = series), linewidth = 1.35, show.legend = TRUE) +
    # Binned mean ± CI: error bar + point (no connecting line)
    geom_errorbar(data = df_bins,
                  aes(x = p_mean, ymin = ci_l, ymax = ci_u, colour = "Binned mean ± 95% CI"),
                  width = 0.018, linewidth = 0.55, alpha = 0.95, show.legend = TRUE) +
    geom_point(data = df_bins,
               aes(x = p_mean, y = y_mean, colour = "Binned mean ± 95% CI"),
               size = 2.4, show.legend = TRUE) +
    # Legend formatting
    scale_color_manual(values = col_map, breaks = breaks_lines_points, name = NULL) +
    scale_linetype_manual(values = lty_map, name = NULL, guide = "none") +
    scale_fill_manual(values = fill_map,  name = NULL) +
    guides(
      color = guide_legend(
        order = 1,
        override.aes = list(
          shape    = c(NA, NA, NA, 16),             # only the binned key shows a dot
          linetype = c("22", "solid", "solid", "blank"),
          linewidth= c(1.1, 1.35, 1.35, NA),
          fill     = NA
        )
      ),
      fill = guide_legend(order = 2,
                          override.aes = list(
                            linetype = "blank",   # show a pure color tile (no dot)
                            shape = NA, colour = NA, alpha = ci_alpha
                          ))
    ) +
    coord_equal(xlim = c(0,1), ylim = c(0,1), expand = FALSE) +
    labs(
      title = paste0("Calibration: ", setname, " — ", tag),
      x = ifelse(grepl("^Raw", tag), "Predicted probability (raw)", "Predicted probability (calibrated)"),
      y = "Observed frequency"
    ) +
    theme_classic(base_size = 12) +
    theme(
      panel.border  = element_rect(colour = "black", fill = NA, linewidth = 0.7),
      plot.title    = element_text(face = "bold"),
      legend.key    = element_rect(fill = NA, colour = NA),
      legend.background = element_rect(fill = scales::alpha("white", 0.85), colour = NA),
      legend.key.width = unit(2.1, "lines")
    )
  
  if (metrics_pos != "none"){
    spec <- list(tl=list(x=0.03,y=0.97,h=0,v=1), tr=list(x=0.97,y=0.97,h=1,v=1),
                 bl=list(x=0.03,y=0.03,h=0,v=0), br=list(x=0.97,y=0.03,h=1,v=0))[[metrics_pos]]
    plt <- plt +
      annotate("label", x = spec$x, y = spec$y, hjust = spec$h, vjust = spec$v,
               label = metrics_txt, size = metrics_size,
               fill = scales::alpha("white", metrics_alpha), label.size = 0)
  }
  
  f_pdf <- file.path(out_dir, paste0(file_stub, ".pdf"))
  ggsave(f_pdf, plt, width = 7.2, height = 6.2, dpi = dpi,
         device = if (capabilities("cairo")) cairo_pdf else pdf)
  message("Saved: ", f_pdf)
  invisible(list(file = f_pdf, plot = plt))
}

# ============== Train-only bootstrap aggregation band ==============
plot_train_bootstrap_band <- function(y, p, out_dir = "calib_plots",
                                      file_stub = "Train_calibration_bootstrap_band",
                                      B = 200, bins = 10, span = 0.75, dpi = 300,
                                      add_auc = TRUE, metrics_pos = "tr",
                                      metrics_alpha = 0.85, metrics_size = 3.6,
                                      freq_bins = 40, freq_height = 0.10,
                                      ci_alpha = 0.40){
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  y <- as.integer(y); p <- .safe_prob(p)
  xgrid <- seq(max(min(p),1e-6), min(max(p),1-1e-6), length.out = 200)
  
  n <- length(y)
  Yboot <- replicate(B, {
    id <- sample.int(n, n, replace = TRUE)
    .nonparam_curve(y[id], p[id], xgrid = xgrid, span = span)
  })
  
  qlo <- apply(Yboot, 1, quantile, probs = 0.05, na.rm = TRUE)
  qmd <- apply(Yboot, 1, quantile, probs = 0.50, na.rm = TRUE)
  qhi <- apply(Yboot, 1, quantile, probs = 0.95, na.rm = TRUE)
  
  df_band <- tibble(x = xgrid, y = qmd, y_l = qlo, y_u = qhi,
                    series = "Bootstrap median", band = "Bootstrap 90% CI")
  
  df_bins  <- .bin_summary(y, p, bins)
  lc       <- logistic_cal_line(y, p)
  df_ideal <- tibble(x = c(0,1), y = c(0,1), series = "Ideal")
  df_freq  <- .freq_spikes(p, bins = freq_bins, max_height = freq_height)
  
  auc_val <- if (add_auc && length(unique(y)) == 2L) as.numeric(pROC::roc(y, p, quiet = TRUE)$auc) else NA_real_
  brier   <- mean((p - y)^2)
  ab_fit  <- glm(y ~ qlogis(p), family = binomial())
  a_hat   <- unname(coef(ab_fit)[1]); b_hat <- unname(coef(ab_fit)[2])
  eavg <- mean(abs(df_bins$y_mean - df_bins$p_mean)); emax <- max(abs(df_bins$y_mean - df_bins$p_mean))
  metrics_txt <- paste0(
    if (!is.na(auc_val)) paste0("AUC   ", sprintf('%.3f', auc_val), "\n") else "",
    "Brier   ", sprintf('%.3f', brier), "\n",
    "Intercept   ", sprintf('%.3f', a_hat), "\n",
    "Slope   ", sprintf('%.3f', b_hat), "\n",
    "Eavg/Emax   ", sprintf('%.3f', eavg), "/", sprintf('%.3f', emax)
  )
  
  col_map <- c("Ideal" = "grey40",
               "Logistic calibration" = "#2C7FB8",
               "Bootstrap median" = "#1F9E89",
               "Binned mean ± 95% CI" = "#E68613")
  lty_map  <- c("Ideal"="22","Logistic calibration"="solid","Bootstrap median"="solid")
  fill_map <- c("Bootstrap 90% CI" = "#C8E6DE")
  breaks_lines_points <- c("Ideal","Logistic calibration","Bootstrap median","Binned mean ± 95% CI")
  
  plt <- ggplot() +
    geom_segment(data = df_freq, aes(x = x, xend = x, y = 0, yend = h),
                 inherit.aes = FALSE, linewidth = 0.7, alpha = 0.55, colour = "grey45", show.legend = FALSE) +
    geom_ribbon(data = df_band, aes(x = x, ymin = y_l, ymax = y_u, fill = band),
                alpha = ci_alpha, colour = NA, show.legend = TRUE) +
    geom_line(data = df_ideal, aes(x = x, y = y, colour = series, linetype = series), linewidth = 1.1, show.legend = TRUE) +
    geom_line(data = lc$df,   aes(x = x, y = y, colour = series, linetype = series), linewidth = 1.35, show.legend = TRUE) +
    geom_line(data = df_band, aes(x = x, y = y, colour = series, linetype = series), linewidth = 1.35, show.legend = TRUE) +
    geom_errorbar(data = df_bins, aes(x = p_mean, ymin = ci_l, ymax = ci_u, colour = "Binned mean ± 95% CI"),
                  width = 0.018, linewidth = 0.55, alpha = 0.95, show.legend = TRUE) +
    geom_point(data = df_bins, aes(x = p_mean, y = y_mean, colour = "Binned mean ± 95% CI"),
               size = 2.4, show.legend = TRUE) +
    scale_color_manual(values = col_map, breaks = breaks_lines_points, name = NULL) +
    scale_linetype_manual(values = lty_map, name = NULL, guide = "none") +
    scale_fill_manual(values = fill_map, name = NULL) +
    guides(
      color = guide_legend(
        order = 1,
        override.aes = list(
          shape    = c(NA, NA, NA, 16),
          linetype = c("22", "solid", "solid", "blank"),
          linewidth= c(1.1, 1.35, 1.35, NA),
          fill     = NA
        )
      ),
      fill = guide_legend(order = 2,
                          override.aes = list(
                            linetype = "blank", shape = NA, colour = NA, alpha = ci_alpha
                          ))
    ) +
    coord_equal(xlim = c(0,1), ylim = c(0,1), expand = FALSE) +
    labs(title = "Calibration: Train — Bootstrap aggregated curve",
         x = "Predicted probability", y = "Observed frequency") +
    theme_classic(base_size = 12) +
    theme(
      panel.border  = element_rect(colour = "black", fill = NA, linewidth = 0.7),
      plot.title    = element_text(face = "bold"),
      legend.key    = element_rect(fill = NA, colour = NA),
      legend.background = element_rect(fill = scales::alpha("white", 0.85), colour = NA),
      legend.key.width = unit(2.1, "lines")
    ) +
    annotate("label", x = 0.97, y = 0.97, hjust = 1, vjust = 1,
             label = metrics_txt, size = metrics_size,
             fill = scales::alpha("white", metrics_alpha), label.size = 0)
  
  f_pdf <- file.path(out_dir, paste0(file_stub, ".pdf"))
  ggsave(f_pdf, plt, width = 7.2, height = 6.2, dpi = dpi,
         device = if (capabilities("cairo")) cairo_pdf else pdf)
  message("Saved: ", f_pdf)
  invisible(list(file = f_pdf, plot = plt))
}

# ============== Orchestrator: fit calibrator + plot Train/Test/Valid ==============
calibrate_and_plot_from_preds <- function(
    y_train, p_train,
    y_test,  p_test,
    y_valid, p_valid,
    cal_method = c("isotonic","platt"),
    iso_interp = c("linear","constant"),
    out_dir = "calib_plots",
    bins = 10, span = 0.75, dpi = 300, add_auc = TRUE,
    n_boot_train = 0, freq_bins = 40, freq_height = 0.10,
    ci_alpha = 0.40
){
  cal_method <- match.arg(cal_method)
  iso_interp <- match.arg(iso_interp)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # Fit calibrator on Train predictions (prefer OOF/OOB if available)
  if (cal_method == "platt") {
    cal <- platt_fit(y_train, p_train)
    apply_cal <- function(p) platt_apply(p, cal)
    tag_cal <- "Calibrated (Platt from Train)"
  } else {
    cal <- isotonic_fit(y_train, p_train)
    apply_cal <- function(p) isotonic_apply(p, cal, iso_interp)
    tag_cal <- "Calibrated (Isotonic from Train)"
  }
  
  # Apply to each split
  p_test_cal  <- .safe_prob(apply_cal(p_test))
  p_valid_cal <- .safe_prob(apply_cal(p_valid))
  p_train_cal <- .safe_prob(apply_cal(p_train))
  
  # Train: raw + calibrated
  calib_plot_one(y_train, .safe_prob(p_train), "Train", "Raw",
                 out_dir = out_dir, file_stub = "Train_calibration_raw",
                 bins = bins, span = span, dpi = dpi, add_auc = add_auc,
                 freq_bins = freq_bins, freq_height = freq_height, ci_alpha = ci_alpha)
  calib_plot_one(y_train, p_train_cal, "Train", tag_cal,
                 out_dir = out_dir, file_stub = "Train_calibration_calibrated",
                 bins = bins, span = span, dpi = dpi, add_auc = add_auc,
                 freq_bins = freq_bins, freq_height = freq_height, ci_alpha = ci_alpha)
  
  # Optional: bootstrap-aggregated nonparametric band on Train
  if (n_boot_train > 0){
    plot_train_bootstrap_band(y_train, .safe_prob(p_train),
                              out_dir = out_dir,
                              file_stub = "Train_calibration_bootstrap_band",
                              B = n_boot_train, bins = bins, span = span, dpi = dpi,
                              add_auc = add_auc,
                              freq_bins = freq_bins, freq_height = freq_height,
                              ci_alpha = ci_alpha)
  }
  
  # Test: raw + calibrated
  calib_plot_one(y_test,  .safe_prob(p_test),  "Test", "Raw",
                 out_dir = out_dir, file_stub = "Test_calibration_raw",
                 bins = bins, span = span, dpi = dpi, add_auc = add_auc,
                 freq_bins = freq_bins, freq_height = freq_height, ci_alpha = ci_alpha)
  calib_plot_one(y_test,  p_test_cal,          "Test", tag_cal,
                 out_dir = out_dir, file_stub = "Test_calibration_calibrated",
                 bins = bins, span = span, dpi = dpi, add_auc = add_auc,
                 freq_bins = freq_bins, freq_height = freq_height, ci_alpha = ci_alpha)
  
  # Validation: raw + calibrated
  calib_plot_one(y_valid, .safe_prob(p_valid), "Validation", "Raw",
                 out_dir = out_dir, file_stub = "Valid_calibration_raw",
                 bins = bins, span = span, dpi = dpi, add_auc = add_auc,
                 freq_bins = freq_bins, freq_height = freq_height, ci_alpha = ci_alpha)
  calib_plot_one(y_valid, p_valid_cal,         "Validation", tag_cal,
                 out_dir = out_dir, file_stub = "Valid_calibration_calibrated",
                 bins = bins, span = span, dpi = dpi, add_auc = add_auc,
                 freq_bins = freq_bins, freq_height = freq_height, ci_alpha = ci_alpha)
  
  # Export metrics table
  metric_one <- function(y, p, set){
    auc_val <- if (length(unique(y))==2L) as.numeric(pROC::roc(y, .safe_prob(p), quiet=TRUE)$auc) else NA_real_
    brier   <- mean((.safe_prob(p) - y)^2)
    ab      <- glm(as.integer(y) ~ qlogis(.safe_prob(p)), family = binomial())
    ee      <- .eavg_emax(y, .safe_prob(p), bins = bins)
    tibble(Set=set, AUC=auc_val, Brier=brier,
           Intercept=unname(coef(ab)[1]), Slope=unname(coef(ab)[2]),
           Eavg=ee["Eavg"], Emax=ee["Emax"])
  }
  metrics <- bind_rows(
    metric_one(y_train, p_train,     "Train Raw"),
    metric_one(y_train, p_train_cal, "Train Calibrated"),
    metric_one(y_test,  p_test,      "Test Raw"),
    metric_one(y_test,  p_test_cal,  "Test Calibrated"),
    metric_one(y_valid, p_valid,     "Validation Raw"),
    metric_one(y_valid, p_valid_cal, "Validation Calibrated")
  )
  write_csv(metrics, file.path(out_dir, "metrics_train_test_valid.csv"))
  
  invisible(list(calibrator = cal, metrics = metrics))
}

# ======= RF / ranger adapters (optional) =======
.pick_prob_col <- function(prob_mat, pos_level){
  if (is.null(colnames(prob_mat))) {
    if (is.numeric(pos_level)) prob_mat[, pos_level]
    else stop("Probability matrix has no column names; please pass numeric pos_level (1/2).")
  } else {
    if (is.numeric(pos_level)) prob_mat[, pos_level] else prob_mat[, pos_level]
  }
}
rf_oob_prob_randomForest <- function(rf_model, pos_level){
  .pick_prob_col(rf_model$votes, pos_level)
}
rf_predict_prob_randomForest <- function(rf_model, newdata, pos_level){
  pr <- predict(rf_model, newdata = newdata, type = "prob")
  .pick_prob_col(pr, pos_level)
}
rgr_oob_prob_ranger <- function(rgr_model, pos_level){
  .pick_prob_col(rgr_model$predictions, pos_level)
}
rgr_predict_prob_ranger <- function(rgr_model, newdata, pos_level){
  pr <- predict(rgr_model, data = newdata)$predictions
  .pick_prob_col(pr, pos_level)
}

# ======= LASSO (cv.glmnet) adapters (optional) =======
lasso_oof_prob_cvglmnet <- function(cvfit, s = c("lambda.1se","lambda.min")){
  s <- match.arg(s)
  if (is.null(cvfit$fit.preval))
    stop("cv.glmnet must be trained with keep=TRUE to expose OOF probabilities (fit.preval).")
  lam <- if (s == "lambda.min") cvfit$lambda.min else cvfit$lambda.1se
  j <- which.min(abs(cvfit$lambda - lam))
  p <- cvfit$fit.preval[, j]
  if (any(p < 0 | p > 1, na.rm = TRUE)) p <- plogis(p)
  .safe_prob(as.numeric(p))
}
lasso_predict_prob_cvglmnet <- function(cvfit, newx, s = c("lambda.1se","lambda.min")){
  s <- match.arg(s)
  as.numeric(.safe_prob(predict(cvfit, newx = newx, s = s, type = "response")))
}
lasso_predict_prob_glmnet <- function(fit, newx, lambda){
  as.numeric(.safe_prob(predict(fit, newx = newx, s = lambda, type = "response")))
}

p_train_lasso <- as.numeric(predictions_train_re)
p_test_lasso  <- as.numeric(predictions_test_re)
p_valid_lasso <- as.numeric(predictions_valid_re)

# Isotonic
calibrate_and_plot_from_preds(
  y_train = y_train_re, p_train = p_train_lasso,
  y_test  = y_test_re,  p_test  = p_test_lasso,
  y_valid = y_valid_re, p_valid = p_valid_lasso,
  cal_method   = "isotonic", iso_interp = "linear",
  out_dir      = "calibration_output_lasso_iso",
  bins = 10, span = 0.8, dpi = 300, add_auc = TRUE,
  n_boot_train = 1000, ci_alpha = 0.40
)

# Platt
calibrate_and_plot_from_preds(
  y_train = y_train_re, p_train = p_train_lasso,
  y_test  = y_test_re,  p_test  = p_test_lasso,
  y_valid = y_valid_re, p_valid = p_valid_lasso,
  cal_method   = "platt",
  out_dir      = "calibration_output_lasso_platt",
  bins = 10, span = 0.8, dpi = 300, add_auc = TRUE,
  n_boot_train = 1000, ci_alpha = 0.40
)