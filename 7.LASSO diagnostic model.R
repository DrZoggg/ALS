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
