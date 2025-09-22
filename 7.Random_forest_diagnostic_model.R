setwd('Diag_models/RF')

# Packages
library(vegan)        # Ecological and multivariate statistics
library(ggplot2)      # Visualization
library(RColorBrewer) # Color palettes
library(tidyverse)    # Data wrangling & visualization

# Datasets
load("../../Data_preprocessing/train_cohort.rda")  # train_data
load("../../Data_preprocessing/test_cohort.rda")   # test_data
load("../../Data_preprocessing/valid_cohort.rda")  # validation_data
genes <- readRDS("../../ML_features_select/common_genes.rds")

# -------------------------------
# Build RF on the training set
# -------------------------------
diagnosis <- train_cohort_clin$diagnosis
datExpr   <- as.data.frame(cbind(diagnosis, t(train_cohort_expr)))
datExpr   <- datExpr[, c('diagnosis', genes)]
str(datExpr)

# Rename the first column to 'group'
colnames(datExpr)[1] <- 'group'
datExpr$group <- as.factor(datExpr$group)

# Convert expression columns (except the first) to numeric
datExpr[, -1] <- lapply(datExpr[, -1], function(x) as.numeric(as.character(x)))

# Verify structure after conversion
str(datExpr)

# Format for RF input:
#   rows = samples, columns = genes
#   first column = group label; remaining columns = expression matrix
x <- as.matrix(datExpr[, -1]) %>% scale()
y <- datExpr[, 1]

library(randomForest)
rf_model <- randomForest(
  x,
  y,
  ntree = 500,
  mtry = 3,
  importance = TRUE,   # compute variable importance
  type = "classification"
)
rf_model
saveRDS(rf_model, file = 'rf_model.rds')

# -------------------------------
# Error rate vs. number of trees
# -------------------------------
library(ggplot2)

# Prepare long-format error-rate data for ggplot
err_rate_long <- data.frame(
  trees      = rep(1:nrow(rf_model$err.rate), 3),
  error_rate = c(rf_model$err.rate[, 1],
                 rf_model$err.rate[, 2],
                 rf_model$err.rate[, 3]),
  type       = rep(c("OOB", "ALS", "CON"), each = nrow(rf_model$err.rate))
)

# Line plot of error rate over trees
p <- ggplot(err_rate_long, aes(x = trees, y = error_rate, color = type)) +
  geom_line(alpha = 1) +  # full opacity
  labs(title = "Error Rate vs. Number of Trees",
       x = "Number of Trees", y = "Error Rate") +
  theme_minimal() +
  theme(
    legend.position = "topright",                         # place legend (kept as-is)
    panel.border    = element_rect(colour = "black", fill = NA, size = 1),
    legend.title    = element_blank(),
    legend.key      = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks      = element_line(color = "black")
  ) +
  scale_color_manual(values = c("#FF6A6A", "#4682B4", "black")) +  # custom colors
  scale_x_continuous(breaks = seq(0, 500, by = 50)) +              # keep x breaks
  scale_y_continuous(limits = c(0, 0.45),
                     breaks = seq(0, 0.45, by = 0.05)) +
  guides(color = guide_legend(
    title = "Type",
    override.aes = list(size = 1.5, linetype = "solid")
  )) +
  geom_text(
    data = err_rate_long[which(err_rate_long$trees == max(err_rate_long$trees)), ],
    aes(x = trees, y = error_rate, label = type),
    color = c("black", "#FF6A6A", "#4682B4"),
    hjust = 2, vjust = -2
  )
p
ggsave(p, file = 'RFmodel_trees数目与错误率关系折线图.pdf', height = 5, width = 5)

# -------------------------------
# Variable importance (base plot)
# -------------------------------
varImpPlot(rf_model)

# -------------------------------
# Variable importance (ggplot + patchwork)
# -------------------------------
library(randomForest)
library(ggplot2)
library(dplyr)
library(patchwork)

# 1) Extract importance and convert to data frame
imp <- importance(rf_model)
imp_df <- data.frame(
  Feature             = rownames(imp),
  MeanDecreaseAccuracy = imp[, "MeanDecreaseAccuracy"],
  MeanDecreaseGini     = imp[, "MeanDecreaseGini"]
)

# 2) Sort for plotting
imp_df_mda <- imp_df %>% arrange(MeanDecreaseAccuracy)  # ascending by MDA
imp_df_mdg <- imp_df %>% arrange(MeanDecreaseGini)      # ascending by MDG

# 3) Two bar charts
p1 <- ggplot(imp_df_mda,
             aes(x = reorder(Feature, MeanDecreaseAccuracy),
                 y = MeanDecreaseAccuracy)) +
  geom_col(fill = "#9AC0CD") +
  coord_flip() +
  labs(title = NULL, x = NULL, y = "Mean Decrease Accuracy") +
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    panel.grid   = element_blank(),
    axis.text    = element_text(color = "black", face = "bold"),
    axis.title   = element_text(color = "black", face = "bold")
  )

p2 <- ggplot(imp_df_mdg,
             aes(x = reorder(Feature, MeanDecreaseGini),
                 y = MeanDecreaseGini)) +
  geom_col(fill = '#FA8072') +
  coord_flip() +
  labs(title = NULL, x = NULL, y = "Mean Decrease Gini") +
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    panel.grid   = element_blank(),
    axis.text    = element_text(color = "black", face = "bold"),
    axis.title   = element_text(color = "black", face = "bold")
  )

# 4) Compose with a global title
final_plot <- p1 + p2 +
  plot_annotation(
    title = "Random Forest Model",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  )

print(final_plot)
ggsave(final_plot, file = '随机森林模型特征重要性.pdf', height = 5, width = 5)

# -------------------------------
# ROC/AUC — Train set
# -------------------------------
# Predict class probabilities (ALS vs CON)
predictions_prob <- predict(rf_model, x, type = "prob")

library(pROC)
roc_curve <- roc(y, predictions_prob[, 'ALS'], levels = c("CON", "ALS"))
plot(roc_curve)
auc(roc_curve)
train_ROC <- roc_curve  # save train ROC

# -------------------------------
# ROC/AUC — Test set
# -------------------------------
diagnosis <- test_cohort_clin$diagnosis
datExpr   <- as.data.frame(cbind(diagnosis, t(test_cohort_exp)))
datExpr   <- datExpr[, c('diagnosis', genes)]
str(datExpr)

# Rename the first column to 'group'
colnames(datExpr)[1] <- 'group'
datExpr$group <- as.factor(datExpr$group)

# Convert expression columns to numeric
datExpr[, -1] <- lapply(datExpr[, -1], function(x) as.numeric(as.character(x)))

# Verify structure
str(datExpr)

# Reuse training scaling parameters
x_test <- scale(datExpr[-1],
                center = attr(x, "scaled:center"),
                scale  = attr(x, "scaled:scale"))
y_test <- datExpr[, 1]

# Predict class probabilities and compute ROC/AUC
predictions_prob <- predict(rf_model, x_test, type = "prob")
library(pROC)
roc_curve <- roc(y_test, predictions_prob[, 'ALS'], levels = c("CON", "ALS"))
plot(roc_curve)
auc(roc_curve)
test_ROC <- roc_curve

# -------------------------------
# ROC/AUC — External validation
# -------------------------------
diagnosis <- valid_cohort_clin$diagnosis
datExpr   <- as.data.frame(cbind(diagnosis, t(valid_cohort_expr)))
datExpr   <- datExpr[, c('diagnosis', genes)]
str(datExpr)

# Rename the first column to 'group'
colnames(datExpr)[1] <- 'group'
datExpr$group <- as.factor(datExpr$group)

# Convert expression columns to numeric
datExpr[, -1] <- lapply(datExpr[, -1], function(x) as.numeric(as.character(x)))

# Keep ALS and CON only
datExpr_filtered <- datExpr[datExpr$group %in% c("ALS", "CON"), ]
str(datExpr_filtered)

# Reuse training scaling parameters
x_val <- scale(datExpr_filtered[-1],
               center = attr(x, "scaled:center"),
               scale  = attr(x, "scaled:scale"))
y_val <- datExpr_filtered[, 1]

# Predict class probabilities and compute ROC/AUC
predictions_prob <- predict(rf_model, x_val, type = "prob")
library(pROC)
roc_curve <- roc(y_val, predictions_prob[, 'ALS'], levels = c("CON", "ALS"))
plot(roc_curve)
auc(roc_curve)
validation_ROC <- roc_curve

# -------------------------------
# Plot combined ROC curves
# -------------------------------
library(pROC)
library(ggplot2)
library(dplyr)

# Extract ROC data
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

# Combine
roc_data <- bind_rows(train_data, test_data, validation_data)
options(repr.plot.width = 5, repr.plot.height = 5)

# AUCs
train_auc      <- auc(train_ROC)
test_auc       <- auc(test_ROC)
validation_auc <- auc(validation_ROC)

# Legend labels with AUC
labels <- c(
  paste("Train (ROC =", round(train_auc, 3), ")"),
  paste("Test (ROC =", round(test_auc, 3), ")"),
  paste("Validation (ROC =", round(validation_auc, 3), ")")
)

# Keep legend order
roc_data$dataset <- factor(roc_data$dataset, levels = c("Train", "Test", "Validation"))

# ROC plot (specificity on x, sensitivity on y)
p <- ggplot(roc_data, aes(x = specificity, y = sensitivity, color = dataset)) +
  geom_line(size = 1.2) +
  geom_segment(aes(x = 0, y = 1, xend = 1, yend = 0),
               linetype = "dashed", color = "gray") +
  labs(
    title    = "Receiver Operating Characteristic (ROC)",
    subtitle = "ALS vs CON",
    x = "Specificity",
    y = "Sensitivity",
    color = "Random Forest Model:"
  ) +
  scale_x_reverse(limits = c(1, 0), breaks = seq(0, 1, 0.25), expand = c(0.02, 0.02)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0.02, 0.02)) +
  scale_color_manual(values = c("Train" = "#CDB38B", "Test" = "#8DB6CD", "Validation" = "#CD5555"),
                     labels = labels) +
  theme_minimal() +
  theme(
    plot.title       = element_text(hjust = 0.5, face = "bold", size = 10),
    plot.subtitle    = element_text(hjust = 0.5, face = "bold", size = 8, color = "black"),
    axis.title       = element_text(face = "bold", size = 8),
    axis.text        = element_text(color = "black", size = 6),
    axis.ticks       = element_line(size = 1),
    axis.ticks.length = unit(-0.2, "cm"),
    panel.border     = element_rect(color = "black", fill = NA, size = 1),
    panel.grid       = element_blank(),
    legend.title     = element_text(face = "bold", size = 8, color = "black"),
    legend.text      = element_text(face = "bold", size = 8, color = "black"),
    legend.position  = c(0.80, 0.15)
  )

print(p)
ggsave(p, file = 'RF_ROC.Plot.pdf', height = 5, width = 5)

# -------------------------------
# Summarize ROC metrics
# -------------------------------
library(pROC)
library(dplyr)

# Extract metrics from a pROC::roc object
metrics_from_roc <- function(roc_obj, dataset) {
  auc_val <- as.numeric(auc(roc_obj))
  best <- coords(
    roc_obj, x = "best", best.method = "youden",
    ret = c("threshold", "sensitivity", "specificity", "accuracy"),
    transpose = FALSE
  )
  
  n_cases    <- length(roc_obj$cases)
  n_controls <- length(roc_obj$controls)
  
  tibble(
    dataset            = dataset,
    N                  = n_cases + n_controls,
    `Case proportion`  = n_cases / (n_cases + n_controls),  # renamed
    AUC                = auc_val,
    Accuracy           = as.numeric(best["accuracy"]),
    Sensitivity        = as.numeric(best["sensitivity"]),
    Specificity        = as.numeric(best["specificity"]),
    Threshold          = as.numeric(best["threshold"]),
    Direction          = roc_obj$direction
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
readr::write_csv(metrics_tbl, "RF_metrics_summary.csv")
