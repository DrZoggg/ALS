setwd("ML_features_select/RF_features_filter")

library(vegan)        # Ecological and multivariate statistical analyses
library(ggplot2)      # Data visualization
library(RColorBrewer) # Color palettes
library(tidyverse)    # Data wrangling & visualization toolbox

### Up- and down-regulated module genes
load("../DEA&EA/up&down_regulated_module_DEGs.rda")

### Load WGCNA input matrix (only one outlier sample was filtered out)
load("../WGCNA/Step01-WGCNA_input.Rda")

### Clinical information
load("../Data_preprocessing/disc_cohort.rda")

all_dif_module <- c(down_dif_express_module, up_dif_express_module)
diagnosis <- disc_cohort_clin[rownames(datExpr), ]$diagnosis
datExpr <- as.data.frame(cbind(diagnosis, datExpr))
datExpr <- datExpr[, c('diagnosis', all_dif_module)]
str(datExpr)

# Set the first column name to 'group'
colnames(datExpr)[1] <- 'group'
datExpr$group <- as.factor(datExpr$group)

# Convert all expression columns (except the first) to numeric
datExpr[, -1] <- lapply(datExpr[, -1], function(x) as.numeric(as.character(x)))
datExpr[, -1] <- as.matrix(datExpr[, -1]) %>% scale()

# Inspect structure to confirm the conversion
str(datExpr)

data_pca <- datExpr

# Split dataset
set.seed(123)  # Set seed for reproducibility
sample_indices <- sample(1:nrow(data_pca), nrow(data_pca) * 0.3)  # Randomly sample 30% rows as test set
test_data  <- data_pca[sample_indices, ]        # Extract test set by indices
train_data <- data_pca[-sample_indices, ]       # Remaining rows as training set

library(randomForest)
rf_model <- randomForest(
  x = train_data[, -1],
  y = train_data[, 1],
  ntree = 500,
  mtry = 21,
  importance = TRUE,    # Compute variable importance
  type = "classification"
)
rf_model

library(ggplot2)

# Assume you already have data frame `err_rate`
# Convert to long format for ggplot
err_rate_long <- data.frame(
  trees = rep(1:nrow(rf_model$err.rate), 3),
  error_rate = c(rf_model$err.rate[, 1], rf_model$err.rate[, 2], rf_model$err.rate[, 3]),
  type = rep(c("OOB", "ALS", "CON"), each = nrow(rf_model$err.rate))
)

# Plot error rate vs. number of trees
p <- ggplot(err_rate_long, aes(x = trees, y = error_rate, color = type)) +
  geom_line(alpha = 1) +  # Full opacity
  labs(title = "Error Rate vs. Number of Trees", x = "Number of Trees", y = "Error Rate") +
  theme_minimal() +
  theme(
    legend.position = "topright",                           # Place legend at top-right
    panel.border = element_rect(colour = "black", fill = NA, size = 1),  # Add border
    legend.title = element_blank(),                         # Remove legend title
    legend.key = element_rect(fill = "white"),              # White legend background
    panel.grid.major = element_blank(),                     # Remove major grid
    panel.grid.minor = element_blank(),                     # Remove minor grid
    axis.ticks = element_line(color = "black")              # Keep axis ticks
  ) +
  scale_color_manual(values = c("#FF6A6A", "#4682B4", "black")) +  # Manual colors
  scale_x_continuous(breaks = seq(0, 500, by = 50)) +              # Keep x-axis breaks
  scale_y_continuous(limits = c(0, 0.45), breaks = seq(0, 0.45, by = 0.05)) +  # y-axis range & breaks
  guides(color = guide_legend(
    title = "Type",  # Legend title
    override.aes = list(size = 1.5, linetype = "solid")  # Legend lines
  )) +
  geom_text(
    data = err_rate_long[which(err_rate_long$trees == max(err_rate_long$trees)), ],
    aes(x = trees, y = error_rate, label = type),
    color = c("black", "#FF6A6A", "#4682B4"),  # Match colors to each type
    hjust = 2,  # Adjust text position (horizontal)
    vjust = -2  # Adjust text position (vertical)
  )
p
ggsave(p, file = 'RFmodel_trees数目与错误率关系折线图.pdf', height = 5, width = 5)

# Plot variable importance of the Random Forest model
varImpPlot(rf_model)

# Predict class probabilities on the test set (type = "prob")
predictions_prob <- predict(rf_model, test_data, type = "prob")

library(pROC)

# Assume "ALS" is the positive class and "CON" is the negative class
roc_curve <- roc(test_data$group, predictions_prob[, 'ALS'], levels = c("CON", "ALS"))

# Find the best cutoff
best_cutoff <- coords(roc_curve, "best", ret = c("threshold", "sensitivity", "specificity"))

pdf(file = 'RF.ROC曲线.pdf', height = 5, width = 5)
# Draw ROC curve
plot(roc_curve, main = "ROC Curve of Random Forest")

# Add dashed lines at the best threshold
abline(h = best_cutoff[2], col = "red",  lty = 2)  # Sensitivity reference
abline(v = best_cutoff[3], col = "blue", lty = 2)  # Specificity reference

# Add annotations for Sens/Spec next to dashed lines
text(0.9, 0.75, paste("Sens:", round(best_cutoff[2], 2)), col = "red",  pos = 4, cex = 1)
text(0.9, 0.70, paste("Spec:", round(best_cutoff[3], 2)), col = "blue", pos = 4, cex = 1)

# Compute and display AUC
auc_value <- auc(roc_curve)
text(0.9, 0.65, paste("AUC:", round(auc_value, 3)), pos = 4, cex = 1)
dev.off()

# Draw ROC curve again (on screen)
plot(roc_curve, main = "ROC Curve of Random Forest")
abline(h = best_cutoff[2], col = "red",  lty = 2)
abline(v = best_cutoff[3], col = "blue", lty = 2)
text(0.9, 0.75, paste("Sens:", round(best_cutoff[2], 2)), col = "red",  pos = 4, cex = 1)
text(0.9, 0.70, paste("Spec:", round(best_cutoff[3], 2)), col = "blue", pos = 4, cex = 1)
auc_value <- auc(roc_curve)
text(0.9, 0.65, paste("AUC:", round(auc_value, 3)), pos = 4, cex = 1)

library(rfPermute)  # Evaluate significance of RF variable importance

# Run rfPermute to compute variable-importance p-values
set.seed(123)  # Random seed
factor_rfP <- rfPermute(
  x = train_data[, -1],
  y = train_data[, 1],
  ntree = 500,
  mtry = 21,
  importance = TRUE,    # Compute variable importance
  type = "classification",
  nrep = 1000,          # Permutation times for p-values
  num.cores = 24        # Parallel CPU cores
)
factor_rfP  # Print rfPermute results

# Extract and scale variable importance
importance_factor.scale <- data.frame(
  importance(factor_rfP, scale = TRUE),  # Scaled importance
  check.names = FALSE                    # Keep column names unchanged
)
head(importance_factor.scale)  # Preview scores

# Build importance table (Mean Decrease Accuracy)
importance_df <- data.frame(
  Feature   = rownames(importance_factor.scale),
  Importance = importance_factor.scale$MeanDecreaseAccuracy,
  p_value    = importance_factor.scale$MeanDecreaseAccuracy.pval
)

library(ggplot2)

# Order by MDA and mark p < 0.05; pick top 20% threshold
importance_df <- importance_df[order(importance_df$Importance, decreasing = TRUE), ]
top_20_percent <- round(0.2 * nrow(importance_df))
importance_df$highlight <- ifelse(importance_df$p_value < 0.05, "#FFC1C1", "gray")  # Bar color

# Create threshold line value
threshold_value <- importance_df$Importance[top_20_percent]

# Visualize MDA importance
p1 <- ggplot(importance_df, aes(x = reorder(Feature, Importance), y = Importance, fill = highlight)) +
  geom_bar(stat = "identity") +
  coord_flip() +  # Flip coordinates
  labs(
    title = "Feature Importance (MDA)",
    y = "Mean Decrease in Accuracy",
    x = paste("Features (n=", nrow(importance_df), ")", sep = "")
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),  # Hide y-axis labels
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Solid border
    legend.position = "none"  # No legend
  ) +
  geom_hline(yintercept = threshold_value, linetype = "dashed", color = "grey", size = 1) +  # Threshold line
  scale_fill_identity() +  # Use fill colors as-is
  scale_y_continuous(breaks = seq(0, max(importance_df$Importance), by = 1))  # y-axis ticks
p1
ggsave(p1, file = 'RFmodel_449featuresMDA.pdf', height = 5, width = 5)

# Build importance table (Gini decrease)
gini_df <- data.frame(
  Feature = rownames(importance_factor.scale),
  Gini    = importance_factor.scale$MeanDecreaseGini,
  p_value = importance_factor.scale$MeanDecreaseGini.pval
)

library(ggplot2)

# Order by Gini and mark p < 0.05; pick top 20% threshold
gini_df <- gini_df[order(gini_df$Gini, decreasing = TRUE), ]
top_20_percent <- round(0.2 * nrow(gini_df))
gini_df$highlight <- ifelse(gini_df$p_value < 0.05, "#FFC1C1", "gray")  # Bar color

# Create threshold line value
threshold_value_gini <- gini_df$Gini[top_20_percent]

# Visualize Gini importance
p1 <- ggplot(gini_df, aes(x = reorder(Feature, Gini), y = Gini, fill = highlight)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(
    title = "Feature Importance (Gini Index)",
    y = "Gini Index",
    x = paste("Features (n=", nrow(gini_df), ")", sep = "")
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    legend.position = "none"
  ) +
  geom_hline(yintercept = threshold_value_gini, linetype = "dashed", color = "grey", size = 1) +
  scale_fill_identity() +
  scale_y_continuous(breaks = seq(0, max(gini_df$Gini), by = 1))
p1
ggsave(p1, file = 'RFmodel_449featuresGiniIndex.pdf', height = 5, width = 5)


################## Extract genes that meet criteria (top 20% by importance, p-value < 0.05)

# Sort by MDA in descending order
importance_df_sorted <- importance_df[order(importance_df$Importance, decreasing = TRUE), ]

# Compute number of genes in the top 20%
top_20_percent_importance <- round(0.2 * nrow(importance_df_sorted))

# Keep top 20% and filter p < 0.05
importance_df_filtered <- importance_df_sorted[1:top_20_percent_importance, ]
importance_df_filtered <- importance_df_filtered[importance_df_filtered$p_value < 0.05, ]

# Sort by Gini in descending order
gini_df_sorted <- gini_df[order(gini_df$Gini, decreasing = TRUE), ]

# Compute number of genes in the top 20%
top_20_percent_gini <- round(0.2 * nrow(gini_df_sorted))

# Keep top 20% and filter p < 0.05
gini_df_filtered <- gini_df_sorted[1:top_20_percent_gini, ]
gini_df_filtered <- gini_df_filtered[gini_df_filtered$p_value < 0.05, ]

# Extract candidate gene names
importance_genes <- importance_df_filtered$Feature
gini_genes <- gini_df_filtered$Feature

# Intersection of candidates from both metrics
common_genes <- intersect(importance_genes, gini_genes)

# Optionally, save as a data.frame
common_genes_df <- data.frame(Feature = common_genes)

# Preview final feature list
print(common_genes_df)
saveRDS(common_genes_df, file = 'features_RFmodel.rds')

library(VennDiagram)
library(grid)

# Unique elements for the Venn diagram
importance_genes <- unique(importance_df_filtered$Feature)
gini_genes <- unique(gini_df_filtered$Feature)

# Save Venn diagram as a 5x5 inch PDF
pdf("gene_intersection_venn.pdf", width = 5, height = 5)

# Venn diagram of overlap between importance_df and gini_df features
venn.plot <- venn.diagram(
  x = list(MDA = importance_genes, Gini = gini_genes),   # Category labels: MDA and Gini
  category.names = c("MDA", "Gini"),
  filename = NULL,       # Do not write an image file here (we draw to the PDF device)
  output = TRUE,
  imagetype = "png",     # Image type (kept as-is)
  col = "black",         # Circle color
  fill = c("#FF6347", "#4682B4"),  # Fill colors
  alpha = 0.5,           # Transparency
  label.col = "white",   # Label color
  cex = 1.5,             # Label font size
  cat.cex = 1.5,         # Category label size
  cat.pos = c(-115, 115),# Category label positions
  cat.dist = c(0.1, 0.1),# Distance between category labels and circles
  margin = 0.1
)

# Draw the Venn diagram to the PDF device
grid.draw(venn.plot)

# Close the PDF device to save the file
dev.off()

# Optionally draw again to the current graphics device
grid.draw(venn.plot)
