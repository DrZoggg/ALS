setwd('ML_IntersectGene_Eval')
options(stringsAsFactors = F)
library(data.table)
library(tibble)
library(dplyr)
library(GEOquery)
library(clusterProfiler)
library(org.Hs.eg.db)
load("../Data_preprocessing/valid_cohort.rda")
ls()
table(colnames(valid_cohort_expr) == valid_cohort_clin$Sample)
valid_cohort_expr <- valid_cohort_expr[,valid_cohort_clin$diagnosis %in% c('ALS','CON')]
valid_cohort_clin <- valid_cohort_clin[valid_cohort_clin$diagnosis %in% c('ALS','CON'),]
my_mRNA <- valid_cohort_expr
print(paste0('# Dimensions of the entire DISC cohort:',dim(my_mRNA)))
# Split into male and female cohorts for stratified analysis
valid_cohort_clin_male <- valid_cohort_clin[valid_cohort_clin$Sex == 'm',]
valid_cohort_expr_male <- valid_cohort_expr[,valid_cohort_clin$Sex == 'm',]

valid_cohort_clin_female <- valid_cohort_clin[valid_cohort_clin$Sex == 'v',]
valid_cohort_expr_female <- valid_cohort_expr[,valid_cohort_clin$Sex == 'v',]


###########################ABCA1
ABCA1_exp <- valid_cohort_expr_male['ABCA1',] 
colnames(ABCA1_exp) <- valid_cohort_clin_male$diagnosis
ABCA1_exp <- t(ABCA1_exp)
str(ABCA1_exp)
# Convert matrix to data frame
ABCA1_exp_df <- data.frame(Diagnosis = rownames(ABCA1_exp), ABCA1_Expression = as.numeric(ABCA1_exp[,1]))
# Ensure Diagnosis is a factor variable
ABCA1_exp_df$Diagnosis <- as.factor(ABCA1_exp_df$Diagnosis)
box_data <- subset(ABCA1_exp_df, Diagnosis %in% c("ALS", "CON"))
library(ggplot2)
library(ggpubr)

group_colors <- c("ALS" = "#9AC5CD", "CON" = "#FFDAB9")

p <- ggplot(box_data, aes(x = Diagnosis, y = ABCA1_Expression, fill = Diagnosis)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.1, alpha = 0.8, color = "black", size = 0.5) +
  labs(title = "ABCA1 Expression across Different Groups", 
       subtitle = "Male participants", 
       x = "Diagnosis", 
       y = "ABCA1 Expression") +
  scale_fill_manual(values = group_colors) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        legend.position = "right") +
  stat_compare_means(method = "wilcox.test",
                     comparisons = list(c("ALS", "CON")),
                     label = "p.signif") # ⭐ Automatically add significance stars

print(p)
ggsave(p,file='ABCA1Boxplot_male.pdf',height=5,width=5)
 
library(pROC)    
library(ggplot2)   
library(grid)      

 # ==== Data preparation ====
 
roc_df <- droplevels(subset(box_data, Diagnosis %in% c("ALS","CON")))
roc_df$Diagnosis <- factor(roc_df$Diagnosis, levels = c("CON","ALS"))  

 # Compute ROC curve object
roc_abca1 <- roc(response = roc_df$Diagnosis,
                 predictor = roc_df$ABCA1_Expression,   
                 levels = c("CON","ALS"),
                 direction = "auto", quiet = TRUE)

# Extract AUC value
auc_val <- as.numeric(auc(roc_abca1))

 # Find optimal threshold (Youden index) and extract common performance metrics
best <- coords(roc_abca1, x = "best", best.method = "youden",
               ret = c("threshold","sensitivity","specificity","accuracy"),
               transpose = FALSE)
sen <- as.numeric(best["sensitivity"])   
spe <- as.numeric(best["specificity"])   
acc <- as.numeric(best["accuracy"])     

# Extract points along the ROC curve (Specificity vs Sensitivity)
cd  <- coords(roc_abca1, x = "all", ret = c("specificity","sensitivity"), transpose = FALSE)
df_roc <- data.frame(Specificity = cd[, "specificity"], Sensitivity = cd[, "sensitivity"])

 # ==== Customize info box position (after x-axis reversal, bottom-right ≈ (0.05, 0.08)) ====
label_x <- 0.05  
label_y <- 0.08   

# ==== Plotting ====
options(repr.plot.width=6, repr.plot.height=6)  # Set plot window size

p <- ggplot(df_roc, aes(x = Specificity, y = Sensitivity)) +
  geom_line(linewidth = 1.4, color = "#2E86DE") +  # Plot ROC curve (blue)
  
# Add 45° gray dashed line (random classification reference)
  annotate("segment", x = 1, y = 0, xend = 0, yend = 1,
           linetype = "dashed", linewidth = 0.9, color = "grey50") +
  
# Add red dashed lines for Sensitivity and Specificity at optimal threshold
  geom_hline(yintercept = sen, linetype = "dashed", linewidth = 0.9, color = "red") +
  geom_vline(xintercept = spe, linetype = "dashed", linewidth = 0.9, color = "red") +
  
# Add info box (AUC, Accuracy, Sensitivity, Specificity) at bottom right of plot
  annotate("label", x = label_x, y = label_y, hjust = 1, vjust = 0,
           label = sprintf("AUC = %.3f\nAccuracy = %.3f\nSensitivity = %.3f\nSpecificity = %.3f",
                           auc_val, acc, sen, spe),
           size = 4.1, fill = "#F7F9FB", color = "black", label.size = 0.3) +
  
  # Add title, subtitle, and axis labels
  labs(
    title = "ABCA1 (Single-Gene ROC)",
    subtitle = "Male (ALS vs CON)",
    x = "Specificity",
    y = "Sensitivity"
  ) +
  
  # Reverse x-axis: Specificity 1 → 0, fixed ticks with padding at both ends
  scale_x_reverse(breaks = seq(0, 1, 0.25),
                  limits = c(1, 0),
                  expand = expansion(mult = c(0.02, 0.02))) +
  # Set y-axis range to 0–1 with 0.25 intervals
  scale_y_continuous(breaks = seq(0, 1, 0.25),
                     limits = c(0, 1),
                     expand = expansion(mult = c(0.02, 0.02))) +
  
# Set theme style
  theme_classic(base_size = 12) +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold", size = 15), # Center and bold main title
    plot.subtitle = element_text(hjust = 0.5, face = "bold", size = 11),  # Center and bold subtitle
    axis.title    = element_text(face = "bold", size = 13),               # Bold axis titles
    axis.text     = element_text(color = "black", size = 11),             # Customize axis tick labels
    axis.ticks    = element_line(linewidth = 0.9, color = "black"),       # Customize axis ticks
    axis.ticks.length = unit(0.22, "cm"),                                 # Adjust tick length
    panel.border  = element_rect(color = "black", fill = NA, linewidth = 1.2),  # Remove black border around the plot
    axis.line     = element_blank(),                                      # Remove default x/y axes to avoid double lines
    plot.margin   = margin(10, 12, 10, 12)                                # Margin
  ) +
  coord_fixed()  # Keep x:y ratio = 1:1

print(p)

# Save the figure
ggsave("ABCA1_ROC_male.pdf", p, width = 6, height = 6)

ABCA1_exp <- valid_cohort_expr_female['ABCA1',] 
colnames(ABCA1_exp) <- valid_cohort_clin_female$diagnosis
ABCA1_exp <- t(ABCA1_exp)
str(ABCA1_exp)
# Convert matrix to data frame
ABCA1_exp_df <- data.frame(Diagnosis = rownames(ABCA1_exp), ABCA1_Expression = as.numeric(ABCA1_exp[,1]))
# Ensure Diagnosis is a factor variable
ABCA1_exp_df$Diagnosis <- as.factor(ABCA1_exp_df$Diagnosis)
box_data <- subset(ABCA1_exp_df, Diagnosis %in% c("ALS", "CON"))

library(ggplot2)
library(ggpubr)

group_colors <- c("ALS" = "#9AC5CD", "CON" = "#FFDAB9")

p <- ggplot(box_data, aes(x = Diagnosis, y = ABCA1_Expression, fill = Diagnosis)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.1, alpha = 0.8, color = "black", size = 0.5) +
  labs(title = "ABCA1 Expression across Different Groups", 
       subtitle = "Female participants", 
       x = "Diagnosis", 
       y = "ABCA1 Expression") +
  scale_fill_manual(values = group_colors) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        legend.position = "right") +
  stat_compare_means(method = "wilcox.test",
                     comparisons = list(c("ALS", "CON")),
                     label = "p.signif") # ⭐ Automatically add significance stars

print(p)
ggsave(p,file='ABCA1Boxplot_female.pdf',height=5,width=5)
 
library(pROC)    
library(ggplot2)   
library(grid)      

 # ==== Data preparation ====
 
roc_df <- droplevels(subset(box_data, Diagnosis %in% c("ALS","CON")))
roc_df$Diagnosis <- factor(roc_df$Diagnosis, levels = c("CON","ALS"))  

 # Compute ROC curve object
roc_abca1 <- roc(response = roc_df$Diagnosis,
                 predictor = roc_df$ABCA1_Expression,   
                 levels = c("CON","ALS"),
                 direction = "auto", quiet = TRUE)

# Extract AUC value
auc_val <- as.numeric(auc(roc_abca1))

 # Find optimal threshold (Youden index) and extract common performance metrics
best <- coords(roc_abca1, x = "best", best.method = "youden",
               ret = c("threshold","sensitivity","specificity","accuracy"),
               transpose = FALSE)
sen <- as.numeric(best["sensitivity"])   
spe <- as.numeric(best["specificity"])   
acc <- as.numeric(best["accuracy"])     

# Extract points along the ROC curve (Specificity vs Sensitivity)
cd  <- coords(roc_abca1, x = "all", ret = c("specificity","sensitivity"), transpose = FALSE)
df_roc <- data.frame(Specificity = cd[, "specificity"], Sensitivity = cd[, "sensitivity"])

 # ==== Customize info box position (after x-axis reversal, bottom-right ≈ (0.05, 0.08)) ====
label_x <- 0.05   
label_y <- 0.08   

# ==== Plotting ====
options(repr.plot.width=6, repr.plot.height=6)  # Set plot window size

p <- ggplot(df_roc, aes(x = Specificity, y = Sensitivity)) +
  geom_line(linewidth = 1.4, color = "#2E86DE") +  # Plot ROC curve (blue)
  
# Add 45° gray dashed line (random classification reference)
  annotate("segment", x = 1, y = 0, xend = 0, yend = 1,
           linetype = "dashed", linewidth = 0.9, color = "grey50") +
  
# Add red dashed lines for Sensitivity and Specificity at optimal threshold
  geom_hline(yintercept = sen, linetype = "dashed", linewidth = 0.9, color = "red") +
  geom_vline(xintercept = spe, linetype = "dashed", linewidth = 0.9, color = "red") +
  
# Add info box (AUC, Accuracy, Sensitivity, Specificity) at bottom right of plot
  annotate("label", x = label_x, y = label_y, hjust = 1, vjust = 0,
           label = sprintf("AUC = %.3f\nAccuracy = %.3f\nSensitivity = %.3f\nSpecificity = %.3f",
                           auc_val, acc, sen, spe),
           size = 4.1, fill = "#F7F9FB", color = "black", label.size = 0.3) +
  
  # Add title, subtitle, and axis labels
  labs(
    title = "ABCA1 (Single-Gene ROC)",
    subtitle = "Female (ALS vs CON)",
    x = "Specificity",
    y = "Sensitivity"
  ) +
  
  # Reverse x-axis: Specificity 1 → 0, fixed ticks with padding at both ends
  scale_x_reverse(breaks = seq(0, 1, 0.25),
                  limits = c(1, 0),
                  expand = expansion(mult = c(0.02, 0.02))) +
  # Set y-axis range to 0–1 with 0.25 intervals
  scale_y_continuous(breaks = seq(0, 1, 0.25),
                     limits = c(0, 1),
                     expand = expansion(mult = c(0.02, 0.02))) +
  
# Set theme style
  theme_classic(base_size = 12) +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold", size = 15), # Center and bold main title
    plot.subtitle = element_text(hjust = 0.5, face = "bold", size = 11),  # Center and bold subtitle
    axis.title    = element_text(face = "bold", size = 13),               # Bold axis titles
    axis.text     = element_text(color = "black", size = 11),             # Customize axis tick labels
    axis.ticks    = element_line(linewidth = 0.9, color = "black"),       # Customize axis ticks
    axis.ticks.length = unit(0.22, "cm"),                                 # Adjust tick length
    panel.border  = element_rect(color = "black", fill = NA, linewidth = 1.2),  # Remove black border around the plot
    axis.line     = element_blank(),                                      # Remove default x/y axes to avoid double lines
    plot.margin   = margin(10, 12, 10, 12)                                # Margin
  ) +
  coord_fixed()  # Keep x:y ratio = 1:1

print(p)

# Save the figure
ggsave("ABCA1_ROC_female.pdf", p, width = 6, height = 6)

#############################DDX51
ABCA1_exp <- valid_cohort_expr_male['DDX51',] 
colnames(ABCA1_exp) <- valid_cohort_clin_male$diagnosis
ABCA1_exp <- t(ABCA1_exp)
str(ABCA1_exp)
# Convert matrix to data frame
ABCA1_exp_df <- data.frame(Diagnosis = rownames(ABCA1_exp), ABCA1_Expression = as.numeric(ABCA1_exp[,1]))
# Ensure Diagnosis is a factor variable
ABCA1_exp_df$Diagnosis <- as.factor(ABCA1_exp_df$Diagnosis)
box_data <- subset(ABCA1_exp_df, Diagnosis %in% c("ALS", "CON"))

library(ggplot2)
library(ggpubr)

group_colors <- c("ALS" = "#9AC5CD", "CON" = "#FFDAB9")

p <- ggplot(box_data, aes(x = Diagnosis, y = ABCA1_Expression, fill = Diagnosis)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.1, alpha = 0.8, color = "black", size = 0.5) +
  labs(title = "DDX51 Expression across Different Groups", 
       subtitle = "Male participants",
       x = "Diagnosis", 
       y = "DDX51 Expression") +
  scale_fill_manual(values = group_colors) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        legend.position = "right") +
  stat_compare_means(method = "wilcox.test",
                     comparisons = list(c("ALS", "CON")),
                     label = "p.signif") # ⭐ Automatically add significance stars

print(p)
ggsave(p,file='DDX51Boxplot_male.pdf',height=5,width=5)
 
library(pROC)    
library(ggplot2)   
library(grid)      

 # ==== Data preparation ====
 
roc_df <- droplevels(subset(box_data, Diagnosis %in% c("ALS","CON")))
roc_df$Diagnosis <- factor(roc_df$Diagnosis, levels = c("CON","ALS"))  

 # Compute ROC curve object
roc_abca1 <- roc(response = roc_df$Diagnosis,
                 predictor = roc_df$ABCA1_Expression,   
                 levels = c("CON","ALS"),
                 direction = "auto", quiet = TRUE)

# Extract AUC value
auc_val <- as.numeric(auc(roc_abca1))

 # Find optimal threshold (Youden index) and extract common performance metrics
best <- coords(roc_abca1, x = "best", best.method = "youden",
               ret = c("threshold","sensitivity","specificity","accuracy"),
               transpose = FALSE)
sen <- as.numeric(best["sensitivity"])   
spe <- as.numeric(best["specificity"])   
acc <- as.numeric(best["accuracy"])     

# Extract points along the ROC curve (Specificity vs Sensitivity)
cd  <- coords(roc_abca1, x = "all", ret = c("specificity","sensitivity"), transpose = FALSE)
df_roc <- data.frame(Specificity = cd[, "specificity"], Sensitivity = cd[, "sensitivity"])

 # ==== Customize info box position (after x-axis reversal, bottom-right ≈ (0.05, 0.08)) ====
label_x <- 0.05   
label_y <- 0.08   

# ==== Plotting ====
options(repr.plot.width=6, repr.plot.height=6)  # Set plot window size

p <- ggplot(df_roc, aes(x = Specificity, y = Sensitivity)) +
  geom_line(linewidth = 1.4, color = "#2E86DE") +  # Plot ROC curve (blue)
  
# Add 45° gray dashed line (random classification reference)
  annotate("segment", x = 1, y = 0, xend = 0, yend = 1,
           linetype = "dashed", linewidth = 0.9, color = "grey50") +
  
# Add red dashed lines for Sensitivity and Specificity at optimal threshold
  geom_hline(yintercept = sen, linetype = "dashed", linewidth = 0.9, color = "red") +
  geom_vline(xintercept = spe, linetype = "dashed", linewidth = 0.9, color = "red") +
  
# Add info box (AUC, Accuracy, Sensitivity, Specificity) at bottom right of plot
  annotate("label", x = label_x, y = label_y, hjust = 1, vjust = 0,
           label = sprintf("AUC = %.3f\nAccuracy = %.3f\nSensitivity = %.3f\nSpecificity = %.3f",
                           auc_val, acc, sen, spe),
           size = 4.1, fill = "#F7F9FB", color = "black", label.size = 0.3) +
  
  # Add title, subtitle, and axis labels
  labs(
    title = "DDX51 (Single-Gene ROC)",
    subtitle = "Male (ALS vs CON)",
    x = "Specificity",
    y = "Sensitivity"
  ) +
  
  # Reverse x-axis: Specificity 1 → 0, fixed ticks with padding at both ends
  scale_x_reverse(breaks = seq(0, 1, 0.25),
                  limits = c(1, 0),
                  expand = expansion(mult = c(0.02, 0.02))) +
  # Set y-axis range to 0–1 with 0.25 intervals
  scale_y_continuous(breaks = seq(0, 1, 0.25),
                     limits = c(0, 1),
                     expand = expansion(mult = c(0.02, 0.02))) +
  
# Set theme style
  theme_classic(base_size = 12) +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold", size = 15), # Center and bold main title
    plot.subtitle = element_text(hjust = 0.5, face = "bold", size = 11),  # Center and bold subtitle
    axis.title    = element_text(face = "bold", size = 13),               # Bold axis titles
    axis.text     = element_text(color = "black", size = 11),             # Customize axis tick labels
    axis.ticks    = element_line(linewidth = 0.9, color = "black"),       # Customize axis ticks
    axis.ticks.length = unit(0.22, "cm"),                                 # Adjust tick length
    panel.border  = element_rect(color = "black", fill = NA, linewidth = 1.2),  # Remove black border around the plot
    axis.line     = element_blank(),                                      # Remove default x/y axes to avoid double lines
    plot.margin   = margin(10, 12, 10, 12)                                # Margin
  ) +
  coord_fixed()  # Keep x:y ratio = 1:1

print(p)

# Save the figure
ggsave("DDX51_ROC_male.pdf", p, width = 6, height = 6)
ABCA1_exp <- valid_cohort_expr_female['DDX51',] 
colnames(ABCA1_exp) <- valid_cohort_clin_female$diagnosis
ABCA1_exp <- t(ABCA1_exp)
str(ABCA1_exp)
# Convert matrix to data frame
ABCA1_exp_df <- data.frame(Diagnosis = rownames(ABCA1_exp), ABCA1_Expression = as.numeric(ABCA1_exp[,1]))
# Ensure Diagnosis is a factor variable
ABCA1_exp_df$Diagnosis <- as.factor(ABCA1_exp_df$Diagnosis)
box_data <- subset(ABCA1_exp_df, Diagnosis %in% c("ALS", "CON"))
library(ggplot2)
library(ggpubr)

group_colors <- c("ALS" = "#9AC5CD", "CON" = "#FFDAB9")

p <- ggplot(box_data, aes(x = Diagnosis, y = ABCA1_Expression, fill = Diagnosis)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.1, alpha = 0.8, color = "black", size = 0.5) +
  labs(title = "DDX51 Expression across Different Groups", 
       subtitle = "Female participants",
       x = "Diagnosis", 
       y = "DDX51 Expression") +
  scale_fill_manual(values = group_colors) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        legend.position = "right") +
  stat_compare_means(method = "wilcox.test",
                     comparisons = list(c("ALS", "CON")),
                     label = "p.signif") # ⭐ Automatically add significance stars

print(p)
ggsave(p,file='DDX51Boxplot_female.pdf',height=5,width=5)

 
library(pROC)    
library(ggplot2)   
library(grid)      

 # ==== Data preparation ====
 
roc_df <- droplevels(subset(box_data, Diagnosis %in% c("ALS","CON")))
roc_df$Diagnosis <- factor(roc_df$Diagnosis, levels = c("CON","ALS"))  

 # Compute ROC curve object
roc_abca1 <- roc(response = roc_df$Diagnosis,
                 predictor = roc_df$ABCA1_Expression,   
                 levels = c("CON","ALS"),
                 direction = "auto", quiet = TRUE)

# Extract AUC value
auc_val <- as.numeric(auc(roc_abca1))

 # Find optimal threshold (Youden index) and extract common performance metrics
best <- coords(roc_abca1, x = "best", best.method = "youden",
               ret = c("threshold","sensitivity","specificity","accuracy"),
               transpose = FALSE)
sen <- as.numeric(best["sensitivity"])   
spe <- as.numeric(best["specificity"])   
acc <- as.numeric(best["accuracy"])     

# Extract points along the ROC curve (Specificity vs Sensitivity)
cd  <- coords(roc_abca1, x = "all", ret = c("specificity","sensitivity"), transpose = FALSE)
df_roc <- data.frame(Specificity = cd[, "specificity"], Sensitivity = cd[, "sensitivity"])

 # ==== Customize info box position (after x-axis reversal, bottom-right ≈ (0.05, 0.08)) ====
label_x <- 0.05   
label_y <- 0.08   

# ==== Plotting ====
options(repr.plot.width=6, repr.plot.height=6)  # Set plot window size

p <- ggplot(df_roc, aes(x = Specificity, y = Sensitivity)) +
  geom_line(linewidth = 1.4, color = "#2E86DE") +  # Plot ROC curve (blue)
  
# Add 45° gray dashed line (random classification reference)
  annotate("segment", x = 1, y = 0, xend = 0, yend = 1,
           linetype = "dashed", linewidth = 0.9, color = "grey50") +
  
# Add red dashed lines for Sensitivity and Specificity at optimal threshold
  geom_hline(yintercept = sen, linetype = "dashed", linewidth = 0.9, color = "red") +
  geom_vline(xintercept = spe, linetype = "dashed", linewidth = 0.9, color = "red") +
  
# Add info box (AUC, Accuracy, Sensitivity, Specificity) at bottom right of plot
  annotate("label", x = label_x, y = label_y, hjust = 1, vjust = 0,
           label = sprintf("AUC = %.3f\nAccuracy = %.3f\nSensitivity = %.3f\nSpecificity = %.3f",
                           auc_val, acc, sen, spe),
           size = 4.1, fill = "#F7F9FB", color = "black", label.size = 0.3) +
  
  # Add title, subtitle, and axis labels
  labs(
    title = "DDX51 (Single-Gene ROC)",
    subtitle = "Female (ALS vs CON)",
    x = "Specificity",
    y = "Sensitivity"
  ) +
  
  # Reverse x-axis: Specificity 1 → 0, fixed ticks with padding at both ends
  scale_x_reverse(breaks = seq(0, 1, 0.25),
                  limits = c(1, 0),
                  expand = expansion(mult = c(0.02, 0.02))) +
  # Set y-axis range to 0–1 with 0.25 intervals
  scale_y_continuous(breaks = seq(0, 1, 0.25),
                     limits = c(0, 1),
                     expand = expansion(mult = c(0.02, 0.02))) +
  
# Set theme style
  theme_classic(base_size = 12) +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold", size = 15), # Center and bold main title
    plot.subtitle = element_text(hjust = 0.5, face = "bold", size = 11),  # Center and bold subtitle
    axis.title    = element_text(face = "bold", size = 13),               # Bold axis titles
    axis.text     = element_text(color = "black", size = 11),             # Customize axis tick labels
    axis.ticks    = element_line(linewidth = 0.9, color = "black"),       # Customize axis ticks
    axis.ticks.length = unit(0.22, "cm"),                                 # Adjust tick length
    panel.border  = element_rect(color = "black", fill = NA, linewidth = 1.2),  # Remove black border around the plot
    axis.line     = element_blank(),                                      # Remove default x/y axes to avoid double lines
    plot.margin   = margin(10, 12, 10, 12)                                # Margin
  ) +
  coord_fixed()  # Keep x:y ratio = 1:1

print(p)

# Save the figure
ggsave("DDX51_ROC_female.pdf", p, width = 6, height = 6)

########################TMEM71
ABCA1_exp <- valid_cohort_expr_male['TMEM71',] 
colnames(ABCA1_exp) <- valid_cohort_clin_male$diagnosis
ABCA1_exp <- t(ABCA1_exp)
str(ABCA1_exp)
# Convert matrix to data frame
ABCA1_exp_df <- data.frame(Diagnosis = rownames(ABCA1_exp), ABCA1_Expression = as.numeric(ABCA1_exp[,1]))
# Ensure Diagnosis is a factor variable
ABCA1_exp_df$Diagnosis <- as.factor(ABCA1_exp_df$Diagnosis)
box_data <- subset(ABCA1_exp_df, Diagnosis %in% c("ALS", "CON"))
library(ggplot2)
library(ggpubr)

group_colors <- c("ALS" = "#9AC5CD", "CON" = "#FFDAB9")

p <- ggplot(box_data, aes(x = Diagnosis, y = ABCA1_Expression, fill = Diagnosis)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.1, alpha = 0.8, color = "black", size = 0.5) +
  labs(title = "TMEM71 Expression across Different Groups", 
       subtitle = "Male participants",
       x = "Diagnosis", 
       y = "TMEM71 Expression") +
  scale_fill_manual(values = group_colors) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        legend.position = "right") +
  stat_compare_means(method = "wilcox.test",
                     comparisons = list(c("ALS", "CON")),
                     label = "p.signif") # ⭐ Automatically add significance stars

print(p)
ggsave(p,file='TMEM71Boxplot_male.pdf',height=5,width=5)
 
library(pROC)    
library(ggplot2)   
library(grid)      

 # ==== Data preparation ====
 
roc_df <- droplevels(subset(box_data, Diagnosis %in% c("ALS","CON")))
roc_df$Diagnosis <- factor(roc_df$Diagnosis, levels = c("CON","ALS"))  

 # Compute ROC curve object
roc_abca1 <- roc(response = roc_df$Diagnosis,
                 predictor = roc_df$ABCA1_Expression,   
                 levels = c("CON","ALS"),
                 direction = "auto", quiet = TRUE)

# Extract AUC value
auc_val <- as.numeric(auc(roc_abca1))

 # Find optimal threshold (Youden index) and extract common performance metrics
best <- coords(roc_abca1, x = "best", best.method = "youden",
               ret = c("threshold","sensitivity","specificity","accuracy"),
               transpose = FALSE)
sen <- as.numeric(best["sensitivity"])   
spe <- as.numeric(best["specificity"])   
acc <- as.numeric(best["accuracy"])     

# Extract points along the ROC curve (Specificity vs Sensitivity)
cd  <- coords(roc_abca1, x = "all", ret = c("specificity","sensitivity"), transpose = FALSE)
df_roc <- data.frame(Specificity = cd[, "specificity"], Sensitivity = cd[, "sensitivity"])

 # ==== Customize info box position (after x-axis reversal, bottom-right ≈ (0.05, 0.08)) ====
label_x <- 0.05   
label_y <- 0.08   

# ==== Plotting ====
options(repr.plot.width=6, repr.plot.height=6)  # Set plot window size

p <- ggplot(df_roc, aes(x = Specificity, y = Sensitivity)) +
  geom_line(linewidth = 1.4, color = "#2E86DE") +  # Plot ROC curve (blue)
  
# Add 45° gray dashed line (random classification reference)
  annotate("segment", x = 1, y = 0, xend = 0, yend = 1,
           linetype = "dashed", linewidth = 0.9, color = "grey50") +
  
# Add red dashed lines for Sensitivity and Specificity at optimal threshold
  geom_hline(yintercept = sen, linetype = "dashed", linewidth = 0.9, color = "red") +
  geom_vline(xintercept = spe, linetype = "dashed", linewidth = 0.9, color = "red") +
  
# Add info box (AUC, Accuracy, Sensitivity, Specificity) at bottom right of plot
  annotate("label", x = label_x, y = label_y, hjust = 1, vjust = 0,
           label = sprintf("AUC = %.3f\nAccuracy = %.3f\nSensitivity = %.3f\nSpecificity = %.3f",
                           auc_val, acc, sen, spe),
           size = 4.1, fill = "#F7F9FB", color = "black", label.size = 0.3) +
  
  # Add title, subtitle, and axis labels
  labs(
    title = "TMEM71 (Single-Gene ROC)",
    subtitle = "Male (ALS vs CON)",
    x = "Specificity",
    y = "Sensitivity"
  ) +
  
  # Reverse x-axis: Specificity 1 → 0, fixed ticks with padding at both ends
  scale_x_reverse(breaks = seq(0, 1, 0.25),
                  limits = c(1, 0),
                  expand = expansion(mult = c(0.02, 0.02))) +
  # Set y-axis range to 0–1 with 0.25 intervals
  scale_y_continuous(breaks = seq(0, 1, 0.25),
                     limits = c(0, 1),
                     expand = expansion(mult = c(0.02, 0.02))) +
  
# Set theme style
  theme_classic(base_size = 12) +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold", size = 15), # Center and bold main title
    plot.subtitle = element_text(hjust = 0.5, face = "bold", size = 11),  # Center and bold subtitle
    axis.title    = element_text(face = "bold", size = 13),               # Bold axis titles
    axis.text     = element_text(color = "black", size = 11),             # Customize axis tick labels
    axis.ticks    = element_line(linewidth = 0.9, color = "black"),       # Customize axis ticks
    axis.ticks.length = unit(0.22, "cm"),                                 # Adjust tick length
    panel.border  = element_rect(color = "black", fill = NA, linewidth = 1.2),  # Remove black border around the plot
    axis.line     = element_blank(),                                      # Remove default x/y axes to avoid double lines
    plot.margin   = margin(10, 12, 10, 12)                                # Margin
  ) +
  coord_fixed()  # Keep x:y ratio = 1:1

print(p)

# Save the figure
ggsave("TMEM71_ROC_male.pdf", p, width = 6, height = 6)
ABCA1_exp <- valid_cohort_expr_female['TMEM71',] 
colnames(ABCA1_exp) <- valid_cohort_clin_female$diagnosis
ABCA1_exp <- t(ABCA1_exp)
str(ABCA1_exp)
# Convert matrix to data frame
ABCA1_exp_df <- data.frame(Diagnosis = rownames(ABCA1_exp), ABCA1_Expression = as.numeric(ABCA1_exp[,1]))
# Ensure Diagnosis is a factor variable
ABCA1_exp_df$Diagnosis <- as.factor(ABCA1_exp_df$Diagnosis)
box_data <- subset(ABCA1_exp_df, Diagnosis %in% c("ALS", "CON"))
library(ggplot2)
library(ggpubr)

group_colors <- c("ALS" = "#9AC5CD", "CON" = "#FFDAB9")

p <- ggplot(box_data, aes(x = Diagnosis, y = ABCA1_Expression, fill = Diagnosis)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.1, alpha = 0.8, color = "black", size = 0.5) +
  labs(title = "TMEM71 Expression across Different Groups", 
       subtitle = "Female participants",
       x = "Diagnosis", 
       y = "TMEM71 Expression") +
  scale_fill_manual(values = group_colors) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        legend.position = "right") +
  stat_compare_means(method = "wilcox.test",
                     comparisons = list(c("ALS", "CON")),
                     label = "p.signif") # ⭐ Automatically add significance stars

print(p)
ggsave(p,file='TMEM71Boxplot_female.pdf',height=5,width=5)

 
library(pROC)    
library(ggplot2)   
library(grid)      

 # ==== Data preparation ====
 
roc_df <- droplevels(subset(box_data, Diagnosis %in% c("ALS","CON")))
roc_df$Diagnosis <- factor(roc_df$Diagnosis, levels = c("CON","ALS"))  

 # Compute ROC curve object
roc_abca1 <- roc(response = roc_df$Diagnosis,
                 predictor = roc_df$ABCA1_Expression,   
                 levels = c("CON","ALS"),
                 direction = "auto", quiet = TRUE)

# Extract AUC value
auc_val <- as.numeric(auc(roc_abca1))

 # Find optimal threshold (Youden index) and extract common performance metrics
best <- coords(roc_abca1, x = "best", best.method = "youden",
               ret = c("threshold","sensitivity","specificity","accuracy"),
               transpose = FALSE)
sen <- as.numeric(best["sensitivity"])   
spe <- as.numeric(best["specificity"])   
acc <- as.numeric(best["accuracy"])     

# Extract points along the ROC curve (Specificity vs Sensitivity)
cd  <- coords(roc_abca1, x = "all", ret = c("specificity","sensitivity"), transpose = FALSE)
df_roc <- data.frame(Specificity = cd[, "specificity"], Sensitivity = cd[, "sensitivity"])

 # ==== Customize info box position (after x-axis reversal, bottom-right ≈ (0.05, 0.08)) ====
label_x <- 0.05   
label_y <- 0.08   

# ==== Plotting ====
options(repr.plot.width=6, repr.plot.height=6)  # Set plot window size

p <- ggplot(df_roc, aes(x = Specificity, y = Sensitivity)) +
  geom_line(linewidth = 1.4, color = "#2E86DE") +  # Plot ROC curve (blue)
  
# Add 45° gray dashed line (random classification reference)
  annotate("segment", x = 1, y = 0, xend = 0, yend = 1,
           linetype = "dashed", linewidth = 0.9, color = "grey50") +
  
# Add red dashed lines for Sensitivity and Specificity at optimal threshold
  geom_hline(yintercept = sen, linetype = "dashed", linewidth = 0.9, color = "red") +
  geom_vline(xintercept = spe, linetype = "dashed", linewidth = 0.9, color = "red") +
  
# Add info box (AUC, Accuracy, Sensitivity, Specificity) at bottom right of plot
  annotate("label", x = label_x, y = label_y, hjust = 1, vjust = 0,
           label = sprintf("AUC = %.3f\nAccuracy = %.3f\nSensitivity = %.3f\nSpecificity = %.3f",
                           auc_val, acc, sen, spe),
           size = 4.1, fill = "#F7F9FB", color = "black", label.size = 0.3) +
  
  # Add title, subtitle, and axis labels
  labs(
    title = "TMEM71 (Single-Gene ROC)",
    subtitle = "Female (ALS vs CON)",
    x = "Specificity",
    y = "Sensitivity"
  ) +
  
  # Reverse x-axis: Specificity 1 → 0, fixed ticks with padding at both ends
  scale_x_reverse(breaks = seq(0, 1, 0.25),
                  limits = c(1, 0),
                  expand = expansion(mult = c(0.02, 0.02))) +
  # Set y-axis range to 0–1 with 0.25 intervals
  scale_y_continuous(breaks = seq(0, 1, 0.25),
                     limits = c(0, 1),
                     expand = expansion(mult = c(0.02, 0.02))) +
  
# Set theme style
  theme_classic(base_size = 12) +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold", size = 15), # Center and bold main title
    plot.subtitle = element_text(hjust = 0.5, face = "bold", size = 11),  # Center and bold subtitle
    axis.title    = element_text(face = "bold", size = 13),               # Bold axis titles
    axis.text     = element_text(color = "black", size = 11),             # Customize axis tick labels
    axis.ticks    = element_line(linewidth = 0.9, color = "black"),       # Customize axis ticks
    axis.ticks.length = unit(0.22, "cm"),                                 # Adjust tick length
    panel.border  = element_rect(color = "black", fill = NA, linewidth = 1.2),  # Remove black border around the plot
    axis.line     = element_blank(),                                      # Remove default x/y axes to avoid double lines
    plot.margin   = margin(10, 12, 10, 12)                                # Margin
  ) +
  coord_fixed()  # Keep x:y ratio = 1:1

print(p)

# Save the figure
ggsave("TMEM71_ROC_female.pdf", p, width = 6, height = 6)

#############QPCT
ABCA1_exp <- valid_cohort_expr_male['QPCT',] 
colnames(ABCA1_exp) <- valid_cohort_clin_male$diagnosis
ABCA1_exp <- t(ABCA1_exp)
str(ABCA1_exp)
# Convert matrix to data frame
ABCA1_exp_df <- data.frame(Diagnosis = rownames(ABCA1_exp), ABCA1_Expression = as.numeric(ABCA1_exp[,1]))
# Ensure Diagnosis is a factor variable
ABCA1_exp_df$Diagnosis <- as.factor(ABCA1_exp_df$Diagnosis)
box_data <- subset(ABCA1_exp_df, Diagnosis %in% c("ALS", "CON"))

library(ggplot2)
library(ggpubr)

group_colors <- c("ALS" = "#9AC5CD", "CON" = "#FFDAB9")

p <- ggplot(box_data, aes(x = Diagnosis, y = ABCA1_Expression, fill = Diagnosis)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.1, alpha = 0.8, color = "black", size = 0.5) +
  labs(title = "QPCT Expression across Different Groups",
       subtitle = "Male participants",
       x = "Diagnosis", 
       y = "QPCT Expression") +
  scale_fill_manual(values = group_colors) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        legend.position = "right") +
  stat_compare_means(method = "wilcox.test",
                     comparisons = list(c("ALS", "CON")),
                     label = "p.signif") # ⭐ Automatically add significance stars

print(p)

ggsave(p,file='QPCTBoxplot_male.pdf',height=5,width=5)
 
library(pROC)    
library(ggplot2)   
library(grid)      

 # ==== Data preparation ====
 
roc_df <- droplevels(subset(box_data, Diagnosis %in% c("ALS","CON")))
roc_df$Diagnosis <- factor(roc_df$Diagnosis, levels = c("CON","ALS"))  

 # Compute ROC curve object
roc_abca1 <- roc(response = roc_df$Diagnosis,
                 predictor = roc_df$ABCA1_Expression,   
                 levels = c("CON","ALS"),
                 direction = "auto", quiet = TRUE)

# Extract AUC value
auc_val <- as.numeric(auc(roc_abca1))

 # Find optimal threshold (Youden index) and extract common performance metrics
best <- coords(roc_abca1, x = "best", best.method = "youden",
               ret = c("threshold","sensitivity","specificity","accuracy"),
               transpose = FALSE)
sen <- as.numeric(best["sensitivity"])   
spe <- as.numeric(best["specificity"])   
acc <- as.numeric(best["accuracy"])     

# Extract points along the ROC curve (Specificity vs Sensitivity)
cd  <- coords(roc_abca1, x = "all", ret = c("specificity","sensitivity"), transpose = FALSE)
df_roc <- data.frame(Specificity = cd[, "specificity"], Sensitivity = cd[, "sensitivity"])

 # ==== Customize info box position (after x-axis reversal, bottom-right ≈ (0.05, 0.08)) ====
label_x <- 0.05   
label_y <- 0.08   

# ==== Plotting ====
options(repr.plot.width=6, repr.plot.height=6)  # Set plot window size

p <- ggplot(df_roc, aes(x = Specificity, y = Sensitivity)) +
  geom_line(linewidth = 1.4, color = "#2E86DE") +  # Plot ROC curve (blue)
  
# Add 45° gray dashed line (random classification reference)
  annotate("segment", x = 1, y = 0, xend = 0, yend = 1,
           linetype = "dashed", linewidth = 0.9, color = "grey50") +
  
# Add red dashed lines for Sensitivity and Specificity at optimal threshold
  geom_hline(yintercept = sen, linetype = "dashed", linewidth = 0.9, color = "red") +
  geom_vline(xintercept = spe, linetype = "dashed", linewidth = 0.9, color = "red") +
  
# Add info box (AUC, Accuracy, Sensitivity, Specificity) at bottom right of plot
  annotate("label", x = label_x, y = label_y, hjust = 1, vjust = 0,
           label = sprintf("AUC = %.3f\nAccuracy = %.3f\nSensitivity = %.3f\nSpecificity = %.3f",
                           auc_val, acc, sen, spe),
           size = 4.1, fill = "#F7F9FB", color = "black", label.size = 0.3) +
  
  # Add title, subtitle, and axis labels
  labs(
    title = "QPCT (Single-Gene ROC)",
    subtitle = "Male (ALS vs CON)",
    x = "Specificity",
    y = "Sensitivity"
  ) +
  
  # Reverse x-axis: Specificity 1 → 0, fixed ticks with padding at both ends
  scale_x_reverse(breaks = seq(0, 1, 0.25),
                  limits = c(1, 0),
                  expand = expansion(mult = c(0.02, 0.02))) +
  # Set y-axis range to 0–1 with 0.25 intervals
  scale_y_continuous(breaks = seq(0, 1, 0.25),
                     limits = c(0, 1),
                     expand = expansion(mult = c(0.02, 0.02))) +
  
# Set theme style
  theme_classic(base_size = 12) +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold", size = 15), # Center and bold main title
    plot.subtitle = element_text(hjust = 0.5, face = "bold", size = 11),  # Center and bold subtitle
    axis.title    = element_text(face = "bold", size = 13),               # Bold axis titles
    axis.text     = element_text(color = "black", size = 11),             # Customize axis tick labels
    axis.ticks    = element_line(linewidth = 0.9, color = "black"),       # Customize axis ticks
    axis.ticks.length = unit(0.22, "cm"),                                 # Adjust tick length
    panel.border  = element_rect(color = "black", fill = NA, linewidth = 1.2),  # Remove black border around the plot
    axis.line     = element_blank(),                                      # Remove default x/y axes to avoid double lines
    plot.margin   = margin(10, 12, 10, 12)                                # Margin
  ) +
  coord_fixed()  # Keep x:y ratio = 1:1

print(p)

# Save the figure
ggsave("QPCT_ROC_male.pdf", p, width = 6, height = 6)
ABCA1_exp <- valid_cohort_expr_female['QPCT',] 
colnames(ABCA1_exp) <- valid_cohort_clin_female$diagnosis
ABCA1_exp <- t(ABCA1_exp)
str(ABCA1_exp)
# Convert matrix to data frame
ABCA1_exp_df <- data.frame(Diagnosis = rownames(ABCA1_exp), ABCA1_Expression = as.numeric(ABCA1_exp[,1]))
# Ensure Diagnosis is a factor variable
ABCA1_exp_df$Diagnosis <- as.factor(ABCA1_exp_df$Diagnosis)
box_data <- subset(ABCA1_exp_df, Diagnosis %in% c("ALS", "CON"))
library(ggplot2)
library(ggpubr)

group_colors <- c("ALS" = "#9AC5CD", "CON" = "#FFDAB9")

p <- ggplot(box_data, aes(x = Diagnosis, y = ABCA1_Expression, fill = Diagnosis)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.1, alpha = 0.8, color = "black", size = 0.5) +
  labs(title = "QPCT Expression across Different Groups", 
       subtitle = "Female participants",
       x = "Female Diagnosis", 
       y = "QPCT Expression") +
  scale_fill_manual(values = group_colors) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        legend.position = "right") +
  stat_compare_means(method = "wilcox.test",
                     comparisons = list(c("ALS", "CON")),
                     label = "p.signif") # ⭐ Automatically add significance stars

print(p)

ggsave(p,file='QPCTBoxplot_female.pdf',height=5,width=5)
 
library(pROC)    
library(ggplot2)   
library(grid)      

 # ==== Data preparation ====
 
roc_df <- droplevels(subset(box_data, Diagnosis %in% c("ALS","CON")))
roc_df$Diagnosis <- factor(roc_df$Diagnosis, levels = c("CON","ALS"))  

 # Compute ROC curve object
roc_abca1 <- roc(response = roc_df$Diagnosis,
                 predictor = roc_df$ABCA1_Expression,   
                 levels = c("CON","ALS"),
                 direction = "auto", quiet = TRUE)

# Extract AUC value
auc_val <- as.numeric(auc(roc_abca1))

 # Find optimal threshold (Youden index) and extract common performance metrics
best <- coords(roc_abca1, x = "best", best.method = "youden",
               ret = c("threshold","sensitivity","specificity","accuracy"),
               transpose = FALSE)
sen <- as.numeric(best["sensitivity"])   
spe <- as.numeric(best["specificity"])   
acc <- as.numeric(best["accuracy"])     

# Extract points along the ROC curve (Specificity vs Sensitivity)
cd  <- coords(roc_abca1, x = "all", ret = c("specificity","sensitivity"), transpose = FALSE)
df_roc <- data.frame(Specificity = cd[, "specificity"], Sensitivity = cd[, "sensitivity"])

 # ==== Customize info box position (after x-axis reversal, bottom-right ≈ (0.05, 0.08)) ====
label_x <- 0.05   
label_y <- 0.08   

# ==== Plotting ====
options(repr.plot.width=6, repr.plot.height=6)  # Set plot window size

p <- ggplot(df_roc, aes(x = Specificity, y = Sensitivity)) +
  geom_line(linewidth = 1.4, color = "#2E86DE") +  # Plot ROC curve (blue)
  
# Add 45° gray dashed line (random classification reference)
  annotate("segment", x = 1, y = 0, xend = 0, yend = 1,
           linetype = "dashed", linewidth = 0.9, color = "grey50") +
  
# Add red dashed lines for Sensitivity and Specificity at optimal threshold
  geom_hline(yintercept = sen, linetype = "dashed", linewidth = 0.9, color = "red") +
  geom_vline(xintercept = spe, linetype = "dashed", linewidth = 0.9, color = "red") +
  
# Add info box (AUC, Accuracy, Sensitivity, Specificity) at bottom right of plot
  annotate("label", x = label_x, y = label_y, hjust = 1, vjust = 0,
           label = sprintf("AUC = %.3f\nAccuracy = %.3f\nSensitivity = %.3f\nSpecificity = %.3f",
                           auc_val, acc, sen, spe),
           size = 4.1, fill = "#F7F9FB", color = "black", label.size = 0.3) +
  
  # Add title, subtitle, and axis labels
  labs(
    title = "QPCT (Single-Gene ROC)",
    subtitle = "Female (ALS vs CON)",
    x = "Specificity",
    y = "Sensitivity"
  ) +
  
  # Reverse x-axis: Specificity 1 → 0, fixed ticks with padding at both ends
  scale_x_reverse(breaks = seq(0, 1, 0.25),
                  limits = c(1, 0),
                  expand = expansion(mult = c(0.02, 0.02))) +
  # Set y-axis range to 0–1 with 0.25 intervals
  scale_y_continuous(breaks = seq(0, 1, 0.25),
                     limits = c(0, 1),
                     expand = expansion(mult = c(0.02, 0.02))) +
  
# Set theme style
  theme_classic(base_size = 12) +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold", size = 15), # Center and bold main title
    plot.subtitle = element_text(hjust = 0.5, face = "bold", size = 11),  # Center and bold subtitle
    axis.title    = element_text(face = "bold", size = 13),               # Bold axis titles
    axis.text     = element_text(color = "black", size = 11),             # Customize axis tick labels
    axis.ticks    = element_line(linewidth = 0.9, color = "black"),       # Customize axis ticks
    axis.ticks.length = unit(0.22, "cm"),                                 # Adjust tick length
    panel.border  = element_rect(color = "black", fill = NA, linewidth = 1.2),  # Remove black border around the plot
    axis.line     = element_blank(),                                      # Remove default x/y axes to avoid double lines
    plot.margin   = margin(10, 12, 10, 12)                                # Margin
  ) +
  coord_fixed()  # Keep x:y ratio = 1:1

print(p)

# Save the figure
ggsave("QPCT_ROC_female.pdf", p, width = 6, height = 6)

############SRPK1
ABCA1_exp <- valid_cohort_expr_male['SRPK1',] 
colnames(ABCA1_exp) <- valid_cohort_clin_male$diagnosis
ABCA1_exp <- t(ABCA1_exp)
str(ABCA1_exp)
# Convert matrix to data frame
ABCA1_exp_df <- data.frame(Diagnosis = rownames(ABCA1_exp), ABCA1_Expression = as.numeric(ABCA1_exp[,1]))
# Ensure Diagnosis is a factor variable
ABCA1_exp_df$Diagnosis <- as.factor(ABCA1_exp_df$Diagnosis)
box_data <- subset(ABCA1_exp_df, Diagnosis %in% c("ALS", "CON"))
library(ggplot2)
library(ggpubr)

group_colors <- c("ALS" = "#9AC5CD", "CON" = "#FFDAB9")

p <- ggplot(box_data, aes(x = Diagnosis, y = ABCA1_Expression, fill = Diagnosis)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.1, alpha = 0.8, color = "black", size = 0.5) +
  labs(title = "SRPK1 Expression across Different Groups", 
       subtitle = "Male participants",
       x = "Diagnosis", 
       y = "SRPK1 Expression") +
  scale_fill_manual(values = group_colors) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        legend.position = "right") +
  stat_compare_means(method = "wilcox.test",
                     comparisons = list(c("ALS", "CON")),
                     label = "p.signif") # ⭐ Automatically add significance stars

print(p)
ggsave(p,file='SRPK1Boxplot_male.pdf',height=5,width=5)
 
library(pROC)    
library(ggplot2)   
library(grid)      

 # ==== Data preparation ====
 
roc_df <- droplevels(subset(box_data, Diagnosis %in% c("ALS","CON")))
roc_df$Diagnosis <- factor(roc_df$Diagnosis, levels = c("CON","ALS"))  

 # Compute ROC curve object
roc_abca1 <- roc(response = roc_df$Diagnosis,
                 predictor = roc_df$ABCA1_Expression,   
                 levels = c("CON","ALS"),
                 direction = "auto", quiet = TRUE)

# Extract AUC value
auc_val <- as.numeric(auc(roc_abca1))

 # Find optimal threshold (Youden index) and extract common performance metrics
best <- coords(roc_abca1, x = "best", best.method = "youden",
               ret = c("threshold","sensitivity","specificity","accuracy"),
               transpose = FALSE)
sen <- as.numeric(best["sensitivity"])   
spe <- as.numeric(best["specificity"])   
acc <- as.numeric(best["accuracy"])     

# Extract points along the ROC curve (Specificity vs Sensitivity)
cd  <- coords(roc_abca1, x = "all", ret = c("specificity","sensitivity"), transpose = FALSE)
df_roc <- data.frame(Specificity = cd[, "specificity"], Sensitivity = cd[, "sensitivity"])

 # ==== Customize info box position (after x-axis reversal, bottom-right ≈ (0.05, 0.08)) ====
label_x <- 0.05   
label_y <- 0.08   

# ==== Plotting ====
options(repr.plot.width=6, repr.plot.height=6)  # Set plot window size

p <- ggplot(df_roc, aes(x = Specificity, y = Sensitivity)) +
  geom_line(linewidth = 1.4, color = "#2E86DE") +  # Plot ROC curve (blue)
  
# Add 45° gray dashed line (random classification reference)
  annotate("segment", x = 1, y = 0, xend = 0, yend = 1,
           linetype = "dashed", linewidth = 0.9, color = "grey50") +
  
# Add red dashed lines for Sensitivity and Specificity at optimal threshold
  geom_hline(yintercept = sen, linetype = "dashed", linewidth = 0.9, color = "red") +
  geom_vline(xintercept = spe, linetype = "dashed", linewidth = 0.9, color = "red") +
  
# Add info box (AUC, Accuracy, Sensitivity, Specificity) at bottom right of plot
  annotate("label", x = label_x, y = label_y, hjust = 1, vjust = 0,
           label = sprintf("AUC = %.3f\nAccuracy = %.3f\nSensitivity = %.3f\nSpecificity = %.3f",
                           auc_val, acc, sen, spe),
           size = 4.1, fill = "#F7F9FB", color = "black", label.size = 0.3) +
  
  # Add title, subtitle, and axis labels
  labs(
    title = "SRPK1 (Single-Gene ROC)",
    subtitle = "Male (ALS vs CON)",
    x = "Specificity",
    y = "Sensitivity"
  ) +
  
  # Reverse x-axis: Specificity 1 → 0, fixed ticks with padding at both ends
  scale_x_reverse(breaks = seq(0, 1, 0.25),
                  limits = c(1, 0),
                  expand = expansion(mult = c(0.02, 0.02))) +
  # Set y-axis range to 0–1 with 0.25 intervals
  scale_y_continuous(breaks = seq(0, 1, 0.25),
                     limits = c(0, 1),
                     expand = expansion(mult = c(0.02, 0.02))) +
  
# Set theme style
  theme_classic(base_size = 12) +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold", size = 15), # Center and bold main title
    plot.subtitle = element_text(hjust = 0.5, face = "bold", size = 11),  # Center and bold subtitle
    axis.title    = element_text(face = "bold", size = 13),               # Bold axis titles
    axis.text     = element_text(color = "black", size = 11),             # Customize axis tick labels
    axis.ticks    = element_line(linewidth = 0.9, color = "black"),       # Customize axis ticks
    axis.ticks.length = unit(0.22, "cm"),                                 # Adjust tick length
    panel.border  = element_rect(color = "black", fill = NA, linewidth = 1.2),  # Remove black border around the plot
    axis.line     = element_blank(),                                      # Remove default x/y axes to avoid double lines
    plot.margin   = margin(10, 12, 10, 12)                                # Margin
  ) +
  coord_fixed()  # Keep x:y ratio = 1:1

print(p)

# Save the figure
ggsave("SRPK1_ROC_male.pdf", p, width = 6, height = 6)
ABCA1_exp <- valid_cohort_expr_female['SRPK1',] 
colnames(ABCA1_exp) <- valid_cohort_clin_female$diagnosis
ABCA1_exp <- t(ABCA1_exp)
str(ABCA1_exp)
# Convert matrix to data frame
ABCA1_exp_df <- data.frame(Diagnosis = rownames(ABCA1_exp), ABCA1_Expression = as.numeric(ABCA1_exp[,1]))
# Ensure Diagnosis is a factor variable
ABCA1_exp_df$Diagnosis <- as.factor(ABCA1_exp_df$Diagnosis)
box_data <- subset(ABCA1_exp_df, Diagnosis %in% c("ALS", "CON"))
library(ggplot2)
library(ggpubr)

group_colors <- c("ALS" = "#9AC5CD", "CON" = "#FFDAB9")

p <- ggplot(box_data, aes(x = Diagnosis, y = ABCA1_Expression, fill = Diagnosis)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.1, alpha = 0.8, color = "black", size = 0.5) +
  labs(title = "SRPK1 Expression across Different Groups", 
       subtitle = "Female participants",
       x = "Diagnosis", 
       y = "SRPK1 Expression") +
  scale_fill_manual(values = group_colors) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        legend.position = "right") +
  stat_compare_means(method = "wilcox.test",
                     comparisons = list(c("ALS", "CON")),
                     label = "p.signif") # ⭐ Automatically add significance stars

print(p)
ggsave(p,file='SRPK1Boxplot_female.pdf',height=5,width=5)
 
library(pROC)    
library(ggplot2)   
library(grid)      

 # ==== Data preparation ====
 
roc_df <- droplevels(subset(box_data, Diagnosis %in% c("ALS","CON")))
roc_df$Diagnosis <- factor(roc_df$Diagnosis, levels = c("CON","ALS"))  

 # Compute ROC curve object
roc_abca1 <- roc(response = roc_df$Diagnosis,
                 predictor = roc_df$ABCA1_Expression,   
                 levels = c("CON","ALS"),
                 direction = "auto", quiet = TRUE)

# Extract AUC value
auc_val <- as.numeric(auc(roc_abca1))

 # Find optimal threshold (Youden index) and extract common performance metrics
best <- coords(roc_abca1, x = "best", best.method = "youden",
               ret = c("threshold","sensitivity","specificity","accuracy"),
               transpose = FALSE)
sen <- as.numeric(best["sensitivity"])   
spe <- as.numeric(best["specificity"])   
acc <- as.numeric(best["accuracy"])     

# Extract points along the ROC curve (Specificity vs Sensitivity)
cd  <- coords(roc_abca1, x = "all", ret = c("specificity","sensitivity"), transpose = FALSE)
df_roc <- data.frame(Specificity = cd[, "specificity"], Sensitivity = cd[, "sensitivity"])

 # ==== Customize info box position (after x-axis reversal, bottom-right ≈ (0.05, 0.08)) ====
label_x <- 0.05   
label_y <- 0.08   

# ==== Plotting ====
options(repr.plot.width=6, repr.plot.height=6)  # Set plot window size

p <- ggplot(df_roc, aes(x = Specificity, y = Sensitivity)) +
  geom_line(linewidth = 1.4, color = "#2E86DE") +  # Plot ROC curve (blue)
  
# Add 45° gray dashed line (random classification reference)
  annotate("segment", x = 1, y = 0, xend = 0, yend = 1,
           linetype = "dashed", linewidth = 0.9, color = "grey50") +
  
# Add red dashed lines for Sensitivity and Specificity at optimal threshold
  geom_hline(yintercept = sen, linetype = "dashed", linewidth = 0.9, color = "red") +
  geom_vline(xintercept = spe, linetype = "dashed", linewidth = 0.9, color = "red") +
  
# Add info box (AUC, Accuracy, Sensitivity, Specificity) at bottom right of plot
  annotate("label", x = label_x, y = label_y, hjust = 1, vjust = 0,
           label = sprintf("AUC = %.3f\nAccuracy = %.3f\nSensitivity = %.3f\nSpecificity = %.3f",
                           auc_val, acc, sen, spe),
           size = 4.1, fill = "#F7F9FB", color = "black", label.size = 0.3) +
  
  # Add title, subtitle, and axis labels
  labs(
    title = "SRPK1 (Single-Gene ROC)",
    subtitle = "Female (ALS vs CON)",
    x = "Specificity",
    y = "Sensitivity"
  ) +
  
  # Reverse x-axis: Specificity 1 → 0, fixed ticks with padding at both ends
  scale_x_reverse(breaks = seq(0, 1, 0.25),
                  limits = c(1, 0),
                  expand = expansion(mult = c(0.02, 0.02))) +
  # Set y-axis range to 0–1 with 0.25 intervals
  scale_y_continuous(breaks = seq(0, 1, 0.25),
                     limits = c(0, 1),
                     expand = expansion(mult = c(0.02, 0.02))) +
  
# Set theme style
  theme_classic(base_size = 12) +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold", size = 15), # Center and bold main title
    plot.subtitle = element_text(hjust = 0.5, face = "bold", size = 11),  # Center and bold subtitle
    axis.title    = element_text(face = "bold", size = 13),               # Bold axis titles
    axis.text     = element_text(color = "black", size = 11),             # Customize axis tick labels
    axis.ticks    = element_line(linewidth = 0.9, color = "black"),       # Customize axis ticks
    axis.ticks.length = unit(0.22, "cm"),                                 # Adjust tick length
    panel.border  = element_rect(color = "black", fill = NA, linewidth = 1.2),  # Remove black border around the plot
    axis.line     = element_blank(),                                      # Remove default x/y axes to avoid double lines
    plot.margin   = margin(10, 12, 10, 12)                                # Margin
  ) +
  coord_fixed()  # Keep x:y ratio = 1:1

print(p)

# Save the figure
ggsave("SRPK1_ROC_female.pdf", p, width = 6, height = 6)


##################RPS6KA5
ABCA1_exp <- valid_cohort_expr_male['RPS6KA5',] 
colnames(ABCA1_exp) <- valid_cohort_clin_male$diagnosis
ABCA1_exp <- t(ABCA1_exp)
str(ABCA1_exp)
# Convert matrix to data frame
ABCA1_exp_df <- data.frame(Diagnosis = rownames(ABCA1_exp), ABCA1_Expression = as.numeric(ABCA1_exp[,1]))
# Ensure Diagnosis is a factor variable
ABCA1_exp_df$Diagnosis <- as.factor(ABCA1_exp_df$Diagnosis)
box_data <- subset(ABCA1_exp_df, Diagnosis %in% c("ALS", "CON"))
library(ggplot2)
library(ggpubr)

group_colors <- c("ALS" = "#9AC5CD", "CON" = "#FFDAB9")

p <- ggplot(box_data, aes(x = Diagnosis, y = ABCA1_Expression, fill = Diagnosis)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.1, alpha = 0.8, color = "black", size = 0.5) +
  labs(title = "RPS6KA5 Expression across Different Groups", 
       subtitle = "Male participants",
       x = "Diagnosis", 
       y = "RPS6KA5 Expression") +
  scale_fill_manual(values = group_colors) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        legend.position = "right") +
  stat_compare_means(method = "wilcox.test",
                     comparisons = list(c("ALS", "CON")),
                     label = "p.signif") # ⭐ Automatically add significance stars

print(p)
ggsave(p,file='RPS6KA5Boxplot_male.pdf',height=5,width=5)

 
library(pROC)    
library(ggplot2)   
library(grid)      

 # ==== Data preparation ====
 
roc_df <- droplevels(subset(box_data, Diagnosis %in% c("ALS","CON")))
roc_df$Diagnosis <- factor(roc_df$Diagnosis, levels = c("CON","ALS"))  

 # Compute ROC curve object
roc_abca1 <- roc(response = roc_df$Diagnosis,
                 predictor = roc_df$ABCA1_Expression,   
                 levels = c("CON","ALS"),
                 direction = "auto", quiet = TRUE)

# Extract AUC value
auc_val <- as.numeric(auc(roc_abca1))

 # Find optimal threshold (Youden index) and extract common performance metrics
best <- coords(roc_abca1, x = "best", best.method = "youden",
               ret = c("threshold","sensitivity","specificity","accuracy"),
               transpose = FALSE)
sen <- as.numeric(best["sensitivity"])   
spe <- as.numeric(best["specificity"])   
acc <- as.numeric(best["accuracy"])     

# Extract points along the ROC curve (Specificity vs Sensitivity)
cd  <- coords(roc_abca1, x = "all", ret = c("specificity","sensitivity"), transpose = FALSE)
df_roc <- data.frame(Specificity = cd[, "specificity"], Sensitivity = cd[, "sensitivity"])

 # ==== Customize info box position (after x-axis reversal, bottom-right ≈ (0.05, 0.08)) ====
label_x <- 0.05   
label_y <- 0.08   

# ==== Plotting ====
options(repr.plot.width=6, repr.plot.height=6)  # Set plot window size

p <- ggplot(df_roc, aes(x = Specificity, y = Sensitivity)) +
  geom_line(linewidth = 1.4, color = "#2E86DE") +  # Plot ROC curve (blue)
  
# Add 45° gray dashed line (random classification reference)
  annotate("segment", x = 1, y = 0, xend = 0, yend = 1,
           linetype = "dashed", linewidth = 0.9, color = "grey50") +
  
# Add red dashed lines for Sensitivity and Specificity at optimal threshold
  geom_hline(yintercept = sen, linetype = "dashed", linewidth = 0.9, color = "red") +
  geom_vline(xintercept = spe, linetype = "dashed", linewidth = 0.9, color = "red") +
  
# Add info box (AUC, Accuracy, Sensitivity, Specificity) at bottom right of plot
  annotate("label", x = label_x, y = label_y, hjust = 1, vjust = 0,
           label = sprintf("AUC = %.3f\nAccuracy = %.3f\nSensitivity = %.3f\nSpecificity = %.3f",
                           auc_val, acc, sen, spe),
           size = 4.1, fill = "#F7F9FB", color = "black", label.size = 0.3) +
  
  # Add title, subtitle, and axis labels
  labs(
    title = "RPS6KA5 (Single-Gene ROC)",
    subtitle = "Male (ALS vs CON)",
    x = "Specificity",
    y = "Sensitivity"
  ) +
  
  # Reverse x-axis: Specificity 1 → 0, fixed ticks with padding at both ends
  scale_x_reverse(breaks = seq(0, 1, 0.25),
                  limits = c(1, 0),
                  expand = expansion(mult = c(0.02, 0.02))) +
  # Set y-axis range to 0–1 with 0.25 intervals
  scale_y_continuous(breaks = seq(0, 1, 0.25),
                     limits = c(0, 1),
                     expand = expansion(mult = c(0.02, 0.02))) +
  
# Set theme style
  theme_classic(base_size = 12) +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold", size = 15), # Center and bold main title
    plot.subtitle = element_text(hjust = 0.5, face = "bold", size = 11),  # Center and bold subtitle
    axis.title    = element_text(face = "bold", size = 13),               # Bold axis titles
    axis.text     = element_text(color = "black", size = 11),             # Customize axis tick labels
    axis.ticks    = element_line(linewidth = 0.9, color = "black"),       # Customize axis ticks
    axis.ticks.length = unit(0.22, "cm"),                                 # Adjust tick length
    panel.border  = element_rect(color = "black", fill = NA, linewidth = 1.2),  # Remove black border around the plot
    axis.line     = element_blank(),                                      # Remove default x/y axes to avoid double lines
    plot.margin   = margin(10, 12, 10, 12)                                # Margin
  ) +
  coord_fixed()  # Keep x:y ratio = 1:1

print(p)

# Save the figure
ggsave("RPS6KA5_ROC_male.pdf", p, width = 6, height = 6)
ABCA1_exp <- valid_cohort_expr_female['RPS6KA5',] 
colnames(ABCA1_exp) <- valid_cohort_clin_female$diagnosis
ABCA1_exp <- t(ABCA1_exp)
str(ABCA1_exp)
# Convert matrix to data frame
ABCA1_exp_df <- data.frame(Diagnosis = rownames(ABCA1_exp), ABCA1_Expression = as.numeric(ABCA1_exp[,1]))
# Ensure Diagnosis is a factor variable
ABCA1_exp_df$Diagnosis <- as.factor(ABCA1_exp_df$Diagnosis)
box_data <- subset(ABCA1_exp_df, Diagnosis %in% c("ALS", "CON"))
library(ggplot2)
library(ggpubr)

group_colors <- c("ALS" = "#9AC5CD", "CON" = "#FFDAB9")

p <- ggplot(box_data, aes(x = Diagnosis, y = ABCA1_Expression, fill = Diagnosis)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.1, alpha = 0.8, color = "black", size = 0.5) +
  labs(title = "RPS6KA5 Expression across Different Groups", 
       subtitle = "Female participants",
       x = "Diagnosis", 
       y = "RPS6KA5 Expression") +
  scale_fill_manual(values = group_colors) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        legend.position = "right") +
  stat_compare_means(method = "wilcox.test",
                     comparisons = list(c("ALS", "CON")),
                     label = "p.signif") # ⭐ Automatically add significance stars

print(p)
ggsave(p,file='RPS6KA5Boxplot_female.pdf',height=5,width=5)

 
library(pROC)    
library(ggplot2)   
library(grid)      

 # ==== Data preparation ====
 
roc_df <- droplevels(subset(box_data, Diagnosis %in% c("ALS","CON")))
roc_df$Diagnosis <- factor(roc_df$Diagnosis, levels = c("CON","ALS"))  

 # Compute ROC curve object
roc_abca1 <- roc(response = roc_df$Diagnosis,
                 predictor = roc_df$ABCA1_Expression,   
                 levels = c("CON","ALS"),
                 direction = "auto", quiet = TRUE)

# Extract AUC value
auc_val <- as.numeric(auc(roc_abca1))

 # Find optimal threshold (Youden index) and extract common performance metrics
best <- coords(roc_abca1, x = "best", best.method = "youden",
               ret = c("threshold","sensitivity","specificity","accuracy"),
               transpose = FALSE)
sen <- as.numeric(best["sensitivity"])   
spe <- as.numeric(best["specificity"])   
acc <- as.numeric(best["accuracy"])     

# Extract points along the ROC curve (Specificity vs Sensitivity)
cd  <- coords(roc_abca1, x = "all", ret = c("specificity","sensitivity"), transpose = FALSE)
df_roc <- data.frame(Specificity = cd[, "specificity"], Sensitivity = cd[, "sensitivity"])

 # ==== Customize info box position (after x-axis reversal, bottom-right ≈ (0.05, 0.08)) ====
label_x <- 0.05   
label_y <- 0.08   

# ==== Plotting ====
options(repr.plot.width=6, repr.plot.height=6)  # Set plot window size

p <- ggplot(df_roc, aes(x = Specificity, y = Sensitivity)) +
  geom_line(linewidth = 1.4, color = "#2E86DE") +  # Plot ROC curve (blue)
  
# Add 45° gray dashed line (random classification reference)
  annotate("segment", x = 1, y = 0, xend = 0, yend = 1,
           linetype = "dashed", linewidth = 0.9, color = "grey50") +
  
# Add red dashed lines for Sensitivity and Specificity at optimal threshold
  geom_hline(yintercept = sen, linetype = "dashed", linewidth = 0.9, color = "red") +
  geom_vline(xintercept = spe, linetype = "dashed", linewidth = 0.9, color = "red") +
  
# Add info box (AUC, Accuracy, Sensitivity, Specificity) at bottom right of plot
  annotate("label", x = label_x, y = label_y, hjust = 1, vjust = 0,
           label = sprintf("AUC = %.3f\nAccuracy = %.3f\nSensitivity = %.3f\nSpecificity = %.3f",
                           auc_val, acc, sen, spe),
           size = 4.1, fill = "#F7F9FB", color = "black", label.size = 0.3) +
  
  # Add title, subtitle, and axis labels
  labs(
    title = "RPS6KA5 (Single-Gene ROC)",
    subtitle = "Female (ALS vs CON)",
    x = "Specificity",
    y = "Sensitivity"
  ) +
  
  # Reverse x-axis: Specificity 1 → 0, fixed ticks with padding at both ends
  scale_x_reverse(breaks = seq(0, 1, 0.25),
                  limits = c(1, 0),
                  expand = expansion(mult = c(0.02, 0.02))) +
  # Set y-axis range to 0–1 with 0.25 intervals
  scale_y_continuous(breaks = seq(0, 1, 0.25),
                     limits = c(0, 1),
                     expand = expansion(mult = c(0.02, 0.02))) +
  
# Set theme style
  theme_classic(base_size = 12) +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold", size = 15), # Center and bold main title
    plot.subtitle = element_text(hjust = 0.5, face = "bold", size = 11),  # Center and bold subtitle
    axis.title    = element_text(face = "bold", size = 13),               # Bold axis titles
    axis.text     = element_text(color = "black", size = 11),             # Customize axis tick labels
    axis.ticks    = element_line(linewidth = 0.9, color = "black"),       # Customize axis ticks
    axis.ticks.length = unit(0.22, "cm"),                                 # Adjust tick length
    panel.border  = element_rect(color = "black", fill = NA, linewidth = 1.2),  # Remove black border around the plot
    axis.line     = element_blank(),                                      # Remove default x/y axes to avoid double lines
    plot.margin   = margin(10, 12, 10, 12)                                # Margin
  ) +
  coord_fixed()  # Keep x:y ratio = 1:1

print(p)

# Save the figure
ggsave("RPS6KA5_ROC_female.pdf", p, width = 6, height = 6)



#################LYRM5
ABCA1_exp <- valid_cohort_expr_male['LYRM5',] 
colnames(ABCA1_exp) <- valid_cohort_clin_male$diagnosis
ABCA1_exp <- t(ABCA1_exp)
str(ABCA1_exp)
# Convert matrix to data frame
ABCA1_exp_df <- data.frame(Diagnosis = rownames(ABCA1_exp), ABCA1_Expression = as.numeric(ABCA1_exp[,1]))
# Ensure Diagnosis is a factor variable
ABCA1_exp_df$Diagnosis <- as.factor(ABCA1_exp_df$Diagnosis)
box_data <- subset(ABCA1_exp_df, Diagnosis %in% c("ALS", "CON"))
library(ggplot2)
library(ggpubr)

group_colors <- c("ALS" = "#9AC5CD", "CON" = "#FFDAB9")

p <- ggplot(box_data, aes(x = Diagnosis, y = ABCA1_Expression, fill = Diagnosis)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.1, alpha = 0.8, color = "black", size = 0.5) +
  labs(title = "LYRM5 Expression across Different Groups", 
       subtitle = "Male participants",
       x = "Diagnosis", 
       y = "LYRM5 Expression") +
  scale_fill_manual(values = group_colors) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        legend.position = "right") +
  stat_compare_means(method = "wilcox.test",
                     comparisons = list(c("ALS", "CON")),
                     label = "p.signif") # ⭐ Automatically add significance stars

print(p)
ggsave(p,file='LYRM5Boxplot_male.pdf',height=5,width=5)
 
library(pROC)    
library(ggplot2)   
library(grid)      

 # ==== Data preparation ====
 
roc_df <- droplevels(subset(box_data, Diagnosis %in% c("ALS","CON")))
roc_df$Diagnosis <- factor(roc_df$Diagnosis, levels = c("CON","ALS"))  

 # Compute ROC curve object
roc_abca1 <- roc(response = roc_df$Diagnosis,
                 predictor = roc_df$ABCA1_Expression,   
                 levels = c("CON","ALS"),
                 direction = "auto", quiet = TRUE)

# Extract AUC value
auc_val <- as.numeric(auc(roc_abca1))

 # Find optimal threshold (Youden index) and extract common performance metrics
best <- coords(roc_abca1, x = "best", best.method = "youden",
               ret = c("threshold","sensitivity","specificity","accuracy"),
               transpose = FALSE)
sen <- as.numeric(best["sensitivity"])   
spe <- as.numeric(best["specificity"])   
acc <- as.numeric(best["accuracy"])     

# Extract points along the ROC curve (Specificity vs Sensitivity)
cd  <- coords(roc_abca1, x = "all", ret = c("specificity","sensitivity"), transpose = FALSE)
df_roc <- data.frame(Specificity = cd[, "specificity"], Sensitivity = cd[, "sensitivity"])

 # ==== Customize info box position (after x-axis reversal, bottom-right ≈ (0.05, 0.08)) ====
label_x <- 0.05   
label_y <- 0.08   

# ==== Plotting ====
options(repr.plot.width=6, repr.plot.height=6)  # Set plot window size

p <- ggplot(df_roc, aes(x = Specificity, y = Sensitivity)) +
  geom_line(linewidth = 1.4, color = "#2E86DE") +  # Plot ROC curve (blue)
  
# Add 45° gray dashed line (random classification reference)
  annotate("segment", x = 1, y = 0, xend = 0, yend = 1,
           linetype = "dashed", linewidth = 0.9, color = "grey50") +
  
# Add red dashed lines for Sensitivity and Specificity at optimal threshold
  geom_hline(yintercept = sen, linetype = "dashed", linewidth = 0.9, color = "red") +
  geom_vline(xintercept = spe, linetype = "dashed", linewidth = 0.9, color = "red") +
  
# Add info box (AUC, Accuracy, Sensitivity, Specificity) at bottom right of plot
  annotate("label", x = label_x, y = label_y, hjust = 1, vjust = 0,
           label = sprintf("AUC = %.3f\nAccuracy = %.3f\nSensitivity = %.3f\nSpecificity = %.3f",
                           auc_val, acc, sen, spe),
           size = 4.1, fill = "#F7F9FB", color = "black", label.size = 0.3) +
  
  # Add title, subtitle, and axis labels
  labs(
    title = "LYRM5 (Single-Gene ROC)",
    subtitle = "Male (ALS vs CON)",
    x = "Specificity",
    y = "Sensitivity"
  ) +
  
  # Reverse x-axis: Specificity 1 → 0, fixed ticks with padding at both ends
  scale_x_reverse(breaks = seq(0, 1, 0.25),
                  limits = c(1, 0),
                  expand = expansion(mult = c(0.02, 0.02))) +
  # Set y-axis range to 0–1 with 0.25 intervals
  scale_y_continuous(breaks = seq(0, 1, 0.25),
                     limits = c(0, 1),
                     expand = expansion(mult = c(0.02, 0.02))) +
  
# Set theme style
  theme_classic(base_size = 12) +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold", size = 15), # Center and bold main title
    plot.subtitle = element_text(hjust = 0.5, face = "bold", size = 11),  # Center and bold subtitle
    axis.title    = element_text(face = "bold", size = 13),               # Bold axis titles
    axis.text     = element_text(color = "black", size = 11),             # Customize axis tick labels
    axis.ticks    = element_line(linewidth = 0.9, color = "black"),       # Customize axis ticks
    axis.ticks.length = unit(0.22, "cm"),                                 # Adjust tick length
    panel.border  = element_rect(color = "black", fill = NA, linewidth = 1.2),  # Remove black border around the plot
    axis.line     = element_blank(),                                      # Remove default x/y axes to avoid double lines
    plot.margin   = margin(10, 12, 10, 12)                                # Margin
  ) +
  coord_fixed()  # Keep x:y ratio = 1:1

print(p)

# Save the figure
ggsave("LYRM5_ROC_male.pdf", p, width = 6, height = 6)
ABCA1_exp <- valid_cohort_expr_female['LYRM5',] 
colnames(ABCA1_exp) <- valid_cohort_clin_female$diagnosis
ABCA1_exp <- t(ABCA1_exp)
str(ABCA1_exp)
# Convert matrix to data frame
ABCA1_exp_df <- data.frame(Diagnosis = rownames(ABCA1_exp), ABCA1_Expression = as.numeric(ABCA1_exp[,1]))
# Ensure Diagnosis is a factor variable
ABCA1_exp_df$Diagnosis <- as.factor(ABCA1_exp_df$Diagnosis)
box_data <- subset(ABCA1_exp_df, Diagnosis %in% c("ALS", "CON"))
library(ggplot2)
library(ggpubr)

group_colors <- c("ALS" = "#9AC5CD", "CON" = "#FFDAB9")

p <- ggplot(box_data, aes(x = Diagnosis, y = ABCA1_Expression, fill = Diagnosis)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.1, alpha = 0.8, color = "black", size = 0.5) +
  labs(title = "LYRM5 Expression across Different Groups", 
       subtitle = "Female participants",
       x = "Diagnosis", 
       y = "LYRM5 Expression") +
  scale_fill_manual(values = group_colors) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        legend.position = "right") +
  stat_compare_means(method = "wilcox.test",
                     comparisons = list(c("ALS", "CON")),
                     label = "p.signif") # ⭐ Automatically add significance stars

print(p)
ggsave(p,file='LYRM5Boxplot_female.pdf',height=5,width=5)
 
library(pROC)    
library(ggplot2)   
library(grid)      

 # ==== Data preparation ====
 
roc_df <- droplevels(subset(box_data, Diagnosis %in% c("ALS","CON")))
roc_df$Diagnosis <- factor(roc_df$Diagnosis, levels = c("CON","ALS"))  

 # Compute ROC curve object
roc_abca1 <- roc(response = roc_df$Diagnosis,
                 predictor = roc_df$ABCA1_Expression,   
                 levels = c("CON","ALS"),
                 direction = "auto", quiet = TRUE)

# Extract AUC value
auc_val <- as.numeric(auc(roc_abca1))

 # Find optimal threshold (Youden index) and extract common performance metrics
best <- coords(roc_abca1, x = "best", best.method = "youden",
               ret = c("threshold","sensitivity","specificity","accuracy"),
               transpose = FALSE)
sen <- as.numeric(best["sensitivity"])   
spe <- as.numeric(best["specificity"])   
acc <- as.numeric(best["accuracy"])     

# Extract points along the ROC curve (Specificity vs Sensitivity)
cd  <- coords(roc_abca1, x = "all", ret = c("specificity","sensitivity"), transpose = FALSE)
df_roc <- data.frame(Specificity = cd[, "specificity"], Sensitivity = cd[, "sensitivity"])

 # ==== Customize info box position (after x-axis reversal, bottom-right ≈ (0.05, 0.08)) ====
label_x <- 0.05   
label_y <- 0.08   

# ==== Plotting ====
options(repr.plot.width=6, repr.plot.height=6)  # Set plot window size

p <- ggplot(df_roc, aes(x = Specificity, y = Sensitivity)) +
  geom_line(linewidth = 1.4, color = "#2E86DE") +  # Plot ROC curve (blue)
  
# Add 45° gray dashed line (random classification reference)
  annotate("segment", x = 1, y = 0, xend = 0, yend = 1,
           linetype = "dashed", linewidth = 0.9, color = "grey50") +
  
# Add red dashed lines for Sensitivity and Specificity at optimal threshold
  geom_hline(yintercept = sen, linetype = "dashed", linewidth = 0.9, color = "red") +
  geom_vline(xintercept = spe, linetype = "dashed", linewidth = 0.9, color = "red") +
  
# Add info box (AUC, Accuracy, Sensitivity, Specificity) at bottom right of plot
  annotate("label", x = label_x, y = label_y, hjust = 1, vjust = 0,
           label = sprintf("AUC = %.3f\nAccuracy = %.3f\nSensitivity = %.3f\nSpecificity = %.3f",
                           auc_val, acc, sen, spe),
           size = 4.1, fill = "#F7F9FB", color = "black", label.size = 0.3) +
  
  # Add title, subtitle, and axis labels
  labs(
    title = "LYRM5 (Single-Gene ROC)",
    subtitle = "Female (ALS vs CON)",
    x = "Specificity",
    y = "Sensitivity"
  ) +
  
  # Reverse x-axis: Specificity 1 → 0, fixed ticks with padding at both ends
  scale_x_reverse(breaks = seq(0, 1, 0.25),
                  limits = c(1, 0),
                  expand = expansion(mult = c(0.02, 0.02))) +
  # Set y-axis range to 0–1 with 0.25 intervals
  scale_y_continuous(breaks = seq(0, 1, 0.25),
                     limits = c(0, 1),
                     expand = expansion(mult = c(0.02, 0.02))) +
  
# Set theme style
  theme_classic(base_size = 12) +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold", size = 15), # Center and bold main title
    plot.subtitle = element_text(hjust = 0.5, face = "bold", size = 11),  # Center and bold subtitle
    axis.title    = element_text(face = "bold", size = 13),               # Bold axis titles
    axis.text     = element_text(color = "black", size = 11),             # Customize axis tick labels
    axis.ticks    = element_line(linewidth = 0.9, color = "black"),       # Customize axis ticks
    axis.ticks.length = unit(0.22, "cm"),                                 # Adjust tick length
    panel.border  = element_rect(color = "black", fill = NA, linewidth = 1.2),  # Remove black border around the plot
    axis.line     = element_blank(),                                      # Remove default x/y axes to avoid double lines
    plot.margin   = margin(10, 12, 10, 12)                                # Margin
  ) +
  coord_fixed()  # Keep x:y ratio = 1:1

print(p)

# Save the figure
ggsave("LYRM5_ROC_female.pdf", p, width = 6, height = 6)




###################SLC25A20
ABCA1_exp <- valid_cohort_expr_male['SLC25A20',] 
colnames(ABCA1_exp) <- valid_cohort_clin_male$diagnosis
ABCA1_exp <- t(ABCA1_exp)
str(ABCA1_exp)
# Convert matrix to data frame
ABCA1_exp_df <- data.frame(Diagnosis = rownames(ABCA1_exp), ABCA1_Expression = as.numeric(ABCA1_exp[,1]))
# Ensure Diagnosis is a factor variable
ABCA1_exp_df$Diagnosis <- as.factor(ABCA1_exp_df$Diagnosis)
box_data <- subset(ABCA1_exp_df, Diagnosis %in% c("ALS", "CON"))
library(ggplot2)
library(ggpubr)

group_colors <- c("ALS" = "#9AC5CD", "CON" = "#FFDAB9")

p <- ggplot(box_data, aes(x = Diagnosis, y = ABCA1_Expression, fill = Diagnosis)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.1, alpha = 0.8, color = "black", size = 0.5) +
  labs(title = "SLC25A20 Expression across Different Groups", 
       subtitle = "Male participants",
       x = "Diagnosis", 
       y = "SLC25A20 Expression") +
  scale_fill_manual(values = group_colors) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        legend.position = "right") +
  stat_compare_means(method = "wilcox.test",
                     comparisons = list(c("ALS", "CON")),
                     label = "p.signif") # ⭐ Automatically add significance stars

print(p)
ggsave(p,file='SLC25A20Boxplot_male.pdf',height=5,width=5)
 
library(pROC)    
library(ggplot2)   
library(grid)      

 # ==== Data preparation ====
 
roc_df <- droplevels(subset(box_data, Diagnosis %in% c("ALS","CON")))
roc_df$Diagnosis <- factor(roc_df$Diagnosis, levels = c("CON","ALS"))  

 # Compute ROC curve object
roc_abca1 <- roc(response = roc_df$Diagnosis,
                 predictor = roc_df$ABCA1_Expression,   
                 levels = c("CON","ALS"),
                 direction = "auto", quiet = TRUE)

# Extract AUC value
auc_val <- as.numeric(auc(roc_abca1))

 # Find optimal threshold (Youden index) and extract common performance metrics
best <- coords(roc_abca1, x = "best", best.method = "youden",
               ret = c("threshold","sensitivity","specificity","accuracy"),
               transpose = FALSE)
sen <- as.numeric(best["sensitivity"])   
spe <- as.numeric(best["specificity"])   
acc <- as.numeric(best["accuracy"])     

# Extract points along the ROC curve (Specificity vs Sensitivity)
cd  <- coords(roc_abca1, x = "all", ret = c("specificity","sensitivity"), transpose = FALSE)
df_roc <- data.frame(Specificity = cd[, "specificity"], Sensitivity = cd[, "sensitivity"])

 # ==== Customize info box position (after x-axis reversal, bottom-right ≈ (0.05, 0.08)) ====
label_x <- 0.05   
label_y <- 0.08   

# ==== Plotting ====
options(repr.plot.width=6, repr.plot.height=6)  # Set plot window size

p <- ggplot(df_roc, aes(x = Specificity, y = Sensitivity)) +
  geom_line(linewidth = 1.4, color = "#2E86DE") +  # Plot ROC curve (blue)
  
# Add 45° gray dashed line (random classification reference)
  annotate("segment", x = 1, y = 0, xend = 0, yend = 1,
           linetype = "dashed", linewidth = 0.9, color = "grey50") +
  
# Add red dashed lines for Sensitivity and Specificity at optimal threshold
  geom_hline(yintercept = sen, linetype = "dashed", linewidth = 0.9, color = "red") +
  geom_vline(xintercept = spe, linetype = "dashed", linewidth = 0.9, color = "red") +
  
# Add info box (AUC, Accuracy, Sensitivity, Specificity) at bottom right of plot
  annotate("label", x = label_x, y = label_y, hjust = 1, vjust = 0,
           label = sprintf("AUC = %.3f\nAccuracy = %.3f\nSensitivity = %.3f\nSpecificity = %.3f",
                           auc_val, acc, sen, spe),
           size = 4.1, fill = "#F7F9FB", color = "black", label.size = 0.3) +
  
  # Add title, subtitle, and axis labels
  labs(
    title = "SLC25A20 (Single-Gene ROC)",
    subtitle = "Male (ALS vs CON)",
    x = "Specificity",
    y = "Sensitivity"
  ) +
  
  # Reverse x-axis: Specificity 1 → 0, fixed ticks with padding at both ends
  scale_x_reverse(breaks = seq(0, 1, 0.25),
                  limits = c(1, 0),
                  expand = expansion(mult = c(0.02, 0.02))) +
  # Set y-axis range to 0–1 with 0.25 intervals
  scale_y_continuous(breaks = seq(0, 1, 0.25),
                     limits = c(0, 1),
                     expand = expansion(mult = c(0.02, 0.02))) +
  
# Set theme style
  theme_classic(base_size = 12) +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold", size = 15), # Center and bold main title
    plot.subtitle = element_text(hjust = 0.5, face = "bold", size = 11),  # Center and bold subtitle
    axis.title    = element_text(face = "bold", size = 13),               # Bold axis titles
    axis.text     = element_text(color = "black", size = 11),             # Customize axis tick labels
    axis.ticks    = element_line(linewidth = 0.9, color = "black"),       # Customize axis ticks
    axis.ticks.length = unit(0.22, "cm"),                                 # Adjust tick length
    panel.border  = element_rect(color = "black", fill = NA, linewidth = 1.2),  # Remove black border around the plot
    axis.line     = element_blank(),                                      # Remove default x/y axes to avoid double lines
    plot.margin   = margin(10, 12, 10, 12)                                # Margin
  ) +
  coord_fixed()  # Keep x:y ratio = 1:1

print(p)

# Save the figure
ggsave("SLC25A20_ROC_male.pdf", p, width = 6, height = 6)
ABCA1_exp <- valid_cohort_expr_female['SLC25A20',] 
colnames(ABCA1_exp) <- valid_cohort_clin_female$diagnosis
ABCA1_exp <- t(ABCA1_exp)
str(ABCA1_exp)
# Convert matrix to data frame
ABCA1_exp_df <- data.frame(Diagnosis = rownames(ABCA1_exp), ABCA1_Expression = as.numeric(ABCA1_exp[,1]))
# Ensure Diagnosis is a factor variable
ABCA1_exp_df$Diagnosis <- as.factor(ABCA1_exp_df$Diagnosis)
box_data <- subset(ABCA1_exp_df, Diagnosis %in% c("ALS", "CON"))
library(ggplot2)
library(ggpubr)

group_colors <- c("ALS" = "#9AC5CD", "CON" = "#FFDAB9")

p <- ggplot(box_data, aes(x = Diagnosis, y = ABCA1_Expression, fill = Diagnosis)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.1, alpha = 0.8, color = "black", size = 0.5) +
  labs(title = "SLC25A20 Expression across Different Groups", 
       subtitle = "Female participants",
       x = "Diagnosis", 
       y = "SLC25A20 Expression") +
  scale_fill_manual(values = group_colors) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        legend.position = "right") +
  stat_compare_means(method = "wilcox.test",
                     comparisons = list(c("ALS", "CON")),
                     label = "p.signif") # ⭐ Automatically add significance stars

print(p)
ggsave(p,file='SLC25A20Boxplot_female.pdf',height=5,width=5)
 
library(pROC)    
library(ggplot2)   
library(grid)      

 # ==== Data preparation ====
 
roc_df <- droplevels(subset(box_data, Diagnosis %in% c("ALS","CON")))
roc_df$Diagnosis <- factor(roc_df$Diagnosis, levels = c("CON","ALS"))  

 # Compute ROC curve object
roc_abca1 <- roc(response = roc_df$Diagnosis,
                 predictor = roc_df$ABCA1_Expression,   
                 levels = c("CON","ALS"),
                 direction = "auto", quiet = TRUE)

# Extract AUC value
auc_val <- as.numeric(auc(roc_abca1))

 # Find optimal threshold (Youden index) and extract common performance metrics
best <- coords(roc_abca1, x = "best", best.method = "youden",
               ret = c("threshold","sensitivity","specificity","accuracy"),
               transpose = FALSE)
sen <- as.numeric(best["sensitivity"])   
spe <- as.numeric(best["specificity"])   
acc <- as.numeric(best["accuracy"])     

# Extract points along the ROC curve (Specificity vs Sensitivity)
cd  <- coords(roc_abca1, x = "all", ret = c("specificity","sensitivity"), transpose = FALSE)
df_roc <- data.frame(Specificity = cd[, "specificity"], Sensitivity = cd[, "sensitivity"])

 # ==== Customize info box position (after x-axis reversal, bottom-right ≈ (0.05, 0.08)) ====
label_x <- 0.05   
label_y <- 0.08   

# ==== Plotting ====
options(repr.plot.width=6, repr.plot.height=6)  # Set plot window size

p <- ggplot(df_roc, aes(x = Specificity, y = Sensitivity)) +
  geom_line(linewidth = 1.4, color = "#2E86DE") +  # Plot ROC curve (blue)
  
# Add 45° gray dashed line (random classification reference)
  annotate("segment", x = 1, y = 0, xend = 0, yend = 1,
           linetype = "dashed", linewidth = 0.9, color = "grey50") +
  
# Add red dashed lines for Sensitivity and Specificity at optimal threshold
  geom_hline(yintercept = sen, linetype = "dashed", linewidth = 0.9, color = "red") +
  geom_vline(xintercept = spe, linetype = "dashed", linewidth = 0.9, color = "red") +
  
# Add info box (AUC, Accuracy, Sensitivity, Specificity) at bottom right of plot
  annotate("label", x = label_x, y = label_y, hjust = 1, vjust = 0,
           label = sprintf("AUC = %.3f\nAccuracy = %.3f\nSensitivity = %.3f\nSpecificity = %.3f",
                           auc_val, acc, sen, spe),
           size = 4.1, fill = "#F7F9FB", color = "black", label.size = 0.3) +
  
  # Add title, subtitle, and axis labels
  labs(
    title = "SLC25A20 (Single-Gene ROC)",
    subtitle = "Female (ALS vs CON)",
    x = "Specificity",
    y = "Sensitivity"
  ) +
  
  # Reverse x-axis: Specificity 1 → 0, fixed ticks with padding at both ends
  scale_x_reverse(breaks = seq(0, 1, 0.25),
                  limits = c(1, 0),
                  expand = expansion(mult = c(0.02, 0.02))) +
  # Set y-axis range to 0–1 with 0.25 intervals
  scale_y_continuous(breaks = seq(0, 1, 0.25),
                     limits = c(0, 1),
                     expand = expansion(mult = c(0.02, 0.02))) +
  
# Set theme style
  theme_classic(base_size = 12) +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold", size = 15), # Center and bold main title
    plot.subtitle = element_text(hjust = 0.5, face = "bold", size = 11),  # Center and bold subtitle
    axis.title    = element_text(face = "bold", size = 13),               # Bold axis titles
    axis.text     = element_text(color = "black", size = 11),             # Customize axis tick labels
    axis.ticks    = element_line(linewidth = 0.9, color = "black"),       # Customize axis ticks
    axis.ticks.length = unit(0.22, "cm"),                                 # Adjust tick length
    panel.border  = element_rect(color = "black", fill = NA, linewidth = 1.2),  # Remove black border around the plot
    axis.line     = element_blank(),                                      # Remove default x/y axes to avoid double lines
    plot.margin   = margin(10, 12, 10, 12)                                # Margin
  ) +
  coord_fixed()  # Keep x:y ratio = 1:1

print(p)

# Save the figure
ggsave("SLC25A20_ROC_female.pdf", p, width = 6, height = 6)




###################NT5DC1
ABCA1_exp <- valid_cohort_expr_male['NT5DC1',] 
colnames(ABCA1_exp) <- valid_cohort_clin_male$diagnosis
ABCA1_exp <- t(ABCA1_exp)
str(ABCA1_exp)
# Convert matrix to data frame
ABCA1_exp_df <- data.frame(Diagnosis = rownames(ABCA1_exp), ABCA1_Expression = as.numeric(ABCA1_exp[,1]))
# Ensure Diagnosis is a factor variable
ABCA1_exp_df$Diagnosis <- as.factor(ABCA1_exp_df$Diagnosis)
box_data <- subset(ABCA1_exp_df, Diagnosis %in% c("ALS", "CON"))
library(ggplot2)
library(ggpubr)

group_colors <- c("ALS" = "#9AC5CD", "CON" = "#FFDAB9")

p <- ggplot(box_data, aes(x = Diagnosis, y = ABCA1_Expression, fill = Diagnosis)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.1, alpha = 0.8, color = "black", size = 0.5) +
  labs(title = "NT5DC1 Expression across Different Groups", 
       subtitle = "Male participants",
       x = "Diagnosis", 
       y = "NT5DC1 Expression") +
  scale_fill_manual(values = group_colors) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        legend.position = "right") +
  stat_compare_means(method = "wilcox.test",
                     comparisons = list(c("ALS", "CON")),
                     label = "p.signif") # ⭐ Automatically add significance stars

print(p)
ggsave(p,file='NT5DC1Boxplot_male.pdf',height=5,width=5)
 
library(pROC)    
library(ggplot2)   
library(grid)      

 # ==== Data preparation ====
 
roc_df <- droplevels(subset(box_data, Diagnosis %in% c("ALS","CON")))
roc_df$Diagnosis <- factor(roc_df$Diagnosis, levels = c("CON","ALS"))  

 # Compute ROC curve object
roc_abca1 <- roc(response = roc_df$Diagnosis,
                 predictor = roc_df$ABCA1_Expression,   
                 levels = c("CON","ALS"),
                 direction = "auto", quiet = TRUE)

# Extract AUC value
auc_val <- as.numeric(auc(roc_abca1))

 # Find optimal threshold (Youden index) and extract common performance metrics
best <- coords(roc_abca1, x = "best", best.method = "youden",
               ret = c("threshold","sensitivity","specificity","accuracy"),
               transpose = FALSE)
sen <- as.numeric(best["sensitivity"])   
spe <- as.numeric(best["specificity"])   
acc <- as.numeric(best["accuracy"])     

# Extract points along the ROC curve (Specificity vs Sensitivity)
cd  <- coords(roc_abca1, x = "all", ret = c("specificity","sensitivity"), transpose = FALSE)
df_roc <- data.frame(Specificity = cd[, "specificity"], Sensitivity = cd[, "sensitivity"])

 # ==== Customize info box position (after x-axis reversal, bottom-right ≈ (0.05, 0.08)) ====
label_x <- 0.05   
label_y <- 0.08   

# ==== Plotting ====
options(repr.plot.width=6, repr.plot.height=6)  # Set plot window size

p <- ggplot(df_roc, aes(x = Specificity, y = Sensitivity)) +
  geom_line(linewidth = 1.4, color = "#2E86DE") +  # Plot ROC curve (blue)
  
# Add 45° gray dashed line (random classification reference)
  annotate("segment", x = 1, y = 0, xend = 0, yend = 1,
           linetype = "dashed", linewidth = 0.9, color = "grey50") +
  
# Add red dashed lines for Sensitivity and Specificity at optimal threshold
  geom_hline(yintercept = sen, linetype = "dashed", linewidth = 0.9, color = "red") +
  geom_vline(xintercept = spe, linetype = "dashed", linewidth = 0.9, color = "red") +
  
# Add info box (AUC, Accuracy, Sensitivity, Specificity) at bottom right of plot
  annotate("label", x = label_x, y = label_y, hjust = 1, vjust = 0,
           label = sprintf("AUC = %.3f\nAccuracy = %.3f\nSensitivity = %.3f\nSpecificity = %.3f",
                           auc_val, acc, sen, spe),
           size = 4.1, fill = "#F7F9FB", color = "black", label.size = 0.3) +
  
  # Add title, subtitle, and axis labels
  labs(
    title = "NT5DC1 (Single-Gene ROC)",
    subtitle = "Male (ALS vs CON)",
    x = "Specificity",
    y = "Sensitivity"
  ) +
  
  # Reverse x-axis: Specificity 1 → 0, fixed ticks with padding at both ends
  scale_x_reverse(breaks = seq(0, 1, 0.25),
                  limits = c(1, 0),
                  expand = expansion(mult = c(0.02, 0.02))) +
  # Set y-axis range to 0–1 with 0.25 intervals
  scale_y_continuous(breaks = seq(0, 1, 0.25),
                     limits = c(0, 1),
                     expand = expansion(mult = c(0.02, 0.02))) +
  
# Set theme style
  theme_classic(base_size = 12) +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold", size = 15), # Center and bold main title
    plot.subtitle = element_text(hjust = 0.5, face = "bold", size = 11),  # Center and bold subtitle
    axis.title    = element_text(face = "bold", size = 13),               # Bold axis titles
    axis.text     = element_text(color = "black", size = 11),             # Customize axis tick labels
    axis.ticks    = element_line(linewidth = 0.9, color = "black"),       # Customize axis ticks
    axis.ticks.length = unit(0.22, "cm"),                                 # Adjust tick length
    panel.border  = element_rect(color = "black", fill = NA, linewidth = 1.2),  # Remove black border around the plot
    axis.line     = element_blank(),                                      # Remove default x/y axes to avoid double lines
    plot.margin   = margin(10, 12, 10, 12)                                # Margin
  ) +
  coord_fixed()  # Keep x:y ratio = 1:1

print(p)

# Save the figure
ggsave("NT5DC1_ROC_male.pdf", p, width = 6, height = 6)
ABCA1_exp <- valid_cohort_expr_female['NT5DC1',] 
colnames(ABCA1_exp) <- valid_cohort_clin_female$diagnosis
ABCA1_exp <- t(ABCA1_exp)
str(ABCA1_exp)
# Convert matrix to data frame
ABCA1_exp_df <- data.frame(Diagnosis = rownames(ABCA1_exp), ABCA1_Expression = as.numeric(ABCA1_exp[,1]))
# Ensure Diagnosis is a factor variable
ABCA1_exp_df$Diagnosis <- as.factor(ABCA1_exp_df$Diagnosis)
box_data <- subset(ABCA1_exp_df, Diagnosis %in% c("ALS", "CON"))
library(ggplot2)
library(ggpubr)

group_colors <- c("ALS" = "#9AC5CD", "CON" = "#FFDAB9")

p <- ggplot(box_data, aes(x = Diagnosis, y = ABCA1_Expression, fill = Diagnosis)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.1, alpha = 0.8, color = "black", size = 0.5) +
  labs(title = "NT5DC1 Expression across Different Groups", 
       subtitle = "Female participants",
       x = "Diagnosis", 
       y = "NT5DC1 Expression") +
  scale_fill_manual(values = group_colors) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        legend.position = "right") +
  stat_compare_means(method = "wilcox.test",
                     comparisons = list(c("ALS", "CON")),
                     label = "p.signif") # ⭐ Automatically add significance stars

print(p)
ggsave(p,file='NT5DC1Boxplot_female.pdf',height=5,width=5)
 
library(pROC)    
library(ggplot2)   
library(grid)      

 # ==== Data preparation ====
 
roc_df <- droplevels(subset(box_data, Diagnosis %in% c("ALS","CON")))
roc_df$Diagnosis <- factor(roc_df$Diagnosis, levels = c("CON","ALS"))  

 # Compute ROC curve object
roc_abca1 <- roc(response = roc_df$Diagnosis,
                 predictor = roc_df$ABCA1_Expression,   
                 levels = c("CON","ALS"),
                 direction = "auto", quiet = TRUE)

# Extract AUC value
auc_val <- as.numeric(auc(roc_abca1))

 # Find optimal threshold (Youden index) and extract common performance metrics
best <- coords(roc_abca1, x = "best", best.method = "youden",
               ret = c("threshold","sensitivity","specificity","accuracy"),
               transpose = FALSE)
sen <- as.numeric(best["sensitivity"])   
spe <- as.numeric(best["specificity"])   
acc <- as.numeric(best["accuracy"])     

# Extract points along the ROC curve (Specificity vs Sensitivity)
cd  <- coords(roc_abca1, x = "all", ret = c("specificity","sensitivity"), transpose = FALSE)
df_roc <- data.frame(Specificity = cd[, "specificity"], Sensitivity = cd[, "sensitivity"])

 # ==== Customize info box position (after x-axis reversal, bottom-right ≈ (0.05, 0.08)) ====
label_x <- 0.05   
label_y <- 0.08   

# ==== Plotting ====
options(repr.plot.width=6, repr.plot.height=6)  # Set plot window size

p <- ggplot(df_roc, aes(x = Specificity, y = Sensitivity)) +
  geom_line(linewidth = 1.4, color = "#2E86DE") +  # Plot ROC curve (blue)
  
# Add 45° gray dashed line (random classification reference)
  annotate("segment", x = 1, y = 0, xend = 0, yend = 1,
           linetype = "dashed", linewidth = 0.9, color = "grey50") +
  
# Add red dashed lines for Sensitivity and Specificity at optimal threshold
  geom_hline(yintercept = sen, linetype = "dashed", linewidth = 0.9, color = "red") +
  geom_vline(xintercept = spe, linetype = "dashed", linewidth = 0.9, color = "red") +
  
# Add info box (AUC, Accuracy, Sensitivity, Specificity) at bottom right of plot
  annotate("label", x = label_x, y = label_y, hjust = 1, vjust = 0,
           label = sprintf("AUC = %.3f\nAccuracy = %.3f\nSensitivity = %.3f\nSpecificity = %.3f",
                           auc_val, acc, sen, spe),
           size = 4.1, fill = "#F7F9FB", color = "black", label.size = 0.3) +
  
  # Add title, subtitle, and axis labels
  labs(
    title = "NT5DC1 (Single-Gene ROC)",
    subtitle = "Female (ALS vs CON)",
    x = "Specificity",
    y = "Sensitivity"
  ) +
  
  # Reverse x-axis: Specificity 1 → 0, fixed ticks with padding at both ends
  scale_x_reverse(breaks = seq(0, 1, 0.25),
                  limits = c(1, 0),
                  expand = expansion(mult = c(0.02, 0.02))) +
  # Set y-axis range to 0–1 with 0.25 intervals
  scale_y_continuous(breaks = seq(0, 1, 0.25),
                     limits = c(0, 1),
                     expand = expansion(mult = c(0.02, 0.02))) +
  
# Set theme style
  theme_classic(base_size = 12) +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold", size = 15), # Center and bold main title
    plot.subtitle = element_text(hjust = 0.5, face = "bold", size = 11),  # Center and bold subtitle
    axis.title    = element_text(face = "bold", size = 13),               # Bold axis titles
    axis.text     = element_text(color = "black", size = 11),             # Customize axis tick labels
    axis.ticks    = element_line(linewidth = 0.9, color = "black"),       # Customize axis ticks
    axis.ticks.length = unit(0.22, "cm"),                                 # Adjust tick length
    panel.border  = element_rect(color = "black", fill = NA, linewidth = 1.2),  # Remove black border around the plot
    axis.line     = element_blank(),                                      # Remove default x/y axes to avoid double lines
    plot.margin   = margin(10, 12, 10, 12)                                # Margin
  ) +
  coord_fixed()  # Keep x:y ratio = 1:1

print(p)

# Save the figure
ggsave("NT5DC1_ROC_female.pdf", p, width = 6, height = 6)
