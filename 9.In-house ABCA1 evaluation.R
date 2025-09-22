############################################
# 0) Libraries and Working Directory
############################################
library(mice)
library(dplyr)
library(broom)
library(readr)
library(ggplot2)
library(ggpubr)
library(tidyr)
setwd('In-house ABCA1 evaluation')


############################################
# 1) Data Loading and Basic Preprocessing
############################################
### Clean and preprocess the 'data' object
data <- read_excel("Serum_ABCA1_and_ClinicalData_ALS_vs_Control.xlsx")

# Convert Group to numeric 0/1 (CON → 0, ALS → 1) then factor with readable labels
data$Group_num <- ifelse(data$Group == "ALS", 1, 0)
data$Group <- factor(data$Group_num, levels = c(0, 1), labels = c("CON", "ALS"))
data$Group_num <- NULL

# Convert Sex to numeric 0/1 (Female → 0, Male → 1) then factor with readable labels
data$Sex_num <- ifelse(data$Sex == "male", 1, 0)
data$Sex <- factor(data$Sex_num, levels = c(0, 1), labels = c("Female", "Male"))
data$Sex_num <- NULL


############################################
# 2) Handling Upper Detection Limits (VB12, Folate)
############################################
# For VB12 and Folate columns with upper detection limits, use the upper bound directly
# Clean VB12
data$`VB12(pg/ml)` <- gsub(">", "", data$`VB12(pg/ml)`)
data$`VB12(pg/ml)` <- as.numeric(data$`VB12(pg/ml)`)

# Clean Folate
data$`Folate(ng/mL)` <- gsub(">", "", data$`Folate(ng/mL)`)
data$`Folate(ng/mL)` <- as.numeric(data$`Folate(ng/mL)`)


############################################
# 3) Binary Medication/Comorbidity Variables → Factors
############################################
# Batch convert 0/1 variables to factors; set reference as 'No'
drug_vars <- c("Riluzole Use", "LipidDrug Use", "Antihypertensive Use", "Antidiabetic Use",'Hypertension','Diabetes')
data[drug_vars] <- lapply(data[drug_vars], function(x) {
  factor(x, levels = c(0, 1), labels = c("No", "Yes"))
})


############################################
# 4) Transform Storage Duration (log)
############################################
# Create a log-transformed storage-time column (natural log)
data$`Storage Months` <- log(data$`Storage Months`)
saveRDS(data,file='data.rds')


############################################
# 5) Missingness Summary and Export
############################################
# Check missingness proportion per column
missing_pct <- sapply(data, function(x) mean(is.na(x)))
missing_pct

# Build a tidy missingness summary table
missing_summary <- data %>%
  summarise(across(
    everything(),
    list(
      Missing_Count = ~sum(is.na(.)),
      Missing_Pct   = ~mean(is.na(.))
    )
  )) %>%
  pivot_longer(
    everything(),
    names_to = c("Variable", ".value"),
    names_pattern = "^(.*)_(Missing_(?:Count|Pct))$"
  ) %>%
  mutate(Total = nrow(data)) %>%
  arrange(Missing_Pct)

write_csv(missing_summary, "missing_summary.csv")


############################################
# 6) Filter Variables by Missingness and Clean Column Names
############################################
# Drop columns with >40% missingness
data_filtered <- data[, missing_pct <= 0.4]

# Drop non-informative columns (first 4 columns)
data_filtered <- data_filtered[,-c(1:4)]

# Clean column names: replace spaces, parentheses, slashes
clean_names <- function(names_vector) {
  names_vector %>%
    gsub("\\(", "_", .) %>%
    gsub("\\)", "", .) %>%
    gsub("/", "_", .) %>%
    gsub(" ", "_", .)
}
colnames(data_filtered) <- clean_names(colnames(data_filtered))


############################################
# 7) Multiple Imputation by Chained Equations (MICE)
############################################
# Run MICE (predictive mean matching), m=5
imp <- mice(data_filtered, method = "pmm", m = 5, maxit = 50, seed = 123)
saveRDS(imp,file='data_miced.rds')
data_miced <- imp

# Overview of imputation
summary(imp)
stripplot(imp)

# Export original and imputed datasets to Excel
library(openxlsx)
m <- imp$m 
wb <- createWorkbook()

# Original data (with missingness)
addWorksheet(wb, "Original")
writeData(wb, "Original", imp$data)

# All imputed completed datasets
for (i in 1:m) {
  addWorksheet(wb, paste0("Imputed_", i))
  writeData(wb, paste0("Imputed_", i), complete(imp, i))
}

saveWorkbook(wb, "MICE_imputed_datasets.xlsx", overwrite = TRUE)


############################################
# 8) Univariate Spearman Screening Across Imputations
############################################
# Identify candidate variables (excluding ABCA1 outcome)
vars <- setdiff(colnames(complete(imp, 1)), "ABCA1_ng_ml")

# Extract completed datasets
imputed_list <- lapply(1:data_miced$m, function(i) complete(data_miced, i))

# Compute Spearman rho and p across imputations for each variable
cor_results <- list()
for (var in vars) {
  rhos <- c()
  pvals <- c()
  
  for (i in 1:length(imputed_list)) {
    data_i <- imputed_list[[i]]
    
    # If factor, convert to numeric coding (e.g., 0/1)
    if (is.factor(data_i[[var]])) {
      data_i[[var]] <- as.numeric(data_i[[var]])
    }
    
    # Spearman correlation
    cor_test_result <- cor.test(data_i[[var]], data_i$ABCA1_ng_ml, method = "spearman", exact = FALSE)
    rhos <- c(rhos, cor_test_result$estimate)
    pvals <- c(pvals, cor_test_result$p.value)
  }
  
  # Aggregate across imputations: mean rho, median p, and max p
  cor_results[[var]] <- data.frame(
    Variable = var,
    Mean_Rho = mean(rhos, na.rm = TRUE),
    Median_p = median(pvals, na.rm = TRUE),
    Max_p = max(pvals, na.rm = TRUE)
  )
}

# Combine results and save
final_cor_results <- bind_rows(cor_results)
readr::write_csv(final_cor_results,
                 "ABCA1_MICE_Spearman_cor_results.csv")

# Select significant variables (by Max_p here)
filtered_vars <- final_cor_results %>%filter(Max_p < 0.05)
filtered_vars  <- filtered_vars[2:4,] # exclude 'Group' row


############################################
# 9) Scatterplots + LOESS (per selected variable)
############################################
plot_one_var <- function(var_name, data, rho_val, p_val) {
  p <- ggplot(data, aes_string(x = var_name, y = "ABCA1_ng_ml")) +
    geom_point(color = '#EE6363', size = 1.5) +
    geom_smooth(method = "loess", se = TRUE, color = "#009ACD",
                fill = "#009ACD", alpha = 0.2, size = 0.75, span = 0.75) +
    geom_rug(sides = "bl", color = '#1F77B4', alpha = 0.85) +
    annotate("text",
             x = min(data[[var_name]], na.rm = TRUE) + 0.05,
             y = max(data$ABCA1_ng_ml, na.rm = TRUE) + 0.25,
             label = paste0(
               "italic(R)~'='~", round(rho_val, 3),
               "*\"\\n\"*",
               "italic(p)~'='~", format.pval(p_val, digits = 3, eps = 1e-3)
             ),
             hjust = 0, vjust = 1, size = 5,
             parse = TRUE)+
    labs(
      x = var_name,
      y = "Serum ABCA1 Protein (ng/ml)",
      title = paste0(var_name, " vs. ABCA1")
    ) +
    theme(
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.8),
      axis.text = element_text(color = "black", size = 12),
      axis.title = element_text(size = 14),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt")
    )
  
  # Add marginal histograms
  p_final <- ggMarginal(p,
                        type = "histogram",
                        xparams = list(fill = "#1F77B4", color = "white"),
                        yparams = list(fill = "#1F77B4", color = "white"),
                        size = 10)
  
  return(p_final)
}

# Batch plotting for selected variables
plots_list <- list()
library(ggplot2)
library(ggpubr)
library(ggExtra)
library(mice)
library(dplyr)
library(broom)

for (i in 1:nrow(filtered_vars)) {
  var_name <- filtered_vars$Variable[i]
  rho_val <- filtered_vars$Mean_Rho[i]
  p_val <- filtered_vars$Max_p[i]
  plots_list[[var_name]] <- plot_one_var(var_name,complete(data_miced, 1), rho_val, p_val)  # use the first completed dataset
}

# Example: show one variable plot
print(plots_list[[1]])  # e.g., 'Group' plot (index-based)

# Save example plots
ggsave(plots_list[[1]], file='Weight_Spearman.pdf',height=4,width=5)
ggsave(plots_list[[2]], file='BMI_Spearman.pdf',height=4,width=5)
ggsave(plots_list[[3]], file='LDL_Spearman.pdf',height=4,width=5)


############################################
# 10) Variable Selection and Multivariable Candidates
############################################
# Four variables (Group, BMI, LDL, and Weight) were initially selected; Weight was excluded due to collinearity with BMI.
# Stepwise AIC selection across imputations
library(MASS)
selected_vars <- list()

for (i in 1:5) {
  data_i <- complete(imp, i)
  full_model <- lm(ABCA1_ng_ml ~ Group + LDL_mmol_L + BMI, data = data_i)
  step_model <- stepAIC(full_model, direction = "both", trace = FALSE)
  vars_i <- names(coef(step_model))[-1]
  selected_vars[[i]] <- vars_i
}

# Count selection frequency across imputations
all_selected <- unlist(selected_vars)
var_freq <- sort(table(all_selected), decreasing = TRUE)
print(var_freq)
# Note: All three variables (Group, LDL_mmol_L, BMI) were retained in all five imputations.


############################################
# 11) Final 3-Variable Model (Group, BMI, LDL) + MI Pooling
############################################
# Build the final formula
final_vars <- c("Group",'BMI',"LDL_mmol_L")
final_formula <- as.formula(
  paste("ABCA1_ng_ml ~", paste(final_vars, collapse = " + "))
)

# Fit the same linear model in each of the 5 imputed datasets
models <- lapply(1:5, function(i) {
  data_i <- complete(imp, i)        # extract the i-th completed dataset from the mids object
  lm(final_formula, data = data_i)  # fit the model on that completed dataset
})

# Pool the models via Rubin's rules
pooled_model <- pool(models)

# Show the pooled inference table
summary(pooled_model)


############################################
# 12) Forest Plot (ggplot) for Pooled Results
############################################
library(broom)
library(ggplot2)

# Tidy pooled model with 95% CI
pooled_result <- tidy(pooled_model, conf.int = TRUE)

# Drop intercept
pooled_result <- pooled_result %>%
  filter(term != "(Intercept)")

# Forest plot via ggplot
ggplot(pooled_result, aes(x = estimate, y = reorder(term, estimate))) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_classic() +
  labs(x = "β Coefficient (95% CI)", y = NULL, title = "Forest plot of multivariable linear regression") +
  theme(
    axis.text = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  )


############################################
# 13) Forest Plot (ggforestplot) Alternative
############################################
library(broom)
library(ggforestplot)

# Tidy pooled model with 95% CI
pooled_result <- tidy(pooled_model, conf.int = TRUE)

# Drop intercept
pooled_result <- pooled_result %>%
  filter(term != "(Intercept)")
colnames(pooled_result) <- gsub("\\.", "", colnames(pooled_result))

# Forest plot via ggforestplot
forestplot(
  df = pooled_result,
  estimate = estimate,
  pvalue = pvalue,
  se = stderror,
  name = term,
  logodds = FALSE
) +
  ggtitle("Forest Plot of Multivariable Linear Regression") +
  theme_classic()


############################################
# 14) MI-Recommended Fitting (with.mids) + Model Fit Metrics
############################################
# Ensure reference level for Group
imp$data$Group <- relevel(factor(imp$data$Group), ref = "CON")

# Fit on the mids object directly (recommended)
fit_mi <- with(data = imp, expr = lm(ABCA1_ng_ml ~ Group + BMI + LDL_mmol_L))
pooled  <- pool(fit_mi)

# Summary with 95% CI
tab <- summary(pooled, conf.int = TRUE, conf.level = 0.95)

# Approximate pooled Adjusted R²: average across completed datasets
r2_each <- sapply(1:imp$m, function(i){
  d <- complete(imp, i)
  summary(lm(ABCA1_ng_ml ~ Group + BMI + LDL_mmol_L, data = d))$adj.r.squared
})
adj_r2_mean <- mean(r2_each, na.rm = TRUE)

# Overall F-statistic per completed dataset; take the mean (reporting purpose)
f_each <- sapply(1:imp$m, function(i) {
  d <- complete(imp, i)
  summary(lm(ABCA1_ng_ml ~ Group + BMI + LDL_mmol_L, data = d))$fstatistic[1]  # F value
})
f_mean <- mean(f_each, na.rm = TRUE)
f_each; f_mean


############################################
# 15) Tidy Table for Plotting (ggplot)
############################################
library(dplyr)
plot_df <- tab %>%
  filter(term != "(Intercept)") %>%
  transmute(term,
            estimate = estimate,
            conf.low = `2.5 %`,
            conf.high = `97.5 %`,
            p.value = p.value)

# Forest plot (ggplot) from pooled 'tab'
library(ggplot2)
ggplot(plot_df, aes(x = estimate, y = reorder(term, estimate))) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_classic() +
  labs(x = "β Coefficient (95% CI)", y = NULL,
       title = "Multivariable linear regression (MI pooled)") +
  theme(axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5))

# Multicollinearity check (VIF) on one completed dataset
library(car)
d1 <- complete(imp, 1)
vif(lm(ABCA1_ng_ml ~ Group + BMI + LDL_mmol_L, data = d1))


############################################
# 16) Group Normality / Variance Homogeneity / t-test / Visualization
############################################
# Shapiro-Wilk normality test by group
shapiro_test_con <- shapiro.test(as.numeric(data[data$Group == "CON", "ABCA1(ng/ml)"][[1]]))
shapiro_test_als <- shapiro.test(as.numeric(data[data$Group == "ALS", "ABCA1(ng/ml)"][[1]]))

shapiro_test_con # normality holds
shapiro_test_als # normality holds

# Levene's test for homogeneity of variance
library(car)
leveneTest(`ABCA1(ng/ml)` ~ Group, data = data) # p-value = 0.2823 indicates homoscedasticity

# Welch's t-test for group comparison
t_test_result <- t.test(`ABCA1(ng/ml)` ~ Group, data = data)
t_test_result
# Welch Two Sample t-test
# 
# data:  ABCA1(ng/ml) by Group
# t = -3.0095, df = 27.106, p-value = 0.005599
# alternative hypothesis: true difference in means between group CON and group ALS is not equal to 0
# 95 percent confidence interval:
#   -0.8161011 -0.1544945
# sample estimates:
#   mean in group CON mean in group ALS 
# 1.634338          2.119635 

# Bar + jitter plot for ABCA1 by Group
group_colors <- c("ALS" = "#9AC5CD", "CON" = "#FFDAB9")
p <- ggplot(data, aes(x = Group, y = `ABCA1(ng/ml)`, fill = Group)) +
  stat_summary(fun = mean, geom = "bar", color = "black", width = 0.6) +  
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, size = 0.8) +
  geom_jitter(aes(fill = Group), shape = 21, size = 2.5, width = 0.15, stroke = 0.3, color = "black", alpha = 0.9) +
  scale_fill_manual(values = group_colors) +
  labs(x = NULL, y = "Serum ABCA1 Protein (ng/ml)") +
  theme_classic() +
  theme(axis.title.y = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12, face = "bold"),
        legend.position = "none",
        axis.line = element_line(size = 0.8),
        axis.ticks = element_line(size = 0.8),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5))

# Add significance stars
p <- p + stat_compare_means(method = "t.test",
                            comparisons = list(c("ALS", "CON")),
                            label = "p.signif",
                            tip.length = 0.01)
p 
ggsave(p,file = 'In-house ABCA1 differential analysis.pdf',height = 5 , width = 5)


############################################
# 17) Storage Time vs ABCA1 Plot (by Group)
############################################
p1 <- ggplot(data, aes(x = `Storage Months`, y = `ABCA1(ng/ml)`, color = Group, fill = Group)) +
  geom_point(size = 1.5) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.1) +
  labs(
    x = expression(ln*"(Storage Duration in Months)"),
    y = "ABCA1 Protein (ng/ml)"
  ) +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  theme_classic(base_size = 8) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    axis.line = element_line(color = "black", size = 0.8)
  )
p1 
ggsave(p1,file = 'Storage time vs ABCA1.pdf',height = 3 , width = 5)


############################################
# 18) Simple Linear Model (Group + Storage Time) on Raw Data
############################################
# Fit a simple linear regression including Group and log(Storage Months)
model_simple <- lm(`ABCA1(ng/ml)` ~ Group + `Storage Months`, data = data)
summary(model_simple)




library(dplyr)
library(ggplot2)
library(ggforce)

data <- data %>%
  mutate(
    Age_group = cut(Age,
                    breaks = c(40, 54, 64, 100),
                    labels = c("<55", "55–64", "≥65"),
                    right = TRUE),
    BMI_group = cut(BMI,
                    breaks = c(18, 22, 25, 40),
                    labels = c("<22", "22–25", "≥25"),
                    right = TRUE)
  )

# 2) Determine subject order (same radius = same subject).
#    For example: order by Group, then Sex, then Age ascending.
plot_order <- data %>%
  mutate(.id = row_number()) %>%
  arrange(Group, Sex, Age) %>%
  mutate(
    idx = row_number(),
    n   = n(),
    theta_min = 2*pi*(idx-1)/n,
    theta_max = 2*pi*idx/n
  ) %>%
  dplyr::select(.id, Group, Sex, Age_group, BMI_group, theta_min, theta_max)

# 3) Build data frames for each ring (set inner/outer radii; same subject aligned radially)
mk_ring <- function(df, var, ring_name, r0, r1){
  df %>%
    transmute(
      ring   = ring_name,
      category = as.character({{ var }}),
      theta_min, theta_max,
      r0 = r0, r1 = r1
    )
}

ring_group <- mk_ring(plot_order, Group,     "Group", r0 = 0.18, r1 = 0.30)
ring_sex   <- mk_ring(plot_order, Sex,       "Sex",   r0 = 0.32, r1 = 0.44)
ring_age   <- mk_ring(plot_order, Age_group, "Age",   r0 = 0.46, r1 = 0.62)
ring_bmi   <- mk_ring(plot_order, BMI_group, "BMI",   r0 = 0.64, r1 = 0.82)

plot_df <- bind_rows(ring_group, ring_sex, ring_age, ring_bmi)

# 4) Define color palettes for each ring (names must match actual factor levels)
pal_group <- c(CON ='#4682B4', ALS = "#FA8072")
pal_sex   <- c(Female = "#F4A460", Male = "#1abc9c")
pal_age   <- c("<55" = "#f8c4a0", "55–64" = "#FFA07A", "≥65" = "#f1c40f")
pal_bmi   <- c("<22" = "#FFC0CB", "22–25" = "#96CDCD", "≥25" = "#F5F5DC")

# Map each ring to a unified 'fill' column
plot_df <- plot_df %>%
  mutate(fill = case_when(
    ring == "Group" ~ pal_group[category],
    ring == "Sex"   ~ pal_sex[category],
    ring == "Age"   ~ pal_age[category],
    ring == "BMI"   ~ pal_bmi[category]
  ))

# 5) Add combined label (ring:category) for legend display
plot_df <- plot_df %>%
  mutate(label = paste(ring, category, sep=":"))

# 6) Plot multi-ring radial chart (same subject aligned across rings)
p <- ggplot(plot_df) +
  geom_arc_bar(
    aes(x0 = 0, y0 = 0, r0 = r0, r = r1,
        start = theta_min, end = theta_max,
        fill = label),   # map to label for legend
    color = "black", size = 0.75
  ) +
  scale_fill_manual(values = setNames(plot_df$fill, plot_df$label)) +
  coord_fixed() +
  theme_void() +
  ggtitle("Distribution of Group, Sex, Age, and BMI across subjects") +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold"))

p



# ===== Sex-stratified comparison of ABCA1 (ALS vs CON) with visualization =====

library(dplyr)
library(ggplot2)
library(ggpubr)
library(car)

# Function to run stratified analysis by Sex
sex_stratified_tests <- function(sex_level, data){
  cat("\n============================\n")
  cat("Sex =", sex_level, "\n")
  
  sub_data <- filter(data, Sex == sex_level)
  
  # Normality test (Shapiro-Wilk) within each group
  shapiro_con <- shapiro.test(sub_data$`ABCA1(ng/ml)`[sub_data$Group == "CON"])
  shapiro_als <- shapiro.test(sub_data$`ABCA1(ng/ml)`[sub_data$Group == "ALS"])
  
  cat("Shapiro-Wilk Test (CON): W =", round(shapiro_con$statistic,3), 
      "p =", round(shapiro_con$p.value,4), "\n")
  cat("Shapiro-Wilk Test (ALS): W =", round(shapiro_als$statistic,3), 
      "p =", round(shapiro_als$p.value,4), "\n")
  
  # Homogeneity of variance test (Levene’s)
  levene_res <- leveneTest(`ABCA1(ng/ml)` ~ Group, data = sub_data)
  print(levene_res)
  
  # t-test (Welch’s by default)
  t_res <- t.test(`ABCA1(ng/ml)` ~ Group, data = sub_data)
  print(t_res)
  
  return(list(shapiro_con = shapiro_con,
              shapiro_als = shapiro_als,
              levene = levene_res,
              t_test = t_res))
}

# Run for both Female and Male
results_female <- sex_stratified_tests("Female", data)
results_male   <- sex_stratified_tests("Male", data)

# Ensure no missing values in key columns and set factor orders
data <- data %>%
  filter(!is.na(Sex), !is.na(Group), !is.na(`ABCA1(ng/ml)`))
data$Group <- factor(data$Group, levels = c("CON","ALS"))
data$Sex   <- factor(data$Sex,   levels = c("Female","Male"))

# Run t-tests within each sex (prints to console)
tt_female <- t.test(`ABCA1(ng/ml)` ~ Group, data = subset(data, Sex == "Female"))
tt_male   <- t.test(`ABCA1(ng/ml)` ~ Group, data = subset(data, Sex == "Male"))
print(tt_female)
print(tt_male)

# Colors (same style as your reference)
group_colors <- c("ALS" = "#9AC5CD", "CON" = "#FFDAB9")

# Helper to build the same plot for a given sex
plot_by_sex <- function(sex_label){
  d <- subset(data, Sex == sex_label)
  y_max <- max(d$`ABCA1(ng/ml)`, na.rm = TRUE)
  ggplot(d, aes(x = Group, y = `ABCA1(ng/ml)`, fill = Group)) +
    stat_summary(fun = mean, geom = "bar", color = "black", width = 0.6) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, size = 0.8) +
    geom_jitter(shape = 21, size = 2.5, width = 0.15, stroke = 0.3,
                color = "black", alpha = 0.9) +
    scale_fill_manual(values = group_colors) +
    labs(x = NULL, y = "Serum ABCA1 Protein (ng/ml)",
         title = paste0("Sex: ", sex_label)) +
    theme_classic() +
    theme(axis.title.y = element_text(size = 14, face = "bold"),
          axis.text = element_text(size = 12, face = "bold"),
          legend.position = "none",
          axis.line = element_line(size = 0.8),
          axis.ticks = element_line(size = 0.8),
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) +
    stat_compare_means(method = "t.test",
                       comparisons = list(c("ALS","CON")),
                       label = "p.signif",
                       tip.length = 0.01,
                       label.y = y_max * 1.05)
}

# Build plots for Female and Male
p_female <- plot_by_sex("Female")
p_male   <- plot_by_sex("Male")

# Show individually
p_female
p_male

# Optional: combine into a single 2-panel figure
p_combined <- ggarrange(p_female, p_male, ncol = 2, labels = c("A","B"))
p_combined

# Save figures
ggsave("ABCA1_by_group_Female.pdf", p_female, width = 5, height = 5)
ggsave("ABCA1_by_group_Male.pdf",   p_male,   width = 5, height = 5)
ggsave("ABCA1_by_group_by_sex_combined.pdf", p_combined, width = 10, height = 5)



# ===== Age (<=65 vs >65) stratified comparison of ABCA1 (ALS vs CON) =====

library(dplyr)
library(ggplot2)
library(ggpubr)
library(car)

# Build age stratum (<=65 vs >65)
data <- data %>%
  filter(!is.na(Age), !is.na(Group), !is.na(`ABCA1(ng/ml)`)) %>%
  mutate(
    Group = factor(Group, levels = c("CON","ALS")),
    Age65 = factor(ifelse(Age > 65, ">65", "<=65"), levels = c("<=65", ">65"))
  )

# Function to run stratified analysis by Age65
age_stratified_tests <- function(age_level, data){
  cat("\n============================\n")
  cat("Age =", age_level, "\n")
  
  sub_data <- dplyr::filter(data, Age65 == age_level)
  
  # Shapiro–Wilk (normality) within each group
  shapiro_con <- shapiro.test(sub_data$`ABCA1(ng/ml)`[sub_data$Group == "CON"])
  shapiro_als <- shapiro.test(sub_data$`ABCA1(ng/ml)`[sub_data$Group == "ALS"])
  
  cat("Shapiro-Wilk (CON): W =", round(shapiro_con$statistic,3), 
      "p =", round(shapiro_con$p.value,4), "\n")
  cat("Shapiro-Wilk (ALS): W =", round(shapiro_als$statistic,3), 
      "p =", round(shapiro_als$p.value,4), "\n")
  
  # Levene’s test (homogeneity of variance)
  levene_res <- leveneTest(`ABCA1(ng/ml)` ~ Group, data = sub_data)
  print(levene_res)
  
  # t-test (Welch by default; same as your sex-stratified code)
  t_res <- t.test(`ABCA1(ng/ml)` ~ Group, data = sub_data)
  print(t_res)
  
  return(list(shapiro_con = shapiro_con,
              shapiro_als = shapiro_als,
              levene = levene_res,
              t_test = t_res))
}

# Run for both strata (<=65 and >65)
results_le65 <- age_stratified_tests("<=65", data)
results_gt65 <- age_stratified_tests(">65",  data)

# --- Visualization (same style as your reference) ---
group_colors <- c("ALS" = "#9AC5CD", "CON" = "#FFDAB9")

plot_by_age <- function(age_label){
  d <- dplyr::filter(data, Age65 == age_label)
  y_max <- max(d$`ABCA1(ng/ml)`, na.rm = TRUE)
  ggplot(d, aes(x = Group, y = `ABCA1(ng/ml)`, fill = Group)) +
    stat_summary(fun = mean, geom = "bar", color = "black", width = 0.6) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, size = 0.8) +
    geom_jitter(shape = 21, size = 2.5, width = 0.15, stroke = 0.3,
                color = "black", alpha = 0.9) +
    scale_fill_manual(values = group_colors) +
    labs(x = NULL, y = "Serum ABCA1 Protein (ng/ml)",
         title = paste0("Age ", age_label)) +
    theme_classic() +
    theme(axis.title.y = element_text(size = 14, face = "bold"),
          axis.text = element_text(size = 12, face = "bold"),
          legend.position = "none",
          axis.line = element_line(size = 0.8),
          axis.ticks = element_line(size = 0.8),
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) +
    stat_compare_means(method = "t.test",
                       comparisons = list(c("ALS","CON")),
                       label = "p.signif",
                       tip.length = 0.01,
                       label.y = y_max * 1.05)
}

p_le65 <- plot_by_age("<=65")
p_gt65 <- plot_by_age(">65")

# Show individually
p_le65
p_gt65

# Combine into one figure
p_age_combined <- ggpubr::ggarrange(p_le65, p_gt65, ncol = 2, labels = c("A","B"))
p_age_combined

# Save (optional)
ggsave("ABCA1_by_group_AgeLE65.pdf", p_le65, width = 5, height = 5)
ggsave("ABCA1_by_group_AgeGT65.pdf", p_gt65, width = 5, height = 5)
ggsave("ABCA1_by_group_byAge65_combined.pdf", p_age_combined, width = 10, height = 5)





# ===== Linear model: ABCA1(ng/ml) ~ Group + Age + Sex (raw data, no imputation) =====
library(dplyr)
library(broom)
library(car)
library(lmtest)
library(emmeans)
library(ggplot2)

# 0) Prepare analysis dataset
df <- data %>%
  filter(!is.na(`ABCA1(ng/ml)`), !is.na(Group), !is.na(Age), !is.na(Sex)) %>%
  mutate(
    Group = factor(Group, levels = c("CON","ALS")),   # set reference
    Sex   = factor(Sex,   levels = c("Female","Male"))
  )

# 1) Fit model
fit <- lm(`ABCA1(ng/ml)` ~ Group + Age + Sex, data = df)

# 2) Coefficients with 95% CI (primary reporting table)
coef_tab <- broom::tidy(fit, conf.int = TRUE, conf.level = 0.95)
print(coef_tab)

# 3) Type II ANOVA (adjusted effect for each term; see 'Group' row)
anova_type2 <- car::Anova(fit, type = 2)
print(anova_type2)

# 4) Model summary metrics
summ   <- summary(fit)
adj_r2 <- summ$adj.r.squared
sigma  <- summ$sigma
cat(sprintf("\nAdjusted R^2 = %.3f | Residual SE = %.3f (df = %d)\n",
            adj_r2, sigma, summ$df[2]))

# 5) Diagnostics (residual normality, heteroskedasticity, collinearity)
res_sw <- shapiro.test(residuals(fit))  # residual normality
bp     <- lmtest::bptest(fit)           # Breusch–Pagan for heteroskedasticity
vifs   <- car::vif(fit)                 # VIFs for multicollinearity
print(res_sw); print(bp); print(vifs)

# 6) Adjusted means by Group (controlling Age & Sex) and pairwise contrast
emm_g <- emmeans::emmeans(fit, ~ Group)
print(emm_g)
print(summary(contrast(emm_g, "pairwise"), infer = TRUE))  # adjusted diff ALS - CON

# 7) Visualization: adjusted means (EMMs) with 95% CI
emm_df <- as.data.frame(emm_g)
group_colors <- c("ALS" = "#9AC5CD", "CON" = "#FFDAB9")
p_emm <- ggplot(emm_df, aes(x = Group, y = emmean, fill = Group)) +
  geom_col(width = 0.6, color = "black") +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.18, size = 0.8) +
  scale_fill_manual(values = group_colors) +
  labs(x = NULL, y = "Adjusted mean ABCA1 (ng/ml)",
       title = "Group effect on ABCA1 adjusted for Age and Sex") +
  theme_classic() +
  theme(axis.title.y = element_text(size = 13, face = "bold"),
        axis.text    = element_text(size = 11, face = "bold"),
        legend.position = "none",
        plot.title   = element_text(size = 14, face = "bold", hjust = 0.5))
p_emm
ggsave("ABCA1_Group_AgeSex_adjusted_means.pdf", p_emm, width = 5, height = 5)



# ===== Coefficient forest plot: β (95% CI) for Group, Age, Sex =====
library(broom)
library(dplyr)
library(ggplot2)

# 1) Tidy coefficients (drop intercept; relabel for display)
coef_df <- broom::tidy(fit, conf.int = TRUE) %>%
  dplyr::filter(term != "(Intercept)") %>%
  dplyr::mutate(
    term = dplyr::recode(term,
                         "GroupALS" = "Group (ALS vs CON)",
                         "Age"      = "Age (years)",
                         "SexMale"  = "Sex (Male vs Female)")
  ) %>%
  dplyr::arrange(estimate) %>%                       # sort by β
  dplyr::mutate(term = factor(term, levels = term))  # keep order for plotting

# 3) Plot
ggplot(coef_df, aes(x = estimate, y = term)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.4) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.3, color = "#FFDAB9") +
  geom_point(shape = 22, fill = "#9AC5CD", size = 3.5, stroke = 1.2,color='#5F9EA0') +
  labs(
    x = expression(beta~"(95% CI)"),
    y = NULL,
    title = "Multivariable Linear Regression for Serum ABCA1"
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
