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
imp <- readRDS('data_miced.rds')
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


# Note: In this cohort, LDL_mmol_L and BMI have complete data (0% missing).
d1 <- complete(data_miced, 1)
d1_male   <- droplevels(subset(d1,   Sex == "Male"))
d1_female <- droplevels(subset(d1,   Sex == "Female"))
# ---- safety checks (optional) ----
stopifnot(all(c("ABCA1_ng_ml","LDL_mmol_L","BMI") %in% names(d1_male)))
stopifnot(all(c("ABCA1_ng_ml","LDL_mmol_L","BMI") %in% names(d1_female)))

# ---- helper: compute Spearman (rho & p) on a given subset ----
get_spearman <- function(df, xvar, yvar = "ABCA1_ng_ml") {
  ok <- stats::complete.cases(df[, c(xvar, yvar)])
  n  <- sum(ok)
  if (n < 3) return(list(rho = NA_real_, p = NA_real_, n = n))
  ct <- suppressWarnings(cor.test(df[[xvar]][ok], df[[yvar]][ok],
                                  method = "spearman", exact = FALSE))
  list(rho = unname(ct$estimate), p = ct$p.value, n = n)
}

# ---- loop over variables and sex subsets; build plots & save PDFs ----
plots_list    <- list()
filtered_vars <- c("LDL_mmol_L","BMI")

for (sex in c("male","female")) {
  df_sex <- if (sex == "male") d1_male else d1_female
  
  for (i in seq_along(filtered_vars)) {
    # 1) variable name
    var_name <- filtered_vars[i]
    
    # 2) compute Spearman on this subset
    sp <- get_spearman(df_sex, var_name)  # list(rho=..., p=..., n=...)
    
    # 3) build plot with the computed rho/p
    p <- plot_one_var(var_name, df_sex, sp$rho, sp$p)
    
    # 4) store & save
    key <- paste0(var_name, "_", sex)
    plots_list[[key]] <- p
    
    ggsave(plot = p,
           filename = paste0(var_name, "_", sex, "_Spearman.pdf"),
           height = 4, width = 5, device = "pdf")
  }
}

# Example: show one plot in the viewer
print(plots_list[["LDL_mmol_L_male"]])


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


















############################################################
# Baseline tables (ALS vs CON) — automatic test selection
# - Keep variables only if non-missing n ≥ 24 (missing <20%)
# - Auto-choose tests per variable:
#     * Categorical: Fisher's exact if any cell/expected < 5, else Chi-square
#     * Continuous: Welch t-test if both groups ~normal (Shapiro p≥0.05),
#                   else Wilcoxon rank-sum
# - Exports: Excel (.xlsx) and PDF (.pdf)
# NOTE: PDF export requires a LaTeX installation (e.g., TinyTeX).
############################################################

# ---- 0) Packages ----
# install.packages(c("readxl","dplyr","compareGroups"))
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(compareGroups)
})

# ---- 1) Load data ----
file_path <- "Serum_ABCA1_and_ClinicalData_ALS_vs_Control.xlsx"
dat <- read_excel(file_path)

# ---- 2) Preprocessing (factor coding & 0/1 to No/Yes) ----
to_no_yes_factor <- function(x) {
  if (is.factor(x)) return(x)
  if (is.numeric(x) || is.integer(x)) {
    ux <- unique(na.omit(x))
    if (all(ux %in% c(0, 1))) return(factor(x, levels = c(0,1), labels = c("No","Yes")))
  }
  if (is.character(x)) {
    xl <- tolower(trimws(x))
    if (all(na.omit(xl) %in% c("0","1","no","yes","none","absent","present"))) {
      y <- ifelse(xl %in% c("1","yes","present"), "Yes",
                  ifelse(xl %in% c("0","no","none","absent"), "No", NA))
      return(factor(y, levels = c("No","Yes")))
    }
  }
  x
}

# Group as factor (must be "CON" and "ALS")
if (!"Group" %in% names(dat)) stop("Column `Group` not found (needs 'CON'/'ALS').")
dat <- dat |> mutate(Group = factor(Group, levels = c("CON","ALS")))

# Sex to factor if present (robust recoding)
if ("Sex" %in% names(dat)) {
  dat <- dat |>
    mutate(
      Sex = case_when(
        tolower(as.character(Sex)) %in% c("female","f","0") ~ "Female",
        tolower(as.character(Sex)) %in% c("male","m","1")   ~ "Male",
        TRUE ~ as.character(Sex)
      ),
      Sex = factor(Sex, levels = c("Female","Male"))
    )
}

# Convert 0/1 comorbidity/medication variables to No/Yes factors (if present)
bin_cols <- c("Hypertension","Diabetes","Riluzole Use","LipidDrug Use",
              "Antihypertensive Use","Antidiabetic Use")
for (v in intersect(bin_cols, names(dat))) dat[[v]] <- to_no_yes_factor(dat[[v]])

# ---- 3) Final variable lists (as agreed) ----
# Main (≤15): TC in main; Glucose moved to supplement
main_vars_final <- c(
  "Age","Sex","BMI","Height(cm)","Weight(kg)",
  "Hypertension","Diabetes","Antihypertensive Use","Antidiabetic Use",
  "HbA1c","ALT(U/L)",
  "LDL(mmol/L)","HDL(mmol/L)","TC(mmol/L)","TG(mmol/L)",
  "WBC(e9/L)"
)

supp_vars_final <- c(
  "Riluzole Use","LipidDrug Use",
  "Glucose(mmol/L)","Creatinine(umol/L)","Hcy(umol/L)",
  "Na(mmol/L)","K(mmol/L)","Cl(mmol/L)","Ca(mmol/L)","P(mmol/L)",
  "AST(U/L)","ALP(U/L)","DBIL(umol/L)","IBIL(umol/L)","ALB(g/L)","GLB(g/L)",
  "HB(g/L)","PLT(e9/L)","RBC(e12/L)",
  "NEUT(e9/L)","LYM(e9/L)","MONO(e9/L)","EOS(e9/L)","BASO(e9/L)"
)

# ---- 4) Missingness rule: keep only variables with n ≥ 24 ----
total_n <- nrow(dat)
nn <- sapply(dat, function(x) sum(!is.na(x)))
keep_main <- intersect(main_vars_final, names(dat))
keep_supp <- intersect(supp_vars_final, names(dat))
keep_main <- keep_main[nn[keep_main] >= 24]
keep_supp <- keep_supp[nn[keep_supp] >= 24]

dropped_main <- setdiff(main_vars_final, keep_main)
dropped_supp <- setdiff(supp_vars_final, keep_supp)
if (length(dropped_main)) message("Dropped from MAIN (n<24): ", paste(dropped_main, collapse=", "))
if (length(dropped_supp)) message("Dropped from SUPP (n<24): ", paste(dropped_supp, collapse=", "))

# ---- 5) Auto-select methods per variable (compareGroups::method vector) ----
# method = 1 -> nonparametric (Wilcoxon for continuous; Fisher for categorical)
# method = 2 -> parametric (Welch t for continuous; Chi-square for categorical)
is_categorical <- function(x) is.factor(x) || is.character(x)

choose_method_for_var <- function(varname, data) {
  x <- data[[varname]]
  g <- data$Group
  ok <- !is.na(x) & !is.na(g)
  x <- x[ok]; g <- g[ok]
  if (length(x) == 0) return(NA_integer_)  # let compareGroups decide (shouldn't happen after filtering)
  
  # CATEGORICAL: Fisher if any cell or expected < 5; else Chi-square
  if (is_categorical(x)) {
    f <- if (!is.factor(x)) factor(x) else x
    # If too many levels (sparse multi-level), prefer Fisher
    if (nlevels(f) > 6) return(1L)
    tab <- table(g, f)
    if (any(tab < 5)) return(1L)
    # try expected counts; if error, fallback to Fisher
    exp_ok <- tryCatch({
      suppressWarnings({
        ce <- suppressWarnings(chisq.test(tab, correct = FALSE))
        all(ce$expected >= 5)
      })
    }, error = function(e) FALSE)
    if (!exp_ok) return(1L) else return(2L)
  }
  
  # CONTINUOUS: Welch t if both groups ~normal; else Wilcoxon
  # Need enough per-group size for Shapiro (n>=3)
  x0 <- x[g == "CON"]; x1 <- x[g == "ALS"]
  normal0 <- if (length(x0) >= 3) shapiro.test(x0)$p.value >= 0.05 else FALSE
  normal1 <- if (length(x1) >= 3) shapiro.test(x1)$p.value >= 0.05 else FALSE
  if (normal0 && normal1) return(2L) else return(1L)
}

# Build method vectors
method_vec_main <- rep(NA_integer_, length(keep_main))
names(method_vec_main) <- keep_main
for (v in keep_main) method_vec_main[v] <- choose_method_for_var(v, dat)

method_vec_supp <- NULL
if (length(keep_supp)) {
  method_vec_supp <- rep(NA_integer_, length(keep_supp))
  names(method_vec_supp) <- keep_supp
  for (v in keep_supp) method_vec_supp[v] <- choose_method_for_var(v, dat)
}

# Optional: print a compact summary of selected methods
print_selected_methods <- function(method_vec, title) {
  if (is.null(method_vec)) return()
  cat("\n----", title, "----\n", sep = "")
  for (nm in names(method_vec)) {
    m <- method_vec[[nm]]
    cat(sprintf("%-25s : %s\n", nm,
                ifelse(is.na(m), "auto (NA)",
                       ifelse(m == 1, "nonparametric (Wilcoxon/Fisher)", "parametric (Welch t/Chi-square)"))))
  }
}
print_selected_methods(method_vec_main, "MAIN methods")
print_selected_methods(method_vec_supp, "SUPP methods")

# ---- 6) Fit compareGroups and create tables ----
res_main <- compareGroups(
  Group ~ .,
  data   = dat[, c("Group", keep_main), drop = FALSE],
  method = method_vec_main
)
tab_main <- createTable(res_main, show.p.overall = TRUE, show.n = TRUE, digits = 2)

tab_supp <- NULL
if (length(keep_supp)) {
  res_supp <- compareGroups(
    Group ~ .,
    data   = dat[, c("Group", keep_supp), drop = FALSE],
    method = method_vec_supp
  )
  tab_supp <- createTable(res_supp, show.p.overall = TRUE, show.n = TRUE, digits = 2)
}

# ---- 7) Export: Excel & PDF ----
export2xls(tab_main, file = "Table1_Baseline_main.xlsx")
if (!is.null(tab_supp)) export2xls(tab_supp, file = "TableS_Baseline_supp.xlsx")

safe_export2pdf <- function(tbl, path) {
  tryCatch(
    { export2pdf(tbl, file = path); message("PDF exported: ", path) },
    error = function(e) {
      message("PDF export failed: ", e$message,
              "\nTip: use export2html() + webshot2/pagedown to convert to PDF.")
    }
  )
}
safe_export2pdf(tab_main, "Table1_Baseline_main.pdf")
if (!is.null(tab_supp)) safe_export2pdf(tab_supp, "TableS_Baseline_supp.pdf")

message("\n✅ Done:")
message(" - Main:  Table1_Baseline_main.xlsx / Table1_Baseline_supp.pdf")
if (is.null(tab_supp)) {
  message(" - Supplementary: not created (no variables with n ≥ 24).")
} else {
  message(" - Supp:  TableS_Baseline_supp.xlsx / TableS_Baseline_supp.pdf")
}

############################################################
# METHODS (for manuscript; English)
# We summarized baseline characteristics of ALS vs control cohorts.
# Variables were included only if missingness <20% (i.e., available n ≥ 24 of 30).
# For continuous variables, we assessed normality per group using the Shapiro–Wilk test.
# If both groups were approximately normal (p ≥ 0.05), we used a two-sample Welch t-test
# and reported mean (SD). Otherwise, we used the Wilcoxon rank-sum test and reported
# median [IQR]. For categorical variables, we first examined contingency tables; if any
# observed cell count or expected cell count was < 5, we used Fisher’s exact test; otherwise,
# we used the Chi-square test. Two-sided p-values were reported as overall group comparisons.
# Table 1 (main) contained ≤15 key variables (Age, Sex, BMI, comorbidities/medications,
# and core labs including HbA1c, eGFR, ALT, LDL, HDL, TC, TG, WBC). Additional variables
# meeting the missingness threshold were presented in the supplementary table.
############################################################
















############################################################
# ALSFRS-R (total & 4 domains) vs ABCA1 — Full Pipeline
# - Read & preprocessing
# - Spearman correlations (+ FDR) with scatter/LOESS plots
# - Adjusted linear regression with standardized betas
# - Median split (High vs Low) group tests:
#     If normality & equal variances → Student's t-test + BAR (mean±SE)
#     Else → Wilcoxon rank-sum + BOXPLOT
# - Save ALL figures as PDF ONLY
############################################################

# ---- 0) Packages ----
# install.packages(c("readxl","dplyr","ggplot2","ggpubr","ggExtra","broom","rstatix","pROC","car","splines"))
suppressPackageStartupMessages({
  library(readxl);  library(dplyr);   library(ggplot2); library(ggpubr)
  library(ggExtra); library(broom);   library(rstatix); library(pROC)
  library(car);     library(splines)
})

# ---- 1) Load data ----
dat <- read_excel("ALS_ABCA1_ALSFRSR_Dataset.xlsx")

# ---- 2) Harmonize column names ----
# Handle trailing spaces in names (e.g., "Riluzole Use ")
names(dat) <- trimws(names(dat))

# Map common names to snake_case for coding convenience
rename_map <- c(
  "ABCA1(ng/ml)" = "ABCA1_ng_ml",
  "ALSFRS-R"     = "ALSFRS_total",
  "FINE MOTOR"   = "FINE_MOTOR",
  "GROSS MOTOR"  = "GROSS_MOTOR",
  "eGFR(ml/min)" = "eGFR_mL_min",
  "LDL(mmol/L)"  = "LDL_mmol_L",
  "HDL(mmol/L)"  = "HDL_mmol_L",
  "TC(mmol/L)"   = "TC_mmol_L",
  "TG(mmol/L)"   = "TG_mmol_L",
  "Riluzole Use" = "Riluzole_Use",
  "LipidDrug Use"= "LipidDrug_Use",
  "Antihypertensive Use" = "Antihypertensive_Use",
  "Antidiabetic Use"     = "Antidiabetic_Use",
  "Height(cm)"   = "Height_cm",
  "Weight(kg)"   = "Weight_kg"
)
for (old in intersect(names(rename_map), names(dat))) {
  names(dat)[names(dat) == old] <- rename_map[[old]]
}

# ---- 3) Basic preprocessing & types ----
# Factorize Sex; medication uses to No/Yes if present
dat <- dat %>%
  mutate(
    Sex = if ("Sex" %in% names(.))
      factor(trimws(as.character(Sex)), levels = c("Female","Male")) else NULL,
    Riluzole_Use  = if ("Riluzole_Use"  %in% names(.))
      factor(Riluzole_Use,  levels = c(0,1), labels = c("No","Yes")) else NULL,
    LipidDrug_Use = if ("LipidDrug_Use" %in% names(.))
      factor(LipidDrug_Use, levels = c(0,1), labels = c("No","Yes")) else NULL
  )

# Ensure numeric for ALSFRS-R scores and ABCA1
score_vars    <- c("ALSFRS_total","BULBAR","FINE_MOTOR","GROSS_MOTOR","RESPIRATORY")
scores_to_use <- intersect(score_vars, names(dat))
if (!"ABCA1_ng_ml" %in% names(dat)) stop("ABCA1 column not found.")
dat[scores_to_use] <- lapply(dat[scores_to_use], function(x) suppressWarnings(as.numeric(x)))
dat$ABCA1_ng_ml    <- suppressWarnings(as.numeric(dat$ABCA1_ng_ml))

message("N = ", nrow(dat))
message("Scores available: ", paste(scores_to_use, collapse = ", "))

# ---- helper: drop covariates with only one level/value ----
drop_one_level_covars <- function(df, covars){
  keep <- c(); dropped <- c()
  for (v in covars) {
    if (!v %in% names(df)) next
    x <- df[[v]]; x <- x[!is.na(x)]
    ok <- if (is.factor(x)) nlevels(droplevels(x)) >= 2 else length(unique(x)) >= 2
    if (ok) keep <- c(keep, v) else dropped <- c(dropped, v)
  }
  list(keep = keep, dropped = dropped)
}

# ---------------------------------------------------------
# A) Spearman correlations (ABCA1 vs 5 ALSFRS scores)
# ---------------------------------------------------------
spearman_one <- function(df, score){
  ok <- complete.cases(df[, c("ABCA1_ng_ml", score)])
  if (sum(ok) < 3) return(tibble(term = score, rho = NA_real_, p = NA_real_, n = sum(ok)))
  ct <- suppressWarnings(cor.test(df$ABCA1_ng_ml[ok], df[[score]][ok], method = "spearman", exact = FALSE))
  tibble(term = score, rho = unname(ct$estimate), p = ct$p.value, n = sum(ok))
}
cor_tab <- bind_rows(lapply(scores_to_use, function(s) spearman_one(dat, s))) %>%
  mutate(p_fdr = p.adjust(p, method = "fdr"))
print(cor_tab)

# ---- E1) Spearman-style visualization (scatter + LOESS + margins + annotation) ----
spearman_annot_plot <- function(df, xvar, yvar = "ABCA1_ng_ml") {
  ok <- complete.cases(df[, c(xvar, yvar)])
  rho <- pval <- NA_real_
  if (sum(ok) >= 3) {
    ct <- suppressWarnings(cor.test(df[[xvar]][ok], df[[yvar]][ok], method = "spearman", exact = FALSE))
    rho  <- unname(ct$estimate); pval <- ct$p.value
  }
  p <- ggplot(df, aes(.data[[xvar]], .data[[yvar]])) +
    geom_point(color = '#EE6363', size = 1.8, alpha = 0.9) +
    geom_smooth(method = "loess", se = TRUE,
                color = "#009ACD", fill = "#009ACD", alpha = 0.2,
                size = 0.75, span = 0.75) +
    geom_rug(sides = "bl", color = '#1F77B4', alpha = 0.85) +
    annotate("text",
             x = min(df[[xvar]], na.rm = TRUE),
             y = max(df[[yvar]], na.rm = TRUE),
             hjust = 0, vjust = 1, size = 5, parse = TRUE,
             label = paste0(
               "italic(R)[s]*' = ", ifelse(is.na(rho), "NA", sprintf('%.3f', rho)), "'",
               "*\"\\n\"*",
               "italic(p)*' = ", ifelse(is.na(pval), "NA", format.pval(pval, digits = 3, eps = 1e-3)), "'"
             )) +
    labs(x = xvar, y = "Serum ABCA1 Protein (ng/mL)",
         title = paste0(xvar, " vs. ABCA1")) +
    theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.8),
      axis.text = element_text(color = "black", size = 12),
      axis.title = element_text(size = 14),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.margin = margin(5, 5, 5, 5)
    )
  ggMarginal(p, type = "histogram",
             xparams = list(fill = "#1F77B4", color = "white"),
             yparams = list(fill = "#1F77B4", color = "white"),
             size = 10)
}

# Save one PDF per score (PDF only)
for (s in scores_to_use) {
  p_s <- spearman_annot_plot(dat, s, "ABCA1_ng_ml")
  ggsave(filename = paste0("ABCA1_vs_", s, "_Spearman.pdf"),
         plot = p_s, width = 5, height = 4, device = "pdf")
}

# ---------------------------------------------------------
# B) Adjusted linear regression (std beta for score)
# ---------------------------------------------------------
conf_cands <- c("Age","Sex","BMI","eGFR_mL_min","LDL_mmol_L","LipidDrug_Use","Riluzole_Use")
conf_ok    <- drop_one_level_covars(dat, intersect(conf_cands, names(dat)))
if (length(conf_ok$dropped)) message("Dropped covariates (one level/value): ", paste(conf_ok$dropped, collapse = ", "))
covars     <- conf_ok$keep

fit_list <- lapply(scores_to_use, function(s){
  f <- as.formula(paste("ABCA1_ng_ml ~", paste(c(s, covars), collapse = " + ")))
  m <- lm(f, data = dat)
  
  # Standardized beta for the score (z-score Y and X)
  X  <- model.matrix(m)
  sy <- as.numeric(scale(dat$ABCA1_ng_ml))
  sx <- scale(X[,-1, drop = FALSE])
  m_std <- lm(sy ~ sx)
  
  score_pos <- which(colnames(sx) == s)
  beta <- conf.low <- conf.high <- pval <- NA_real_
  if (length(score_pos) == 1) {
    co <- summary(m_std)$coefficients
    ci <- confint(m_std)
    beta     <- co[score_pos + 1, "Estimate"]
    pval     <- co[score_pos + 1, "Pr(>|t|)"]
    conf.low <- ci[score_pos + 1, 1]
    conf.high<- ci[score_pos + 1, 2]
  } else {
    if (s %in% rownames(summary(m)$coefficients)) pval <- summary(m)$coefficients[s,4]
  }
  
  tibble(
    term     = s,
    beta_std = beta,
    conf.low = conf.low,
    conf.high= conf.high,
    p_raw    = pval,
    p_fdr    = p.adjust(pval, "fdr"),
    adjR2    = summary(m)$adj.r.squared,
    n_model  = length(residuals(m))
  )
}) %>% bind_rows()
print(fit_list)

# ---------------------------------------------------------
# C) Median split (High vs Low) + auto test selection
#     If (both groups normal) & (Levene p>=.05) → Student t-test (var.equal=TRUE) + BAR
#     Else → Wilcoxon + BOXPLOT
# ---------------------------------------------------------
mk_group <- function(x) ifelse(x >= median(x, na.rm = TRUE), "High","Low")
for (s in scores_to_use) dat[[paste0(s, "_grp")]] <- mk_group(dat[[s]])

# Helper to generate significance label from p
p_to_label <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 0.001) "p < 0.001" else sprintf("p = %.3f", p)
}

# Group test + effect size for one grouped variable
group_test_one <- function(gvar){
  if (!gvar %in% names(dat)) return(NULL)
  ok <- complete.cases(dat[, c("ABCA1_ng_ml", gvar)])
  if (sum(ok) < 3) return(tibble(group = gvar, test = NA, stat = NA, p = NA, effect = NA, n = sum(ok)))
  
  df_sub <- dat[ok, c("ABCA1_ng_ml", gvar)]
  names(df_sub) <- c("ABCA1_ng_ml","grp")
  by_g <- split(df_sub$ABCA1_ng_ml, df_sub$grp)
  if (length(by_g) < 2 || min(lengths(by_g)) < 2)
    return(tibble(group = gvar, test = NA, stat = NA, p = NA, effect = NA, n = sum(ok)))
  
  # Normality + homoscedasticity
  p0 <- if (length(by_g[[1]]) >= 3) shapiro.test(by_g[[1]])$p.value else 0
  p1 <- if (length(by_g[[2]]) >= 3) shapiro.test(by_g[[2]])$p.value else 0
  lev_p <- tryCatch(car::leveneTest(df_sub$ABCA1_ng_ml ~ df_sub$grp)$`Pr(>F)`[1], error = function(e) NA_real_)
  
  if (!is.na(p0) && !is.na(p1) && p0 >= .05 && p1 >= .05 && !is.na(lev_p) && lev_p >= .05) {
    # Student's t-test (equal variances)
    t <- t.test(ABCA1_ng_ml ~ grp, data = df_sub, var.equal = TRUE)
    eff <- rstatix::cohens_d(df_sub, ABCA1_ng_ml ~ grp, var.equal = TRUE)
    tibble(group = gvar, test = "Student t", stat = unname(t$statistic), p = t$p.value,
           effect = eff$effsize[1], n = sum(ok), lev_p = lev_p)
  } else {
    # Wilcoxon
    w <- wilcox.test(ABCA1_ng_ml ~ grp, data = df_sub, exact = FALSE)
    eff <- rstatix::wilcox_effsize(df_sub, ABCA1_ng_ml ~ grp)
    tibble(group = gvar, test = "Wilcoxon rank-sum", stat = unname(w$statistic), p = w$p.value,
           effect = eff$effsize[1], n = sum(ok), lev_p = lev_p)
  }
}

grp_vars    <- paste0(scores_to_use, "_grp")
grp_results <- bind_rows(lapply(grp_vars, group_test_one)) %>%
  mutate(p_fdr = p.adjust(p, method = "fdr"))
print(grp_results)

# ---- E2) Visualization for grouped comparisons ----
group_colors <- c("High" = "#9AC5CD", "Low" = "#FFDAB9")

plot_grouped <- function(df, score, y = "ABCA1_ng_ml") {
  gvar <- paste0(score, "_grp")
  if (!gvar %in% names(df)) return(NULL)
  d <- df %>% select(all_of(c(gvar, y))) %>% rename(grp = !!gvar) %>% filter(complete.cases(.))
  if (nrow(d) < 3 || length(unique(d$grp)) < 2) return(NULL)
  
  # Normality & Levene for choosing plot/test
  by_g <- split(d[[y]], d$grp)
  p0 <- if (length(by_g[[1]]) >= 3) shapiro.test(by_g[[1]])$p.value else 0
  p1 <- if (length(by_g[[2]]) >= 3) shapiro.test(by_g[[2]])$p.value else 0
  lev_p <- tryCatch(car::leveneTest(d[[y]] ~ d$grp)$`Pr(>F)`[1], error = function(e) NA_real_)
  y_max <- max(d[[y]], na.rm = TRUE)
  
  if (!is.na(p0) && !is.na(p1) && p0 >= .05 && p1 >= .05 && !is.na(lev_p) && lev_p >= .05) {
    # Student t-test p-value for annotation
    pval <- t.test(ABCA1_ng_ml ~ grp, data = d, var.equal = TRUE)$p.value
    p <- ggplot(d, aes(x = grp, y = .data[[y]], fill = grp)) +
      stat_summary(fun = mean, geom = "bar", color = "black", width = 0.6) +
      stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, size = 0.8) +
      geom_jitter(shape = 21, size = 2.5, width = 0.15, stroke = 0.3,
                  color = "black", alpha = 0.9) +
      scale_fill_manual(values = group_colors) +
      labs(x = NULL, y = "Serum ABCA1 Protein (ng/mL)",
           title = paste0(score, " (High vs Low)")) +
      theme_classic() +
      theme(axis.title.y = element_text(size = 14, face = "bold"),
            axis.text = element_text(size = 12, face = "bold"),
            legend.position = "none",
            axis.line = element_line(size = 0.8),
            axis.ticks = element_line(size = 0.8),
            plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) +
      annotate("text", x = 1.5, y = y_max * 1.05, label = p_to_label(pval),
               size = 5, fontface = "bold")
  } else {
    # Wilcoxon p-value for annotation
    pval <- wilcox.test(ABCA1_ng_ml ~ grp, data = d, exact = FALSE)$p.value
    p <- ggplot(d, aes(x = grp, y = .data[[y]], fill = grp)) +
      geom_boxplot(width = 0.55, color = "black", outlier.shape = NA, alpha = 0.9) +
      geom_jitter(shape = 21, size = 2.5, width = 0.15, stroke = 0.3,
                  color = "black", alpha = 0.9) +
      scale_fill_manual(values = group_colors) +
      labs(x = NULL, y = "Serum ABCA1 Protein (ng/mL)",
           title = paste0(score, " (High vs Low)")) +
      theme_classic() +
      theme(axis.title.y = element_text(size = 14, face = "bold"),
            axis.text = element_text(size = 12, face = "bold"),
            legend.position = "none",
            axis.line = element_line(size = 0.8),
            axis.ticks = element_line(size = 0.8),
            plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) +
      annotate("text", x = 1.5, y = y_max * 1.05, label = p_to_label(pval),
               size = 5, fontface = "bold")
  }
  p
}

for (s in scores_to_use) {
  p_g <- plot_grouped(dat, s, "ABCA1_ng_ml")
  if (!is.null(p_g)) {
    ggsave(filename = paste0("ABCA1_", s, "_HighLow.pdf"),
           plot = p_g, width = 5, height = 5, device = "pdf")
  }
}

# ---------------------------------------------------------
# D) Optional: Early vs Non-early ROC (AUC)
# ---------------------------------------------------------
if ("ALSFRS_total" %in% scores_to_use) {
  q3 <- quantile(dat$ALSFRS_total, 0.75, na.rm = TRUE)
  dat$early <- factor(ifelse(dat$ALSFRS_total >= q3, "Early","NonEarly"))
  if (nlevels(droplevels(dat$early)) == 2) {
    roc1 <- roc(dat$early, dat$ABCA1_ng_ml, levels = c("NonEarly","Early"), direction = ">")
    message(sprintf("AUC (ABCA1 for Early vs NonEarly): %.3f", as.numeric(roc1$auc)))
    # If you want a ROC PDF: uncomment next two lines
    # pdf("ABCA1_ROC_Early_vs_NonEarly.pdf", width = 5, height = 5)
    # plot.roc(roc1, print.auc = TRUE); dev.off()
  } else {
    message("Skipping ROC: early/non-early collapsed to one class in this sample.")
  }
}

readr::write_csv(cor_tab,   "ALSFRSR_ABCA1_spearman_results.csv")
readr::write_csv(fit_list,  "ALSFRSR_ABCA1_adjusted_lm_results.csv")
readr::write_csv(grp_results,"ALSFRSR_ABCA1_group_tests.csv")
############################################################
# METHODS (commented for manuscript)
#
# Study question.
#   We examined whether serum ABCA1 levels relate to ALS functional status at
#   assessment using the ALSFRS-R total score and its four domains
#   (BULBAR, FINE MOTOR, GROSS MOTOR, RESPIRATORY).
#
# Data handling.
#   Column names were harmonized by trimming leading/trailing spaces and
#   standardizing key variables to snake_case (e.g., “Riluzole Use ” → Riluzole_Use).
#   Sex was encoded as a binary factor (Female/Male). Medication use variables
#   (e.g., lipid-lowering drugs, riluzole) were encoded as factors (No/Yes).
#   ALSFRS-R scores and ABCA1 (ng/mL) were analyzed on their native scales.
#
# Primary association (correlation).
#   Spearman rank correlations were computed between ABCA1 and each ALSFRS-R
#   score (total and four domains). Two-sided p-values were reported and
#   false-discovery rate (FDR) correction was applied across the five tests.
#   For visualization, scatterplots with LOESS smoothing, marginal histograms,
#   and in-panel annotations of Spearman’s rho and unadjusted p-values were used.
#
# Adjusted analyses (regression).
#   For each ALSFRS-R score, a separate linear regression modeled ABCA1 as the
#   outcome with the score as the main predictor, adjusting for available
#   covariates (Age, Sex, BMI, eGFR, LDL, lipid-lowering drug use, riluzole use).
#   Prior to modeling, covariates with only one observed level or value in the
#   sample were removed a priori to avoid contrasts errors. To aid interpretation,
#   standardized beta coefficients (and 95% CIs) for the ALSFRS-R score were
#   obtained by refitting models with z-scored outcome and predictors.
#   Nonlinearity of the ALSFRS-R total score was explored using a natural cubic
#   spline (df = 3) specification (results descriptive).
#
# Grouped comparisons.
#   Each ALSFRS-R score was dichotomized at the within-sample median (High vs Low).
#   For ABCA1 differences between groups, normality within each group was assessed
#   using the Shapiro–Wilk test and variance homogeneity was assessed using
#   Levene’s test. If both groups were approximately normal (p ≥ 0.05) and
#   variances were homogeneous (p ≥ 0.05), Student’s two-sample t-test
#   (var.equal = TRUE) was used and results were visualized using bar charts
#   of mean ± SE with jittered individual points. Otherwise, the Wilcoxon
#   rank-sum test was applied and visualized using boxplots with jitter.
#   Effect sizes (Cohen’s d for t-tests; rank-biserial for Wilcoxon) accompanied
#   p-values; FDR correction across the five grouped comparisons was reported.
#
# Exploratory discrimination.
#   “Early ALS” was defined a priori as the upper quartile of the ALSFRS-R total
#   score. The ability of ABCA1 to distinguish Early vs Non-early was assessed
#   using ROC analysis (AUC). This exploratory step is sample-size-limited and
#   intended for hypothesis generation.
#
# Statistical environment.
#   Analyses were performed in R using packages: readxl, dplyr, ggplot2, ggpubr,
#   ggExtra, broom, rstatix, car, pROC, and splines. All tests were two-sided with
#   α = 0.05. Figures were saved as PDF only.
############################################################
