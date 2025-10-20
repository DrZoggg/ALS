#############################################
######### GSE250455 Analysis Script 
######### (aggregate duplicate genes by max)
#############################################

## =========================
## 0) Load required packages
## =========================
invisible(lapply(c("tidyverse","DESeq2","ggplot2","ggrepel","pheatmap"),
                 library, character.only = TRUE))
library(readr)

## =========================
## 1) Load raw count data
## =========================
setwd("External Validation Cohort")
expr_df <- read_csv("GSE250455_counts_raw.csv")

## =========================
## 2) Build count matrix 
##    - Aggregate duplicated gene names using max value
## =========================
counts_df <- expr_df %>%
  select(gene_name, everything(), -gene_id) %>%
  group_by(gene_name) %>%
  summarise(across(where(is.numeric), ~ max(.x, na.rm = TRUE)), .groups = "drop")

counts <- counts_df %>%
  column_to_rownames("gene_name") %>%
  as.matrix()
storage.mode(counts) <- "integer"

dup_n <- nrow(expr_df) - nrow(counts)
cat("Number of duplicate gene rows removed after aggregation: ", dup_n, "\n")

## =========================
## 3) Construct sample metadata
##    - Parse subject IDs and Pre/Post condition
## =========================
sample_table <- tibble(sample = colnames(counts)) %>%
  tidyr::separate(sample, into = c("subject","condition"), sep = "-", remove = FALSE) %>%
  mutate(
    condition = factor(condition, levels = c("Pre","Post")),
    subject   = factor(subject)
  ) %>% as.data.frame()
rownames(sample_table) <- sample_table$sample

## =========================
## 4) Filter out lowly expressed genes
##    - Keep genes with count ≥10 in at least 3 samples
## =========================
keep <- rowSums(counts >= 10) >= 3
counts_f <- counts[keep, , drop = FALSE]
cat("Number of genes retained after filtering: ", nrow(counts_f), "\n")

## =========================
## 5) Differential expression with DESeq2
##    - Paired design: ~ subject + condition
## =========================
dds <- DESeqDataSetFromMatrix(countData = counts_f,
                              colData   = sample_table,
                              design    = ~ subject + condition)
dds <- DESeq(dds)

# Test Post vs Pre
res <- results(dds, contrast = c("condition","Post","Pre"))

# List available coefficients
resultsNames(dds)

# Log fold-change shrinkage (DESeq2 recommended approach)
res_shrunk <- lfcShrink(dds, coef = "condition_Post_vs_Pre", type = "apeglm")

res_df <- as.data.frame(res_shrunk) %>%
  rownames_to_column("gene")

# Export all DE results
res_df %>%
  arrange(padj) %>%
  write.csv("GSE250455_DESeq2_Post_vs_Pre_all.csv", row.names = FALSE)

## =========================
## 6) Extract ABCA1 result
## =========================
find_gene_row <- function(df, symbol) {
  df %>% filter(toupper(gene) == toupper(symbol))
}
abca1_row <- find_gene_row(res_df, "ABCA1")

if (nrow(abca1_row) == 0) {
  message("ABCA1 not found in differential expression results. Check gene names.")
} else {
  cat(sprintf("\nDESeq2 (Post vs Pre) ABCA1: log2FC = %.3f, p = %.4g, padj = %.4g\n",
              abca1_row$log2FoldChange, abca1_row$pvalue, abca1_row$padj))
  write.csv(abca1_row, "GSE250455_DESeq2_ABCA1.csv", row.names = FALSE)
}

## =========================
## 7) VST normalization (for visualization)
## =========================
vsd <- vst(dds, blind = FALSE)
vsd_mat <- assay(vsd)

get_gene_vsd <- function(mat, symbol) {
  r <- which(toupper(rownames(mat)) == toupper(symbol))
  if (length(r) == 0) return(NULL)
  v <- mat[r[1], ]
  names(v) <- colnames(mat)
  v
}
abca1_v <- get_gene_vsd(vsd_mat, "ABCA1")

## =========================
## 8) Visualization: ABCA1 Pre vs Post
##    - Bar + jitter plot
##    - Dashed lines for paired samples
##    - Annotated with DESeq2 p-value
## =========================
# Build per-sample dataframe
plot_df <- data.frame(
  sample   = names(abca1_v),
  vsd_expr = as.numeric(abca1_v)
) %>%
  left_join(sample_table %>% select(sample, subject, condition), by = "sample")

# Summarize mean ± SE per condition
bar_df <- plot_df %>%
  group_by(condition) %>%
  summarise(
    mean = mean(vsd_expr, na.rm = TRUE),
    se   = sd(vsd_expr, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# Colors for groups
group_colors <- c("Pre" = "#9AC5CD", "Post" = "#FFDAB9")

# Plot
p <- ggplot(plot_df, aes(x = condition, y = vsd_expr, fill = condition)) +
  # Mean bars + error bars
  stat_summary(fun = mean, geom = "bar", color = "black", width = 0.6) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, size = 0.8) +
  
  # Paired dashed lines per subject
  geom_line(aes(group = subject), color = "gray20", linetype = "dashed", alpha = 0.9) +
  
  # Individual sample points
  geom_jitter(shape = 21, size = 3, width = 0.12, stroke = 0.6,
              color = "black", alpha = 0.9) +
  
  scale_fill_manual(values = group_colors) +
  labs(
    title = "GSE250455: ALS Pre vs Post Exercise — ABCA1",
    x = NULL,
    y = "VST-normalized ABCA1 expression"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.y    = element_text(size = 14, face = "bold"),
    axis.text       = element_text(size = 12, face = "bold"),
    axis.line       = element_line(size = 0.8),
    axis.ticks      = element_line(size = 0.8),
    plot.title      = element_text(size = 16, face = "bold", hjust = 0.5)
  )

# Significance annotation
pval <- if (nrow(abca1_row)) abca1_row$pvalue else NA_real_
y_max <- max(bar_df$mean + bar_df$se, na.rm = TRUE)
y_brkt <- y_max * 1.12
y_label <- y_brkt * 1.03

p <- p +
  annotate("segment", x = 1, xend = 2, y = y_brkt, yend = y_brkt, linewidth = 0.8) +
  annotate("text", x = 1.5, y = y_label,
           label = sprintf("p = %.3g (DESeq2)", pval), size = 5)

# Print plot
p
ggsave("GSE250455_ABCA1_plot.pdf", p, width = 6, height = 5)


rm(list=ls())




#############################################
######### GSE212131 Analysis Script
######### (RMA-normalized microarray; Short vs Long ALS)
#############################################

## =========================
## 0) Load required packages
## =========================
invisible(lapply(c("tidyverse","limma","ggplot2"), library, character.only = TRUE))
library(readr)

## =========================
## 1) Load phenotype and expression data
## =========================
pheno <- read_csv("GSE212131_pheno.csv")           # Must contain 'geo_accession' and 'title'
expr  <- read_csv("GSE212131_RMA_normalized.csv")  # Rows: probes/genes; Columns: GSM sample IDs

## =========================
## 2) Build expression matrix 
##    - First column is gene symbol
##    - Keep only GSM columns present in pheno$geo_accession
##    - Aggregate duplicate gene symbols by column-wise max
## =========================
gsm_cols <- intersect(names(expr)[-1], pheno$geo_accession)  # Exclude first column (gene names)

# Select gene name column + GSM columns that appear in pheno$geo_accession
expr_sub <- expr %>%
  dplyr::select(1, all_of(gsm_cols)) %>%
  dplyr::rename(gene_symbol = 1)

# Aggregate duplicate genes by column-wise max
expr_by_gene <- expr_sub %>%
  dplyr::group_by(gene_symbol) %>%
  dplyr::summarise(dplyr::across(where(is.numeric), ~ max(.x, na.rm = TRUE)), .groups = "drop")

# Convert to matrix: rows = gene, cols = GSM
emat <- expr_by_gene %>%
  tibble::column_to_rownames("gene_symbol") %>%
  as.matrix()

# Report number of rows collapsed
dup_removed <- nrow(expr) - nrow(emat)
cat("Number of duplicate gene rows collapsed by max: ", dup_removed, "\n")

## =========================
## 3) Build design matrix: Short vs Long ALS
##    - Parse group from pheno$title ("Short-1" -> "Short"; "Long-2" -> "Long")
##    - Align sample order between metadata and expression matrix
## =========================
pheno2 <- pheno %>%
  dplyr::mutate(
    group = sub("^([^-]+).*", "\\1", title)   # Extract token before the first hyphen
  ) %>%
  dplyr::filter(geo_accession %in% colnames(emat)) %>%
  dplyr::mutate(group = factor(group, levels = c("Short","Long"))) %>%
  dplyr::arrange(match(geo_accession, colnames(emat)))

# Reorder columns of emat to match pheno2
emat <- emat[, pheno2$geo_accession, drop = FALSE]

# Quick sanity check
print(table(pheno2$group))

## =========================
## 4) Differential expression with limma (RMA data)
##    - Model: ~ 0 + group  (Short vs Long)
##    - Contrast: Long - Short (you can invert if needed)
## =========================
design <- model.matrix(~ 0 + pheno2$group)
colnames(design) <- levels(pheno2$group)  # "Short","Long"

fit <- lmFit(emat, design)
contrast.mat <- makeContrasts(Long_vs_Short = Long - Short, levels = design)
fit2 <- contrasts.fit(fit, contrast.mat)
fit2 <- eBayes(fit2)

# Full table
tt <- topTable(fit2, coef = "Long_vs_Short", n = Inf, sort.by = "P")
tt <- tibble::rownames_to_column(tt, var = "gene")

# Export full results
readr::write_csv(tt, "GSE212131_limma_Long_vs_Short_all.csv")

## =========================
## 5) Extract ABCA1 row & report stats
## =========================
abca1_row <- tt %>% dplyr::filter(toupper(gene) == "ABCA1")
if (nrow(abca1_row) == 0) {
  warning("ABCA1 not found in limma results. Check gene symbol column.")
} else {
  cat(sprintf("limma (Long vs Short) ABCA1: logFC = %.3f, p = %.4g, adj.P.Val = %.4g\n",
              abca1_row$logFC, abca1_row$P.Value, abca1_row$adj.P.Val))
  readr::write_csv(abca1_row, "GSE212131_limma_ABCA1.csv")
}

## =========================
## 6) Build per-sample dataframe for ABCA1 (for plotting)
##    - Use the already-normalized RMA expression
## =========================
get_gene_expr <- function(mat, symbol) {
  r <- which(toupper(rownames(mat)) == toupper(symbol))
  if (length(r) == 0) return(NULL)
  v <- mat[r[1], ]
  names(v) <- colnames(mat)
  v
}
abca1_expr <- get_gene_expr(emat, "ABCA1")
if (is.null(abca1_expr)) stop("ABCA1 not found in expression matrix.")

plot_df <- data.frame(
  sample = names(abca1_expr),
  expr   = as.numeric(abca1_expr)
) %>%
  dplyr::left_join(pheno2 %>% dplyr::select(geo_accession, group), 
                   by = c("sample" = "geo_accession"))

## =========================
## 7) Summary stats for bars (mean ± SE per group)
## =========================
bar_df <- plot_df %>%
  dplyr::group_by(group) %>%
  dplyr::summarise(
    mean = mean(expr, na.rm = TRUE),
    se   = stats::sd(expr, na.rm = TRUE) / sqrt(dplyr::n()),
    .groups = "drop"
  )

## =========================
## 8) Visualization: ABCA1 (Short vs Long)
##    - Bar + jitter (no paired lines; cross-sectional groups)
##    - Annotate exact p-value from limma
## =========================

# Colors for groups (use more vibrant colors)
group_colors <- c("Short" = "#9AC5CD", "Long" =  "#FFDAB9")

# Build plot
p <- ggplot(plot_df, aes(x = group, y = expr, fill = group)) +
  # Mean bars + error bars (mean ± SE)
  stat_summary(fun = mean, geom = "bar", color = "black", width = 0.7, size = 0.8) +  # Wider bars, black border
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.15, size = 0.8) +  # Error bars
  
  # Individual points (adjust jitter size, transparency, and border)
  geom_jitter(shape = 21, size = 4, width = 0.12, stroke = 0.8, color = "black", alpha = 0.7) +  # Bigger, more visible
  
  scale_fill_manual(values = group_colors) +
  labs(
    title = "GSE212131: ABCA1 Expression in Short vs Long ALS",
    x = NULL,
    y = "ABCA1 Expression (RMA-normalized)"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.y    = element_text(size = 14, face = "bold"),
    axis.text       = element_text(size = 12, face = "bold"),
    axis.line       = element_line(size = 0.8),
    axis.ticks      = element_line(size = 0.8),
    plot.title      = element_text(size = 16, face = "bold", hjust = 0.5),  # Center title
    plot.margin     = margin(10, 10, 10, 10)
  )

# Significance annotation (limma p-value)
pval <- if (nrow(abca1_row)) abca1_row$P.Value else NA_real_
y_max   <- max(bar_df$mean + bar_df$se, na.rm = TRUE)
y_brkt  <- y_max * 1.12  # Place bracket slightly above the highest bar
y_label <- y_brkt * 1.05  # Position p-value label just above the bracket

p <- p +
  annotate("segment", x = 1, xend = 2, y = y_brkt, yend = y_brkt, linewidth = 0.8) +  # Bracket line between bars
  annotate("text",    x = 1.5, y = y_label,
           label = sprintf("p = %.3g (limma)", pval), size = 5, hjust = 0.5, fontface = "bold")  # Centered p-value

# Print plot
p

## =========================
## 9) (Optional) Save figure files
ggsave("GSE212131_ABCA1_Short_vs_Long.pdf",  p, width = 6, height = 5)
rm(list=ls())






#############################################
########## GSE234297 Analysis Script
########## (Bulk RNA-seq; ALS vs Healthy Control)
#############################################

## =========================
## 0) Load required packages
## =========================
invisible(lapply(c("tidyverse","DESeq2","ggplot2","ggrepel"), library, character.only = TRUE))
library(readr)

## =========================
## 1) Load counts and phenotype data
## =========================
GSE234297_counts <- read_csv("GSE_expr_and_phen/GSE234297__counts_raw.csv")  # Counts matrix
GSE234297_pheno  <- read_csv("GSE_expr_and_phen/GSE234297_pheno.csv")         # Phenotype data

## =========================
## 2) Prepare expression matrix
##    - First column is gene names
##    - Keep only ALS (Case) and CON (Control) columns
## =========================
# Ensure the expression matrix has genes in the first column
expr_df <- GSE234297_counts %>%
  dplyr::rename(gene_name = 1)  # Rename first column to 'gene_name'

# Extract ALS and CON samples (columns)
sample_cols <- colnames(expr_df)[2:ncol(expr_df)]  # Skip first column (gene names)

# Create a group vector based on column names (Case and Control)
group <- ifelse(grepl("Case", sample_cols), "ALS", "Control")

# Extract expression matrix for ALS and CON samples
expr_matrix <- expr_df %>%
  dplyr::select(gene_name, all_of(sample_cols)) %>%
  column_to_rownames("gene_name") %>%
  as.matrix()

# Sanity check on data
cat("Expression matrix dimensions: ", dim(expr_matrix), "\n")
cat("First few columns (samples): ", colnames(expr_matrix)[1:5], "\n")

## =========================
## 3) Prepare phenotype data (group information)
## =========================
# Create phenotype dataframe with group information
pheno <- data.frame(
  geo_accession = sample_cols,
  group = factor(group, levels = c("Control", "ALS"))
)

# Make sure the sample order in phenotype matches the expression matrix columns
pheno <- pheno[match(colnames(expr_matrix), pheno$geo_accession), ]
expr_matrix <- round(expr_matrix)  # Convert to integers by rounding
## =========================
## 4) DESeq2 Differential Expression Analysis
##    - Use design: ~ group
##    - Compare ALS vs Control
## =========================
# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = expr_matrix,
                              colData = pheno,
                              design = ~ group)

# Run DESeq2 analysis
dds <- DESeq(dds)

# Extract results (ALS vs Control)
res <- results(dds, contrast = c("group", "ALS", "Control"))

# Shrink logFC for more stable results
res_shrunk <- lfcShrink(dds, coef = "group_ALS_vs_Control", type = "apeglm")

# Convert to dataframe and reorder
res_df <- as.data.frame(res_shrunk) %>%
  rownames_to_column("gene")

# Export full results
write.csv(res_df, "GSE234297_DESeq2_Results.csv", row.names = FALSE)

## =========================
## 5) Extract ABCA1 results (DESeq2)
## =========================
abca1_row <- res_df %>% filter(toupper(gene) == "ABCA1")

if (nrow(abca1_row) == 0) {
  message("ABCA1 not found in DESeq2 results. Please check gene symbol.")
} else {
  cat(sprintf("DESeq2 (ALS vs Control) ABCA1: log2FC = %.3f, p = %.4g, padj = %.4g\n",
              abca1_row$log2FoldChange, abca1_row$pvalue, abca1_row$padj))
  write.csv(abca1_row, "GSE234297_DESeq2_ABCA1.csv", row.names = FALSE)
}

## =========================
## 6) VST normalization (for visualization)
## =========================
# VST normalization for visualization
vsd <- vst(dds, blind = FALSE)
vsd_mat <- assay(vsd)

# Extract ABCA1 VST expression
abca1_vsd <- vsd_mat["ABCA1", ]

# Create plot dataframe (including group info)
plot_df <- data.frame(
  sample = names(abca1_vsd),
  vsd_expr = as.numeric(abca1_vsd)
) %>%
  left_join(pheno %>% select(geo_accession, group), by = c("sample" = "geo_accession"))

## =========================
## 7) Summary statistics for plotting (mean ± SE)
## =========================
bar_df <- plot_df %>%
  group_by(group) %>%
  summarise(
    mean = mean(vsd_expr, na.rm = TRUE),
    se   = sd(vsd_expr, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )


## =========================
## 8) Visualization: ABCA1 (ALS vs Control)
##    - Violin plot with significance annotation and improved boxplot
## =========================

# Define colors for groups (Control vs ALS)
group_colors <- c("Control" = "#FFDAB9", "ALS" = "#9AC5CD")

# Build violin plot with improved boxplot
p <- ggplot(plot_df, aes(x = group, y = vsd_expr, fill = group)) + 
  # Violin plot to show data distribution
  geom_violin(trim = FALSE, color = "black", alpha = 1) +  # Adjust transparency of the violin
  
  # Add individual points
  geom_jitter(shape = 21, size = 0.5, width = 0.12, stroke = 0.6, 
              color = "black", alpha = 0.7) +
  
  # Add boxplot within violin plot with custom colors and thick black median line
  geom_boxplot(width = 0.15, color = "black", fill = "white", 
               alpha = 0.8, outlier.shape = NA, 
               median.props = list(size = 2, color = "black")) +  # Custom median line style
  
  scale_fill_manual(values = group_colors) +  # Custom colors
  labs(
    title = "GSE234297: ABCA1 Expression in ALS vs Control",
    x = NULL,
    y = "VST-normalized ABCA1 Expression"  # Correct y-axis label for VST-normalized data
  ) +
  theme_classic() +
  theme(
    legend.position = "none",  # Remove legend
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, face = "bold"),
    axis.line = element_line(size = 0.8),
    axis.ticks = element_line(size = 0.8),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)  # Centered title
  )

# Add significance annotation (from DESeq2 results)
pval <- if (nrow(abca1_row)) abca1_row$pvalue else NA_real_

# Set the threshold for significance
significance <- ifelse(pval <= 0.001, "***", ifelse(pval <= 0.01, "**", ifelse(pval <= 0.05, "*", "n.s.")))

y_max   <- max(plot_df$vsd_expr, na.rm = TRUE)
y_brkt  <- y_max * 1.1  # Position bracket above the max value
y_label <- y_brkt * 1.03  # Position p-value label above the bracket

p <- p +
  # Add bracket between violins (horizontal line)
  annotate("segment", x = 1, xend = 2, y = y_brkt, yend = y_brkt, linewidth = 0.8) +  # Bracket between violins
  annotate("text", x = 1.5, y = y_label,
           label = sprintf("%s\np = %.3g (DESeq2)", significance, pval), size = 5)  # P-value annotation

# Print the plot
p
ggsave("GSE234297_ABCA1_ALS_vs_CON.pdf",  p, width = 6, height = 5)
rm(list=ls())











#############################################
# GSE153960 (postmortem brain bulk RNA-seq)
# Per-tissue DE analysis for ABCA1:
#   - ALS vs Control (always when both exist)
#   - ALS vs OND (Other Neurological Disorders) when OND exists
# Pipeline:
#   counts -> collapse duplicates -> filter & align samples
#   -> per-tissue DESeq2 (~ group_simple) -> extract ABCA1
#   -> VST visualization (violin + box) faceted and per-tissue PDFs
#############################################

## =========================
## 0) Dependencies
## =========================
suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(readr)
  library(stringr)
})

## =========================
## 1) I/O paths
## =========================
counts_path <- "GSE_expr_and_phen/GSE153960_counts_raw.csv"   # first column = gene
pheno_path  <- "GSE_expr_and_phen/GSE153960_pheno.csv"
out_dir     <- "GSE153960_ABCA1_output"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

## =========================
## 2) Load counts and phenotype
## =========================
pheno_raw  <- read_csv(pheno_path, show_col_types = FALSE)
counts_raw <- read_csv(counts_path, show_col_types = FALSE)  # first column = gene

## =========================
## 3) Build counts matrix
##    - Collapse duplicated genes by row-wise maximum
##    - Keep integer storage mode (for DESeq2)
## =========================
counts_df <- as.data.frame(counts_raw)
stopifnot(ncol(counts_df) >= 2)
gene_col <- colnames(counts_df)[1]

counts_collapsed <- counts_df %>%
  dplyr::rename(gene = !!gene_col) %>%
  group_by(gene) %>%
  summarise(across(where(is.numeric), ~ max(.x, na.rm = TRUE)), .groups = "drop")

counts_mat <- counts_collapsed %>%
  column_to_rownames("gene") %>%
  as.matrix()
storage.mode(counts_mat) <- "integer"

## =========================
## 4) Prepare phenotype (keep ALS / Control / OND) and align samples
## =========================
keep_groups <- c("ALS Spectrum MND",
                 "Non-Neurological Control",
                 "Other Neurological Disorders")

pheno <- pheno_raw %>%
  filter(group %in% keep_groups) %>%
  select(GSM, group, source_name) %>%
  mutate(
    group_simple = case_when(
      str_detect(group, regex("^ALS Spectrum MND$", ignore_case = TRUE))              ~ "ALS",
      str_detect(group, regex("^Non-Neurological Control$", ignore_case = TRUE))      ~ "Control",
      str_detect(group, regex("^Other Neurological D[iI]sorders$", ignore_case = TRUE)) ~ "OND",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(group_simple))

# Align by GSM intersection and consistent order
common_gsm <- intersect(colnames(counts_mat), pheno$GSM)
counts_mat <- counts_mat[, common_gsm, drop = FALSE]
pheno      <- pheno[match(common_gsm, pheno$GSM), , drop = FALSE]

# Factors (reference = Control; order: Control, OND, ALS)
pheno$group_simple <- factor(pheno$group_simple, levels = c("Control","OND","ALS"))
pheno$source_name  <- droplevels(factor(pheno$source_name))

## 5) Utility functions
## =========================

# Keep genes with sufficient counts in this subset
keep_genes_subset <- function(counts_subset, min_count = 10, min_samples = 3) {
  rowSums(counts_subset >= min_count) >= min_samples
}

# Safe filename for tissue
safe_tissue <- function(x) {
  x %>%
    stringr::str_replace_all("[^A-Za-z0-9]+", "_") %>%
    stringr::str_replace("^_+|_+$", "")
}

# Consistent colors for groups
group_colors <- c("Control" = "#FFDAB9", "OND" = "#A3BE8C", "ALS" = "#9AC5CD")

# Normalize a DESeq2 results object to a data.frame with stable columns
normalize_res_df <- function(res_obj) {
  df <- as.data.frame(res_obj) %>% tibble::rownames_to_column("gene")
  needed <- c("gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")
  for (nm in needed) {
    if (!nm %in% names(df)) df[[nm]] <- NA_real_
  }
  df
}

# Safe extractor for a pairwise contrast A vs B within a single-factor DESeq2 model
# Priority:
#   1) If 'name' exists for A_vs_B, use it and try LFC shrinkage (apeglm)
#   2) Else if 'name' exists for B_vs_A, use it and flip signs
#   3) Else try contrast=c(factor, A, B)
#   4) Else return NULL
safe_results_contrast <- function(dds, factor_name = "group_simple", A = "ALS", B = "Control") {
  nms <- resultsNames(dds)
  nm1 <- paste0(factor_name, "_", A, "_vs_", B)
  nm2 <- paste0(factor_name, "_", B, "_vs_", A)
  
  if (nm1 %in% nms) {
    res <- results(dds, name = nm1)
    sh  <- try(lfcShrink(dds, coef = nm1, type = "apeglm"), silent = TRUE)
    if (!inherits(sh, "try-error")) res <- sh
    return(res)
  }
  
  if (nm2 %in% nms) {
    res2 <- results(dds, name = nm2)
    # flip to get A_vs_B
    res2$log2FoldChange <- -res2$log2FoldChange
    if (!is.null(res2$stat)) res2$stat <- -res2$stat
    return(res2)
  }
  
  res_try <- try(results(dds, contrast = c(factor_name, A, B)), silent = TRUE)
  if (!inherits(res_try, "try-error")) return(res_try)
  
  return(NULL)
}

# --- Plot function with points + p-value brackets (DESeq2) ---
# Build one violin+box plot (per tissue) on VST-normalized ABCA1
# p_AC = p-value for ALS vs Control (numeric or NA)
# p_AO = p-value for ALS vs OND     (numeric or NA)
plot_abca1_vst_tissue <- function(vsd, pheno_t, tissue_label, p_AC = NA_real_, p_AO = NA_real_) {
  mat <- assay(vsd)
  if (!"ABCA1" %in% rownames(mat)) return(NULL)
  
  # Prepare plotting dataframe
  v <- as.numeric(mat["ABCA1", ])
  df <- data.frame(
    sample   = colnames(mat),
    vsd_expr = v
  ) %>%
    dplyr::left_join(pheno_t %>% dplyr::select(GSM, group_simple), by = c("sample" = "GSM"))
  
  # Force order: ALS first, Control second, then OND if present
  desired_order <- c("ALS","Control","OND")
  present <- desired_order[desired_order %in% unique(as.character(df$group_simple))]
  df$group_simple <- factor(df$group_simple, levels = present)
  
  # Base plot: violin + white box + sample points (size = 70% of previous)
  p <- ggplot(df, aes(x = group_simple, y = vsd_expr, fill = group_simple)) +
    geom_violin(trim = FALSE, color = "black", alpha = 0.7) +
    geom_boxplot(width = 0.18, fill = "white", color = "black",
                 alpha = 1, outlier.shape = NA, linewidth = 0.7) +
    geom_jitter(shape = 21, size = 1.26, stroke = 0.5, color = "black",  # 70% of 1.8
                width = 0.12, alpha = 0.5) +
    scale_fill_manual(values = group_colors, drop = FALSE) +
    labs(
      title = paste0("GSE153960 — ", tissue_label, ": ABCA1 (VST)"),
      x = NULL,
      y = "VST-normalized ABCA1 expression"
    ) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.title.y    = element_text(size = 12, face = "bold"),
      axis.text       = element_text(size = 11, face = "bold"),
      axis.line       = element_line(size = 0.6),
      axis.ticks      = element_line(size = 0.6),
      plot.title      = element_text(size = 14, face = "bold", hjust = 0.5)
    )
  
  # ---- P-value annotations (DESeq2) ----
  add_bracket <- function(pv, g1, g2, y, lbl_size = 3.8) {
    if (is.na(pv)) return(NULL)
    levs <- levels(df$group_simple)
    x1 <- match(g1, levs); x2 <- match(g2, levs)
    if (any(is.na(c(x1, x2)))) return(NULL)
    xmid <- (x1 + x2) / 2
    list(
      annotate("segment", x = x1, xend = x2, y = y, yend = y, linewidth = 0.7),
      annotate("text", x = xmid, y = y * 1.02,
               label = sprintf("p = %.3g (DESeq2)", pv), size = lbl_size)
    )
  }
  
  y_max <- max(df$vsd_expr, na.rm = TRUE)
  y1 <- y_max * 1.08
  y2 <- y_max * 1.16
  
  # ALS vs Control
  if (all(c("ALS","Control") %in% levels(df$group_simple))) {
    p <- p + add_bracket(p_AC, "Control", "ALS", y1)
  }
  # ALS vs OND
  if (all(c("ALS","OND") %in% levels(df$group_simple))) {
    yy <- if (all(c("ALS","Control") %in% levels(df$group_simple))) y2 else y1
    p  <- p + add_bracket(p_AO, "OND", "ALS", yy)
  }
  
  p
}

## =========================
## 6) Per-tissue DESeq2 and ABCA1 extraction
## =========================
tissues <- levels(pheno$source_name)

abca1_summary <- list()
per_tissue_plots <- list()

for (tissue in tissues) {
  idx <- pheno$source_name == tissue
  counts_t <- counts_mat[, idx, drop = FALSE]
  pheno_t  <- pheno[idx, , drop = FALSE]
  
  # Skip if less than 2 groups in this tissue
  present_groups <- sort(unique(as.character(pheno_t$group_simple)))
  if (length(present_groups) < 2) next
  
  # Drop unused levels and set a convenient reference
  pheno_t$group_simple <- droplevels(pheno_t$group_simple)
  if ("Control" %in% levels(pheno_t$group_simple)) {
    pheno_t$group_simple <- relevel(pheno_t$group_simple, ref = "Control")
  } else if ("OND" %in% levels(pheno_t$group_simple)) {
    pheno_t$group_simple <- relevel(pheno_t$group_simple, ref = "OND")
  } else {
    pheno_t$group_simple <- relevel(pheno_t$group_simple, ref = levels(pheno_t$group_simple)[1])
  }
  
  # Low-count filter
  keep_t <- keep_genes_subset(counts_t, min_count = 10, min_samples = 3)
  counts_t_f <- counts_t[keep_t, , drop = FALSE]
  
  # Build DESeq2
  col_df <- pheno_t %>% as.data.frame() %>% column_to_rownames("GSM")
  dds_t <- DESeqDataSetFromMatrix(countData = counts_t_f,
                                  colData   = col_df,
                                  design    = ~ group_simple)
  dds_t <- DESeq(dds_t)
  
  # Collect ABCA1 results
  res_list <- list()
  
  # ALS vs Control (if both exist)
  if (all(c("ALS","Control") %in% levels(pheno_t$group_simple))) {
    res_AC <- safe_results_contrast(dds_t, factor_name = "group_simple", A = "ALS", B = "Control")
    if (!is.null(res_AC)) {
      df <- normalize_res_df(res_AC)
      row <- df %>% dplyr::filter(toupper(gene) == "ABCA1")
      if (nrow(row) == 1) res_list[["ALS_vs_Control"]] <- row
    }
  }
  
  # ALS vs OND (if both exist)
  if (all(c("ALS","OND") %in% levels(pheno_t$group_simple))) {
    res_AO <- safe_results_contrast(dds_t, factor_name = "group_simple", A = "ALS", B = "OND")
    if (!is.null(res_AO)) {
      df <- normalize_res_df(res_AO)
      row <- df %>% dplyr::filter(toupper(gene) == "ABCA1")
      if (nrow(row) == 1) res_list[["ALS_vs_OND"]] <- row
    }
  }
  
  # Save per-tissue ABCA1 results CSV (tolerate missing 'stat')
  if (length(res_list) > 0) {
    out_tbl <- dplyr::bind_rows(res_list, .id = "contrast")
    wanted  <- c("contrast","gene","log2FoldChange","lfcSE","stat","pvalue","padj")
    out_tbl <- out_tbl[, intersect(wanted, names(out_tbl)), drop = FALSE]
    readr::write_csv(out_tbl, file.path(out_dir, paste0("ABCA1_DESeq2_", safe_tissue(tissue), ".csv")))
    abca1_summary[[tissue]] <- out_tbl %>% dplyr::mutate(tissue = tissue, .before = 1)
  }
  
  # Extract p-values for plot annotations
  p_AC <- if (!is.null(res_list[["ALS_vs_Control"]])) res_list[["ALS_vs_Control"]]$pvalue[1] else NA_real_
  p_AO <- if (!is.null(res_list[["ALS_vs_OND"]]))     res_list[["ALS_vs_OND"]]$pvalue[1]     else NA_real_
  
  # VST visualization with points and p-value brackets (DESeq2)
  vsd_t <- vst(dds_t, blind = FALSE)
  p_t   <- plot_abca1_vst_tissue(vsd_t, pheno_t, tissue_label = tissue, p_AC = p_AC, p_AO = p_AO)
  if (!is.null(p_t)) {
    per_tissue_plots[[tissue]] <- p_t
    ggsave(filename = file.path(out_dir, paste0("ABCA1_VST_", safe_tissue(tissue), ".pdf")),
           plot = p_t, width = 6.5, height = 5.0, device = cairo_pdf)
  }
}

## =========================
## 7) Write overall summary and a faceted PDF
## =========================
if (length(abca1_summary) > 0) {
  summary_tbl <- dplyr::bind_rows(abca1_summary) %>%
    dplyr::relocate(tissue, .before = 1)
  readr::write_csv(summary_tbl, file.path(out_dir, "ABCA1_DESeq2_summary_all_tissues.csv"))
} else {
  message("No ABCA1 rows found in DE results (check gene symbol case).")
}

# Optional faceted overview PDF (VST), skip if nothing collected
if (length(per_tissue_plots) > 0) {
  facet_df_list <- list()
  for (tissue in names(per_tissue_plots)) {
    idx <- pheno$source_name == tissue
    counts_t <- counts_mat[, idx, drop = FALSE]
    pheno_t  <- pheno[idx, , drop = FALSE]
    
    pheno_t$group_simple <- droplevels(pheno_t$group_simple)
    if ("Control" %in% levels(pheno_t$group_simple)) {
      pheno_t$group_simple <- relevel(pheno_t$group_simple, ref = "Control")
    } else if ("OND" %in% levels(pheno_t$group_simple)) {
      pheno_t$group_simple <- relevel(pheno_t$group_simple, ref = "OND")
    }
    
    keep_t <- keep_genes_subset(counts_t, min_count = 10, min_samples = 3)
    counts_t_f <- counts_t[keep_t, , drop = FALSE]
    col_df <- pheno_t %>% as.data.frame() %>% column_to_rownames("GSM")
    dds_t <- DESeqDataSetFromMatrix(countData = counts_t_f,
                                    colData   = col_df,
                                    design    = ~ group_simple)
    dds_t <- DESeq(dds_t)
    vsd_t <- vst(dds_t, blind = FALSE)
    
    mat <- assay(vsd_t)
    if (!"ABCA1" %in% rownames(mat)) next
    v <- as.numeric(mat["ABCA1", ])
    df <- data.frame(
      sample   = colnames(mat),
      vsd_expr = v
    ) %>%
      dplyr::left_join(pheno_t %>% dplyr::select(GSM, group_simple), by = c("sample" = "GSM")) %>%
      dplyr::mutate(tissue = tissue)
    facet_df_list[[tissue]] <- df
  }
  
  if (length(facet_df_list) > 0) {
    facet_df <- dplyr::bind_rows(facet_df_list)
    facet_df$group_simple <- factor(facet_df$group_simple, levels = c("Control","OND","ALS"))
    
    p_facet <- ggplot(facet_df, aes(x = group_simple, y = vsd_expr, fill = group_simple)) +
      geom_violin(trim = FALSE, color = "black", alpha = 0.7) +
      geom_boxplot(width = 0.15, fill = "white", color = "black",
                   alpha = 1, outlier.shape = NA, linewidth = 0.7) +
      # optional: show points in overview as well
      geom_jitter(shape = 21, size = 0.9, stroke = 0.4, color = "black",
                  width = 0.12, alpha = 0.4) +
      scale_fill_manual(values = group_colors, drop = FALSE) +
      facet_wrap(~ tissue, scales = "free_y") +
      labs(
        title = "GSE153960 — ABCA1 (VST) across tissues",
        x = NULL,
        y = "VST-normalized ABCA1 expression"
      ) +
      theme_classic() +
      theme(
        legend.position = "none",
        axis.title.y    = element_text(size = 12, face = "bold"),
        axis.text       = element_text(size = 10, face = "bold"),
        strip.text      = element_text(size = 10, face = "bold"),
        plot.title      = element_text(size = 14, face = "bold", hjust = 0.5)
      )
    
    ggsave(filename = file.path(out_dir, "ABCA1_VST_faceted_all_tissues.pdf"),
           plot = p_facet, width = 12, height = 8, device = cairo_pdf)
  }
}

message("Done. Outputs saved under: ", normalizePath(out_dir))
rm(list=ls())



#############################################
# GSE124439 (postmortem brain bulk RNA-seq)
# Per-tissue ABCA1 analysis:
#   - ALS vs Control
#   - ALS vs OND (Other Neurological Disorders)
# Pipeline:
#   counts -> collapse duplicates -> align samples
#   -> per-tissue DESeq2 (~ group_simple) -> extract ABCA1
#   -> VST visualization (violin + box + points) with DESeq2 p-values
#############################################

## =========================
## 0) Dependencies
## =========================
suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(readr)
  library(stringr)
})

## =========================
## 1) I/O paths
## =========================
counts_path <- "GSE_expr_and_phen/GSE124439_counts_raw.csv"  # first col = gene
pheno_path  <- "GSE_expr_and_phen/GSE124439_pheno.csv"
out_dir     <- "GSE124439_ABCA1_output"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

## =========================
## 2) Load counts and phenotype
## =========================
pheno_raw  <- read_csv(pheno_path,  show_col_types = FALSE)
counts_raw <- read_csv(counts_path, show_col_types = FALSE)

## =========================
## 3) Build counts matrix
##    - Collapse duplicated genes by row-wise maximum
##    - Keep integer storage mode (for DESeq2)
## =========================
counts_df <- as.data.frame(counts_raw)
stopifnot(ncol(counts_df) >= 2)
gene_col <- colnames(counts_df)[1]

counts_collapsed <- counts_df %>%
  dplyr::rename(gene = !!gene_col) %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(across(where(is.numeric), ~ max(.x, na.rm = TRUE)), .groups = "drop")

counts_mat <- counts_collapsed %>%
  column_to_rownames("gene") %>%
  as.matrix()
storage.mode(counts_mat) <- "integer"

## =========================
## 4) Prepare phenotype (ALS/Control/OND) and align samples
## =========================
# sample_group has 3 classes:
#   "ALS Spectrum MND", "Non-Neurological Control", "Other Neurological Disorders"
# source_name_ch1 has 2 tissues: "Frontal Cortex", "Motor Cortex"

keep_groups <- c("ALS Spectrum MND",
                 "Non-Neurological Control",
                 "Other Neurological Disorders")

pheno <- pheno_raw %>%
  dplyr::filter(sample_group %in% keep_groups) %>%
  dplyr::select(GSM, sample_group, source_name_ch1) %>%
  dplyr::mutate(
    group_simple = dplyr::case_when(
      str_detect(sample_group, regex("^ALS Spectrum MND$", ignore_case = TRUE))            ~ "ALS",
      str_detect(sample_group, regex("^Non-Neurological Control$", ignore_case = TRUE))    ~ "Control",
      str_detect(sample_group, regex("^Other Neurological D[iI]sorders$", ignore_case = TRUE)) ~ "OND",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(group_simple)) %>%
  dplyr::rename(tissue = source_name_ch1)

# Align by GSM intersection and consistent order
common_gsm <- intersect(colnames(counts_mat), pheno$GSM)
counts_mat <- counts_mat[, common_gsm, drop = FALSE]
pheno      <- pheno[match(common_gsm, pheno$GSM), , drop = FALSE]

# Factors (order for plotting later will be overridden to ALS -> Control -> OND)
pheno$group_simple <- factor(pheno$group_simple, levels = c("Control","OND","ALS"))
pheno$tissue       <- droplevels(factor(pheno$tissue))   # "Frontal Cortex", "Motor Cortex"

## =========================
## 5) Utility functions
## =========================

# Keep genes with sufficient counts in this subset
keep_genes_subset <- function(counts_subset, min_count = 10, min_samples = 3) {
  rowSums(counts_subset >= min_count) >= min_samples
}

# Safe filename for tissue
safe_tissue <- function(x) {
  x %>%
    stringr::str_replace_all("[^A-Za-z0-9]+", "_") %>%
    stringr::str_replace("^_+|_+$", "")
}

# Consistent colors for groups
group_colors <- c("Control" = "#FFDAB9", "OND" = "#A3BE8C", "ALS" = "#9AC5CD")

# Normalize a DESeq2 results object to a data.frame with stable columns
normalize_res_df <- function(res_obj) {
  df <- as.data.frame(res_obj) %>% tibble::rownames_to_column("gene")
  needed <- c("gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")
  for (nm in needed) {
    if (!nm %in% names(df)) df[[nm]] <- NA_real_
  }
  df
}

# Safe extractor for a pairwise contrast A vs B within a single-factor DESeq2 model
# Priority:
#   1) If 'name' exists for A_vs_B, use it and try LFC shrinkage (apeglm)
#   2) Else if 'name' exists for B_vs_A, use it and flip signs
#   3) Else try contrast=c(factor, A, B)
#   4) Else return NULL
safe_results_contrast <- function(dds, factor_name = "group_simple", A = "ALS", B = "Control") {
  nms <- resultsNames(dds)
  nm1 <- paste0(factor_name, "_", A, "_vs_", B)
  nm2 <- paste0(factor_name, "_", B, "_vs_", A)
  
  if (nm1 %in% nms) {
    res <- results(dds, name = nm1)
    sh  <- try(lfcShrink(dds, coef = nm1, type = "apeglm"), silent = TRUE)
    if (!inherits(sh, "try-error")) res <- sh
    return(res)
  }
  
  if (nm2 %in% nms) {
    res2 <- results(dds, name = nm2)
    # flip to get A_vs_B
    res2$log2FoldChange <- -res2$log2FoldChange
    if (!is.null(res2$stat)) res2$stat <- -res2$stat
    return(res2)
  }
  
  res_try <- try(results(dds, contrast = c(factor_name, A, B)), silent = TRUE)
  if (!inherits(res_try, "try-error")) return(res_try)
  
  return(NULL)
}

# --- Plot function with points + p-value brackets (DESeq2) ---
# Build one violin+box plot (per tissue) on VST-normalized ABCA1
# p_AC = p-value for ALS vs Control (numeric or NA)
# p_AO = p-value for ALS vs OND     (numeric or NA)
plot_abca1_vst_tissue <- function(vsd, pheno_t, tissue_label, p_AC = NA_real_, p_AO = NA_real_) {
  mat <- assay(vsd)
  if (!"ABCA1" %in% rownames(mat)) return(NULL)
  
  # Prepare plotting dataframe
  v <- as.numeric(mat["ABCA1", ])
  df <- data.frame(
    sample   = colnames(mat),
    vsd_expr = v
  ) %>%
    dplyr::left_join(pheno_t %>% dplyr::select(GSM, group_simple), by = c("sample" = "GSM"))
  
  # Order groups in plot: ALS (1st), Control (2nd), OND (3rd if present)
  desired_order <- c("ALS","Control","OND")
  present <- desired_order[desired_order %in% unique(as.character(df$group_simple))]
  df$group_simple <- factor(df$group_simple, levels = present)
  
  # Base plot: violin + white box + sample points (size = 1.26, alpha = 0.5)
  p <- ggplot(df, aes(x = group_simple, y = vsd_expr, fill = group_simple)) +
    geom_violin(trim = FALSE, color = "black", alpha = 0.7) +
    geom_boxplot(width = 0.18, fill = "white", color = "black",
                 alpha = 1, outlier.shape = NA, linewidth = 0.7) +
    geom_jitter(shape = 21, size = 1.26, stroke = 0.5, color = "black",
                width = 0.12, alpha = 0.5) +
    scale_fill_manual(values = group_colors, drop = FALSE) +
    labs(
      title = paste0("GSE124439 — ", tissue_label, ": ABCA1 (VST)"),
      x = NULL,
      y = "VST-normalized ABCA1 expression"
    ) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.title.y    = element_text(size = 12, face = "bold"),
      axis.text       = element_text(size = 11, face = "bold"),
      axis.line       = element_line(size = 0.6),
      axis.ticks      = element_line(size = 0.6),
      plot.title      = element_text(size = 14, face = "bold", hjust = 0.5)
    )
  
  # ---- P-value annotations (DESeq2) ----
  add_bracket <- function(pv, g1, g2, y, lbl_size = 3.8) {
    if (is.na(pv)) return(NULL)
    levs <- levels(df$group_simple)
    x1 <- match(g1, levs); x2 <- match(g2, levs)
    if (any(is.na(c(x1, x2)))) return(NULL)
    xmid <- (x1 + x2) / 2
    list(
      annotate("segment", x = x1, xend = x2, y = y, yend = y, linewidth = 0.7),
      annotate("text", x = xmid, y = y * 1.02,
               label = sprintf("p = %.3g (DESeq2)", pv), size = lbl_size)
    )
  }
  
  y_max <- max(df$vsd_expr, na.rm = TRUE)
  y1 <- y_max * 1.08
  y2 <- y_max * 1.16
  
  # ALS vs Control
  if (all(c("ALS","Control") %in% levels(df$group_simple))) {
    p <- p + add_bracket(p_AC, "Control", "ALS", y1)
  }
  # ALS vs OND
  if (all(c("ALS","OND") %in% levels(df$group_simple))) {
    yy <- if (all(c("ALS","Control") %in% levels(df$group_simple))) y2 else y1
    p  <- p + add_bracket(p_AO, "OND", "ALS", yy)
  }
  
  p
}

## =========================
## 6) Per-tissue DESeq2 and ABCA1 extraction
## =========================
tissues <- levels(pheno$tissue)

abca1_summary <- list()
per_tissue_plots <- list()

for (tissue in tissues) {
  idx <- pheno$tissue == tissue
  counts_t <- counts_mat[, idx, drop = FALSE]
  pheno_t  <- pheno[idx, , drop = FALSE]
  
  # Skip if less than 2 groups in this tissue
  present_groups <- sort(unique(as.character(pheno_t$group_simple)))
  if (length(present_groups) < 2) next
  
  # Drop unused levels and set a convenient reference
  pheno_t$group_simple <- droplevels(pheno_t$group_simple)
  if ("Control" %in% levels(pheno_t$group_simple)) {
    pheno_t$group_simple <- relevel(pheno_t$group_simple, ref = "Control")
  } else if ("OND" %in% levels(pheno_t$group_simple)) {
    pheno_t$group_simple <- relevel(pheno_t$group_simple, ref = "OND")
  } else {
    pheno_t$group_simple <- relevel(pheno_t$group_simple, ref = levels(pheno_t$group_simple)[1])
  }
  
  # Low-count filter
  keep_t <- keep_genes_subset(counts_t, min_count = 10, min_samples = 3)
  counts_t_f <- counts_t[keep_t, , drop = FALSE]
  
  # Build DESeq2
  col_df <- pheno_t %>% as.data.frame() %>% column_to_rownames("GSM")
  dds_t <- DESeqDataSetFromMatrix(countData = counts_t_f,
                                  colData   = col_df,
                                  design    = ~ group_simple)
  dds_t <- DESeq(dds_t)
  
  # Collect ABCA1 results
  res_list <- list()
  
  # ALS vs Control (if both exist)
  if (all(c("ALS","Control") %in% levels(pheno_t$group_simple))) {
    res_AC <- safe_results_contrast(dds_t, factor_name = "group_simple", A = "ALS", B = "Control")
    if (!is.null(res_AC)) {
      df <- normalize_res_df(res_AC)
      row <- df %>% dplyr::filter(toupper(gene) == "ABCA1")
      if (nrow(row) == 1) res_list[["ALS_vs_Control"]] <- row
    }
  }
  
  # ALS vs OND (if both exist)
  if (all(c("ALS","OND") %in% levels(pheno_t$group_simple))) {
    res_AO <- safe_results_contrast(dds_t, factor_name = "group_simple", A = "ALS", B = "OND")
    if (!is.null(res_AO)) {
      df <- normalize_res_df(res_AO)
      row <- df %>% dplyr::filter(toupper(gene) == "ABCA1")
      if (nrow(row) == 1) res_list[["ALS_vs_OND"]] <- row
    }
  }
  
  # Save per-tissue ABCA1 results CSV
  if (length(res_list) > 0) {
    out_tbl <- dplyr::bind_rows(res_list, .id = "contrast")
    wanted  <- c("contrast","gene","log2FoldChange","lfcSE","stat","pvalue","padj")
    out_tbl <- out_tbl[, intersect(wanted, names(out_tbl)), drop = FALSE]
    readr::write_csv(out_tbl, file.path(out_dir, paste0("ABCA1_DESeq2_", safe_tissue(tissue), ".csv")))
    abca1_summary[[tissue]] <- out_tbl %>% dplyr::mutate(tissue = tissue, .before = 1)
  }
  
  # Extract p-values for plot annotations
  p_AC <- if (!is.null(res_list[["ALS_vs_Control"]])) res_list[["ALS_vs_Control"]]$pvalue[1] else NA_real_
  p_AO <- if (!is.null(res_list[["ALS_vs_OND"]]))     res_list[["ALS_vs_OND"]]$pvalue[1]     else NA_real_
  
  # VST visualization with points and p-value brackets (DESeq2)
  vsd_t <- vst(dds_t, blind = FALSE)
  p_t   <- plot_abca1_vst_tissue(vsd_t, pheno_t, tissue_label = tissue, p_AC = p_AC, p_AO = p_AO)
  if (!is.null(p_t)) {
    per_tissue_plots[[tissue]] <- p_t
    ggsave(filename = file.path(out_dir, paste0("ABCA1_VST_", safe_tissue(tissue), ".pdf")),
           plot = p_t, width = 6.5, height = 5.0, device = cairo_pdf)
  }
}

## =========================
## 7) Write overall summary and a faceted PDF
## =========================
if (length(abca1_summary) > 0) {
  summary_tbl <- dplyr::bind_rows(abca1_summary) %>%
    dplyr::relocate(tissue, .before = 1)
  readr::write_csv(summary_tbl, file.path(out_dir, "ABCA1_DESeq2_summary_all_tissues.csv"))
} else {
  message("No ABCA1 rows found in DE results (check gene symbol case).")
}

# Optional faceted overview PDF (VST), skip if nothing collected
if (length(per_tissue_plots) > 0) {
  facet_df_list <- list()
  for (tissue in names(per_tissue_plots)) {
    idx <- pheno$tissue == tissue
    counts_t <- counts_mat[, idx, drop = FALSE]
    pheno_t  <- pheno[idx, , drop = FALSE]
    
    pheno_t$group_simple <- droplevels(pheno_t$group_simple)
    if ("Control" %in% levels(pheno_t$group_simple)) {
      pheno_t$group_simple <- relevel(pheno_t$group_simple, ref = "Control")
    } else if ("OND" %in% levels(pheno_t$group_simple)) {
      pheno_t$group_simple <- relevel(pheno_t$group_simple, ref = "OND")
    }
    
    keep_t <- keep_genes_subset(counts_t, min_count = 10, min_samples = 3)
    counts_t_f <- counts_t[keep_t, , drop = FALSE]
    col_df <- pheno_t %>% as.data.frame() %>% column_to_rownames("GSM")
    dds_t <- DESeqDataSetFromMatrix(countData = counts_t_f,
                                    colData   = col_df,
                                    design    = ~ group_simple)
    dds_t <- DESeq(dds_t)
    vsd_t <- vst(dds_t, blind = FALSE)
    
    mat <- assay(vsd_t)
    if (!"ABCA1" %in% rownames(mat)) next
    v <- as.numeric(mat["ABCA1", ])
    df <- data.frame(
      sample   = colnames(mat),
      vsd_expr = v
    ) %>%
      dplyr::left_join(pheno_t %>% dplyr::select(GSM, group_simple), by = c("sample" = "GSM")) %>%
      dplyr::mutate(tissue = tissue)
    facet_df_list[[tissue]] <- df
  }
  
  if (length(facet_df_list) > 0) {
    facet_df <- dplyr::bind_rows(facet_df_list)
    
    # Same left-to-right order in facets would vary; fine to keep default here
    p_facet <- ggplot(facet_df, aes(x = group_simple, y = vsd_expr, fill = group_simple)) +
      geom_violin(trim = FALSE, color = "black", alpha = 0.7) +
      geom_boxplot(width = 0.15, fill = "white", color = "black",
                   alpha = 1, outlier.shape = NA, linewidth = 0.7) +
      geom_jitter(shape = 21, size = 0.9, stroke = 0.4, color = "black",
                  width = 0.12, alpha = 0.4) +
      scale_fill_manual(values = group_colors, drop = FALSE) +
      facet_wrap(~ tissue, scales = "free_y") +
      labs(
        title = "GSE124439 — ABCA1 (VST) across tissues",
        x = NULL,
        y = "VST-normalized ABCA1 expression"
      ) +
      theme_classic() +
      theme(
        legend.position = "none",
        axis.title.y    = element_text(size = 12, face = "bold"),
        axis.text       = element_text(size = 10, face = "bold"),
        strip.text      = element_text(size = 10, face = "bold"),
        plot.title      = element_text(size = 14, face = "bold", hjust = 0.5)
      )
    
    ggsave(filename = file.path(out_dir, "ABCA1_VST_faceted_all_tissues.pdf"),
           plot = p_facet, width = 10.5, height = 7.5, device = cairo_pdf)
  }
}

message("Done. Outputs saved under: ", normalizePath(out_dir))
rm(list=ls())











#############################################
## GSE120374 — mouse spinal cord bulk RNA-seq
## Aim: G93A (ALS) vs WT within each age (p30/p70/p100/p120)
## Gene of interest: Abca1
#############################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(readr)
  library(stringr)
})

## =========================
## 1) I/O
## =========================
counts_path <- "GSE_expr_and_phen/GSE120374_counts.csv"
pheno_path  <- "GSE_expr_and_phen/GSE120374_pheno.csv"
out_dir     <- "GSE120374_Abca1_output"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

## =========================
## 2) Load
## =========================
counts_raw <- read_csv(counts_path, show_col_types = FALSE)
pheno_raw  <- read_csv(pheno_path,  show_col_types = FALSE)

## =========================
## 3) Build counts matrix (clean column names / gene column)
## =========================
colnames(counts_raw)[1] <- "gene"
colnames(counts_raw) <- str_trim(colnames(counts_raw))

counts_mat <- counts_raw %>%
  column_to_rownames("gene") %>%
  as.matrix()
storage.mode(counts_mat) <- "integer"

## =========================
## 4) Parse GSM + suffix (WT/G93A/UNK) from column names and merge with pheno
## =========================
# Extract GSM from column names; suffix allows WT/G93A/UNK (case-insensitive); if missing, keep as NA
key <- tibble(sample_col = colnames(counts_mat)) %>%
  transmute(
    sample_col,
    GSM    = str_extract(sample_col, "GSM\\d+"),
    suffix = str_to_upper(str_extract(sample_col, "(WT|G93A|UNK)$"))
  ) %>% 
  filter(!is.na(GSM))

# Normalize age to p30/p70/p100/p120
normalize_age <- function(x) {
  d <- str_extract(as.character(x), "\\d+")
  ifelse(is.na(d), NA_character_, paste0("p", d))
}

# Keep only G93A and WT (remove ATG7* and mixed, etc.)
pheno_f <- pheno_raw %>%
  mutate(
    genotype = as.character(genotype),
    keep_geno = case_when(
      str_detect(genotype, fixed("B6SJLSOD1-G93A", ignore_case = TRUE)) ~ TRUE,
      str_detect(genotype, fixed("B6SJLSOD1-WT",   ignore_case = TRUE)) ~ TRUE,
      TRUE ~ FALSE
    ),
    age = normalize_age(age),
    mouse_id = as.character(mouse_id)
  ) %>%
  filter(keep_geno)

# Diagnostic: check sample sizes per age × genotype (from pheno)
message("\n[PHENO overview: age × genotype]")
print(table(pheno_f$age, pheno_f$genotype))

# Prioritize group by column suffix; if UNK/missing, fall back to genotype in pheno
sample_df <- key %>%
  left_join(pheno_f, by = "GSM") %>%
  mutate(
    suffix = ifelse(is.na(suffix), "UNK", suffix),
    group_simple = case_when(
      suffix == "G93A" ~ "ALS",
      suffix == "WT"   ~ "Control",
      # UNK: fall back to genotype in pheno
      str_detect(genotype, fixed("B6SJLSOD1-G93A", ignore_case = TRUE)) ~ "ALS",
      str_detect(genotype, fixed("B6SJLSOD1-WT",   ignore_case = TRUE)) ~ "Control",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(group_simple), !is.na(age), !is.na(mouse_id))

# Diagnostic: distribution of suffixes in count columns
message("\n[Suffix distribution in counts columns]")
print(table(sample_df$suffix, useNA = "ifany"))

# Keep only columns in counts that are matched to pheno
counts_mat <- counts_mat[, sample_df$sample_col, drop = FALSE]

# Diagnostic: sample counts per age × group before collapsing (from intersected count columns)
message("\n[Diagnostics BEFORE collapsing (per age × group)]")
print(table(sample_df$age, sample_df$group_simple))

## =========================
## 5) Collapse duplicates from the same mouse_id: average by gene × mouse_id × group × age
## =========================
SAFE_SEP <- "__SEP__"

expr_long <- counts_mat %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "sample_col", values_to = "count") %>%
  left_join(sample_df, by = "sample_col") %>%
  group_by(gene, mouse_id, group_simple, age) %>%
  summarise(count = round(mean(count)), .groups = "drop") %>%
  unite("sample_id", mouse_id, group_simple, age, sep = SAFE_SEP, remove = FALSE)

# Collapsed pheno (safely split back from sample_id)
pheno <- expr_long %>%
  distinct(sample_id) %>%
  tidyr::separate(sample_id,
                  into   = c("mouse_id","group_simple","age"),
                  sep    = SAFE_SEP,
                  remove = FALSE)

# Diagnostic: sample counts per age × group after collapsing
message("\n[Diagnostics AFTER collapsing (per age × group)]")
print(table(pheno$age, pheno$group_simple))

# Cast to wide format to get the collapsed expression matrix
counts_collapsed <- expr_long %>%
  select(gene, sample_id, count) %>%
  pivot_wider(names_from = sample_id, values_from = count) %>%
  column_to_rownames("gene") %>%
  as.matrix()
storage.mode(counts_collapsed) <- "integer"

## =========================
## 6) DESeq2: per age (ALS vs Control); and plot Abca1
## =========================
group_colors <- c("ALS" = "#9AC5CD", "Control" = "#FFDAB9")
ages_ordered <- c("p30","p70","p100","p120")
ages <- intersect(ages_ordered, sort(unique(pheno$age)))

abca1_results <- list()

for (ag in ages) {
  message("\nProcessing age = ", ag)
  pheno_ag  <- pheno %>% filter(age == ag)
  
  # At least two groups required
  if (length(unique(pheno_ag$group_simple)) < 2) {
    message("  Skip (only one group present).")
    next
  }
  
  counts_ag <- counts_collapsed[, pheno_ag$sample_id, drop = FALSE]
  
  # DE: Control as reference
  pheno_ag$group_simple <- factor(pheno_ag$group_simple, levels = c("Control","ALS"))
  if (nlevels(pheno_ag$group_simple) < 2) {
    message("  Skip (factor dropped).")
    next
  }
  
  col_df <- pheno_ag %>% as.data.frame() %>% column_to_rownames("sample_id")
  dds <- DESeqDataSetFromMatrix(countData = counts_ag, colData = col_df, design = ~ group_simple)
  dds <- DESeq(dds)
  
  res <- results(dds, contrast = c("group_simple","ALS","Control"))  # Wald test
  res_df <- as.data.frame(res) %>% tibble::rownames_to_column("gene")
  
  # Abca1 (mouse)
  ab_row <- res_df %>% filter(tolower(gene) == "abca1")
  if (nrow(ab_row) > 0) {
    abca1_results[[ag]] <- ab_row %>%
      select(gene, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj) %>%
      mutate(age = ag, .before = 1)
  }
  
  ## Visualization (VST) — ALS on the left, Control on the right; point size 1.26 (~70%)
  vsd <- vst(dds, blind = FALSE)
  mat <- assay(vsd)
  rn_lc <- tolower(rownames(mat))
  if (!"abca1" %in% rn_lc) {
    message("  Abca1 not found at ", ag)
    next
  }
  gname <- rownames(mat)[match("abca1", rn_lc)]
  v <- as.numeric(mat[gname, ])
  
  dfp <- tibble(
    sample   = colnames(mat),
    vsd_expr = v
  ) %>%
    left_join(pheno_ag, by = c("sample" = "sample_id"))
  
  # Plot order: ALS on the left, Control on the right
  dfp$group_simple <- factor(dfp$group_simple, levels = c("ALS","Control"))
  
  pval_txt <- if (nrow(ab_row) > 0 && !is.na(ab_row$pvalue[1])) sprintf("p = %.3g (DESeq2 Wald)", ab_row$pvalue[1]) else NA
  
  p <- ggplot(dfp, aes(x = group_simple, y = vsd_expr, fill = group_simple)) +
    geom_violin(trim = FALSE, color = "black", alpha = 0.7) +
    geom_boxplot(width = 0.18, fill = "white", color = "black",
                 alpha = 1, outlier.shape = NA, linewidth = 0.7) +
    geom_jitter(shape = 21, size = 1.26, stroke = 0.4, color = "black",
                width = 0.12, alpha = 0.5) +
    scale_fill_manual(values = group_colors, drop = FALSE) +
    labs(
      title = paste0("GSE120374 — Abca1 (VST) at ", ag),
      x = NULL,
      y = "VST-normalized Abca1"
    ) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.title.y    = element_text(size = 12, face = "bold"),
      axis.text       = element_text(size = 11, face = "bold"),
      axis.line       = element_line(size = 0.6),
      axis.ticks      = element_line(size = 0.6),
      plot.title      = element_text(size = 14, face = "bold", hjust = 0.5)
    )
  
  # p-value line and label
  if (!is.na(pval_txt) && length(unique(dfp$group_simple)) == 2) {
    y_max <- max(dfp$vsd_expr, na.rm = TRUE)
    y_brk <- y_max * 1.08
    p <- p +
      annotate("segment", x = 1, xend = 2, y = y_brk, yend = y_brk, linewidth = 0.7) +
      annotate("text", x = 1.5, y = y_brk * 1.02, label = pval_txt, size = 3.8)
  }
  
  ggsave(file.path(out_dir, paste0("Abca1_VST_age_", ag, ".pdf")),
         plot = p, width = 6.5, height = 5.0, device = cairo_pdf)
}

## =========================
## 7) Export summary
## =========================
if (length(abca1_results) > 0) {
  bind_rows(abca1_results) %>%
    arrange(factor(age, levels = c("p30","p70","p100","p120"))) %>%
    write_csv(file.path(out_dir, "Abca1_DESeq2_per_age.csv"))
}

message("\nDone. Outputs saved under: ", normalizePath(out_dir), "\n")




#############################################
## GSE120374 — mouse spinal cord bulk RNA-seq
## Aim: G93A (ALS) vs WT within age (p30/p70/p100/p120)
## Gene of interest: Abca1
## Test: DESeq2 Wald (per-age two-group model)
#############################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(readr)
  library(stringr)
})

## =========================
## 1) I/O
## =========================
counts_path <- "GSE_expr_and_phen/GSE120374_counts.csv"  # first column = gene name
pheno_path  <- "GSE_expr_and_phen/GSE120374_pheno.csv"   # contains GSM, mouse_id, age, genotype, etc.
out_dir     <- "GSE120374_Abca1_output"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

## =========================
## 2) Load
## =========================
counts_raw <- read_csv(counts_path, show_col_types = FALSE)
pheno_raw  <- read_csv(pheno_path,  show_col_types = FALSE)

## =========================
## 3) Build counts matrix
## =========================
colnames(counts_raw)[1] <- "gene"
counts_mat <- counts_raw %>%
  column_to_rownames("gene") %>%
  as.matrix()
storage.mode(counts_mat) <- "integer"

## =========================
## 4) Parse column names (GSM and suffix); combine with pheno to determine group & age
## =========================
# Column names may be GSMxxxx_WT / GSMxxxx_G93A / GSMxxxx_UNK
key <- tibble(sample_col = colnames(counts_mat)) %>%
  mutate(
    GSM    = str_extract(sample_col, "^GSM\\d+"),
    suffix = str_extract(sample_col, "(WT|G93A|UNK)$")
  ) %>%
  filter(!is.na(GSM), !is.na(suffix))

# Normalize age format
normalize_age <- function(x) {
  d <- str_extract(as.character(x), "\\d+")
  ifelse(is.na(d), NA_character_, paste0("p", d))
}

# Keep only WT and G93A (remove ATG7* and mixed)
pheno_f <- pheno_raw %>%
  mutate(
    genotype = as.character(genotype),
    keep_geno = case_when(
      str_detect(genotype, "B6SJLSOD1-G93A") ~ TRUE,
      str_detect(genotype, "B6SJLSOD1-WT")   ~ TRUE,
      TRUE ~ FALSE
    )
  ) %>% filter(keep_geno) %>%
  mutate(age = normalize_age(age))

# Group assignment:
# - If column suffix is WT/G93A, use directly
# - If suffix is UNK, fallback to pheno genotype
sample_df <- key %>%
  left_join(pheno_f, by = "GSM") %>%
  mutate(
    group_from_suffix = case_when(
      suffix == "G93A" ~ "ALS",
      suffix == "WT"   ~ "Control",
      TRUE             ~ NA_character_
    ),
    group_from_pheno  = case_when(
      str_detect(genotype, "G93A") ~ "ALS",
      str_detect(genotype, "WT")   ~ "Control",
      TRUE                         ~ NA_character_
    ),
    group_simple = coalesce(group_from_suffix, group_from_pheno)
  ) %>%
  transmute(
    sample_col,
    GSM,
    mouse_id = as.character(mouse_id),
    age      = as.character(age),
    group_simple
  ) %>%
  filter(!is.na(mouse_id), !is.na(age), group_simple %in% c("ALS","Control"))

# Keep only columns matched to phenotype
counts_mat <- counts_mat[, sample_df$sample_col, drop = FALSE]

message("\n[Diagnostics BEFORE collapsing (per age × group)]")
print(table(sample_df$age, sample_df$group_simple))

## =========================
## 5) Collapse replicates of the same mouse_id × age × group
## =========================
SAFE_SEP <- "__SEP__"

expr_long <- counts_mat %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "sample_col", values_to = "count") %>%
  left_join(sample_df, by = "sample_col") %>%
  group_by(gene, mouse_id, group_simple, age) %>%
  summarise(count = round(mean(count)), .groups = "drop") %>%
  unite("sample_id", mouse_id, group_simple, age, sep = SAFE_SEP, remove = FALSE)

# Collapsed phenotype
pheno <- expr_long %>%
  distinct(sample_id) %>%
  tidyr::separate(
    sample_id,
    into   = c("mouse_id","group_simple","age"),
    sep    = SAFE_SEP,
    remove = FALSE
  )

message("\n[Diagnostics AFTER  collapsing (per age × group)]")
print(table(pheno$age, pheno$group_simple))

# Collapsed expression matrix
counts_collapsed <- expr_long %>%
  select(gene, sample_id, count) %>%
  pivot_wider(names_from = sample_id, values_from = count) %>%
  column_to_rownames("gene") %>%
  as.matrix()
storage.mode(counts_collapsed) <- "integer"

## =========================
## 6) DESeq2: per-age Wald test; collect Abca1 p-values & VST values
## =========================
group_colors <- c("ALS" = "#9AC5CD", "Control" = "#FFDAB9")
ages_ordered <- c("p30","p70","p100","p120")  # fixed order
ages <- intersect(ages_ordered, sort(unique(pheno$age)))

abca1_results <- list()
plot_df_list  <- list()
ann_df_list   <- list()

for (ag in ages) {
  message("\nProcessing age = ", ag)
  pheno_ag  <- pheno %>% filter(age == ag)
  
  # Require at least two groups
  if (length(unique(pheno_ag$group_simple)) < 2) {
    message("  Skip (only one group).")
    next
  }
  counts_ag <- counts_collapsed[, pheno_ag$sample_id, drop = FALSE]
  
  # DE: Control as reference
  pheno_ag$group_simple <- factor(pheno_ag$group_simple, levels = c("Control","ALS"))
  if (nlevels(pheno_ag$group_simple) < 2) {
    message("  Skip (factor dropped).")
    next
  }
  
  col_df <- pheno_ag %>% as.data.frame() %>% column_to_rownames("sample_id")
  dds <- DESeqDataSetFromMatrix(countData = counts_ag, colData = col_df, design = ~ group_simple)
  dds <- DESeq(dds)  # Wald
  
  res <- results(dds, contrast = c("group_simple","ALS","Control"))  # Wald test
  res_df <- as.data.frame(res) %>% tibble::rownames_to_column("gene")
  
  # Abca1 (case-insensitive)
  ab_row <- res_df %>% filter(tolower(gene) == "abca1")
  if (nrow(ab_row) > 0) {
    abca1_results[[ag]] <- ab_row %>%
      select(gene, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj) %>%
      mutate(age = ag, .before = 1)
  }
  
  ## VST & extract Abca1
  vsd <- vst(dds, blind = FALSE)
  mat <- assay(vsd)
  rn_lc <- tolower(rownames(mat))
  if (!"abca1" %in% rn_lc) {
    message("  Abca1 not found at ", ag)
    next
  }
  gname <- rownames(mat)[match("abca1", rn_lc)]
  v <- as.numeric(mat[gname, ])
  
  # —— Fix: do not rownames_to_column("sample_id"), directly join by sample_id —— #
  dfp <- tibble(
    sample   = colnames(mat),
    vsd_expr = v
  ) %>%
    left_join(pheno_ag, by = c("sample" = "sample_id")) %>%   # ✅ fixed here
    mutate(
      age          = factor(age, levels = ages_ordered),
      group_simple = factor(group_simple, levels = c("ALS","Control"))  # ALS left, Control right within panel
    )
  
  plot_df_list[[ag]] <- dfp
  
  # Panel p-value (Wald)
  pval_txt <- if (nrow(ab_row) > 0 && !is.na(ab_row$pvalue[1])) sprintf("p = %.3g (DESeq2 Wald)", ab_row$pvalue[1]) else NA
  if (!is.na(pval_txt)) {
    y_max <- max(dfp$vsd_expr, na.rm = TRUE)
    ann_df_list[[ag]] <- tibble(
      age   = factor(ag, levels = ages_ordered),
      x1    = 1, x2 = 2,
      y     = y_max * 1.08,
      y_txt = y_max * 1.10,
      label = pval_txt
    )
  }
  
  # Optional: save individual age plots
  p_single <- ggplot(dfp, aes(x = group_simple, y = vsd_expr, fill = group_simple)) +
    geom_violin(trim = FALSE, color = "black", alpha = 0.7) +
    geom_boxplot(width = 0.18, fill = "white", color = "black",
                 alpha = 1, outlier.shape = NA, linewidth = 0.7) +
    geom_jitter(shape = 21, size = 1.26, stroke = 0.4, color = "black",
                width = 0.12, alpha = 0.5) +
    scale_fill_manual(values = group_colors, drop = FALSE) +
    labs(
      title = paste0("GSE120374 — Abca1 (VST) at ", ag),
      x = NULL,
      y = "VST-normalized Abca1"
    ) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.title.y    = element_text(size = 12, face = "bold"),
      axis.text       = element_text(size = 11, face = "bold"),
      axis.line       = element_line(size = 0.6),
      axis.ticks      = element_line(size = 0.6),
      plot.title      = element_text(size = 14, face = "bold", hjust = 0.5)
    )
  
  if (!is.na(pval_txt)) {
    p_single <- p_single +
      annotate("segment", x = 1, xend = 2, y = y_max * 1.08, yend = y_max * 1.08, linewidth = 0.7) +
      annotate("text", x = 1.5, y = y_max * 1.10, label = pval_txt, size = 3.8)
  }
  
  ggsave(file.path(out_dir, paste0("Abca1_VST_age_", ag, ".pdf")),
         plot = p_single, width = 6.5, height = 5.0, device = cairo_pdf)
}

## Save numerical results per age
if (length(abca1_results) > 0) {
  ab_sum <- bind_rows(abca1_results) %>%
    relocate(age, .before = 1)
  write_csv(ab_sum, file.path(out_dir, "Abca1_DESeq2_Wald_per_age.csv"))
}

## =========================
## 7) Combined plot: 4 ages together (facet in one row, y-axis fixed)
## =========================
if (length(plot_df_list) > 0) {
  plot_df <- bind_rows(plot_df_list)
  ann_df  <- if (length(ann_df_list) > 0) bind_rows(ann_df_list) else NULL
  
  p_all <- ggplot(plot_df, aes(x = group_simple, y = vsd_expr, fill = group_simple)) +
    geom_violin(trim = FALSE, color = "black", alpha = 0.7) +
    geom_boxplot(width = 0.18, fill = "white", color = "black",
                 alpha = 1, outlier.shape = NA, linewidth = 0.7) +
    geom_jitter(shape = 21, size = 1.26, stroke = 0.4, color = "black",
                width = 0.12, alpha = 0.5) +
    scale_fill_manual(values = group_colors, drop = FALSE) +
    facet_wrap(~ age, nrow = 1, scales = "fixed") +  # y-axis fixed
    labs(
      title = "GSE120374 — Abca1 (VST) across ages",
      x = NULL,
      y = "VST-normalized Abca1"
    ) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.title.y    = element_text(size = 12, face = "bold"),
      axis.text       = element_text(size = 10, face = "bold"),
      strip.text      = element_text(size = 11, face = "bold"),
      plot.title      = element_text(size = 14, face = "bold", hjust = 0.5)
    )
  
  # Add p-value annotation per panel
  if (!is.null(ann_df) && nrow(ann_df) > 0) {
    p_all <- p_all +
      geom_segment(data = ann_df,
                   aes(x = x1, xend = x2, y = y, yend = y),
                   inherit.aes = FALSE, linewidth = 0.7) +
      geom_text(data = ann_df,
                aes(x = 1.5, y = y_txt, label = label),
                inherit.aes = FALSE, size = 3.8)
  }
  
  ggsave(file.path(out_dir, "Abca1_VST_all_ages_faceted.pdf"),
         plot = p_all, width = 12, height = 4.5, device = cairo_pdf)
}

message("\nDone. Outputs saved under: ", normalizePath(out_dir))






#############################################
## GSE272626 — human spinal cord bulk RNA-seq
## Focus: ALS vs CTL within tissue (cervical / lumbar)
## Gene of interest: ABCA1
## Test: DESeq2 Wald (two-group, per tissue)
#############################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(readr)
  library(stringr)
})

## =========================
## 0) Paths & output dir
## =========================
counts_path <- "GSE_expr_and_phen/GSE272626_counts.csv"  # 1st column = gene, others = GSM
pheno_path  <- "GSE_expr_and_phen/GSE272626_pheno.csv"   # contains GSM, disease_state, tissue
out_dir     <- "GSE272626_ABCA1_output"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

## =========================
## 1) Load
## =========================
counts_raw <- read_csv(counts_path, show_col_types = FALSE)
pheno_raw  <- read_csv(pheno_path,  show_col_types = FALSE)

## =========================
## 2) Build counts matrix
##    - collapse duplicated genes by row-wise max
## =========================
stopifnot(ncol(counts_raw) >= 2)
colnames(counts_raw)[1] <- "gene"

counts_collapsed <- counts_raw %>%
  group_by(gene) %>%
  summarise(across(where(is.numeric), ~ max(.x, na.rm = TRUE)), .groups = "drop")

counts_mat <- counts_collapsed %>%
  column_to_rownames("gene") %>%
  as.matrix()
storage.mode(counts_mat) <- "integer"

## =========================
## 3) Prepare phenotype
##    - keep only ALS / CTL
##    - keep only cervical / lumbar spinal cord
##    - align by GSM
## =========================
pheno <- pheno_raw %>%
  transmute(
    GSM = as.character(GSM),
    disease_state = as.character(disease_state),
    tissue = as.character(tissue)
  ) %>%
  filter(disease_state %in% c("ALS","CTL"),
         tissue %in% c("cervical spinal cord","lumbar spinal cord")) %>%
  mutate(group_simple = recode(disease_state, "ALS" = "ALS", "CTL" = "Control"))

common_gsm <- intersect(colnames(counts_mat), pheno$GSM)
counts_mat <- counts_mat[, common_gsm, drop = FALSE]
pheno      <- pheno[match(common_gsm, pheno$GSM), , drop = FALSE]

# Factor setting: Control as reference in DE analysis; ALS on the left in plots
pheno$group_simple <- factor(pheno$group_simple, levels = c("Control","ALS"))
pheno$tissue       <- droplevels(factor(pheno$tissue))

## quick checks
message("\n[Sample counts per tissue × group]")
print(with(pheno, table(tissue, group_simple)))

## =========================
## 4) Utilities
## =========================
keep_genes_subset <- function(counts_subset, min_count = 10, min_samples = 3) {
  rowSums(counts_subset >= min_count) >= min_samples
}
safe_tissue <- function(x) {
  x %>% str_replace_all("[^A-Za-z0-9]+", "_") %>% str_replace("^_+|_+$", "")
}
group_colors <- c("ALS" = "#9AC5CD", "Control" = "#FFDAB9")

## =========================
## 5) Per-tissue DE & plotting (ABCA1)
## =========================
tissues <- levels(pheno$tissue)
abca1_results <- list()
plot_list     <- list()
ann_list      <- list()   # p-value annotation data

for (tissue in tissues) {
  message("\nProcessing tissue = ", tissue)
  idx       <- pheno$tissue == tissue
  pheno_t   <- droplevels(pheno[idx, , drop = FALSE])
  counts_t  <- counts_mat[, pheno_t$GSM, drop = FALSE]
  
  # At least two groups required
  if (length(unique(pheno_t$group_simple)) < 2) {
    message("  Skip: only one group present.")
    next
  }
  
  # Low-expression filtering (within this tissue)
  keep_t <- keep_genes_subset(counts_t, min_count = 10, min_samples = 3)
  counts_t <- counts_t[keep_t, , drop = FALSE]
  
  # DESeq2 (Control as reference)
  col_df <- pheno_t %>% as.data.frame() %>% column_to_rownames("GSM")
  dds <- DESeqDataSetFromMatrix(countData = counts_t,
                                colData   = col_df,
                                design    = ~ group_simple)
  dds <- DESeq(dds)  # Wald
  
  res <- results(dds, contrast = c("group_simple","ALS","Control"))  # Wald test
  res_df <- as.data.frame(res) %>% tibble::rownames_to_column("gene")
  
  # ABCA1 (case-insensitive)
  ab_row <- res_df %>% filter(tolower(gene) == "abca1")
  if (nrow(ab_row) > 0) {
    abca1_results[[tissue]] <- ab_row %>%
      select(gene, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj) %>%
      mutate(tissue = tissue, .before = 1)
  }
  
  # VST & plot
  vsd <- vst(dds, blind = FALSE)
  mat <- assay(vsd)
  rn_lc <- tolower(rownames(mat))
  if (!"abca1" %in% rn_lc) {
    message("  ABCA1 not found in VST matrix (skipped plot).")
    next
  }
  gname <- rownames(mat)[match("abca1", rn_lc)]
  v <- as.numeric(mat[gname, ])
  
  dfp <- tibble(
    sample   = colnames(mat),
    vsd_expr = v
  ) %>%
    left_join(pheno_t, by = c("sample" = "GSM")) %>%
    mutate(
      # For plotting: ALS on left, Control on right
      group_simple = factor(ifelse(group_simple == "ALS", "ALS", "Control"),
                            levels = c("ALS","Control")),
      tissue = tissue
    )
  
  pval_txt <- if (nrow(ab_row) > 0 && !is.na(ab_row$pvalue[1]))
    sprintf("p = %.3g (DESeq2 Wald)", ab_row$pvalue[1]) else NA
  
  p <- ggplot(dfp, aes(x = group_simple, y = vsd_expr, fill = group_simple)) +
    geom_violin(trim = FALSE, color = "black", alpha = 0.7) +
    geom_boxplot(width = 0.18, fill = "white", color = "black",
                 alpha = 1, outlier.shape = NA, linewidth = 0.7) +
    geom_jitter(shape = 21, size = 1.26, stroke = 0.4, color = "black",  # point size = 70% of before
                width = 0.12, alpha = 0.5) +
    scale_fill_manual(values = group_colors, drop = FALSE) +
    labs(
      title = paste0("GSE272626 — ", tissue, ": ABCA1 (VST)"),
      x = NULL,
      y = "VST-normalized ABCA1"
    ) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.title.y    = element_text(size = 12, face = "bold"),
      axis.text       = element_text(size = 11, face = "bold"),
      axis.line       = element_line(size = 0.6),
      axis.ticks      = element_line(size = 0.6),
      plot.title      = element_text(size = 14, face = "bold", hjust = 0.5)
    )
  
  # p-value line annotation
  if (!is.na(pval_txt) && length(unique(dfp$group_simple)) == 2) {
    y_max <- max(dfp$vsd_expr, na.rm = TRUE)
    p <- p +
      annotate("segment", x = 1, xend = 2, y = y_max * 1.08, yend = y_max * 1.08, linewidth = 0.7) +
      annotate("text", x = 1.5, y = y_max * 1.10, label = pval_txt, size = 3.8)
    ann_list[[tissue]] <- tibble(
      tissue = tissue, x1 = 1, x2 = 2,
      y = y_max * 1.08, y_txt = y_max * 1.10,
      label = pval_txt
    )
  }
  
  plot_list[[tissue]] <- dfp
  
  # Save single-tissue PDF
  ggsave(
    filename = file.path(out_dir, paste0("ABCA1_VST_", safe_tissue(tissue), ".pdf")),
    plot = p, width = 6.5, height = 5.0, device = cairo_pdf
  )
}

## =========================
## 6) Save summary & combined facet PDF
## =========================
# Summary CSV
if (length(abca1_results) > 0) {
  summary_tbl <- bind_rows(abca1_results) %>%
    relocate(tissue, .before = 1)
  write_csv(summary_tbl, file.path(out_dir, "ABCA1_DESeq2_Wald_per_tissue.csv"))
}

# Combined plot: two tissues in the same figure, y-axis consistent
if (length(plot_list) > 0) {
  facet_df <- bind_rows(plot_list)
  facet_df$tissue <- factor(facet_df$tissue, levels = c("cervical spinal cord","lumbar spinal cord"))
  
  ann_df <- if (length(ann_list) > 0) bind_rows(ann_list) else NULL
  
  p_facet <- ggplot(facet_df, aes(x = group_simple, y = vsd_expr, fill = group_simple)) +
    geom_violin(trim = FALSE, color = "black", alpha = 0.7) +
    geom_boxplot(width = 0.18, fill = "white", color = "black",
                 alpha = 1, outlier.shape = NA, linewidth = 0.7) +
    geom_jitter(shape = 21, size = 1.26, stroke = 0.4, color = "black",
                width = 0.12, alpha = 0.5) +
    scale_fill_manual(values = group_colors, drop = FALSE) +
    facet_wrap(~ tissue, nrow = 1, scales = "fixed") +
    labs(
      title = "GSE272626 — ABCA1 (VST) in cervical vs lumbar spinal cord",
      x = NULL,
      y = "VST-normalized ABCA1"
    ) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.title.y    = element_text(size = 12, face = "bold"),
      axis.text       = element_text(size = 10, face = "bold"),
      strip.text      = element_text(size = 11, face = "bold"),
      plot.title      = element_text(size = 14, face = "bold", hjust = 0.5)
    )
  
  if (!is.null(ann_df) && nrow(ann_df) > 0) {
    p_facet <- p_facet +
      geom_segment(data = ann_df, aes(x = x1, xend = x2, y = y, yend = y),
                   inherit.aes = FALSE, linewidth = 0.7) +
      geom_text(data = ann_df, aes(x = 1.5, y = y_txt, label = label),
                inherit.aes = FALSE, size = 3.8)
  }
  
  ggsave(
    filename = file.path(out_dir, "ABCA1_VST_cervical_lumbar_faceted.pdf"),
    plot = p_facet, width = 12, height = 4.5, device = cairo_pdf
  )
}

message("\nDone. Outputs saved under: ", normalizePath(out_dir))






