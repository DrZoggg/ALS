# ---------------------------
# Load packages
# ---------------------------
options(stringsAsFactors = FALSE)
library(data.table)
library(tibble)
library(dplyr)
library(GEOquery)
library(clusterProfiler)
library(org.Hs.eg.db)

# ---------------------------
# Discovery / Validation cohorts (GEO)
# ---------------------------
disc_cohort  <- getGEO(filename = "Raw_data/GSE112681/GSE112681-GPL6947_series_matrix.txt.gz", getGPL = FALSE)
valid_cohort <- getGEO(filename = "Raw_data/GSE112681/GSE112681-GPL10558_series_matrix.txt.gz", getGPL = FALSE)
str(disc_cohort)

# ---------------------------
# Load expression matrix (discovery cohort)
# ---------------------------
expr <- exprs(disc_cohort) %>% as.data.frame()
expr <- expr[apply(expr, 1, function(x) sum(x != 0) > 0), ]  # Remove genes with all-zero expression across samples

# Annotation table (contains Symbol, ID, Entrez_Gene_ID for mapping)
ann <- fread("F:/坚果云/文件/MND生信分析/Raw_data/Gene_annotation_platform/GPL6947-13512.txt", data.table = FALSE)
ann <- ann[, c("ID", "Symbol", "Entrez_Gene_ID", "Definition")]

# Merge annotation and expression (drop the first ID column after merge)
expr_all <- merge(ann, expr, by.x = 1, by.y = 0)[, -1]
print(paste0("Raw row count: ", nrow(expr_all)))

# Collapse duplicates by gene Symbol using the columnwise maximum; keep non-empty symbols
expr_all <- expr_all %>%
  group_by(Symbol) %>%
  summarise_all(max) %>%
  filter(Symbol != "") %>%
  mutate(SYMBOL = Symbol) %>%            # Duplicate Symbol into a 'SYMBOL' column
  column_to_rownames("Symbol") %>%       # Set Symbol as rownames
  relocate(SYMBOL)                        # Move SYMBOL to be the first column

print(paste0("Filtered row count: ", nrow(expr_all)))

# ---------------------------
# Map Entrez IDs to Ensembl IDs
# ---------------------------
ID <- expr_all$Entrez_Gene_ID
ID_mapping <- bitr(
  ID,
  fromType = "ENTREZID",
  toType   = "ENSEMBL",
  OrgDb    = org.Hs.eg.db
)

# Attach ENSEMBL to the table (matched by Entrez ID)
expr_all$ENSEMBL <- ID_mapping$ENSEMBL[match(expr_all$Entrez_Gene_ID, ID_mapping$ENTREZID)]

# Place ENSEMBL as the 4th column
expr_all <- expr_all[, c(1:3, ncol(expr_all), 4:(ncol(expr_all) - 1))]

# Keep a small mapping frame (SYMBOL/Entrez/Definition/ENSEMBL)
head(expr_all, 3)
mapping <- expr_all[, 1:4]

# ---------------------------
# Clinical data (discovery cohort)
# ---------------------------
clin <- pData(disc_cohort)

# Inspect processing descriptions (for record)
print("Data processing details (discovery cohort):")
unique(clin$data_processing)
unique(clin$data_processing.1)
unique(clin$data_processing.2)
unique(clin$data_processing.3)
unique(clin$data_processing.4)
print("End of processing details.")

# Keep useful fields
clin <- clin[, c("diagnosis:ch1", "Sex:ch1", "site_onset:ch1", "survival_yr:ch1", "age_onset:ch1", "censored:ch1")]
colnames(clin) <- gsub(":ch1", "", colnames(clin))

# Add Sample column and reorder columns
clin$Sample <- rownames(clin)
clin <- clin[, c(7, 1:6)]

# Sort by diagnosis (descending)
clin <- arrange(clin, desc(diagnosis))
head(clin, 3)

# ---------------------------
# Normalize expression across arrays (discovery)
# ---------------------------
expr_all <- as.data.frame(limma::normalizeBetweenArrays(expr_all[, clin$Sample]))
head(expr_all, 3)

# Save discovery cohort objects
disc_cohort_clin    <- clin
disc_cohort_expr    <- expr_all
disc_cohort_mapping <- mapping

# ---------------------------
# Validation cohort
# ---------------------------
str(valid_cohort)

# Load expression matrix (validation cohort)
expr <- exprs(valid_cohort) %>% as.data.frame()
expr <- expr[apply(expr, 1, function(x) sum(x != 0) > 0), ]  # Remove genes with all-zero expression

# Annotation for the validation platform
ann <- fread("Raw_data/Gene_annotation_platform/GPL10558-50081.txt", data.table = FALSE)
ann <- ann[, c("ID", "Symbol", "Entrez_Gene_ID", "Definition")]

expr_all <- merge(ann, expr, by.x = 1, by.y = 0)[, -1]
print(paste0("Raw row count: ", nrow(expr_all)))

expr_all <- expr_all %>%
  group_by(Symbol) %>%
  summarise_all(max) %>%
  filter(Symbol != "") %>%
  mutate(SYMBOL = Symbol) %>%
  column_to_rownames("Symbol") %>%
  relocate(SYMBOL)

print(paste0("Filtered row count: ", nrow(expr_all)))

# Map Entrez -> Ensembl (validation)
ID <- expr_all$Entrez_Gene_ID
ID_mapping <- bitr(
  ID,
  fromType = "ENTREZID",
  toType   = "ENSEMBL",
  OrgDb    = org.Hs.eg.db
)
expr_all$ENSEMBL <- ID_mapping$ENSEMBL[match(expr_all$Entrez_Gene_ID, ID_mapping$ENTREZID)]
expr_all <- expr_all[, c(1:3, ncol(expr_all), 4:(ncol(expr_all) - 1))]

head(expr_all, 3)
mapping <- expr_all[, 1:4]

# Clinical data (validation)
clin <- pData(valid_cohort)

print("Data processing details (validation cohort):")
unique(clin$data_processing)
unique(clin$data_processing.1)
unique(clin$data_processing.2)
unique(clin$data_processing.3)
unique(clin$data_processing.4)
print("End of processing details.")

# Keep useful fields
clin <- clin[, c("diagnosis:ch1", "Sex:ch1", "site_onset:ch1", "survival_yr:ch1", "age_onset:ch1", "censored:ch1")]
colnames(clin) <- gsub(":ch1", "", colnames(clin))

clin$Sample <- rownames(clin)
clin <- clin[, c(7, 1:6)]
clin <- arrange(clin, desc(diagnosis))
head(clin, 3)

# Save validation cohort objects
valid_cohort_clin     <- clin
valid_cohort_expr     <- expr_all
valid_cohort_mapping  <- mapping

# ---------------------------
# Split discovery cohort into train/test (stratified by diagnosis, 70/30)
# ---------------------------
library(caret)

set.seed(1432)  # For reproducibility
train_index <- createDataPartition(disc_cohort_clin$diagnosis, p = 0.7, list = FALSE)

# Split clinical data
disc_cohort_clin_disc <- disc_cohort_clin[train_index, ]   # training (discovery split)
disc_cohort_clin_test <- disc_cohort_clin[-train_index, ]  # testing  (discovery split)

# Get sample IDs for train/test splits
train_samples <- disc_cohort_clin_disc$Sample
test_samples  <- disc_cohort_clin_test$Sample

# Subset expression matrices by train/test sample IDs
# NOTE: Ensure column names of 'disc_cohort_expr' match the values in 'Sample'.
train_cohort_expr <- disc_cohort_expr[, train_samples]
test_cohort_exp   <- disc_cohort_expr[, test_samples]   # (object name follows your original)

# Also subset the clinical frames (if rownames are sample IDs; otherwise match/merge as needed)
train_cohort_clin <- disc_cohort_clin[train_samples, ]
test_cohort_clin  <- disc_cohort_clin[test_samples, ]

# Inspect resulting dimensions
dim(train_cohort_expr)  # e.g., ~ 18971 rows × (≈ 0.7 × 741) columns
dim(test_cohort_exp)    # e.g., ~ 18971 rows × (≈ 0.3 × 741) columns

# ---------------------------
# Save intermediate objects
# ---------------------------
setwd("Data_preprocessing")
save(disc_cohort_clin, disc_cohort_expr, disc_cohort_mapping, file = "disc_cohort.rda")
save(train_cohort_expr, train_cohort_clin,                   file = "train_cohort.rda")
save(test_cohort_exp,  test_cohort_clin,                    file = "test_cohort.rda")
save(valid_cohort_clin, valid_cohort_expr, valid_cohort_mapping, file = "valid_cohort.rda")

# ---------------------------
# Keep shared genes between discovery and test split
# ---------------------------
interID <- Reduce(
  intersect,
  list(
    rownames(disc_cohort_expr),
    rownames(test_cohort_expr)  # NOTE: Your object above is 'test_cohort_exp'; change if intended.
  )
)
print(paste0("Total number of intersecting genes: ", length(interID)))

# Keep only mRNA genes based on 'Definition' (ends with 'mRNA')
disc_cohort_mapping <- disc_cohort_mapping[disc_cohort_mapping$SYMBOL %in% interID, ]
disc_cohort_mapping <- disc_cohort_mapping %>%
  filter(grepl("mRNA.$", Definition))  # retain rows whose Definition ends with "mRNA"

interID <- disc_cohort_mapping$SYMBOL
print(paste0("Number of intersecting mRNA genes: ", length(interID)))

# Subset expression matrix to mRNA gene set
my_mRNA <- disc_cohort_expr[interID, ]
print(paste0("Discovery cohort (mRNA) dimension: ", paste(dim(my_mRNA), collapse = " × ")))

# ---------------------------
# Prepare WGCNA input
# ---------------------------
wgcna_input <- my_mRNA
setwd("../WGCNA")
save(wgcna_input, file = "WGCNA_input.rda")
