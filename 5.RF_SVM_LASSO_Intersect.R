setwd("ML_features_select")

LASSO <- readRDS("./SVM&LASSO_features_filter/LASSO_features.rds")
SVM   <- readRDS("./SVM&LASSO_features_filter/SVM_features.rds")
RF    <- readRDS("RF_features_filter/features_RFmodel.rds")

# -------------------------------
# Venn diagram of LASSO / SVM / RF
# -------------------------------
library(ggvenn)
library(UpSetR)

biocolor <- c('#AFEEEE', "#FFE4C4", '#FF7256')
gene_sets <- list(LASSO = LASSO$lasso_fea,
                  SVM   = SVM,
                  RF    = RF$Feature)

options(repr.plot.width = 5, repr.plot.height = 5)
p1 <- ggvenn(
  gene_sets,
  fill_color    = biocolor,
  fill_alpha    = .70,
  stroke_linetype = "solid",
  stroke_size   = 0.2,
  text_size     = 2.5,
  set_name_size = 3
) +
  theme(
    # Add margins: top, right, bottom, left
    plot.margin = margin(10, 10, 10, 10)
  )
p1
ggsave(p1, file = 'Venn.features.pdf', height = 5, width = 5, dp = 1000)

# -------------------------------
# Get the intersection among LASSO, SVM and RF feature sets
# -------------------------------
gene_sets <- list(LASSO = LASSO$lasso_fea,
                  SVM   = SVM,
                  RF    = RF$Feature)

common_genes <- Reduce(intersect, gene_sets)

# Preview intersection
print(common_genes)
saveRDS(common_genes, file = 'common_genes.rds')

# -------------------------------
# Export sets to an Excel workbook
# -------------------------------
# install.packages("openxlsx")  # if needed
library(openxlsx)

gene_sets <- list(LASSO = LASSO$lasso_fea,
                  SVM   = SVM,
                  RF    = RF$Feature)

common_genes <- Reduce(intersect, gene_sets)

wb <- createWorkbook()

addWorksheet(wb, "Common_Genes")
writeData(wb, sheet = "Common_Genes", common_genes)

addWorksheet(wb, "LASSO_Genes")
writeData(wb, sheet = "LASSO_Genes", LASSO$lasso_fea)

addWorksheet(wb, "RF_Genes")
writeData(wb, sheet = "RF_Genes", RF$Feature)

addWorksheet(wb, "SVM_Genes")
writeData(wb, sheet = "SVM_Genes", SVM)

saveWorkbook(wb, "Gene_Sets.xlsx", overwrite = TRUE)
