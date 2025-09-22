library(stringr)
library(dplyr)
library(limma)
library(data.table)
library(tibble)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pheatmap)
library(ggplot2)

load("WGCNA/Step01-WGCNA_input.Rda")
load("Data_preprocessing/disc_cohort.rda")
ls()

disc_cohort_expr <- disc_cohort_expr[, rownames(datExpr)]
disc_cohort_clin <- disc_cohort_clin[rownames(datExpr), ]
str(disc_cohort_clin)

# Set groups
data_gene <- disc_cohort_expr
group <- factor(disc_cohort_clin$diagnosis)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
contrast <- makeContrasts('ALS - CON', levels = design)

############# Differential expression analysis #############
# data_gene: rows = genes, columns = samples
fit  <- lmFit(data_gene, design)
fit1 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit1)
alldiff <- topTable(fit2, 1, Inf)
alldiff <- tibble::rownames_to_column(alldiff, 'ID')
data_FC <- alldiff[, c(1, 2, 3, 5, 6)]

data_FC$FDR <- p.adjust(data_FC$P.Value, method = "BH")
data_FC <- data_FC[, c(1:4, 6)]
str(data_FC)

# Add a 'Threshold' column for DE status
data_FC$Threshold <- ifelse(
  data_FC$FDR < 0.05 & abs(data_FC$logFC) > 0.25,
  ifelse(data_FC$logFC > 0.25, "Up", "Down"),
  "No"
)

saveRDS(data_FC, file = 'DEA_data.rds')

data_FC_filtered <- data_FC %>% filter(P.Value < 0.05 & FDR < 0.05)

load("/WGCNA/WGCNA_Enrich.rda")
ls()

# Show intersection counts
up_gene_set   <- data_FC_filtered[data_FC_filtered$logFC > 0.25, ]$ID
down_gene_set <- data_FC_filtered[data_FC_filtered$logFC < -0.25, ]$ID
table(green_module %in% up_gene_set)
table(salmon_module %in% down_gene_set)

setwd("DEA&EA")
library(openxlsx)

# Create a new Excel workbook
wb <- createWorkbook()

Up   <- data_FC_filtered[which(data_FC_filtered$logFC >  0.25), ]  # Upregulated genes (logFC threshold adjustable)
Down <- data_FC_filtered[which(data_FC_filtered$logFC < -0.25), ]  # Downregulated genes

# Sheet 1: All
addWorksheet(wb, "All")
writeData(wb, sheet = "All", data_FC)

# Sheet 2: Up
addWorksheet(wb, "Up")
writeData(wb, sheet = "Up", Up)

# Sheet 3: Down
addWorksheet(wb, "Down")
writeData(wb, sheet = "Down", Down)

saveWorkbook(wb, file = 'disc_cohort_DEA.xlsx', overwrite = TRUE)

# Export text files for Cytoscape input
load("../WGCNA/WGCNA_Enrich.rda")
salmon_module <- data.frame(SYMBOL = salmon_module)
green_module  <- data.frame(SYMBOL = green_module)
write.table(salmon_module, file = 'Salmon_PPI_SYMBOL.txt', row.names = F, col.names = F, quote = F, sep = "\t")
write.table(green_module,  file = 'Green_PPI_SYMBOL.txt',  row.names = F, col.names = F, quote = F, sep = "\t")

# ---------------------- Venn plots ----------------------
library(ggvenn)
library(UpSetR)

biocolor <- c('#9BCD9B', "#CD5555")
x <- list(green_module,
          up_gene_set)
names(x) <- c('green module',
              'up gene set')

options(repr.plot.width = 4, repr.plot.height = 4)
p1 <- ggvenn(
  x,
  fill_color    = biocolor,
  fill_alpha    = 0.85,
  stroke_linetype = "solid",
  stroke_size   = 0.2,
  text_size     = 2.5,
  set_name_size = 3
) +
  theme(
    plot.margin = margin(10, 10, 10, 10)  # enlarge all four plot margins: top, right, bottom, left
  )
p1
ggsave(p1, file = 'Venn.上调模块差异基因.pdf', height = 4, width = 4, dp = 1000)

x <- list(
  salmon_module,
  down_gene_set
)
names(x) <- c('salmon_module', 'down gene set')

biocolor <- c('#EE9572', "#8DB6CD")
options(repr.plot.width = 4, repr.plot.height = 4)
p1 <- ggvenn(
  x,
  fill_color    = biocolor,
  fill_alpha    = .70,
  stroke_linetype = "solid",
  stroke_size   = 0.2,
  text_size     = 2.5,
  set_name_size = 3
) +
  theme(
    plot.margin = margin(10, 10, 10, 10)  # enlarge all four plot margins
  )
p1
ggsave(p1, file = 'Venn.下调模块差异基因.pdf', height = 4, width = 4, dp = 1000)

#### Prepare heatmap data ---
# Select upregulated genes: top 25 by FDR
upregulated_genes <- data_FC_filtered %>%
  dplyr::filter(logFC > 0) %>%
  arrange(FDR) %>%
  dplyr::slice_head(n = 25) %>%
  pull(ID)

# Select downregulated genes: top 25 by FDR
downregulated_genes <- data_FC_filtered %>%
  dplyr::filter(logFC < 0) %>%
  arrange(FDR) %>%
  dplyr::slice_head(n = 25) %>%
  pull(ID)

# Merge top up/down gene IDs (total 50)
genes_sig <- c(upregulated_genes, downregulated_genes)

# Subset expression matrix for the 50 genes
heatmap_df <- data_gene[genes_sig, ]

options(repr.plot.width = 7.5, repr.plot.height = 6)

#### Draw heatmap ---
Targets <- disc_cohort_clin %>% arrange(diagnosis)
annotation_col <- Targets %>% dplyr::select(diagnosis)  # choose the group column

library(pheatmap)
library(graphics)

heatmap_df1 <- as.data.frame(lapply(heatmap_df, as.numeric))
rownames(heatmap_df1) <- rownames(heatmap_df)

# Define annotation colors
annotation_colors <- list(
  diagnosis = c("ALS" = "#FB6A4A", "CON" = "#4292C6")
)

# Set breaks to map values to [-2.5, 2.5]
breaks <- seq(-2.5, 2.5, length.out = 101)

p1 <- pheatmap(
  as.matrix(heatmap_df1),
  scale               = "row",   # z-score per gene (row)
  border_color        = NA,
  main                = "disc cohort heatmap",
  annotation_col      = annotation_col,
  annotation_legend   = TRUE,
  annotation_colors   = annotation_colors,
  treeheight_row      = FALSE,
  cluster_cols        = FALSE,
  treeheight_col      = 30,
  fontsize_row        = 6,
  show_colnames       = FALSE,
  legend_breaks       = c(-2, 0, 2),
  legend_labels       = c(-2, 0, 2),
  breaks              = breaks
)

ggsave(p1, file = "disc_top25df_heatmap.pdf", width = 7.5, height = 6)

library(ggplot2)
library(dplyr)

# Data for volcano plot (keep FDR)
data <- data_FC %>%
  dplyr::select(ID, logFC, FDR)

# Volcano thresholds
logFC_cutoff <- 0.25
fdr_cutoff   <- 0.05

# Add DE status labels
data$Threshold <- factor(
  ifelse(data$FDR < fdr_cutoff & abs(data$logFC) >= logFC_cutoff,
         ifelse(data$logFC > logFC_cutoff, "Up", "Down"), "No"),
  levels = c("Down", "No", "Up")
)

range(data$logFC)  # inspect range to decide symmetrical x-limits if needed
data$gene <- data$ID

library(dplyr)
top_up <- data %>%
  filter(Threshold == "Up") %>%
  arrange(FDR) %>%
  slice_head(n = 10)

top_down <- data %>%
  filter(Threshold == "Down") %>%
  arrange(FDR) %>%
  slice_head(n = 10)

num_up   <- sum(data$Threshold == "Up")
num_down <- sum(data$Threshold == "Down")

library(ggrepel)

# Plot params
options(repr.plot.width = 5, repr.plot.height = 5)

# Volcano plot
p1 <- ggplot(data, aes(x = logFC, y = -log10(FDR), color = factor(Threshold, levels = c("Down", "No", "Up")))) +
  geom_point(alpha = 0.7, size = 0.85) +
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), lty = 2, col = "black", lwd = 0.4) +
  geom_hline(yintercept = -log10(fdr_cutoff),              lty = 2, col = "black", lwd = 0.4) +
  scale_color_manual(
    values = c("Down" = "#6CA6CD", "No" = "grey", "Up" = "#EE6363"),
    name   = "Expression"
  ) +
  guides(size = "none") +
  xlab("log2 Fold Change") +
  ylab("-log10 FDR") +
  ggtitle("Differential Expression in ALS vs Control") +
  geom_text_repel(data = top_up,   aes(label = gene), size = 3.2, color = "#EE6363", max.overlaps = Inf) +
  geom_text_repel(data = top_down, aes(label = gene), size = 3.2, color = "#6CA6CD", max.overlaps = Inf) +
  annotate("text", x =  1,  y = 48, label = paste0("Up: ",   num_up),   hjust = 1, size = 4.5, fontface = "bold", color = "#EE6363") +
  annotate("text", x = -1,  y = 48, label = paste0("Down: ", num_down), hjust = 0, size = 4.5, fontface = "bold", color = "#6CA6CD") +
  theme_bw() +
  theme(
    axis.title      = element_text(size = 12),
    axis.text       = element_text(size = 10),
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    legend.title    = element_text(size = 12),
    legend.text     = element_text(size = 10),
    legend.position = c(1, 1),
    legend.justification = c(1, 1),
    legend.background   = element_rect(fill = NULL, color = "black"),
    panel.grid.major    = element_blank(),
    panel.grid.minor    = element_blank()
  ) +
  xlim(-1.2, 1.2) + ylim(0, 50)

p1
ggsave(p1, file = "Volcano_DEG_with_Top10.pdf", width = 5, height = 5)


#########################################
############## Enrichment analysis

# Load files and packages
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)

DEA_data <- readRDS("/DEA_data.rds")
load("../WGCNA/WGCNA_Enrich.rda")
load("../Data_preprocessing/disc_cohort.rda")
ls()

# Map back to 'mapping' to get Entrez_Gene_ID
green_module  <- disc_cohort_mapping[green_module, ]
salmon_module <- disc_cohort_mapping[salmon_module, ]
str(green_module)

# Check presence of Entrez_Gene_ID
table(is.na(green_module$Entrez_Gene_ID))
table(is.na(salmon_module$Entrez_Gene_ID))

green_module_gene <- green_module$Entrez_Gene_ID
green_go <- enrichGO(
  gene          = green_module_gene,
  ont           = 'ALL',
  OrgDb         = org.Hs.eg.db,
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)
green_kk <- enrichKEGG(
  gene          = green_module_gene,
  organism      = "hsa",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)
head(green_kk)
save(green_go, green_kk, file = "GO_KEGG.green_enrich.rda")

Up <- DEA_data[DEA_data$Threshold == 'Up', ]
# Map back to 'mapping' to get Entrez_Gene_ID
up_gene <- disc_cohort_mapping[Up$ID, ]
str(up_gene)

library(dplyr)
up_gene  <- up_gene %>% filter(!is.na(Entrez_Gene_ID))
Pos_gene <- up_gene$Entrez_Gene_ID

Pos_go <- enrichGO(
  gene          = Pos_gene,
  ont           = 'ALL',
  OrgDb         = org.Hs.eg.db,
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)
Pos_kk <- enrichKEGG(
  gene          = Pos_gene,
  organism      = "hsa",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)
head(Pos_kk)

save(Pos_go, Pos_kk, file = "GO_KEGG.UP_enrich.rda")

Pos_go <- setReadable(Pos_go, org.Hs.eg.db, keyType = "ENTREZID")@result %>% filter(pvalue < 0.05 & qvalue < 0.05)
Pos_kk <- setReadable(Pos_kk, org.Hs.eg.db, keyType = "ENTREZID")@result %>% filter(pvalue < 0.05 & qvalue < 0.05)
Pos_kk <- Pos_kk %>% filter(!category == 'Human Diseases')  # do not display disease pathways

# Select top GO terms by |zScore|
go <- Pos_go[] %>%
  arrange(desc(abs(zScore))) %>%
  slice_head(n = 15) %>%
  group_by(ONTOLOGY)

# KEGG: select top by |zScore|
kegg <- Pos_kk[] %>%
  dplyr::arrange(desc(abs(zScore))) %>%
  dplyr::slice(1:5) %>%
  dplyr::select(ID:Count) %>%
  dplyr::mutate(ONTOLOGY = "KEGG") %>%
  dplyr::select(ONTOLOGY, everything())

Enrichement_results <- rbind(go[], kegg[])

# Keep top 20 by |zScore|
Enrichement_results <- Enrichement_results %>%
  dplyr::ungroup() %>%
  arrange(desc(abs(zScore))) %>%
  dplyr::slice(1:20)
str(Enrichement_results)

library(tidyverse)       # data wrangling & viz
library(clusterProfiler) # GO/KEGG enrichment
library(ggforce)         # plotting helpers
library(RColorBrewer)    # palettes
library(grid)            # grid graphics
library(ggh4x)           # ggplot2 extensions
library(ggfun)           # ggplot2 helper functions

# Combine GO and KEGG results
data <- Enrichement_results %>%
  dplyr::mutate(ONTOLOGY = factor(ONTOLOGY, levels = rev(c("BP", "CC", "KEGG")), ordered = TRUE)) %>%
  dplyr::arrange(ONTOLOGY, desc(zScore))

# Tidy geneID strings
data$geneID <- gsub("\\/", ", ", data$geneID)
data$Number <- data$Count / 10  # for point size aesthetics

# Named vector of |logFC| for labeling
logFC_abs <- abs(DEA_data$logFC)
names(logFC_abs) <- DEA_data$ID

# For each row, keep top-15 genes by |logFC|
data$geneID_top15 <- sapply(data$geneID, function(gene_str) {
  genes <- unlist(strsplit(gene_str, ",\\s*"))
  gene_scores <- logFC_abs[genes]
  sorted_genes <- names(sort(gene_scores, decreasing = TRUE, na.last = NA))
  top15 <- head(sorted_genes, 15)
  paste(top15, collapse = ", ")
})

options(repr.plot.width = 7.5, repr.plot.height = 8)

# Colors for bars
color <- rev(c('#FA8072', "#9AC0CD", "#bebbd7"))
p <- ggplot(data) +
  geom_bar(
    data = data,
    aes(x = zScore, y = interaction(Description, ONTOLOGY), fill = ONTOLOGY),
    stat = "identity"
  ) +
  scale_fill_manual(values = color, name = "GO&KEGG") +
  geom_text(
    aes(x = 0.1, y = interaction(Description, ONTOLOGY), label = Description, fontface = 'bold'),
    size = 3, hjust = 0, color = "black"
  ) +
  geom_text(
    aes(x = 0.1, y = interaction(Description, ONTOLOGY), label = geneID_top15),
    size = 2, hjust = 0, vjust = 2.5, color = "black"
  ) +
  geom_point(
    aes(x = -0.75, y = interaction(Description, ONTOLOGY), size = Count, fill = ONTOLOGY),
    shape = 21
  ) +
  geom_text(
    aes(x = -0.75, y = interaction(Description, ONTOLOGY), label = Count),
    size = 2.5, fontface = 'bold'
  ) +
  scale_size(range = c(4, 8), guide = guide_legend(override.aes = list(fill = "black"))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1)), limits = c(-1.5, 12)) +
  labs(x = "zScore", y = "Description") +
  theme(
    legend.title = element_text(face = "bold", color = "#000000", size = 12),
    legend.text  = element_text(color = "#000000", size = 10),
    ggh4x.axis.nestline.y = element_line(size = 3, color = color),
    ggh4x.axis.nesttext.y = element_text(face = "bold", color = color, hjust = 0.5, size = 12, angle = 90),
    legend.background = element_blank(),
    panel.background  = element_blank(),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    axis.title.y      = element_blank(),
    axis.text         = element_text(color = "#000000", size = 12),
    axis.text.y       = element_blank(),
    axis.ticks        = element_blank(),
    axis.title        = element_text(color = "#000000", size = 14),
    axis.ticks.y.left = element_line(color = NA)
  ) +
  scale_y_discrete(labels = function(y) str_wrap(y, width = 60)) +
  guides(
    y     = "axis_nested",
    y.sec = guide_axis_manual(breaks = 1:nrow(data), labels = data$Description),
    fill  = guide_legend(order = 1, reverse = TRUE),
    size  = guide_legend(order = 2, label = TRUE)
  ) +
  guides(fill = guide_legend(reverse = TRUE)) +
  # Add black border around the entire plot
  annotation_custom(grob = rectGrob(gp = gpar(col = "black", fill = NA, lwd = 2)),
                    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)

p
ggsave(p, file = 'UP_GO_KEGG_results.pdf', width = 7.5, height = 8)

green_go <- setReadable(green_go, org.Hs.eg.db, keyType = "ENTREZID")@result %>% filter(pvalue < 0.05 & qvalue < 0.05)
green_kk <- setReadable(green_kk, org.Hs.eg.db, keyType = "ENTREZID")@result %>% filter(pvalue < 0.05 & qvalue < 0.05)

# Select top GO terms by |zScore|
go <- green_go[] %>%
  arrange(desc(abs(zScore))) %>%
  slice_head(n = 15) %>%
  group_by(ONTOLOGY)

# KEGG: select top by |zScore|
kegg <- green_kk[] %>%
  dplyr::arrange(desc(abs(zScore))) %>%
  dplyr::slice(1:5) %>%
  dplyr::select(ID:Count) %>%
  dplyr::mutate(ONTOLOGY = "KEGG") %>%
  dplyr::select(ONTOLOGY, everything())

Enrichement_results <- rbind(go[], kegg[])

# Keep top 20 by |zScore|
Enrichement_results <- Enrichement_results %>%
  dplyr::ungroup() %>%
  arrange(desc(abs(zScore))) %>%
  dplyr::slice(1:20)
str(Enrichement_results)

library(tidyverse)       # data wrangling & viz
library(clusterProfiler) # GO/KEGG enrichment
library(ggforce)         # plotting helpers
library(RColorBrewer)    # palettes
library(grid)            # grid graphics
library(ggh4x)           # ggplot2 extensions
library(ggfun)           # ggplot2 helpers

# Combine GO and KEGG results
data <- Enrichement_results %>%
  dplyr::mutate(ONTOLOGY = factor(ONTOLOGY, levels = rev(c("BP", "CC", "MF", "KEGG")), ordered = TRUE)) %>%
  dplyr::arrange(ONTOLOGY, desc(zScore))

# Tidy geneID strings
data$geneID <- gsub("\\/", ", ", data$geneID)
data$Number <- data$Count / 10  # for point size aesthetics

# Named vector of |logFC| for labeling
logFC_abs <- abs(DEA_data$logFC)
names(logFC_abs) <- DEA_data$ID

# For each row, keep top-15 genes by |logFC|
data$geneID_top15 <- sapply(data$geneID, function(gene_str) {
  genes <- unlist(strsplit(gene_str, ",\\s*"))
  gene_scores <- logFC_abs[genes]
  sorted_genes <- names(sort(gene_scores, decreasing = TRUE, na.last = NA))
  top15 <- head(sorted_genes, 15)
  paste(top15, collapse = ", ")
})

options(repr.plot.width = 7.5, repr.plot.height = 8)

# Colors for bars
color <- rev(c('#FA8072', "#9AC0CD", "#CDB79E", "#bebbd7"))
p <- ggplot(data) +
  geom_bar(
    data = data,
    aes(x = zScore, y = interaction(Description, ONTOLOGY), fill = ONTOLOGY),
    stat = "identity"
  ) +
  scale_fill_manual(values = color, name = "GO&KEGG") +
  geom_text(
    aes(x = 0.1, y = interaction(Description, ONTOLOGY), label = Description, fontface = 'bold'),
    size = 3, hjust = 0, color = "black"
  ) +
  geom_text(
    aes(x = 0.1, y = interaction(Description, ONTOLOGY), label = geneID_top15),
    size = 2, hjust = 0, vjust = 2.5, color = "black"
  ) +
  geom_point(
    aes(x = -0.75, y = interaction(Description, ONTOLOGY), size = Count, fill = ONTOLOGY),
    shape = 21
  ) +
  geom_text(
    aes(x = -0.75, y = interaction(Description, ONTOLOGY), label = Count),
    size = 2.5, fontface = 'bold'
  ) +
  scale_size(range = c(4, 8), guide = guide_legend(override.aes = list(fill = "black"))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1)), limits = c(-1.5, 10)) +
  labs(x = "zScore", y = "Description") +
  theme(
    legend.title = element_text(face = "bold", color = "#000000", size = 12),
    legend.text  = element_text(color = "#000000", size = 10),
    ggh4x.axis.nestline.y = element_line(size = 3, color = color),
    ggh4x.axis.nesttext.y = element_text(face = "bold", color = color, hjust = 0.5, size = 12, angle = 90),
    legend.background = element_blank(),
    panel.background  = element_blank(),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    axis.title.y      = element_blank(),
    axis.text         = element_text(color = "#000000", size = 12),
    axis.text.y       = element_blank(),
    axis.ticks        = element_blank(),
    axis.title        = element_text(color = "#000000", size = 14),
    axis.ticks.y.left = element_line(color = NA)
  ) +
  scale_y_discrete(labels = function(y) str_wrap(y, width = 60)) +
  guides(
    y     = "axis_nested",
    y.sec = guide_axis_manual(breaks = 1:nrow(data), labels = data$Description),
    fill  = guide_legend(order = 1, reverse = TRUE),
    size  = guide_legend(order = 2, label = TRUE)
  ) +
  guides(fill = guide_legend(reverse = TRUE)) +
  # Add black border around the entire plot
  annotation_custom(grob = rectGrob(gp = gpar(col = "black", fill = NA, lwd = 2)),
                    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)

p
ggsave(p, file = 'Green_GO_KEGG_results.pdf', width = 7.5, height = 8)

library(openxlsx)

# Create a new Excel workbook and write results to sheets
wb <- createWorkbook()

addWorksheet(wb, "Up_GO")
writeData(wb, "Up_GO", Pos_go)

addWorksheet(wb, "Up_KEGG")
writeData(wb, "Up_KEGG", Pos_kk)

addWorksheet(wb, "Down_GO")
writeData(wb, "Down_GO", Neg_go)

addWorksheet(wb, "Down_KEGG")
writeData(wb, "Down_KEGG", Neg_kk)

addWorksheet(wb, "Green_GO")
writeData(wb, "Green_GO", green_go)

addWorksheet(wb, "Green_KEGG")
writeData(wb, "Green_KEGG", green_kk)

addWorksheet(wb, "Salmon_GO")
writeData(wb, "Salmon_GO", salmon_go)

addWorksheet(wb, "Salmon_KEGG")
writeData(wb, "Salmon_KEGG", salmon_kk)

# Save workbook
saveWorkbook(wb, file = "MND_GO_KEGG_Enrichment.xlsx", overwrite = TRUE)

############# Save candidate features
up_genes   <- DEA_data[DEA_data$Threshold %in% c('Up'), ]$ID
down_genes <- DEA_data[DEA_data$Threshold %in% c('Down'), ]$ID

up_dif_express_module   <- green_module[green_module %in% up_genes]
down_dif_express_module <- salmon_module[salmon_module %in% down_genes]
length(c(up_dif_express_module, down_dif_express_module))

getwd()
setwd('DEA&EA')
save(up_dif_express_module, down_dif_express_module, file = 'up&down_regulated_module_DEGs.rda')
