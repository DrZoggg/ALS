setwd("./WGCNA/")
options(stringsAsFactors = F)
library(WGCNA)
load("WGCNA_input.rda")

# MAD = median absolute deviation (select top 5000 most variable genes);
# transpose to samples x genes as required by WGCNA
WGCNA_matrix <- t(wgcna_input[order(apply(wgcna_input, 1, mad),
                                    decreasing = T)[1:5000], ])

# Check whether there are samples/genes with too many missing values
datExpr0 <- WGCNA_matrix
gsg <- goodSamplesGenes(datExpr0, verbose = 3)  # check missingness; 'verbose' controls details

# Whether filtering is needed; TRUE = no filtering, FALSE = filtering needed
gsg$allOK

# TRUE means no filtering performed
if (!gsg$allOK) {
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples) > 0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

### Identify outliers via sample clustering and remove them ###
sampleTree <- hclust(dist(datExpr0), method = "average")  # average-linkage hierarchical clustering

options(repr.plot.width = 8, repr.plot.height = 8)
p1 <- plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "", cex.lab = 1.5,
           cex.axis = 1.5, cex.main = 2, label = FALSE
) +  # plot with titles and font sizes
  abline(h = 65, col = "red")
p1

# Remove outlier cluster; 'cutHeight' must match the horizontal line above
clust <- cutreeStatic(sampleTree, cutHeight = 65, minSize = 10)  # static cut; minSize = minimum cluster size
table(clust)  # clustering result
print('过滤的样本是:GSM3078239')

# Keep the cluster that contains the desired samples
keepSamples <- (clust == 1)
datExpr <- datExpr0[keepSamples, ]

# Record numbers of genes and samples for later visualization
nGenes <- ncol(datExpr)   # number of genes
print(paste0('WGCNA分析的基因数目:', nGenes))
nSamples <- nrow(datExpr) # number of samples
print(paste0('WGCNA分析的样本数目:', nSamples))

# Save the final input objects for WGCNA
save(datExpr, nGenes, nSamples, file = "Step01-WGCNA_input.Rda")


#### **Run WGCNA analysis**
rm(list = ls())
library(WGCNA)
library(patchwork)
library(export)
library(ggplot2)

# Enable multithreading
enableWGCNAThreads()

# Load previous step's data
load("Step01-WGCNA_input.Rda")

powers = c(1:30)
# Network topology analysis
# This step uses expression data only (no trait involved yet)
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
test <- sft$fitIndices
# Results are ready; we now visualize diagnostics

# Plot parameters
cex1 <- 0.9
softPowerPick <- 18  # chosen soft-thresholding power

# Fig 1: Scale independence
pdf("soft_threshold_scale_independence.pdf", width = 5, height = 5)
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit (signed R^2)",
     type = "b", pch = 19, col = "steelblue",
     main = "Scale independence")
abline(h = 0.9, col = "red", lty = 2)
points(softPowerPick, 
       -sign(sft$fitIndices[which(sft$fitIndices[, 1] == softPowerPick), 3]) *
         sft$fitIndices[which(sft$fitIndices[, 1] == softPowerPick), 2],
       col = "darkred", pch = 21, bg = "red", cex = 1.5)
text(softPowerPick,
     -sign(sft$fitIndices[which(sft$fitIndices[, 1] == softPowerPick), 3]) *
       sft$fitIndices[which(sft$fitIndices[, 1] == softPowerPick), 2],
     labels = paste0("Power=", softPowerPick), pos = 3, col = "darkred")
dev.off()

# Fig 2: Mean connectivity
pdf("soft_threshold_mean_connectivity.pdf", width = 5, height = 5)
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "b", pch = 19, col = "steelblue",
     main = "Mean connectivity")
points(softPowerPick, 
       sft$fitIndices[which(sft$fitIndices[, 1] == softPowerPick), 5],
       col = "darkred", pch = 21, bg = "red", cex = 1.5)
text(softPowerPick,
     sft$fitIndices[which(sft$fitIndices[, 1] == softPowerPick), 5],
     labels = paste0("Power=", softPowerPick), pos = 3, col = "darkred")
dev.off()

plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit (signed R^2)",
     type = "b", pch = 19, col = "steelblue",
     main = "Scale independence")
abline(h = 0.9, col = "red", lty = 2)
points(softPowerPick, 
       -sign(sft$fitIndices[which(sft$fitIndices[, 1] == softPowerPick), 3]) *
         sft$fitIndices[which(sft$fitIndices[, 1] == softPowerPick), 2],
       col = "darkred", pch = 21, bg = "red", cex = 1.5)
text(softPowerPick,
     -sign(sft$fitIndices[which(sft$fitIndices[, 1] == softPowerPick), 3]) *
       sft$fitIndices[which(sft$fitIndices[, 1] == softPowerPick), 2],
     labels = paste0("Power=", softPowerPick), pos = 3, col = "darkred")

plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "b", pch = 19, col = "steelblue",
     main = "Mean connectivity")
points(softPowerPick, 
       sft$fitIndices[which(sft$fitIndices[, 1] == softPowerPick), 5],
       col = "darkred", pch = 21, bg = "red", cex = 1.5)
text(softPowerPick,
     sft$fitIndices[which(sft$fitIndices[, 1] == softPowerPick), 5],
     labels = paste0("Power=", softPowerPick), pos = 3, col = "darkred")


softPower = 18  # chosen soft-thresholding power
adjacency = adjacency(datExpr, power = softPower)
TOM = TOMsimilarity(adjacency)       # Topological Overlap Matrix
dissTOM = 1 - TOM                    # dissimilarity for gene clustering

# -------------------------------------------------------------------------
geneTree = hclust(as.dist(dissTOM), method = "average")
# Plot the resulting clustering tree (dendrogram)
plot(geneTree, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)

minModuleSize = 30  # user-defined minimum module size

# Module identification using dynamic tree cut
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,  # larger deepSplit → more modules
                            minClusterSize = minModuleSize,
                            cutHeight = 
)
print('dynamicMods:')
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
print('dynamicColors:')
table(dynamicColors)  # 'grey' = genes outside all modules

# Optional: modify module colors ------------------------------------------------
# dynamicColors=gsub('red','DeepSkyBlue',dynamicColors)
# dynamicColors=gsub('green','orange',dynamicColors)
# dynamicColors=gsub('brown','red',dynamicColors)
# dynamicColors=gsub('magenta','green',dynamicColors)
# dynamicColors=gsub('DeepSkyBlue','brown',dynamicColors)
# table(dynamicColors)
# gsub(pattern, replacement, x) where x is 'dynamicColors'
# ------------------------------------------------------------------------------

options(repr.plot.width = 7, repr.plot.height = 5)
plotDendroAndColors(geneTree, dynamicColors, groupLabels = "",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
graph2pdf(file = '1st_WGCNA-Dendrogram.pdf', width = 7, height = 5)

# Topological heatmap ----------------------------------------------------------
nSelect = 2000
set.seed(10)  # reproducibility
nGenes = ncol(datExpr)
select = sample(nGenes, size = nSelect)
selectTOM = dissTOM[select, select]
# Re-cluster the subset (no direct way to restrict the original tree)
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = dynamicColors[select]
# Raise to a power to improve contrast; set diagonal to NA for clarity
plotDiss = selectTOM^softPower
diag(plotDiss) = NA

# png output (high resolution)
png(filename = "mad_Network_heatmap_all_genes.png",
    width = 8, height = 8, units = 'in', res = 1000)
TOMplot(plotDiss,
        selectTree,
        selectColors,
        main = "Network heatmap plot, all genes")
dev.off()

# Module correlations ----------------------------------------------------------
# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes

# Dissimilarity and clustering of module eigengenes
MEDiss = 1 - cor(MEs)
METree = hclust(as.dist(MEDiss), method = "average")

options(repr.plot.width = 7, repr.plot.height = 5)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.15  # merge threshold for modules
abline(h = MEDissThres, col = "red")  # cut line

options(repr.plot.width = 6, repr.plot.height = 9)
MEs <- orderMEs(MEs)
plotEigengeneNetworks(MEs, "Eigengene adjacency heatmap",
                      marDendro = c(3, 3, 2, 6),
                      marHeatmap = c(3, 3, 2, 1), plotDendrograms = T,
                      xLabelsAngle = 90)
pdf(file = 'mad_WGCNA-Eigengene adjacency heatmap.pdf', width = 6, height = 9)
plotEigengeneNetworks(MEs, "Eigengene adjacency heatmap",
                      marDendro = c(0, 3, 2, 4.5),
                      marHeatmap = c(3, 3, 1, 2), plotDendrograms = T,
                      xLabelsAngle = 90)
dev.off()

# Use the chosen soft-thresholding power; type = signed/unsigned controls sign handling
hubs = chooseTopHubInEachModule(datExpr, colorh = dynamicColors, power = 18, type = 'signed')
hubs

# Automatic merging of similar modules
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, quantileSummary = 'mean', verbose = 3)

# Colors of merged modules and their eigengenes
mergedColors = merge$colors
table(mergedColors)
mergedMEs = merge$newMEs

options(repr.plot.width = 9, repr.plot.height = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Module", "Merged Module"),
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE,
                    guideHang = 0.05)
graph2pdf(file = 'mad_WGCNA-Dendrogram_merged.pdf', width = 7, height = 6)

moduleColors = mergedColors
# Construct numerical labels corresponding to colors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder) - 1
MEs = mergedMEs

# Save module colors and labels for subsequent analyses
save(MEs, moduleLabels, moduleColors, geneTree, file = "mad_networkConstruction-stepByStep.Rda")


#### **Module–trait association analysis**

# -------------------------------------------------------------------------
### Link modules to traits (key step, conceptually simple)
### 3. Relating modules to external clinical traits and identifying important genes
#=====================================================================================
#=====================================================================================
rm(list = ls())
setwd("WGCNA")

# Load WGCNA package
library(WGCNA)

# Load input expression and module data
load("Step01-WGCNA_input.Rda")
load("mad_networkConstruction-stepByStep.Rda")

# Define numbers of genes and samples
nGenes = ncol(datExpr)    # number of genes
nSamples = nrow(datExpr)  # number of samples

# Recalculate eigengenes with color labels
test <- moduleEigengenes(datExpr, moduleColors)
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)  # reorder for plotting (cosmetic)

#############################################################################
### Compute correlations; datTraits must be aligned by sample rows
load("WGCNA_input.rda")
load("../Data_preprocessing/disc_cohort.rda")
wgcna_input <- wgcna_input[, rownames(datExpr)]
wgcna_input <- as.data.frame(wgcna_input)
clin <- disc_cohort_clin[rownames(datExpr), ]
clin$diagnosis[clin$diagnosis == "ALS"] <- 1
clin$diagnosis[clin$diagnosis == "CON"] <- 0
clin2 <- data.frame(row.names = rownames(clin), diagnosis = clin$diagnosis)
# datTraits <- as.data.frame(lapply(clin1[,1:2],as.numeric)) # no clinical info
datTraits <- clin2
# rownames(datTraits) <- rownames(clin1)
datTraits$diagnosis <- as.numeric(datTraits$diagnosis)

moduleTraitCor = cor(MEs, datTraits, use = "p")
# p-values for the correlations
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

# Display correlations and p-values together
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)

options(repr.plot.width = 4.5, repr.plot.height = 8)
par(mar = c(3, 8, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = F,
               cex.text = 1,
               zlim = c(-1, 1),
               main = paste("Module-trait relationships"),
               xLabelsAngle = 0)

pdf(file = 'mad_WGCNA-Clin-Cor.pdf', width = 4.5, height = 8)
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = F,
               cex.text = 1,
               zlim = c(-1, 1),
               main = paste("Module-trait relationships"),
               xLabelsAngle = 0)
dev.off()

# Define variable 'Type' containing the trait (diagnosis)
# geneModuleMembership: correlations between each gene and module eigengenes (MM)
# MMPvalue: p-values for MM correlations
Type = as.data.frame(datTraits$diagnosis)
names(Type) = "diagnosis"

# From the 3rd character onward to strip the 'ME' prefix (module color names)
modNames = substring(names(MEs), 3)

# Compute gene-level correlations:
# MM (module membership) and GS (gene significance for the trait)
# use = "p" ignores missing values (NA)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep = "")
names(MMPvalue) = paste("p.MM", modNames, sep = "")

geneTraitSignificance = as.data.frame(cor(datExpr, Type, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(Type), sep = "")
names(GSPvalue) = paste("p.GS.", names(Type), sep = "")

# Negatively associated module: salmon
module = "salmon"
column = match(module, modNames)
table(moduleColors)

# Extract genes in the module
moduleGenes = moduleColors == module

# Essentially a scatter plot
MM <- geneModuleMembership[moduleGenes, column]
GS <- geneTraitSignificance[moduleGenes, 1]
library(ggplot2)
library(ggsci)
dd <- data.frame(MM = MM, GS = GS)
table(dd$MM > 0.6, dd$GS > 0.5)
print(paste0('max(MM)为:', max(MM)))
print(paste0('max(GS)为:', max(GS)))
cor <- cor.test(MM, GS)  # compute correlation between MM and GS for annotation

# Extract p-value, r, and CI
p <- cor$p.value
r <- cor$estimate
r_with_ci <- sprintf("%.3f (%.3f, %.3f)", r, cor$conf.int[1], cor$conf.int[2])
print(paste0('p值为:', p))
print(paste0('r_with_ci:', r_with_ci))

# Load libraries
library(ggplot2)
library(ggpubr)
library(ggExtra)
options(repr.plot.width = 5, repr.plot.height = 5)

# Salmon module color
module_color <- "salmon"
module <- "salmon"

# Compute correlation (again for plotting convenience)
cor_result <- cor.test(dd$MM, dd$GS)
r <- cor_result$estimate
p <- cor_result$p.value

# Main plot (clean background, bold title)
p <- ggplot(dd, aes(x = MM, y = GS)) +
  geom_point(shape = 1, color = '#EE6A50', size = 1.5) +
  geom_smooth(method = "lm", se = TRUE, color = module_color,
              fill = module_color, alpha = 0.15, size = 1.5) +
  geom_rug(sides = "bl", color = module_color, alpha = 0.6) +
  annotate("text",
           x = min(dd$MM) + 0.02,
           y = max(dd$GS),
           label = paste0("R = ", round(r, 3), "\np = ", format.pval(p, digits = 3, eps = 1e-3)),
           hjust = 0, vjust = 1, size = 5, fontface = "italic") +
  labs(
    x = paste("Module membership (", module, " module)", sep = ""),
    y = "Gene significance for ALS",
    title = paste0("MM vs. GS in ", module, " module (ALS)")
  ) +
  theme(
    # No background/grid
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Keep black panel border
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    # Axis text/titles
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(size = 14),
    # Bold centered title
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt")
  )

# Add marginal histograms
p_final <- ggMarginal(p,
                      type = "histogram",
                      xparams = list(fill = module_color, color = "white"),
                      yparams = list(fill = module_color, color = "white"),
                      size = 10)

# Show and save plot
p_final
ggsave("MM_vs_GS_salmon_ALS_clean.pdf", p_final, width = 5, height = 5)

# Save MM and GS per gene
dd$'Gene' <- colnames(datExpr)[moduleGenes]
head(dd)
library(openxlsx)
write.xlsx(dd, file = 'slamon_MM_GS.xlsx')

# Positively associated module: green
module = "green"
column = match(module, modNames)
table(moduleColors)

# Extract genes in the module
moduleGenes = moduleColors == module

# Essentially a scatter plot
MM <- geneModuleMembership[moduleGenes, column]
GS <- geneTraitSignificance[moduleGenes, 1]

library(ggplot2)
library(ggsci)
dd <- data.frame(MM = MM, GS = GS)
table(dd$MM > 0.6, dd$GS > 0.5)
print(paste0('max(MM)为:', max(MM)))
print(paste0('max(GS)为:', max(GS)))
cor <- cor.test(MM, GS)  # correlation between MM and GS

# Extract p-value, r, and CI
p <- cor$p.value
r <- cor$estimate
r_with_ci <- sprintf("%.3f (%.3f, %.3f)", r, cor$conf.int[1], cor$conf.int[2])
print(paste0('p值为:', p))
print(paste0('r_with_ci:', r_with_ci))

# Load libraries
library(ggplot2)
library(ggpubr)
library(ggExtra)
options(repr.plot.width = 5, repr.plot.height = 5)

# Green module color
module_color <- "green"
module <- "green"

# Compute correlation (again for plotting)
cor_result <- cor.test(dd$MM, dd$GS)
r <- cor_result$estimate
p <- cor_result$p.value

# Main plot (clean background, bold title)
p <- ggplot(dd, aes(x = MM, y = GS)) +
  geom_point(shape = 1, color = '#9BCD9B', size = 1.5) +
  geom_smooth(method = "lm", se = TRUE, color = '#2E8B57',
              fill = '#66CDAA', alpha = 0.5, size = 0.1) +
  geom_rug(sides = "bl", color = '#9BCD9B', alpha = 0.6) +
  annotate("text",
           x = min(dd$MM) + 0.02,
           y = max(dd$GS),
           label = paste0("R = ", round(r, 3), "\np = ", format.pval(p, digits = 3, eps = 1e-3)),
           hjust = 0, vjust = 1, size = 5, fontface = "italic") +
  labs(
    x = paste("Module membership (", module, " module)", sep = ""),
    y = "Gene significance for ALS",
    title = paste0("MM vs. GS in ", module, " module (ALS)")
  ) +
  theme(
    # No background/grid
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Keep black panel border
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    # Axis text/titles
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(size = 14),
    # Bold centered title
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt")
  )

# Add marginal histograms
p_final <- ggMarginal(p,
                      type = "histogram",
                      xparams = list(fill = '#9BCD9B', color = "white"),
                      yparams = list(fill = '#9BCD9B', color = "white"),
                      size = 10)

# Show and save plot
p_final
ggsave("MM_vs_GS_green_ALS_clean.pdf", p_final, width = 5, height = 5)

# Save MM and GS per gene
dd$'Gene' <- colnames(datExpr)[moduleGenes]
head(dd)
library(openxlsx)
write.xlsx(dd, file = 'green_MM_GS.xlsx')

# Extract genes in each module
library(dplyr)
gene_name <- rownames(as.data.frame(t(datExpr)))

salmon_module <- gene_name[which(moduleColors == "salmon")]
# pink_module <- as.data.frame(t(exp)) %>% select(pink)
green_module <- gene_name[which(moduleColors == 'green')]
# midnightblue_module <- as.data.frame(t(exp)) %>% select(midnightblue))
print('names of genes cotained in salmon_module:'); head(salmon_module, 20)
print('names of genes cotained in green_module:');  head(green_module, 20)

save(salmon_module, green_module, file = 'WGCNA_Enrich.rda')
