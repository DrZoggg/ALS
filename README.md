# ABCA1 acts as a protective modulator in amyotrophic lateral sclerosis

This repository contains the complete R analysis scripts for the study "**ABCA1 acts as a protective modulator in amyotrophic lateral sclerosis**". 

The project systematically investigates the role of ABCA1 in Amyotrophic Lateral Sclerosis (ALS), beginning with a global burden analysis, followed by biomarker discovery from peripheral blood transcriptomics, diagnostic model construction, and causal validation using Mendelian Randomization. The findings are further validated through an in-house serum ELISA cohort and multiple external transcriptomic datasets from post-mortem tissues, exercise cohorts, and mouse models.

# Analysis Workflow

The scripts are numbered 0-10 to be run in logical order.

0. Global Burden of Disease (Epidemiology)
   Script: 0.GBD_analysis.R


Purpose: Loads, processes, and analyzes GBD 2021 data for motor neuron disease.


Action: Calculates Estimated Annual Percentage Change (EAPC) and generates global trend plots and world maps for DALYs, incidence, prevalence, and mortality.



1. Data Preprocessing & Cohort Generation
   Script: 1.data_proces.R


Purpose: Prepares the primary transcriptomic datasets (GSE112676 and GSE112680) from GEO.


Action:

Loads raw data for the discovery cohort (GSE112676, n=741) and the validation cohort (GSE112680, n=323).


Performs gene annotation (Symbol, Entrez, ENSEMBL) and normalization .

Splits the GSE112676 cohort into a 70% Training set (train_cohort.rda) and a 30% Test set (test_cohort.rda) using stratified sampling.

Saves the independent GSE112680 cohort as the Validation set (valid_cohort.rda).

Prepares and saves the full GSE112676 expression matrix as input for WGCNA (WGCNA_input.rda).



2. WGCNA Co-Expression Analysis
   Script: 2.WGCNA.R


Purpose: Identifies modules of co-expressed genes related to ALS in the full discovery cohort (GSE112676).



Action: Selects top 5000 most variable genes, builds an unsigned network (power=18) , identifies modules, and correlates module eigengenes (MEs) with the ALS diagnosis trait.



Output: Identifies significant modules (e.g., 'salmon' [negative corr.] and 'green' [positive corr.])  and saves their gene lists (WGCNA_Enrich.rda).


3. Differential Expression & Enrichment Analysis (DEA & EA)
   Script: 3.DA&EA.R


Purpose: Identifies differentially expressed genes (DEGs) and performs functional enrichment on the full GSE112676 cohort.

Action:

Runs limma (ALS vs. CON) and identifies 1,079 DEGs (FDR < 0.05, |log₂FC| > 0.25).



Generates volcano plots and heatmaps.


Intersects DEGs with the WGCNA modules (identifying 449 candidate genes).


Performs GO/KEGG enrichment analysis.

Output: Saves the 449 DEGs found within the significant WGCNA modules (up&down_regulated_module_DEGs.rda).

4. Machine Learning Feature Selection
   Scripts: 4.FeatSel_RF.R, 4.FeatSel_SVM&LASSO.R


Purpose: Narrows down the 449 candidate genes from Step 3 using three different ML algorithms.

Action:


..._RF.R: Applies Random Forest (rfPermute) to identify features with high importance (MeanDecreaseAccuracy, Gini) and significance (p < 0.05) .


..._SVM&LASSO.R: Applies LASSO (cv.glmnet with lambda.1se) and SVM-RFE (using caret and msvmRFE.R)  to select features.


Output: Saves feature lists: features_RFmodel.rds, LASSO_features.rds, SVM_features.rds.

5. Biomarker Signature Definition
   Script: 5.RF_SVM_LASSO_Intersect.R


Purpose: Defines the final diagnostic signature by finding the consensus of the three ML methods.

Action: Loads the three feature lists from Step 4 and finds their intersection.


Output: Saves the final 9-gene signature (ABCA1, DDX51, etc.) as common_genes.rds.

6. Single-Gene Signature Validation
   Script: 6.ML_GeneIntersect_Eval.R

Purpose: Evaluates the diagnostic performance of each individual gene in the 9-gene signature.

Action: Loads the independent validation cohort (valid_cohort.rda, GSE112680). For each of the 9 genes, it generates sex-stratified (Male/Female) boxplots and ROC curves (ALS vs. CON).

7. Diagnostic Model Construction & Validation
   Scripts: 7.Random_forest_diagnostic_model.R, 7.LASSO diagnostic model.R


Purpose: Builds and validates a multi-gene diagnostic panel using the 9-gene signature.

Action:

Trains both RF and LASSO models on the train_cohort.rda (GSE112676 split).

Validates the models on the test_cohort.rda (GSE112676 split) and the external valid_cohort.rda (GSE112680).


Output: Generates combined ROC curves (Train/Test/Validation) and model calibration plots (using custom mgcv/glm functions).



8. Causal Inference (Mendelian Randomization)
   Script: 8.MR_ABCA1_to_ALS.R


Purpose: Investigates the causal role of the signature genes, focusing on ABCA1.

Action:

Screens all 9 genes for a causal link to ALS using TwoSampleMR.

Performs a comprehensive, multi-outcome (4 different ALS GWAS) and multi-method sensitivity analysis (MR-PRESSO, Egger, LOO plots, funnel plots, Steiger filtering) for ABCA1 .




9. In-House Cohort Validation (ELISA)
   Script: 9.In-house ABCA1 evaluation.R


Purpose: Validates ABCA1 findings at the protein level using an in-house serum ELISA cohort (15 ALS vs. 15 controls).


Action:

Loads and cleans clinical/ELISA data.

Performs MICE imputation for missing biochemical data.

Generates the main Baseline Table (Table 1) using compareGroups.

Generates the primary figure comparing ABCA1 protein levels (ALS vs. CON) using a t-test.

Performs stratified analyses (by Sex, by Age) and generates plots .

Analyzes correlations between ABCA1 and clinical variables (e.g., BMI, LDL) using pooled Spearman correlations.


Builds a pooled multivariable linear model (ABCA1 ~ Group + BMI + LDL).


Analyzes ABCA1 vs. ALSFRS-R scores (Spearman, LM, median-split plots).


10. External Multi-Cohort Validation
    Script: 10.External Validation Cohort.R


Purpose: Systematically characterizes ABCA1/Abca1 expression in diverse external cohorts (as shown in Figure 8).

Action: Loads and analyzes 7 distinct datasets, performing differential expression analysis for ABCA1/Abca1:

GSE234297: Peripheral Blood (ALS vs. CON); DESeq2. (Fig 8A) 


GSE212131: Peripheral Blood (Short vs. Long duration ALS); limma. (Fig 8B) 



GSE272626: Human Spinal Cord (Cervical & Lumbar; ALS vs. CTL); DESeq2. (Fig 8C) 


GSE120374: Mouse Spinal Cord (SOD1-G93A vs. WT; per-age: p30, p70, p100, p120); DESeq2. (Fig 8D) 



GSE250455: ALS Muscle Biopsy (Pre- vs. Post-exercise); paired DESeq2. (Fig 8E) 


GSE153960: Postmortem Brain/Spinal Cord (Multiple tissues; ALS vs. CON vs. OND); DESeq2. (Fig 8F) 


GSE124439: Postmortem Brain (Frontal & Motor Cortex; ALS vs. CON vs. OND); DESeq2. (Fig 8G) 


Output: Generates standardized violin/box plots for ABCA1 expression for each comparison.

# Citation

If you use this code or data in your research, please cite our paper:

ABCA1 acts as a protective modulator in amyotrophic lateral sclerosis

# Contact

For questions about the code or analysis, please open an issue in this repository.
