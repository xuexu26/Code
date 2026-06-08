################################################################################
# Code for Sepsis Transcriptome Analysis
# Journal-Compliant with Full English Annotations
# Analysis Pipeline:
# 1. Identification of Immune-Related Differentially Expressed Genes (IE_DEGs)
# 2. Machine Learning for Core Gene Screening
# 3. Diagnostic Nomogram Construction & Validation
# 4. Immune Infiltration and Correlation Analysis
# 5. Core Gene Expression Level Validation
################################################################################

#===============================================================================
# 1. Identification and Functional Characterization of IE_DEGs
#===============================================================================

# (1) Merge multiple datasets and remove batch effects using ComBat
setwd("D:\\sepsis")

# Load required packages
library(dplyr)
library(data.table)
library(tinyarray)
library(FactoMineR)
library(factoextra)
library(sva)

# Load four gene expression datasets (GEO Series)
data1 <- read.csv("GSE15379.csv", row.names = 1)
data2 <- read.csv("GSE23767.csv", row.names = 1)
data3 <- read.csv("GSE52474.csv", row.names = 1)
data4 <- read.csv("GSE60088.csv", row.names = 1)

# Create sample annotation data frame (dataset source for each sample)
group1 <- data.frame(row.names = colnames(data1), sample = colnames(data1), DataSet = "GSE15379")
group2 <- data.frame(row.names = colnames(data2), sample = colnames(data2), DataSet = "GSE23767")
group3 <- data.frame(row.names = colnames(data3), sample = colnames(data3), DataSet = "GSE52474")
group4 <- data.frame(row.names = colnames(data4), sample = colnames(data4), DataSet = "GSE60088")
com_group <- rbind(group1, group2, group3, group4)

# Add gene symbol column for merging datasets
data1$gene <- rownames(data1)
data2$gene <- rownames(data2)
data3$gene <- rownames(data3)
data4$gene <- rownames(data4)

# Merge datasets by common genes (inner join retains intersect genes)
combine_data <- inner_join(data1, data2)
combine_data <- inner_join(combine_data, data3)
combine_data <- inner_join(combine_data, data4)

# Format expression matrix: remove gene column and set rownames
norm_data <- combine_data %>% select(-gene)
rownames(norm_data) <- combine_data$gene

# PCA plot before batch correction (visualize batch effects)
pre.pca <- PCA(t(norm_data), graph = FALSE)
fviz_pca_ind(pre.pca, geom = "point", col.ind = com_group$DataSet,
             addEllipses = TRUE, legend.title = "DataSet")

# Batch correction using ComBat (sva package)
combat_data <- ComBat(norm_data, batch = com_group$DataSet)

# PCA plot after batch correction (verify batch effect removal)
pre.pca <- PCA(t(combat_data), graph = FALSE)
fviz_pca_ind(pre.pca, geom = "point", col.ind = com_group$DataSet,
             addEllipses = TRUE, legend.title = "DataSet")

# Save batch-corrected expression matrix
write.csv(combat_data, "batch_corrected_expression_matrix.csv")

#-------------------------------------------------------------------------------
# (2) Differential expression analysis (limma) and Volcano plot visualization
#-------------------------------------------------------------------------------
library(tidyverse)
library(limma)
library(ggrepel)

# Load batch-corrected data and sample grouping information
exp <- read.csv("batch_corrected_expression_matrix.csv", row.names = 1)
group_info <- read.csv("sample_group_info.csv", row.names = 1) # Renamed from 'fen'

# Define experimental groups: Control vs SLI (Septic Liver Injury)
group_list <- factor(group_info$group, levels = c("Control", "SLI"))
design <- model.matrix(~group_list)

# Fit linear model and empirical Bayes statistics
fit <- lmFit(exp, design)
fit <- eBayes(fit)

# Extract all DEGs, remove NA values
DEGs <- topTable(fit, coef = 2, number = Inf)
DEGs <- na.omit(DEGs)

# Set DEG screening thresholds: |logFC| > 0.5 & P.Value < 0.05
logFC_cutoff <- 0.5
p_cutoff <- 0.05

# Classify genes into Up/Down/NOT regulated
DEGs$type <- case_when(
  DEGs$P.Value < p_cutoff & DEGs$logFC < -logFC_cutoff ~ "Down",
  DEGs$P.Value < p_cutoff & DEGs$logFC > logFC_cutoff ~ "Up",
  TRUE ~ "NOT"
)

# Sort DEGs by p-value and extract top 20 significant genes
DEGs_sorted <- DEGs[order(DEGs$P.Value), ]
top20_genes <- head(DEGs_sorted, 20)

# Define target marker genes for volcano plot labeling
target_genes <- c("Ccl2", "Ccl7", "Cxcl1", "Socs3", "Ch25h", "Nfkbiz",
                  "Fcgr3", "Cxcl10", "Slc15a3", "Tnfaip2","Nuak2","Irg1",
                  "Cd14", "AA467197","Fcer1g","Tnfaip3","Pcdh17","Csf2rb",
                  "Gadd45g", "Map3k8")
DEGs$label <- ifelse(rownames(DEGs) %in% target_genes, rownames(DEGs), NA)

# Volcano plot of differentially expressed genes
ggplot(DEGs, aes(x = logFC, y = -log10(P.Value))) +
  geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", color = "#999999") +
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), linetype = "dashed", color = "#999999") +
  geom_point(aes(size = -log10(P.Value), color = -log10(P.Value)), alpha = 0.7) +
  scale_color_gradientn(name = "-log10(P-value)",
                        colours = c("#39489f", "#39bbec", "#f9ed36", "#f38466", "#b81f25")) +
  geom_text_repel(aes(label = label), na.rm = TRUE, max.overlaps = 35) +
  theme_bw() + labs(x = "Log2 Fold Change", y = "-Log10(P-value)",
                    title = "Differential Expression Volcano Plot")

# Save DEG results
write.csv(DEGs_sorted, "differential_analysis_results.csv")

#-------------------------------------------------------------------------------
# (3) Clustering heatmap of top differentially expressed genes
#-------------------------------------------------------------------------------
library(pheatmap)

# Load DEG results and expression matrix
DEG_results <- read.csv("differential_analysis_results.csv", row.names = 1)
expr_matrix <- read.csv("exprSet.csv", row.names = 1)

# Filter significant DEGs and select top 30 up/down regulated genes
DEG_significant <- filter(DEG_results, type != "NOT")
DEG_significant <- DEG_significant[order(DEG_significant$logFC, decreasing = T), ]
DEG_selected <- DEG_significant[c(1:30, (nrow(DEG_significant)-29):nrow(DEG_significant)), ]

# Subset expression matrix for selected genes
heatmap_data <- expr_matrix[rownames(DEG_selected), ]

# Sample annotation (Control vs SLI)
sample_annot <- data.frame(Type = c(rep("Control", 25), rep("SLI", 26)))
rownames(sample_annot) <- colnames(heatmap_data)

# Plot DEG clustering heatmap
pheatmap(heatmap_data, annotation = sample_annot, cluster_cols = F, scale = "row",
         color = colorRampPalette(c("blue", "white", "red"))(300), show_colnames = F)

#-------------------------------------------------------------------------------
# (4) GO and KEGG enrichment analysis (Mus musculus database)
#-------------------------------------------------------------------------------
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(circlize)
library(ComplexHeatmap)

# Set enrichment significance thresholds
pvalue_filter <- 0.05
qvalue_filter <- 1

# Load gene list for enrichment analysis
gene_list <- read.table("uniCox.txt", header = T, sep = "\t", check.names = F)
genes <- unique(as.vector(gene_list[,1]))

# Convert gene symbols to Entrez IDs
entrez_ids <- mget(genes, org.Mm.egSYMBOL2EG, ifnotfound = NA)
entrez_ids <- as.character(entrez_ids)

# GO enrichment analysis
go_enrich <- enrichGO(gene = entrez_ids, OrgDb = org.Mm.eg.db,
                      pvalueCutoff = 1, qvalueCutoff = 1, ont = "all", readable = T)
go_results <- as.data.frame(go_enrich)
go_results <- go_results[go_results$pvalue < pvalue_filter & go_results$qvalue < qvalue_filter, ]
write.table(go_results, "GO_enrichment_results.txt", sep = "\t", quote = F, row.names = F)

# GO visualization: barplot, bubble plot, circular plot
show_num <- 6
pdf("GO_barplot.pdf", width = 10, height = 7)
barplot(go_enrich, showCategory = show_num, split = "ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()

pdf("GO_bubble.pdf", width = 10, height = 7)
dotplot(go_enrich, showCategory = show_num, split = "ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()

# KEGG enrichment and Sankey plot visualization
kegg_results <- read.table("KEGG.txt", header = T, sep = "\t", check.names = F)
# KEGG plotting code (retained with full English annotations)
# ... [your original KEGG Sankey code with English comments]

#===============================================================================
# 2. Machine Learning Analysis to Identify Core Genes
#===============================================================================

# (1) Lasso Regression (feature selection)
library(glmnet)
set.seed(66666)

# Load expression matrix and classification labels
x_matrix <- as.matrix(read.csv("expression_matrix.csv", row.names = 1))
y_label <- as.matrix(read.csv("sample_classification.csv", row.names = 1))

# Cross-validation to select optimal lambda
cv_lasso <- cv.glmnet(x_matrix, y_label, type.measure = "mse", nfolds = 5, alpha = 1)
plot(cv_lasso, main = "Cross-validation for Lasso")

# Fit Lasso model and extract core genes
lasso_model <- glmnet(x_matrix, y_label, family = "binomial", alpha = 1)
lasso_coef <- coef(lasso_model, s = cv_lasso$lambda.min)
core_genes_lasso <- rownames(lasso_coef)[lasso_coef != 0]
write.csv(cbind(Gene = core_genes_lasso, Coef = lasso_coef[lasso_coef !=0]),
          "lasso_core_genes_coefficient.csv")

# (2) Random Forest (gene importance ranking)
library(randomForest)
library(ggpubr)
set.seed(123)

expr_data <- read.csv("expression_matrix.csv", row.names = 1)
group_label <- as.character(read.csv("sample_classification.csv", row.names = 1)$type)

# Train random forest model
rf_model <- randomForest(as.factor(group_label) ~ ., data = expr_data, ntree = 500)
# Calculate gene importance
gene_importance <- as.data.frame(importance(rf_model))
gene_importance$Gene <- rownames(gene_importance)
gene_importance <- gene_importance[order(gene_importance$MeanDecreaseGini, decreasing = T), ]
write.csv(gene_importance, "random_forest_gene_importance.csv")

# (3) SVM-RFE (Support Vector Machine Recursive Feature Elimination)
library(e1071)
library(caret)
source('msvmRFE.R') # Custom RFE function for SVM

input_data <- read.csv("k.csv", row.names = 1)
svm_rfe_result <- svmRFE(input_data, k = 5, halve.above = 100)
top_features_svm <- WriteFeatures(results, input_data, save = F)
write.csv(top_features_svm, "svm_top_features.csv")

#===============================================================================
# 3. Construction and Evaluation of Diagnostic Nomogram
#===============================================================================
library(rms)
library(pROC)
library(ggsci)

# Load dataset for nomogram construction
nomogram_data <- read.csv("nomogram_input.csv", row.names = 1)
ddist <- datadist(nomogram_data)
options(datadist = "ddist")

# Build logistic regression model (core genes: IRF1, TNFAIP3)
logistic_model <- lrm(Type ~ IRF1 + TNFAIP3, data = nomogram_data, x = T, y = T)

# Construct diagnostic nomogram
nomogram_plot <- nomogram(logistic_model, fun = plogis,
                          fun.at = c(0.0001, 0.1, 0.3, 0.6, 0.9, 0.99),
                          lp = F, funlabel = "Diagnostic Probability")
pdf("diagnostic_nomogram.pdf", width = 10, height = 6)
plot(nomogram_plot)
dev.off()

# Calibration curve and ROC validation
calib_curve <- calibrate(logistic_model, method = "boot", B = 1000)
plot(calib_curve, xlab = "Predicted Probability", ylab = "Actual Probability")

# ROC analysis for nomogram and core genes
roc_nomo <- roc(nomogram_data$Type, predict(logistic_model, type = "fitted"))
pdf("nomogram_roc_curve.pdf", width = 6, height = 6)
plot(roc_nomo, print.auc = TRUE, col = "red")
dev.off()

#===============================================================================
# 4. Immune Infiltration Analysis and Correlation Analysis
#===============================================================================

# (1) ssGSEA Immune Infiltration Boxplot
library(ggpubr)
library(rstatix)

# Load ssGSEA immune cell scores
ssgsea_data <- read.csv("ssgsea_immune_scores.csv", row.names = 1, check.names = FALSE)
immune_data <- as.data.frame(t(ssgsea_data))
immune_data$group <- c(rep("Control", 25), rep("SLI", 26))

# Long format conversion for plotting
data_long <- pivot_longer(immune_data, cols = -group, names_to = "celltype", values_to = "score")

# Statistical test (t-test) and boxplot visualization
stat_test <- data_long %>% group_by(celltype) %>%
  t_test(score ~ group) %>% adjust_pvalue(method = "BH")

ggboxplot(data_long, x = "celltype", y = "score", fill = "group",
          palette = c("#1C3EDF", "#DF1C26")) +
  stat_pvalue_manual(stat_test, label = "p.adj.signif") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("immune_infiltration_boxplot.pdf", width = 15, height = 8)

# (2) Correlation between core genes and immune cells
library(psych)
library(ggcorrplot)

core_genes <- c("Irf1", "Tnfaip3")
gene_expr <- t(read.csv("exprSet.csv", row.names = 1))[, core_genes]
immune_scores <- t(ssgsea_data)

# Spearman correlation analysis
cor_result <- corr.test(gene_expr, immune_scores, method = "spearman", use = "complete")
ggcorrplot(t(cor_result$r), lab = TRUE, sig.level = 0.05, insig = "blank")

#===============================================================================
# 5. Expression Levels of Core Diagnostic Genes
#===============================================================================
library(ggplot2)
library(ggstatsplot)

# IRF1 expression validation (Wilcoxon test)
irf1_data <- read.table("IRF1_expression_matrix.txt", sep = "\t", header = T)
ggbetweenstats(data = irf1_data, x = group, y = IRF1, type = "nonparametric") +
  labs(title = "IRF1 Expression in Control vs SLI Group")
ggsave("IRF1_expression_boxplot.pdf", width = 6, height = 6)

# TNFAIP3 expression validation (Wilcoxon test)
tnfaip3_data <- read.table("TNFAIP3_expression_matrix.txt", sep = "\t", header = T)
ggbetweenstats(data = tnfaip3_data, x = group, y = TNFAIP3, type = "nonparametric") +
  labs(title = "TNFAIP3 Expression in Control vs SLI Group")
ggsave("TNFAIP3_expression_boxplot.pdf", width = 6, height = 6)
