#1.Identification and Functional Characterization of IE_DEGs
#(1) merging multiple datasets and removing batch effects
setwd("D:\\Desktop\\data diff")

library(dplyr)
library(data.table)
library(tinyarray)

data1<-read.csv("GSE15379.csv",row.names = 1)
data2<-read.csv("GSE23767.csv",row.names = 1)
data3<-read.csv("GSE52474.csv",row.names = 1)
data4<-read.csv("GSE60088.csv",row.names = 1)

group1=data.frame(row.names = colnames(data1),
                  sample=colnames(data1),
                  DataSet="GSE15379")
group2=data.frame(row.names = colnames(data2),
                  sample=colnames(data2),
                  DataSet="GSE23767")
group3=data.frame(row.names = colnames(data3),
                  sample=colnames(data3),
                  DataSet="GSE52474")
group4=data.frame(row.names = colnames(data4),
                  sample=colnames(data4),
                  DataSet="GSE60088")
com_group=rbind(group1,group2,group3,group4)

data1$gene <- rownames(data1) 
data2$gene <- rownames(data2) 
data3$gene <- rownames(data3) 
data4$gene <- rownames(data4) 

combine_data <- inner_join(data1,data2,)
combine_data <- inner_join(combine_data,data3)
combine_data <- inner_join(combine_data,data4)

norm_data <- combine_data %>% select(-gene) 

rownames(norm_data) <- combine_data$gene 


library(FactoMineR)
library(factoextra)

pre.pca <- PCA(t(norm_data),
               graph = FALSE)
fviz_pca_ind(pre.pca,
             geom= "point",
             col.ind = com_group$DataSet,
             addEllipses = TRUE,
             legend.title="Group"
)

library(sva)
combat_data <- ComBat(norm_data, batch = com_group$DataSet) 
pre.pca <- PCA(t(combat_data),graph = FALSE)
fviz_pca_ind(pre.pca,
             geom= "point",
             col.ind = com_group$DataSet,
             addEllipses = TRUE,
             legend.title="Group"
)
write.csv(combat_data,"ComBat_normalized_expression.csv")



#(2) Differential expression analysis using limma and Volcano plot visualization
library(tidyverse)
library(limma)
exp<-read.csv("ComBat_normalized_expression.csv",row.names = 1)

fen<-read.csv("sample_group_info.csv",row.names = 1)

group_list <- factor(fen$group,levels = c("Control","SLI"))
design <- model.matrix(~group_list)
con <- lmFit(exp,design)
con2 <- eBayes(con)
DEG1 <- topTable(con2, coef = 2, number = Inf)
DEG2 = na.omit(DEG1) 

write.csv(DEG2,"raw_limma_DE_results.csv")
logFC_cut = 0.5
p_cut = 0.05

type1 = (DEG2$P.Value < p_cut)&(DEG2$logFC < -logFC_cut)
type2 = (DEG2$P.Value < p_cut)&(DEG2$logFC > logFC_cut)

DEG2$type = ifelse(type1,"Down",ifelse(type2,"Up","NOT"))

head(DEG2)

write.csv(DEG2,"classified_limma_DE_results.csv")

DEG_sorted <- DEG2[order(DEG2$P.Value), ]

top20_genes <- head(DEG_sorted, 20)
print(top20_genes[, c("logFC", "adj.P.Val", "P.Value", "type")])

top20_simple <- data.frame(
  Gene = rownames(top20_genes),
  logFC = round(top20_genes$logFC, 3),
  P.Value = formatC(top20_genes$adj.P.Val, format = "e", digits = 2),
  Direction = top20_genes$type
)
print(top20_simple)
write.csv(DEG_sorted, "sorted_DEGs_by_Pvalue.csv")
target_genes <- c("Ccl2", "Ccl7", "Cxcl1", "Socs3", "Ch25h", 
                  "Nfkbiz", "Fcgr3", "Cxcl10", "Slc15a3", 
                  "Tnfaip2","Nuak2","Irg1","Cd14", "AA467197","Fcer1g",
                  "Tnfaip3","Pcdh17","Csf2rb",
                  "Gadd45g", "Map3k8"
                  
)
DEG2$label <- ifelse(rownames(DEG2) %in% target_genes, rownames(DEG2), NA)

library(ggrepel)  

ggplot(DEG2, aes(logFC, -log10(P.Value))) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "#999999") + # 更合理的阈值线
  geom_point(aes(size = -log10(P.Value), color = -log10(P.Value)), alpha = 0.7) +
  scale_color_gradientn(
    name = "-log10(P-value)",
    colours = c("#39489f", "#39bbec", "#f9ed36", "#f38466", "#b81f25"),
    limits = c(0, max(-log10(DEG2$P.Value)))
  ) +
  scale_size_continuous(range = c(1, 4)) +
  geom_text_repel(  
    aes(label = label),
    na.rm = TRUE,
    max.overlaps =35,  
    box.padding = 0.5,  
    segment.color = "grey50",  
    min.segment.length = 0.2  
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    axis.text = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    x = "Log2 Fold Change",
    y = "-Log10(P-value)",
    title = "TOP20 Differential Expression Volcano Plot"
  )

#(3) Clustering heatmap of DEGs
library(pheatmap)
library(dplyr)
diff_gene <- read.csv("classified_limma_DE_results.csv",row.names = 1)
exprSet <- read.csv("ComBat_normalized_expression.csv",row.names = 1)
diff_gene <- filter(diff_gene, type != "NOT")
diff_gene <- diff_gene[order(diff_gene$logFC,decreasing = T),]
diff_gene <- diff_gene[c(1:30,(nrow(diff_gene)-29):nrow(diff_gene)),]

rt <- exprSet[rownames(diff_gene),]

ann <- data.frame(name = colnames(rt))
ann$Type <- c(rep("Control",25),rep("SLI",26))
rownames(ann) <- ann$name
ann <- ann[-1]

pheatmap(rt,                        
         annotation = ann,           
         cluster_cols = F,          
         color = colorRampPalette(c("blue", "white", "red"))(300),
         show_colnames = F,          
         show_rownames = T,          
         cluster_rows = T,  
         scale="row",                
         border_color = NA,         
         fontsize = 9,               
         fontsize_row=5,             
         fontsize_col=7,             
         treeheight_row=10,          
         treeheight_col=5,          
         cellwidth = 5,            
         cellheight = 5             
)

#(4)  GO and KEGG enrichment analysis
library(clusterProfiler)
library(org.Mm.eg.db) 
library(enrichplot)
library(ggplot2)
library(circlize)
library(RColorBrewer)
library(dplyr)
library(ComplexHeatmap)
R.utils::setOption("clusterProfiler.download.method","auto")

pvalueFilter=0.05      
qvalueFilter=1       
colorSel="qvalue"
if(qvalueFilter>0.05){
  colorSel="pvalue"
}
ontology.col=c("#00AFBB", "#E7B800", "#90EE90")

rt=read.table("uniCox.txt", header=T, sep="\t", check.names=F)    

genes=unique(as.vector(rt[,1]))
entrezIDs=mget(genes, org.Mm.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)

kk=enrichGO(gene=entrezIDs, OrgDb=org.Mm.eg.db, pvalueCutoff=1, qvalueCutoff=1, ont="all", readable=T)
GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]

write.table(GO, file="GO_enrichment_results.txt", sep="\t", quote=F, row.names = F)

showNum=6
if(nrow(GO)<18){
  showNum=nrow(GO)
}

pdf(file="GObarplot.pdf", width=10, height=7)
bar=barplot(kk, drop=TRUE, showCategory=showNum, label_format=130, split="ONTOLOGY", color=colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bar)
dev.off()

pdf(file="GObubble.pdf", width=10, height=7)
bub=dotplot(kk, showCategory=showNum, orderBy="GeneRatio", label_format=130, split="ONTOLOGY", color=colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bub)
dev.off()

data=GO[order(GO$pvalue),]
datasig=data[data$pvalue<0.05,,drop=F]
BP = datasig[datasig$ONTOLOGY=="BP",,drop=F]
CC = datasig[datasig$ONTOLOGY=="CC",,drop=F]
MF = datasig[datasig$ONTOLOGY=="MF",,drop=F]

BP = head(BP,6)
CC = head(CC,6)
MF = head(MF,6)
data = rbind(BP,CC,MF)
main.col = ontology.col[as.numeric(as.factor(data$ONTOLOGY))]
BgGene = as.numeric(sapply(strsplit(data$BgRatio,"/"),'[',1))
Gene = as.numeric(sapply(strsplit(data$GeneRatio,'/'),'[',1))
ratio = Gene/BgGene
logpvalue = -log(data$pvalue,10)
logpvalue.col = brewer.pal(n = 8, name = "Reds")
f = colorRamp2(breaks = c(0,2,4,6,8,10,15,20), colors = logpvalue.col)
BgGene.col = f(logpvalue)
df = data.frame(GO=data$ID,start=1,end=max(BgGene))
rownames(df) = df$GO
bed2 = data.frame(GO=data$ID,start=1,end=BgGene,BgGene=BgGene,BgGene.col=BgGene.col)
bed3 = data.frame(GO=data$ID,start=1,end=Gene,BgGene=Gene)
bed4 = data.frame(GO=data$ID,start=1,end=max(BgGene),ratio=ratio,col=main.col)
bed4$ratio = bed4$ratio/max(bed4$ratio)*9.5

pdf("GO.circlize.pdf",width=10,height=10)
par(omi=c(0.1,0.1,0.1,1.5))
circos.genomicInitialize(df,plotType="none")
circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
  sector.index = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  circos.text(mean(xlim), mean(ylim), sector.index, cex = 0.8, facing = "bending.inside", niceFacing = TRUE)
}, track.height = 0.08, bg.border = NA,bg.col = main.col)

for(si in get.all.sector.index()) {
  circos.axis(h = "top", labels.cex = 0.6, sector.index = si,track.index = 1,
              major.at=seq(0,max(BgGene),by=100),labels.facing = "clockwise")
}
f = colorRamp2(breaks = c(-1, 0, 1), colors = c("green", "black", "red"))
circos.genomicTrack(bed2, ylim = c(0, 1),track.height = 0.1,bg.border="white",
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      circos.genomicRect(region, value, ytop = 0, ybottom = 1, col = value[,2], 
                                         border = NA, ...)
                      circos.genomicText(region, value, y = 0.4, labels = value[,1], adj=0,cex=0.8,...)
                    })
circos.genomicTrack(bed3, ylim = c(0, 1),track.height = 0.1,bg.border="white",
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      circos.genomicRect(region, value, ytop = 0, ybottom = 1, col = '#BA55D3', 
                                         border = NA, ...)
                      circos.genomicText(region, value, y = 0.4, labels = value[,1], cex=0.9,adj=0,...)
                    })
circos.genomicTrack(bed4, ylim = c(0, 10),track.height = 0.35,bg.border="white",bg.col="grey90",
                    panel.fun = function(region, value, ...) {
                      cell.xlim = get.cell.meta.data("cell.xlim")
                      cell.ylim = get.cell.meta.data("cell.ylim")
                      for(j in 1:9) {
                        y = cell.ylim[1] + (cell.ylim[2]-cell.ylim[1])/10*j
                        circos.lines(cell.xlim, c(y, y), col = "#FFFFFF", lwd = 0.3)
                      }
                      circos.genomicRect(region, value, ytop = 0, ybottom = value[,1], col = value[,2], 
                                         border = NA, ...)
                      #circos.genomicText(region, value, y = 0.3, labels = value[,1], ...)
                    })
circos.clear()
middle.legend = Legend(
  labels = c('Number of Genes','Number of Select','Rich Factor(0-1)'),
  type="points",pch=c(15,15,17),legend_gp = gpar(col=c('pink','#BA55D3',ontology.col[1])),
  title="",nrow=3,size= unit(3, "mm")
)
circle_size = unit(1, "snpc")
draw(middle.legend,x=circle_size*0.42)
main.legend = Legend(
  labels = c("Biological Process","Cellular Component", "Molecular Function"),  type="points",pch=15,
  legend_gp = gpar(col=ontology.col), title_position = "topcenter",
  title = "ONTOLOGY", nrow = 3,size = unit(3, "mm"),grid_height = unit(5, "mm"),
  grid_width = unit(5, "mm")
)
logp.legend = Legend(
  labels=c('(0,2]','(2,4]','(4,6]','(6,8]','(8,10]','(10,15]','(15,20]','>=20'),
  type="points",pch=16,legend_gp=gpar(col=logpvalue.col),title="-log10(Pvalue)",
  title_position = "topcenter",grid_height = unit(5, "mm"),grid_width = unit(5, "mm"),
  size = unit(3, "mm")
)
lgd = packLegend(main.legend,logp.legend)
circle_size = unit(1, "snpc")
print(circle_size)
draw(lgd, x = circle_size*0.85, y=circle_size*0.55,just = "left")
dev.off()


KEGG
library(tidyverse)
library(devtools)
library(ggsankey)
library(ggplot2)
library(cols4all)
library(cowplot)
library(shinyjs)
library(kableExtra)
library(colorblindcheck)

KEGG=read.table("KEGG_enrichment_results.txt", header=T, sep="\t", check.names=F)
KEGG = KEGG[c(1:10),]
colnames(KEGG)

kegg= KEGG[,c("Description","Count","pvalue","GeneRatio")]
kegg$GeneRatio <- sapply(kegg$GeneRatio, function(x) {
  parts <- strsplit(x, "/")[[1]]
  as.numeric(parts[1]) / as.numeric(parts[2])
})

sankey= KEGG[,c("Description","geneID")]
sankey <- sankey %>%
  separate_rows(geneID, sep = "/")
kegg2 <- kegg[length(rownames(kegg)):1,] 
kegg2 <- kegg2 %>%
  mutate(ymax = cumsum(Count)) %>%
  mutate(ymin = ymax -Count) %>%
  mutate(label = (ymin + ymax)/2)

mytheme <- theme(axis.title = element_text(size = 13),
                 axis.text = element_text(size = 11),
                 axis.text.y = element_blank(),
                 axis.ticks.y = element_blank(),
                 legend.title = element_text(size = 13),
                 legend.text = element_text(size = 11))
p1 <- ggplot() +
  geom_point(data = kegg2,
             aes(x = -log10(pvalue),
                 y = label,
                 size = Count,
                 color = GeneRatio)) +
  scale_size_continuous(range=c(2,10)) +
  scale_y_continuous(expand = c(0,0.1),limits = c(0,32)) +
  scale_x_continuous(limits = c(0.1,ceiling(max(-log10(kegg2$pvalue)))+1)) +
  scale_colour_distiller(palette = "Reds", direction = 1) +
  labs(x = "-log10(Pvalue)",
       y = "") +
  theme_bw() +
  mytheme
p1
ggsave("KEGG_bubble_plot.pdf", p1, width = 5, height = 7) 

df <- sankey %>%
  make_long(geneID, Description)

df$node <- factor(df$node,levels = c(sankey$Description %>% unique()%>% rev(),
                                     sankey$geneID %>% unique() %>% rev()))

mycol <- c4a('rainbow_wh_rd',length(unique(df$node)))

p2 <- ggplot(df, aes(x = x,
                     next_x = next_x,
                     node = node,
                     next_node = next_node,
                     fill = node,
                     label = node)) +
  geom_sankey(flow.alpha = 0.5,
              flow.fill = 'grey',
              flow.color = 'grey80',
              node.fill = mycol, 
              smooth = 8,
              width = 0.08) +
  geom_sankey_text(size = 3.2,
                   color = "black")+
  theme_void() +
  theme(legend.position = 'none')+ 
  theme(plot.margin = unit(c(0,8,0,0),units="cm"))
p2
ggsave("KEGG_Sankey_plot.pdf", p2, width = 10, height = 7) 
p3 <- ggdraw() + draw_plot(p2) + draw_plot(p1, scale = 0.5, x = 0.50, y=-0.17, width=0.68, height=1.3)
p3
ggsave("KEGG_combined_bubble_sankey_plot.pdf", p3, width = 10, height = 7) 

#2.Machine Learning Analysis to Identify Core Genes
#(1) Lasso
library(tidyverse)    
library(broom)
library(glmnet)
x<-read.csv("Expression Matrix .csv",row.names = 1)
y<-read.csv("Classification.csv",row.names = 1)
x=as.matrix(x)
y=as.matrix(y)

set.seed(66666)

cvfit=cv.glmnet(x,y,type.measure = "mse",nfolds = 5,alpha=1)

plot(cvfit)

cvfit$lambda.min

c(cvfit$lambda.min, cvfit$lambda.1se)

lasso<-glmnet(x,y,family="binomial",alpha=1,nlambda = 100)

plot(lasso,xvar="lambda",label=F) 

coef <- coef(lasso, s = cvfit$lambda.min)
index <- which(coef != 0)
actCoef <- coef[index]
lassoGene=row.names(coef)[index]
geneCoef=cbind(Gene=lassoGene, Coef=actCoef)
write.csv(geneCoef,"LASSO_Core_Gene_Coefficient.csv")


#(2) RF
library(randomForest)
library(limma)
library(ggpubr)
set.seed(123)
data <-read.csv("Expression Matrix.csv",row.names = 1)
genelist<-read.csv("Sample_Classification_Label.csv",row.names = 1)
group <- as.character(genelist$type)
rf=randomForest(as.factor(group)~., data=data, ntree=500)
pdf(file="RandomForest_Error_Curve.pdf", width=6, height=6)
plot(rf, main="Random Forest Classification Error Curve", lwd=2)
dev.off()

optionTrees=which.min(rf$err.rate[,1])
optionTrees
rf2=randomForest(as.factor(group)~., data=data, ntree=optionTrees)
plot(rf2, main="Random Forest Optimized Error Curve", lwd=2)
importance=importance(x=rf2)
importance=as.data.frame(importance)
importance$size=rownames(importance)
importance=importance[,c(2,1)]

names(importance)=c("Gene","importance")

af=importance[order(importance$importance,decreasing = T),]
af=af[1:11,]

p1=ggplot(af, aes(x=reorder(Gene, importance), y=importance, fill=importance)) +  
  geom_bar(stat="identity", width=0.7) +  
  coord_flip() +  
  scale_fill_gradient(low=ggsci::pal_npg()(2)[1], high=ggsci::pal_npg()(2)[2]) +  
  labs(x="Gene", y="Importance", title="Top 15 Genes by Importance") +  
  theme_bw() +  
  theme(axis.text.x=element_text(angle=0, hjust=1), 
        axis.text.y=element_text(size=12),
        plot.title=element_text(hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())  
p1
pdf(file="RF_Gene_Importance_Barplot.pdf", width=6, height=6)
print(p1)
dev.off()

p2=ggdotchart(af, x = "Gene", y = "importance",
              color = "importance", 
              sorting = "descending",                       
              add = "segments",                             
              add.params = list(color = "lightgray", size = 2), 
              dot.size = 6,                        
              font.label = list(color = "white", size = 9,
                                vjust = 0.5),               
              ggtheme = theme_bw()         ,               
              rotate=TRUE                                       )
p3=p2+ geom_hline(yintercept = 0, linetype = 2, color = "lightgray")+
  gradient_color(palette =c(ggsci::pal_npg()(2)[2],ggsci::pal_npg()(2)[1])      ) 

p3
pdf(file="RF_Gene_Importance_Dotplot.pdf", width=6, height=6)
print(p3)
dev.off()

rfGenes=importance[order(importance[,"importance"], decreasing = TRUE),]
write.csv(rfGenes,"RandomForest_Gene_Feature_Importance.csv")

#(3) SVM
rm(list=ls())
a<-read.csv("machine_learning_input.csv",row.names = 1)
library(tidyverse)
library(glmnet)
source('msvmRFE.R')   
library(VennDiagram)
library(sigFeature)
library(e1071)
library(caret)
library(randomForest)
train<-a
input <- train

svmRFE(input, k = 5, halve.above = 100) 
nfold = 5
nrows = nrow(input)
folds = rep(1:nfold, len=nrows)[sample(nrows)]
folds = lapply(1:nfold, function(x) which(folds == x))
results = lapply(folds, svmRFE.wrap, input, k=5, halve.above=100) 
top.features = WriteFeatures(results, input, save=F) 
head(top.features)

write.csv(top.features,"SVM_RFE_Ranked_Feature_Genes.csv")


featsweep = lapply(1:11, FeatSweep.wrap, results, input) 
featsweep

save(featsweep,file = "SVM_RFE_Feature_Sweep_Result.RData")

no.info = min(prop.table(table(input[,1])))
errors = sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))

PlotErrors(errors, no.info=no.info) 

Plotaccuracy(1-errors,no.info=no.info) 

which.min(errors) 
top<-top.features[1:which.min(errors), "FeatureName"]
write.csv(top,"SVM_RFE_Optimal_Core_Genes.csv")


#3. Construction and Evaluation of a Diagnostic Nomogram
#(1) Nomogram and Calibration curves
library(rms)
library(rmda)
library(Hmisc)
library(glmnet)
library(pROC)
library(ggsci)
rt=read.csv("nomogram_model_input_dataset.csv",row.names = 1)

ddist=datadist(rt)
options(datadist="ddist")

lrmModel=lrm(Type~IRF1+TNFAIP3, data=rt, x=T, y=T)

nomo=nomogram(lrmModel, fun=plogis,
              fun.at=c(0.0001,0.1,0.3,0.6,0.9,0.99),
              lp=F, funlabel="Nomogram")


pdf("Diagnostic_Nomogram.pdf", width=10, height=6)
plot(nomo, cex.axis=0.7)

dev.off()

nomoRisk=predict(lrmModel, type="fitted")
outTab=cbind(rt, Nomogram=nomoRisk)
outTab=rbind(id=colnames(outTab), outTab)
write.table(outTab, file="Nomogram_Risk_Score_Result.txt", sep="\t", quote=F, col.names=F)

cali=calibrate(lrmModel, method="boot", B=1000)

plot(cali,
     xlab="Model Predicted probability",
     ylab="Actual Observed probability", sub=F)

y=outTab[,"Type"]

bioCol=pal_simpsons(palette=c("springfield"), alpha=1)(length(2:ncol(outTab)))
aucText=c()  
k=0  
for(x in colnames(outTab)[2:ncol(outTab)]){
  k=k+1
  roc1=roc(y, as.numeric(outTab[,x]))     
  if(k==1){
    plot(roc1, print.auc=F, col=bioCol[k], legacy.axes=T, main="")
    aucText=c(aucText, paste0(x,", AUC=",sprintf("%.3f",roc1$auc[1])))
  }else{
    
    plot(roc1, print.auc=F, col=bioCol[k], legacy.axes=T, main="ROC Curves of Predictive Indicators", add=TRUE)
   
    aucText=c(aucText, paste0(x,", AUC=",sprintf("%.3f",roc1$auc[1])))
  }
}

legend("bottomright", aucText, lwd=1, bty="n", col=bioCol, cex=0.9)
dev.off()


#(2)ROC curves of nomoscore and the two core genes, respectively.
library(pROC) 
data <- read.csv("Nomogram_Risk_Score.csv")  

true_labels <- data$Type
nomogram_scores <- data$Nomogram

pdf("Nomogram_RiskScore_ROC_Curve.pdf", width = 6, height = 6)

roc_obj <- roc(true_labels, nomogram_scores)
plot(roc_obj, 
     col = "red",                 
     lwd = 2,                     
     xlab = "Specificity",     
     ylab = "Sensitivity",        
     print.auc = TRUE)            
dev.off()




library(pROC)   
rt=read.csv("IRF1_gene_ROC_dataset.csv",row.names = 1)

y=colnames(rt)[1]
x="Irf1" 

ROC=roc(rt[,y], as.vector(rt[,x]))

pdf("IRF1_Single_Gene_ROC.pdf",width=5,height=5)

plot(ROC, print.auc=TRUE, col="red")

dev.off()


library(pROC)   
rt=read.csv("TNFAIP3_gene_ROC_dataset.csv",row.names = 1)

y=colnames(rt)[1]
x="Tnfaip3" 

ROC=roc(rt[,y], as.vector(rt[,x]))

pdf("TNFAIP3_Single_Gene_ROC.pdf",width=5,height=5)

plot(ROC, print.auc=TRUE, col="red")

dev.off()




#(3)  ROC curve of the two core genes.
rt=read.csv("Two_Core_Genes_ROC_Dataset.csv",row.names = 1)  
y=colnames(rt)[1]

bioCol=c("red","blue","green","yellow")
if(ncol(rt)>4){
  bioCol=rainbow(ncol(rt))}

pdf("Two_Genes_Combined_ROC_Curve.pdf",width=5,height=5)

roc1=roc(rt[,y], as.vector(rt[,2]))
geneText=c(colnames(rt)[2])
plot(roc1, col=bioCol[1])

for(i in 3:ncol(rt)){
  roc1=roc(rt[,y], as.vector(rt[,i]))
  lines(roc1, col=bioCol[i-1])
  geneText=c(geneText, colnames(rt)[i])
}

legend("bottomright", geneText,lwd=2,bty="n",col=bioCol[1:(ncol(rt)-1)])
dev.off()


#4. Immune Infiltration Analysis and Correlation Analysis
#(1)Immune Infiltration Analysis

ssgsea_data <- read.csv("ssGSEA_immune_infiltration_score.csv", row.names = 1, check.names = FALSE)

y <- as.data.frame(t(ssgsea_data))

group <- c(rep("Control", 25), rep("SLI", 26))
y$group <- group

library(tidyr)
data_long <- pivot_longer(
  data = y,
  cols = -group,  
  names_to = "celltype",
  values_to = "proportion"
)

library(ggpubr)
library(ggplot2)
library(rstatix)

max_proportion <- max(data_long$proportion, na.rm = TRUE)

stat.test <- data_long %>%
  group_by(celltype) %>%
  t_test(proportion ~ group, ref.group = "Control") %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "celltype", dodge = 0.8)

stat.test$y.position <- max_proportion * 1.15

ggboxplot(
  data = data_long,
  x = "celltype",
  y = "proportion",
  color = "black",
  fill = "group",
  palette = c("#1C3EDF", "#DF1C26"),
  xlab = "ssGSEA",
  ylab = "Score",
  outlier.shape = 20
) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = "top"
  ) +
  
  stat_pvalue_manual(
    stat.test, 
    label = "p.adj.signif",  
 
    bracket.size = 0,  
    tip.length = 0,      
    size = 4,              
    step.increase = 0      
  ) +
 
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
  scale_fill_manual(
    values = c("#1C3EDF", "#DF1C26"),
    name = "group",
    labels = c("Control", "SLI")
  )
ggsave("boxplot.pdf", width = 15, height = 8) 
ggsave("boxplot.png", width = 15, height = 8, dpi = 300)  


#(2)Correlation Analysis
sig_gene <- c( "Irf1","Tnfaip3")

exprSet <- read.csv("ComBat_normalized_expression.csv", header = TRUE, row.names = 1, stringsAsFactors = FALSE)

ssgsea_data <- read.csv(
  "ssGSEA_immune_infiltration_score.csv",
  header = TRUE,    
  row.names = 1,    
  stringsAsFactors = FALSE  
)
x <- t(exprSet)
x <- x[,sig_gene]
y <- t(ssgsea_data)
library(psych)
d <- corr.test(x,y,use="complete",method = 'spearman')
r <- d$r
p <- d$p

library(ggcorrplot)
ggcorrplot(t(d$r), 
           show.legend = T, 
           digits = 2,  sig.level = 0.05,
           insig = 'blank',lab = T)+coord_flip() 

library(pheatmap)
library(reshape2)

if (!is.null(p)){
  ssmt <- p< 0.001
  p[ssmt] <-'***'
  smt <- p >0.001& p < 0.01
  p[ssmt] <-'**'
  smt <- p >0.01& p <0.05
  p[smt] <- '*'
  p[!ssmt&!smt]<- ''
} else {
  p <- F
}
mycol<-colorRampPalette(c("blue","white","tomato"))(100)

pheatmap(r,scale = "none",cluster_row = T, cluster_col = T, border=NA,
         display_numbers = p,fontsize_number = 10, number_color = "white",
         cellwidth = 9, cellheight =13,color=mycol)


#5. Levels of Expression for core Genes
#(1)Levels of Expression for IRF1
library(pcutils)
library(ggstatsplot)
library(PMCMRplus) 
library(ggplot2) 
library(cowplot) 

set.seed(123)

data <- read.table("Expression Matrix IRF1.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T)

{
 
  plist = list()
  for (i in 1:9) {
    plist[[i]] = group_box(tab = data[c("IRF1")], group = "group", metadata = data, 
                           mode = i,p_value1 = "wilcox.test") +
      ggtitle(paste0("mode", i)) +
      theme_classic() + 
      theme(legend.position = "none")
  }
  
  plot_grid(plotlist = plist, ncol = 3)
}

plot(plist[[6]])


#(2)Levels of Expression for TNFAIP3
library(pcutils)
library(ggstatsplot)
library(PMCMRplus) 
library(ggplot2) 
library(cowplot) 

set.seed(123)

data <- read.table("Expression Matrix TNFAIP3.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T)

{
  
  plist = list()
  for (i in 1:9) {
    plist[[i]] = group_box(tab = data[c("TNFAIP3")], group = "group", metadata = data, 
                           mode = i,p_value1 = "wilcox.test") +
      ggtitle(paste0("mode", i)) +
      theme_classic() + 
      theme(legend.position = "none")
  }
  
  plot_grid(plotlist = plist, ncol = 3)
}

plot(plist[[6]])





