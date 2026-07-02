

setwd("D:\\Desktop\\data diff")

library(dplyr)
library(data.table)
library(tinyarray)

########################
# 1. DATA MERGE
########################

data1<-read.csv("GSE15379.csv",row.names=1)
data2<-read.csv("GSE23767.csv",row.names=1)
data3<-read.csv("GSE52474.csv",row.names=1)
data4<-read.csv("GSE60088.csv",row.names=1)

group1=data.frame(row.names=colnames(data1),sample=colnames(data1),DataSet="GSE15379")
group2=data.frame(row.names=colnames(data2),sample=colnames(data2),DataSet="GSE23767")
group3=data.frame(row.names=colnames(data3),sample=colnames(data3),DataSet="GSE52474")
group4=data.frame(row.names=colnames(data4),sample=colnames(data4),DataSet="GSE60088")

com_group=rbind(group1,group2,group3,group4)

data1$gene<-rownames(data1)
data2$gene<-rownames(data2)
data3$gene<-rownames(data3)
data4$gene<-rownames(data4)

combine_data<-Reduce(function(x,y) inner_join(x,y,by="gene"),
                      list(data1,data2,data3,data4))

norm_data<-combine_data %>% select(-gene)
rownames(norm_data)<-combine_data$gene

########################
# 2. PCA + COMBAT
########################

library(FactoMineR)
library(factoextra)
library(sva)

stopifnot(all(colnames(norm_data) %in% rownames(com_group)))

mod <- model.matrix(~ com_group$DataSet)

combat_data <- ComBat(dat=norm_data,
                      batch=com_group$DataSet,
                      mod=mod)

write.csv(combat_data,"ComBat_normalized_expression.csv")

########################
# 3. DEG ANALYSIS
########################

library(limma)

exp<-read.csv("ComBat_normalized_expression.csv",row.names=1)
fen<-read.csv("sample_group_info.csv",row.names=1)

fen <- fen[colnames(exp), , drop=FALSE]
stopifnot(all(rownames(fen)==colnames(exp)))

group_list<-factor(fen$group,levels=c("Control","SLI"))
design<-model.matrix(~group_list)

fit<-lmFit(exp,design)
fit<-eBayes(fit)

DEG<-topTable(fit,coef=2,number=Inf)
DEG<-na.omit(DEG)

DEG$type<-ifelse(DEG$adj.P.Val<0.05 & DEG$logFC>0.5,"Up",
                 ifelse(DEG$adj.P.Val<0.05 & DEG$logFC< -0.5,"Down","NOT"))

write.csv(DEG,"classified_limma_DE_results.csv")

########################
# 4. VOLCANO
########################

library(ggrepel)

DEG$label<-NA

target_genes<-c("Ccl2","Ccl7","Cxcl1","Socs3","Ch25h","Nfkbiz",
"Fcgr3","Cxcl10","Slc15a3","Tnfaip2","Nuak2","Irg1","Cd14",
"AA467197","Fcer1g","Tnfaip3","Pcdh17","Csf2rb","Gadd45g","Map3k8")

DEG$label[rownames(DEG)%in%target_genes]<-rownames(DEG)[rownames(DEG)%in%target_genes]

########################
# 5. HEATMAP
########################

library(pheatmap)

diff_gene<-DEG[DEG$type!="NOT",]
diff_gene<-diff_gene[order(diff_gene$logFC,decreasing=TRUE),]

rt<-exp[rownames(diff_gene),]

ann<-data.frame(Type=fen$group)
rownames(ann)<-rownames(fen)

ann<-ann[colnames(rt),,drop=FALSE]

pheatmap(rt,annotation_col=ann,scale="row")

########################
# 6. GO
########################

library(clusterProfiler)
library(org.Mm.eg.db)

rt<-read.table("uniCox.txt",header=TRUE)
genes<-unique(rt[,1])

entrezIDs<-mget(genes,org.Mm.egSYMBOL2EG,ifnotfound=NA)
entrezIDs<-na.omit(as.character(entrezIDs))

kk<-enrichGO(entrezIDs,OrgDb=org.Mm.eg.db,
             pvalueCutoff=1,qvalueCutoff=1,ont="all",
             readable=TRUE)

GO<-as.data.frame(kk)
GO<-GO[GO$pvalue<0.05,]

########################
# 7. LASSO
########################

library(glmnet)

x<-as.matrix(read.csv("Expression Matrix .csv",row.names=1))
y<-read.csv("Classification.csv",row.names=1)

set.seed(1)

cvfit<-cv.glmnet(x,y$V1,
                 family="binomial",
                 type.measure="class",
                 nfolds=5)

lasso<-glmnet(x,y$V1,family="binomial",alpha=1)

coef<-coef(lasso,s=cvfit$lambda.min)
index<-which(coef!=0)

geneCoef<-data.frame(
  Gene=rownames(coef)[index],
  Coef=as.numeric(coef[index])
)

write.csv(geneCoef,"LASSO_Core_Gene_Coefficient.csv")

########################
# 8. RF
########################

library(randomForest)

data<-read.csv("Expression Matrix.csv",row.names=1)
genelist<-read.csv("Sample_Classification_Label.csv",row.names=1)

group<-genelist$type[colnames(data)]

rf<-randomForest(as.factor(group)~.,data=data,ntree=500)

imp<-importance(rf)
imp<-as.data.frame(imp)
imp$Gene<-rownames(imp)

write.csv(imp,"RandomForest_Gene_Feature_Importance.csv")

########################
# 9. NOMOGRAM
########################

library(rms)

rt<-read.csv("nomogram_model_input_dataset.csv",row.names=1)

ddist<-datadist(rt)
options(datadist="ddist")

lrmModel<-lrm(Type~IRF1+TNFAIP3,data=rt,x=TRUE,y=TRUE)

nomoRisk<-predict(lrmModel,type="fitted")

outTab<-cbind(rt,Nomogram=nomoRisk)

write.table(outTab,"Nomogram_Risk_Score_Result.txt",sep="\t")

########################
# 10. ROC
########################

library(pROC)

roc1<-roc(rt$Type,nomoRisk)
plot(roc1)

########################
# 11. ssGSEA
########################

ssgsea_data<-read.csv("ssGSEA_immune_infiltration_score.csv",row.names=1)

y<-as.data.frame(t(ssgsea_data))
group<-fen$group

y$group<-group

########################
# 12. CORRELATION
########################

library(psych)

exprSet<-read.csv("ComBat_normalized_expression.csv",row.names=1)
ssgsea_data<-read.csv("ssGSEA_immune_infiltration_score.csv",row.names=1)

x<-t(exprSet)[,c("Irf1","Tnfaip3")]
y<-t(ssgsea_data)

common<-intersect(rownames(x),rownames(y))
x<-x[common,,drop=FALSE]
y<-y[common,,drop=FALSE]

d<-corr.test(x,y,use="complete")

r<-d$r
p<-d$p

library(pheatmap)
pheatmap(r)
