
rm(list = ls())
options(stringsAsFactors = F)

# 读取3个软件的差异分析结果
load(file = "data/Step03-limma_voom_nrDEG.Rdata")
load(file = "data/Step03-DESeq2_nrDEG.Rdata")
load(file = "data/Step03-edgeR_nrDEG.Rdata")
ls()

# 提取所有差异表达的基因名
limma_sigGene <- DEG_limma_voom[DEG_limma_voom$regulated!="normal",1]
edgeR_sigGene <- DEG_edgeR[DEG_edgeR$regulated!="normal",1]
DESeq2_sigGene <- DEG_DESeq2[DEG_DESeq2$regulated!="normal",1]

data <- list(limma=limma_sigGene,edgeR=edgeR_sigGene,DESeq2=DESeq2_sigGene)

library(VennDiagram)
# 设置颜色
col <- c('#0099CC','#FF6666','#FFCC99')
venn.diagram(data,lwd=1,lty=1,col=col,fill=col,cat.col=col,cat.cex = 1.8,rotation.degree = 0,cex=1.5,
             alpha = 0.5,reverse=TRUE,width=4000,height = 4000,resolution =600,margin=0.2,
             filename="deg_analysis/DEG_venn.png",imagetype="png")



