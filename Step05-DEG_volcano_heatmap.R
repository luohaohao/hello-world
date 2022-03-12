
rm(list = ls())
options(stringsAsFactors = F)

# 加载原始表达矩阵
load(file = "data/Step01-airwayData.Rdata")

# 读取3个软件的差异分析结果
load(file = "data/Step03-limma_voom_nrDEG.Rdata")
load(file = "data/Step03-DESeq2_nrDEG.Rdata")
load(file = "data/Step03-edgeR_nrDEG.Rdata")
ls()

# 根据需要修改DEG的值
data <- DEG_limma_voom
colnames(data)


# 绘制火山图
library(ggplot2)
colnames(data)
p <- ggplot(data=data, aes(x=logFC, y=-log10(adj.P.Val),color=regulated)) + 
     geom_point(alpha=0.5, size=1.8) + theme_set(theme_set(theme_bw(base_size=20))) + 
     xlab("log2FC") + ylab("-log10(FDR)") +scale_colour_manual(values = c('blue','black','red'))
p

# 提取所有差异表达的基因名
limma_sigGene <- DEG_limma_voom[DEG_limma_voom$regulated!="normal",1]
edgeR_sigGene <- DEG_edgeR[DEG_edgeR$regulated!="normal",1]
DESeq2_sigGene <- DEG_DESeq2[DEG_DESeq2$regulated!="normal",1]

library(pheatmap)
# 绘制热图
dat <- express_cpm[match(limma_sigGene,rownames(express_cpm)),]
dat[1:4,1:4]
group <- data.frame(group=group_list)
rownames(group)=colnames(dat)
p <- pheatmap(dat,scale = "row",show_colnames =T,show_rownames = F, cluster_cols = T,annotation_col=group,main = "limma's DEG",treeheight_row = 0,treeheight_col = 0) 

pdf()
png()
tiff()




