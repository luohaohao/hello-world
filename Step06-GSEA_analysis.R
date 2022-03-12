
# 清空当前环境变量
rm(list = ls())
options(stringsAsFactors = F)

# 加载包
library(GSEABase)
library(clusterProfiler)

# 加载数据
lnames <- load("data/Step03-limma_voom_nrDEG.Rdata")
lnames

DEG <- DEG_limma_voom

## 构造GSEA分析数据
# 去掉没有配对上symbol的行
DEG <- DEG[!DEG$symbol=="NA",]
geneList <- DEG$logFC
names(geneList) <- DEG$symbol
geneList <- sort(geneList,decreasing = T)
head(geneList)


# 选择gmt文件（MigDB中的全部基因集）
geneset <- read.gmt("data/v7.1/c2.cp.kegg.v7.1.symbols.gmt")
egmt <- GSEA(geneList, TERM2GENE=geneset, verbose=T,pvalueCutoff = 1)

kegg_gsea <- as.data.frame(egmt@result)
colnames(kegg_gsea)
write.table(kegg_gsea,"gsea_kegg_fc.xls",row.names = F,sep="\t",quote = F)

library(enrichplot)
gseaplot2(egmt, "KEGG_GAP_JUNCTION",title = "KEGG_GAP_JUNCTION",pvalue_table = T,color = "red")





