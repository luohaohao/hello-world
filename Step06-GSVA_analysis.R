
## 对MigDB中的基因集 做GSVA分析。

rm(list = ls())
options(stringsAsFactors = F)


## 读取基因表达矩阵
lnames <- load(file = "data/Step01-airwayData.Rdata")
lnames

## 将表达矩阵的ensembl ID换成gene symbol
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)

library(clusterProfiler)
id2symbol <- bitr(rownames(express_cpm), fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
head(id2symbol)
symbol <- rep("NA",time=nrow(express_cpm))
symbol[match(id2symbol[,1],rownames(express_cpm))] <- id2symbol[,2]
#express_cpm <- cbind(rownames(express_cpm),symbol,express_cpm)
#colnames(express_cpm)[1] <- "GeneID"
#express_cpm[1:4,1:4]
rownames(express_cpm) <- symbol
exprSet <- express_cpm
# 过滤掉没有匹配上symbol的行
exprSet <- express_cpm[express_cpm[,2]!="NA",]
dim(exprSet)

# 处理多对1的情况
library(limma)
# rownames(exprSet) <- exprSet[,2]
# exprSet <- exprSet[,-c(1:2)]
# dim(exprSet)

exprSet <- limma::avereps(exprSet,ID=rownames(exprSet))
dim(exprSet)
exprSet[1:4,1:4]


## 将表达矩阵转换成通路矩阵
library(GSEABase)
library(GSVA)
geneset <- getGmt("data/v7.1/c2.cp.kegg.v7.1.symbols.gmt")
class(geneset)

es_max <- gsva(exprSet, geneset, mx.diff=F, verbose=T, parallel.sz=8)
dim(es_max)


## 做差异分析
library(limma)

design <- model.matrix(~0+factor(group_list))
colnames(design) <- levels(factor(group_list))
rownames(design) <- colnames(es_max)
design
    
contrast.matrix <- makeContrasts("trt-untrt",levels = design)
contrast.matrix

fit <- lmFit(es_max,design)
fit2 <- contrasts.fit(fit, contrast.matrix) 
fit2 <- eBayes(fit2)

tempOutput <- topTable(fit2, coef=1, n=Inf)
nrDEG <- na.omit(tempOutput) 
head(nrDEG)

# 得到显著通路
nrDEG_sig <- nrDEG[nrDEG$P.Value<0.01 & abs(nrDEG$logFC) > 0.3,]


## 绘制热图
library(pheatmap)
library(stringr)
data <- es_max[match(rownames(nrDEG_sig),rownames(es_max)),]
rownames(data) <- gsub("KEGG_","",rownames(data))

anno <- data.frame(group=group_list)
rownames(anno) <- colnames(data)
p <- pheatmap::pheatmap(data, fontsize_row =10,height = 11,annotation_col = anno,show_colnames = F)




