
rm(list = ls())
options(stringsAsFactors = F)

# 读取基因表达矩阵信息
lname <- load(file = "data/Step01-airwayData.Rdata")
lname 

# 查看分组信息和表达矩阵数据
exprSet <- filter_count
dim(exprSet)
exprSet[1:6,1:6]
table(group_list)

# 加载包
library(DESeq2)

# 第一步，构建DESeq2的DESeq对象
colData <- data.frame(row.names=colnames(exprSet),group_list=group_list)
dds <- DESeqDataSetFromMatrix(countData = exprSet,colData = colData,design = ~ group_list)

# 第二步，进行差异表达分析
dds2 <- DESeq(dds)

# 提取差异分析结果，trt组对untrt组的差异分析结果
tmp <- results(dds2,contrast=c("group_list","trt","untrt"))
DEG_DESeq2 <- as.data.frame(tmp[order(tmp$padj),])
head(DEG_DESeq2)

# 去除差异分析结果中包含NA值的行
DEG_DESeq2 = na.omit(DEG_DESeq2)

# 筛选上下调，设定阈值
fc_cutoff <- 2
fdr <- 0.05

DEG_DESeq2$regulated <- "normal"

loc_up <- intersect(which(DEG_DESeq2$log2FoldChange>log2(fc_cutoff)),which(DEG_DESeq2$padj<fdr))
loc_down <- intersect(which(DEG_DESeq2$log2FoldChange< (-log2(fc_cutoff))),which(DEG_DESeq2$padj<fdr))

DEG_DESeq2$regulated[loc_up] <- "up"
DEG_DESeq2$regulated[loc_down] <- "down"

table(DEG_DESeq2$regulated)

# 添加一列gene symbol
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)

library(clusterProfiler)
id2symbol <- bitr(rownames(DEG_DESeq2), fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Hs.eg.db )
head(id2symbol)

symbol <- rep("NA",time=nrow(DEG_DESeq2))
symbol[match(id2symbol[,1],rownames(DEG_DESeq2))] <- id2symbol[,2]
DEG_DESeq2 <- cbind(rownames(DEG_DESeq2),symbol,DEG_DESeq2)
colnames(DEG_DESeq2)[1] <- "GeneID"
head(DEG_DESeq2)

# 保存
write.table(DEG_DESeq2,"deg_analysis/DEG_DESeq2_all.xls",row.names = F,sep="\t",quote = F)

## 取表达差异倍数和p值,矫正后的pvalue并保存
colnames(DEG_DESeq2)
DEG_DESeq2 <- DEG_DESeq2[,c(1,2,4,7,8,9)]
save(DEG_DESeq2, file = "data/Step03-DESeq2_nrDEG.Rdata")


## 检查是否上下调设置错了
# 挑选一个差异表达基因
head(DEG_DESeq2)

exp <- c(t(express_cpm[match("ENSG00000152583",rownames(express_cpm)),]))
test <- data.frame(value=exp,group=group_list)
library(ggplot2)
ggplot(data=test,aes(x=group,y=value,fill=group)) + geom_boxplot()

