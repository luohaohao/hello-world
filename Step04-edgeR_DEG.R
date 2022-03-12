
# 参考链接：https://www.biostars.org/p/110861/

rm(list = ls())
options(stringsAsFactors = F)

# 读取基因表达矩阵信息并查看分组信息和表达矩阵数据
lname <- load(file = "data/Step01-airwayData.Rdata")
lname

exprSet <- filter_count
dim(exprSet)
exprSet[1:6,1:6]
table(group_list)

# 加载包
library(edgeR)

# 假设数据符合正态分布，构建线性模型。0代表x线性模型的截距为0
design <- model.matrix(~0+factor(group_list))
rownames(design) <- colnames(exprSet)
colnames(design) <- levels(factor(group_list))
design

# 构建edgeR的DGEList对象
DEG <- DGEList(counts=exprSet,group=factor(group_list))

# 增加一列$norm.factors
DEG$samples$lib.size <- colSums(DEG$counts)
DEG$samples

# 归一化基因表达分布
DEG <- calcNormFactors(DEG)

# 计算线性模型的参数
DEG <- estimateGLMCommonDisp(DEG,design)
DEG <- estimateGLMTrendedDisp(DEG, design)
DEG <- estimateGLMTagwiseDisp(DEG, design)

# 拟合线性模型
fit <- glmFit(DEG, design)

# 进行差异分析，1,-1意味着前比后
lrt <- glmLRT(fit, contrast=c(1,-1)) 

# 提取过滤差异分析结果
DEG_edgeR <- as.data.frame(topTags(lrt, n=nrow(DEG)))
head(DEG_edgeR)

# 筛选上下调，设定阈值
fc_cutoff <- 1.5
fdr <- 0.05

DEG_edgeR$regulated <- "normal"

loc_up <- intersect(which(DEG_edgeR$logFC>log2(fc_cutoff)),which(DEG_edgeR$FDR<fdr))
loc_down <- intersect(which(DEG_edgeR$logFC< (-log2(fc_cutoff))),which(DEG_edgeR$FDR<fdr))

DEG_edgeR$regulated[loc_up] <- "up"
DEG_edgeR$regulated[loc_down] <- "down"

table(DEG_edgeR$regulated)

# 添加一列gene symbol
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)

library(clusterProfiler)
id2symbol <- bitr(rownames(DEG_edgeR), fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Hs.eg.db )
nrow(id2symbol)
nrow(DEG_edgeR)
symbol <- "NA"
symbol[match(id2symbol[,1],rownames(DEG_edgeR))] <- id2symbol[,2]
DEG_edgeR <- cbind(rownames(DEG_edgeR),symbol,DEG_edgeR)
colnames(DEG_edgeR)[1] <- "GeneID"

# 保存
DEG_edgeR_up <-  DEG_edgeR[DEG_edgeR$regulated=="up",]
write.table(DEG_edgeR,"deg_analysis/DEG_edgeR_all.xls",row.names = F,sep="\t",quote = F)

## 取表达差异倍数和p值,矫正后的pvalue
colnames(DEG_edgeR)
DEG_edgeR <- DEG_edgeR[,c(1,2,3,6,7,8)]
save(DEG_edgeR, file = "data/Step03-edgeR_nrDEG.Rdata")


## 检查是否上下调设置错了
# 挑选一个差异表达基因
head(DEG_edgeR)

exp <- c(t(express_cpm[match("ENSG00000152583",rownames(express_cpm)),]))
test <- data.frame(value=exp,group=group_list)
library(ggplot2)
ggplot(data=test,aes(x=group,y=value,fill=group)) + geom_boxplot()


