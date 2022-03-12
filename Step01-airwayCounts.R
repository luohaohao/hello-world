
# 魔幻操作，一键清空
rm(list = ls()) 
options(stringsAsFactors = F)

# 加载airway数据集并转换为表达矩阵
library(airway,quietly = T)
data(airway)
class(airway)

rawcount <- assay(airway)
colnames(rawcount)

# 查看表达谱
rawcount[1:4,1:4]

# 去除前的基因表达矩阵情况
dim(rawcount)

# 获取分组信息
group_list <- colData(airway)$dex
group_list

# 过滤在至少在75%的样本中都有表达的基因
keep <- rowSums(rawcount>0) >= floor(0.75*ncol(rawcount))
table(keep)

filter_count <- rawcount[keep,]
filter_count[1:4,1:4]
dim(filter_count)

# 加载edgeR包计算counts per millio(cpm) 表达矩阵,并对结果取log2值
library(edgeR)
express_cpm <- log2(cpm(filter_count)+1)
express_cpm[1:6,1:6]

# 保存表达矩阵和分组结果
save(filter_count,express_cpm,group_list,file = "data/Step01-airwayData.Rdata")


