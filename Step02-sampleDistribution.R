
rm(list = ls())
options(stringsAsFactors = F)

# 加载原始表达的数据
lname <- load(file = "data/Step01-airwayData.Rdata")
lname

exprSet <- express_cpm
exprSet[1:6,1:6]

# 样本表达总体分布-箱式图
library(ggplot2)
# 构造绘图数据
dat <- data.frame(expression=c(exprSet),sample=rep(colnames(exprSet),each=nrow(exprSet)))
head(dat)

p <- ggplot(data = dat,aes(x=sample,y=expression,fill=sample))
p1 <- p + geom_boxplot() + theme(axis.text.x = element_text(angle = 90))  + xlab(NULL) + ylab("log2(CPM+1)")
p1

# 保存图片
pdf(file = "sample_cor/sample_boxplot.pdf",width = 6,height = 8)
print(p1)
dev.off()

# 样本表达总体分布-小提琴图
p2 <- p + geom_violin() + theme(axis.text = element_text(size = 12),axis.text.x = element_text(angle = 90)) + xlab(NULL) + ylab("log2(CPM+1)")
p2

# 保存图片
pdf(file = "sample_cor/sample_violin.pdf",width = 6,height = 8)
print(p2)
dev.off()

# 样本表达总体分布-概率密度分布图
m <- ggplot(data=dat, aes(x=expression))
p3 <- m +  geom_density(aes(fill=sample, colour=sample),alpha = 0.2) + xlab("log2(CPM+1)")
p3

# 保存图片
pdf(file = "sample_cor/sample_density.pdf",width = 7,height = 8)
print(p3)
dev.off()



