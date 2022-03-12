
# 魔幻操作，一键清空
rm(list = ls())  
options(stringsAsFactors = F)

# 加载数据并检查
lname <- load(file = 'data/Step01-airwayData.Rdata')
lname

dat <- express_cpm
dat[1:6,1:6]
dim(dat)


## 1.样本之间的相关性-层次聚类树----
sampleTree <- hclust(dist(t(dat)), method = "average")
plot(sampleTree)



## 2.样本之间的相关性-PCA----
# 第一步，数据预处理
dat <- as.data.frame(t(dat))
dat$group_list <- group_list

# 第二步，绘制PCA图
library(FactoMineR)
library(factoextra)

# 画图仅需要数值型数据，去掉最后一列的分组信息
dat_pca <- PCA(dat[,-ncol(dat)], graph = FALSE)
class(dat_pca)

p <- fviz_pca_ind(dat_pca,
                  geom.ind = "point", # 只显示点，不显示文字
                  col.ind = dat$group_list, # 用不同颜色表示分组
                  palette = c("#00AFBB", "#E7B800"),
                  addEllipses = T ,# 是否圈起来
                  legend.title = "Groups")
p



## 3.样本之间的相关性-cor----
# 利用绝对中位差mad/标准差sd统计学方法进行数据异常值检测
# 将表达量的绝对中位差mad从大到小排列取前500的结果
dat <- express_cpm
tmp <- sort(apply(dat,1, mad),decreasing = T)[1:500]
exprSet <-dat[names(tmp),]

# 使用500个基因的表达量来做相关性图
library(corrplot)
dim(exprSet)

# 计算相关性
M <- cor(exprSet)
g <- corrplot(M,order = "AOE",addCoef.col = "white")

corrplot(M,order = "AOE",type="upper",tl.pos = "d")
corrplot(M,add=TRUE, type="lower", method="number",order="AOE",diag=FALSE,tl.pos="n", cl.pos="n")

# 绘制样本相关性的热图
anno <- data.frame(sampleType=group_list)
rownames(anno) <- colnames(exprSet)
rowanno <- data.frame(sampleType=group_list)
rownames(rowanno) <- colnames(exprSet)
anno
p <- pheatmap::pheatmap(M,display_numbers = T,annotation_col = anno,fontsize = 11,cellheight = 28,cellwidth = 28,annotation_row = rowanno)
p

pdf(file = "data/sample_cor/cor.pdf")
print(p)
dev.off()



