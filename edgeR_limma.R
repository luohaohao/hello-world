Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())

library(dplyr)
library(tidyverse)

###########################################################################################



# 对于基因芯片的差异表达分析而言，由于普遍认为其数据是服从正态分布，
# 因此差异表达分析无非就是用t检验和或者方差分析应用到每一个基因上。高通量一次性找的基因多，
# 于是就需要对多重试验进行矫正，控制假阳性。目前在基因芯片的分析用的最多的就是limma。

# 但是，高通量测序(HTS)的read count普遍认为是服从泊松分布（当然有其他不同意见），
# 不可能直接用正态分布的t检验和方差分析。 当然我们可以简单粗暴的使用对于的非参数检验的方法，
# 但是统计力不够，结果的p值矫正之估计一个差异基因都找不到。老板花了一大笔钱，结果却说没有差异基因，
# 是个负结果，于是好几千经费打了水漂，他肯定是不乐意的。因此，还是得要用参数检验的方法，
# 于是就要说到方差分析和线性模型之间的关系了。
###DESeq2  DEseq2针对有生物学重复的样本。
  
###数据来源https://4va.github.io/biodatasci/r-rnaseq-airway.html#data_needed  
counts <- read.csv("airway_scaledcounts.csv")
metadata <-  read.csv("airway_metadata.csv")
library(DESeq2)

dds <- DESeqDataSetFromMatrix(countData=counts,
                              colData=metadata,
                              design=~dex,
                              tidy=TRUE)
###其中countData存放counts数据，colData存放样本信息的数据，design就是实验设计。
#####数据预处理
dds1 <- dds [ rowSums(counts(dds)) > 1, ]
###虽然DESeq2会自动屏蔽那些低count的基因，但是剔除那些几乎不存在基因的部分能够提高运行速度。

##方差齐性转换
rld <- rlog(dds, blind = FALSE)
head(rld@assays@data@listData[[1]], 3)
vsd <- vst(dds, blind = FALSE)  # 一般大数据集
head(assay(vsd), 3)

##样本距离
sampleDists <- dist(t(assay(rld)))
library("pheatmap")
library("RColorBrewer")
sampleDistMatrix <- as.matrix( sampleDists )

colors <- colorRampPalette( rev(brewer.pal(8, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

sampleTree <- hclust(dist(t(assay(rld))), method = "average")
plot(sampleTree)
##PCA分析
dat <- as.data.frame(t(assay(rld)))
dat$group_list <- metadata$dex

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
##DESeq2提供了专门的方法用于作图
# plotPCA(rld,intgroup=c("dex"))

########差异分析
dds <- DESeq(dds)
#dds <- estimateSizeFactors(dds) 
#计算归一化系数sizeFactor
##dds <- estimateDispersions(dds) 
##估计基因的离散程度
##dds <- nbinomWaldTest(dds) 
##统计检验，差异分析

##获得分析结果
res1 <- results(dds,tidy=TRUE)
# ?results
################################################################################
metadata$dex=c("a","a","a","b","b","c","c","c")

dds <- DESeqDataSetFromMatrix(countData=counts, 
                              colData=metadata, 
                              design=~dex, 
                              tidy=TRUE)

dds <- DESeq(dds)
resultsNames(dds) 
res1 <- results(dds, name="dex_c_vs_a",tidy = T) 
rse1=as.data.frame(res1)
res2 <- results(dds, name="dex_b_vs_a") 
##结果的简单统计
summary(res2) 

res.05 <- results(dds, alpha=0.05) 
#alpha为padj的阈值，默认padj=0.1。
resLFC <- results(dds, lfcThreshold=1) 
#提升log2 fold change threshold，结果中不满足lfc阈值的gene的p值都是1。


##值得一提的是DESeq2软件独有的normlization方法！
dds1 <- rlogTransformation(dds)  ## 得到经过DESeq2软件normlization的表达矩阵！
exprSet_new=assay(dds1)

par(cex = 0.7)#指定符号的大小

n.sample=ncol(counts)#样本数

cols <- rainbow(n.sample*1.2)
?rainbow
colnames(counts)
c1=column_to_rownames(counts,"ensgene")
boxplot(c1, col = cols,main="expression value",las=2)

boxplot(exprSet_new, col = cols,main="expression value",las=2)




##res=arrange(res,padj)
anno <- read.csv("annotables_grch38.csv")
##res_anno = res %>%inner_join(anno,by=c("row"="ensgene"))

res_anno=inner_join(res1,anno,by=c("row"="ensgene"))
res_anno %>% filter(padj<0.05) #筛选padj小于0.05的值

res_anno %>% 
  filter(padj<0.05) %>% 
  write.csv("sigresults.csv")


SPARCL1 <- res_anno %>% filter(symbol=="SPARCL1")
View(SPARCL1)#SPARCL1对的基因ID是ENSG00000103196

plotCounts(dds, gene="ENSG00000152583", intgroup="dex")
###值得注意的是plotCount返回的不仅仅是图，也可以是数据
plotCounts(dds, gene="ENSG00000152583", intgroup="dex", returnData=TRUE)

library(ggplot2)
plotCounts(dds, gene="ENSG00000152583", intgroup="dex", returnData=TRUE) %>% 
  ggplot(aes(dex, count)) + 
  geom_boxplot(aes(fill=dex)) +
  scale_y_log10() + 
  ggtitle("SPARCL1")


##BiocManager::install("apeglm")
library("apeglm")
resLFC <- lfcShrink(dds, coef=2, type="apeglm") 
# 缩小的倍数变化有助于按效应大小对基因进行排序和可视化。log2 FC估计不能解释我们在低read counts下观察到的大的离散程度。
# 与估计离散程度的收缩一样，LFC收缩使用来自所有基因的信息来生成更准确的估计。
# 如果要根据LFC值提取差异基因，需要shrunken values。另外，进行功能分析例如GSEA时，需要提供shrunken values。

plotMA(resLFC, ylim=c(-5,5))  #提高了log2 fold change阈值


# =========================================================================================================


###edgeR  edgeR对于单个样本是比较好的。
##BiocManager::install("edgeR")
  
library(edgeR)
##，使用edgeR中的cpm函数将原始计数转换为CPM和log-CPM值。
colnames(counts)
count1=column_to_rownames(counts,"ensgene")
###排序
count2=count1[,c(1,3,5,7,2,4,6,8)]

group.list=c(rep("control",4),rep("treated",4))
group.list=factor(group.list)
group.list=relevel(group.list,ref = "control")

###第一步： 构建DGEList对象 
dge <- DGEList(counts=count2,group=group.list)

####第二步： 过滤 low counts数据。与DESeq2的预过滤不同，DESeq2的预过滤只是为了改善后续运算性能，
# 简单粗暴的方法   基因至少在某一些文库的count超过10 ~ 15 才被认为是表达。
# keep <- rowSums(dge$counts) > 50
# 利用CPM标准化 0.5(即阈值）等于 10/(最小的文库的 read count数 /1000000)
keep <- rowSums(cpm(dge) > 0.5 ) >=2
table(keep)


###第三步： 根据组成偏好(composition bias)标准化
#利用库的大小来进行标准化TMM
dge <- calcNormFactors(dge) 
?calcNormFactors
##TMM <- calcNormFactors(d, method="TMM")
##TMM <- calcNormFactors(d, method="RLE")
##TMM <- calcNormFactors(d, method="upperquartile")
plotMDS(dge, col = rep(c('red', 'blue'), each = 4), dim = c(1, 2))

#第四步： 实验设计矩阵(Design matrix)  建立分组变量
design <- model.matrix(~0+group.list)
# 多个分组则创建makeContrasts(contrasts=c('group1_group2','group1_group3','group2_group3',levels = design))
# 若coef=2则为group1和group3比较
design1 <- model.matrix(~group.list)
##~0+group :不包括比较矩阵
##~group:包括了比较矩阵

rownames(design)<-colnames(dge)
colnames(design)<-levels(group.list)
#第五步： 估计离散值（Dispersion）。前面已经提到负二项分布（negative binomial，NB)需要均值和离散值两个参数。
#edgeR对每个基因都估测一个经验贝叶斯稳健离散值（mpirical Bayes moderated dispersion），
#还有一个公共离散值（common dispersion，所有基因的经验贝叶斯稳健离散值的均值）以及一个趋势离散值
##install.packages("statmod")
dge <- estimateDisp(dge, design, robust = TRUE)
## estimateDisp()实际上是个组合函数，可以一步得到多个计算结果
##dge <- estimateGLMCommonDisp(dge,design)
##dge <- estimateGLMTrendedDisp(dge, design)
##dge <- estimateGLMTagwiseDisp(dge, design)
plotBCV(dge)

####第六步 差异分析
##负二项式广义对数线性模型


fit <- glmQLFit(dge, design, robust=TRUE)   #拟合模型
cntr.vs.T <- makeContrasts(control-treated, levels=design)
res <- glmQLFTest(fit, contrast=cntr.vs.T)   #统计检验

##用有默认截距的方法
#fit <- glmFit(dge, design1, robust = TRUE)    
#lrt <- glmLRT(fit)  
edger.deg=topTags(res, n = nrow(res$counts))

dge_de <- decideTestsDGE(res, adjust.method = 'fdr', p.value = 0.05) 
#查看默认方法获得的差异基因
summary(dge_de)

plotMD(res, status = dge_de, values = c(1, -1), col = c('blue', 'red'))   
#作图观测
abline(h = c(-1, 1), col = 'gray', lty = 2)


##############类似然负二项式广义对数线性模型



###这里用的是glmQLFTest是因为前面用了glmQLTFit进行拟合，所以需要用QL F-test进行检验。

#ig.edger <- res$table[p.adjust(res$table$PValue, method = "BH") < 0.01, ]
##如果前面用的是glmFit，那么对应的就是glmLRT. 作者称QL F-test更加严格。多重试验矫正用的也是BH方法。

fit2 <- glmFit(dge, design)
#利用似然比检验来检验差异表达基因
fit2 <- glmLRT(fit2, contrast=c(-1,1)) 
#输出靠前的差异表达基因
DEG=topTags(fit2, n=nrow(exp))
DEG=as.data.frame(DEG)


########################################################

##原始数据尺度转换
cpm <- cpm(count1)
lcpm <- cpm(count1, log=TRUE, prior.count=2)
##删除低表达基因
##edgeR包中的filterByExpr函数提供了自动过滤基因的方法，可保留尽可能多的有足够表达计数的基因。
?filterByExpr()
dge <- filterByExpr(dge, group=group.list)

topTags(res,n=10)
is.de <- decideTestsDGE(res)
summary(is.de)
plotMD(res, status=is.de, values=c(1,-1), col=c("red","blue"),
       legend="topright")


# ===========================================================================================
###Linear Models for Microarray Data
library(limma)


##进行差异分析时常用limma。虽然它是针对芯片数据开发的，但也有limma-voom可以分析转录组数据
dge <- DGEList(counts=count2)

group.list=c(rep("control",4),rep("treated",4))
group.list=factor(group.list)
group.list=relevel(group.list,ref = "control")

##分组矩阵 (design) ：就是将表达矩阵的列（各个样本）分成几组（例如最简单的case - control，或者一些时间序列的样本day0, day1, day2 ...）【通过model.matrix()得到】

design <- model.matrix(~0+group.list)
rownames(design)<-colnames(dge)
colnames(design)<-levels(group.list)


### 注意： calcNormFactors并不会标准化数据，只是计算标准化因子
dge <- calcNormFactors(dge)

v <- voom(dge,design, normalize="quantile")
##voom到底做了什么转换？

# 首先原始counts转换成log2的CPM（counts per million reads ），
# 这里的per million reads是根据之前calcNormFactors计算的norm.factors进行规定的；
# 然后根据每个基因的log2CPM制作了线性模型，并计算了残差 ；
# 然后利用了平均表达量（红线）拟合了sqrt(residual standard deviation)；
# 最后得到的平滑曲线可以用来得到每个基因和样本的权重

##lmFit()：线性拟合模型构建
fit <- lmFit(v, design)

constrasts = paste(rev(levels(group.list)),collapse = "-")
##比较矩阵（contrast）：意思就是如何指定函数去进行组间比较【通过makeContrasts()得到】
cont.matrix <- makeContrasts(contrasts=constrasts,levels = design) 
### 比较每个基因
fit2=contrasts.fit(fit,cont.matrix)
###eBayes()使用trend=TRUE对标准误差进行经验贝叶斯平滑，计算每个对比中每个基因的moderated t-statistic和log-odds。
fit2=eBayes(fit2)
##topTable()给出一个最有可能在给定对比下差异表达的基因列表。
###使用coef参数，这里设为1，也就是表示根据上面makeContrasts的第一个（group2-group1）来提取结果
DEG = topTable(fit2, coef=constrasts,sort.by = "P", n=Inf)
##adjust="BH"表示对p值的校正方法，包括了："none", "BH", "BY" and "holm"。
###那么为啥要对P值进行校正呢？
#p值是针对单次检验，假设采用的p值为小于0.05，我们通常认为这个基因在两组样本中的表达是有显著差异的，但是仍旧有5%的概率表示这个基因并不是差异基因。
#但是，当两组样本中有20000个基因采用同样的检验方式进行统计检验时，就会遇到一个问题：单次犯错的概率为0.05， 如果进行20000次检验，那么就会有0.05*20000=1000 个基因在组间的差异被错误估计
# topTable 列出差异显著基因
DEG2=topTable(fit2, coef=1, adjust="BH")
##去掉那些NA值
DEG = na.omit(DEG)


# edgeR 使用经验贝叶斯估计和基于负二项模型的精确检验来确定差异基因。特别地，经验贝叶斯用于通过在基因之间来调节跨基因的过度离散程度。 使用类似于Fisher精确检验但适应过度分散数据的精确检验用于评估每个基因的差异表达。edgeR 在默认情况下，执行TMM标准化程序以考虑样本之间的不同测序深度，通过Benjamini-Hochberg用于控制FDR 。
# Limma包基于线性模型建模。 它最初设计用于分析微阵列数据，但最近已扩展到RNA-seq数据。 根据limma用户指南的当前建议是使用edgeR包的TMM标准化和“voom”转换，其本质上将标准化数据取对数（基数2）并估计它们的均值 - 方差关系以确定在线性建模之前每次观察的权重。 默认情况下，Benjamini-Hochberg程序用于估计FDR 。
# DESeq使用类似于edgeR的负二项式模型，与edgeR类似，执行缩放因子归一化以考虑不同样本的变化的测序深度，并且Benjamini-Hochberg用于控制FDR。 DESeq能够分析具有少量重复的实验。DESeq技术上可以在没有任何生物学重复的情况下进行实验。DESeq2是在DESeq基础上更新的软件。 

# 这两个都属于R包，其相同点在于都是对count data数据进行处理，都是基于负二项分布模型。因此会发现，用两者处理同一组数据，最后在相同阈值下筛选出的大部分基因都是一样的，但是有一部分不同应该是由于其估计离散度的不同方法所导致的。

# =======================================================================
  
  
  