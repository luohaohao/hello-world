##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())
## 加载r包
library(GEOquery)
library(dplyr)
library(tidyverse)
##设置路径
getwd()
# setwd("C:\\shangke\\lession6")

####################################################################################################
### 获取表达矩阵
?getGEO
##下载数据 getGPL = F 这个代表不下载平台注释文件，因为有时候网络不稳定。后面我们会在网页中下载，然后读取。
gset = getGEO('GSE12417', destdir=".",getGPL = F)

e2=gset[[2]]


##当然也可以用环境中的白色括号提取
exp2=e2@assayData[["exprs"]]



#########################################################################################################
### 整理探针转化文件

library(data.table)
###读取测序平台信息，为了进行探针的转化  这里的探针就是exprSe的行名比如：10344614
###平台文件在geo中下载，下载之后要放在指定的路径下面
###用fread函数读取的目的是多线程快速读取
anno=fread("GPL96-57554.txt",sep = "\t",header = T,data.table = F)

colnames(anno)
gene=anno[,c(1,11)]

x=gene$`Gene Symbol`

## "DDR1 /// MIR4640"  这个东西叫字符串
### 字符串切割  字符串的处理公众号专门有一篇推文在讲
a1=strsplit(x,split = " /// ",fixed = T)
##把分割后的字符串的第一个元素提取出来，合并成为一个新的向量
gene.all = sapply(a1,function(x){x[1]})
##其实就是等同于下面

###公众号有一篇图文专门讲apply家族的用法 包括sapply lapply tapply，感兴趣的可以自己去看公众号的推文

### 这时候 我们可以发现anno有22283行，gene.all也是有22283个元素
a3=data.frame(anno$ID,gene.all)
### 这时候的a3创建是否一一对应？  原理：向量是有序的
###向量的特点看r语言基础 第二节课的第6页  向量的特性

###到此为止  我们讲注释的gene symbol文件与表达矩阵都已经准备好了
###########################################################################################
###进行探针的转化
exp=as.data.frame(exp2)###讲矩阵转化为数据框的目的是为了好对表达矩阵进行操作
View(exp) ###行名为探针名 列名为样本名


###用merge函数

exp1=merge(x=a3,y=exp,by.x =1 ,by.y =0 )


###整理表达矩阵
rownames(exp1)=exp1$gene.all
###报错的原因就是有重复的基因名  多个探针对应同一个基因
exp2=distinct(exp1,gene.all,.keep_all = T)
###这时候行名已经从22283变成了13238
rownames(exp2)=exp2$gene.all
###发现还是报错 报错的原因是有NA  缺失值
exp3=na.omit(exp2)
### 从13239变成13237 有一个NA  为啥是一个呢？
rownames(exp3)=exp3$gene.all

###去除第一列和第二列
exp4=exp3[,-c(1,2)]

####读取RNAseq数据
###数据来源https://4va.github.io/biodatasci/r-rnaseq-airway.html#data_needed  
counts <- read.csv("airway_scaledcounts.csv")
metadata <-  read.csv("airway_metadata.csv")

anno <- read.csv("annotables_grch38.csv")


anno=anno[,c(1,3)]
res_anno=inner_join(anno,counts,by=c("ensgene"="ensgene"))
rownames(res_anno)=res_anno$symbol
res_anno=distinct(res_anno,symbol,.keep_all = T)
rownames(res_anno)=res_anno$symbol
res_anno=res_anno[,-c(1,2)]
#########################################################################
blue <- "#33FF33"
yellow <- "#FF6600"


# 提取在三套数据中都出现的基因
comgene <- intersect(rownames(res_anno),rownames(exp4))
# 合并三套数据
combined.expr <- cbind.data.frame(res_anno[comgene,],
                                  exp4[comgene,])

# 绘制PCA散点图，检查批次效应

source("batchPCA.R")
View(batchPCA)
##install.packages("ClassDiscovery")

batchPCA(indata = t(scale(t(combined.expr))),
         batch = rep(c("rnaseq","arry"), times = c(ncol(res_anno),ncol(exp4))),
         fig.dir = ".",
         PCA.fig.title = "before batch",
         cols = c(blue, yellow),
         showID = F,
         cex = 0.7,
         showLegend = T) 

##BiocManager::install("sva")
library(sva)
batch <- data.frame(batch = rep(c("rnaseq","arry"), times = c(ncol(res_anno),ncol(exp4))))
modcombat = model.matrix(~1, data = batch)
combined.expr.combat <- as.data.frame(ComBat(dat=as.matrix(combined.expr), batch=batch$batch, mod=modcombat))


batchPCA(indata = t(scale(t(combined.expr.combat))),
         batch = rep(c("rnaseq","arry"), times = c(ncol(res_anno),ncol(exp4))),
         fig.dir = ".",
         PCA.fig.title = "after Combat ",
         cols = c(blue, yellow),
         showID = F,
         cex = 0.7,
         showLegend = T) # 可以看到三个数据集(批次)混杂在了一次，说明批次效应被基本消除