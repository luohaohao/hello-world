##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())

### 加载r包
library(GEOquery)
library(dplyr)
library(tidyverse)
##设置路径
getwd()
# setwd("C:\\LH\\lession4/deg/")

#########################################################################################
##下载数据 getGPL = F 这个代表不下载平台注释文件，因为有时候网络不稳定。后面我们会在网页中下载，然后读取。
########## 方法一 参数里面改为getGPL =T
gset = getGEO('GSE149507', destdir=".",getGPL =F)
exp = gset[["GSE149507_series_matrix.txt.gz"]]@assayData[["exprs"]]

###############方法二 将网页版的soft文件下载
soft=getGEO(filename ="GSE149507_family.soft")
gpl=soft@gpls[["GPL23270"]]@dataTable@table

##############方法三 直接去网页里面复制

gset=gset[[1]]
pdata=pData(gset)
##得到表达矩阵
exprSet=exprs(gset)

#############################################
###将rntrez id转为gene symbol
##BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
###查看一下有哪些基因ID
keytypes(org.Hs.eg.db)
x1=gpl$ENTREZ_GENE_ID
x1=as.character(x1)
x2=AnnotationDbi::select(org.Hs.eg.db, keys=x1, columns=c("ENTREZID", "SYMBOL"), keytype="ENTREZID")
gpl=gpl[,c(1,2)]
gpl$gene=x2$SYMBOL

###平台和表达矩阵合并
colnames(gpl)
exp=as.data.frame(exprSet)
exp.pl=merge(gpl,exp,by.x=1,by.y=0)
#############################################################################
####将gene名变为行名
# rownames(exp.pl)=exp.pl$gene
# colnames(exp.pl)
exp.pl.1=distinct(exp.pl,gene,.keep_all = T)
###我们把基因名字变为行名
# rownames(exp.pl.1)=exp.pl.1$gene
###发现报错了  这是因为有缺失值的存在
###删除缺失值
exp.pl.2=exp.pl.1[,-4]
exp.pl.2=na.omit(exp.pl.1)
##再次转变为行名
rownames(exp.pl.2)=exp.pl.2$gene
##去除没用信息，只剩下表达矩阵
exp.pl.3=exp.pl.2[,-c(1:3)]

######################################################################################
##查看分组信息
View(pdata)

a=c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35)
a2=a+1
exp.4=exp.pl.3[,c(a2,a)]

##这里我们创建分组信息，意思就是前两个是control，后面两个是处理
group_list=c(rep('control',18),rep('T',18))
##转化为因子  这里可以不用理解 为了后面差异分析的矩阵构建
group_list=factor(group_list)
## 强制限定顺序
group_list <- relevel(group_list, ref="control")

library(limma) 
##除去批次效应
exp.5=normalizeBetweenArrays(exp.4)

boxplot(exp.4,outline=FALSE, notch=T,col=group_list, las=2)
boxplot(exp.5,outline=FALSE, notch=T,col=group_list, las=2)
library(RColorBrewer)
colors<-brewer.pal(18,"Set3")
boxplot(exp.5,col=colors,notch=T,outline=FALSE, las=3,ylim=c(2,10))
boxplot(exp.4,col=colors,notch=T,outline=FALSE, las=3,ylim=c(2,10))


##构建差异分析的矩阵
design=model.matrix(~ group_list)
# 第一列为control列，第二列为1的与control进行比较
colnames(design) <- levels(group_list)
rownames(design) <- colnames(exp.5)


##lmFit()：线性拟合模型构建
fit=lmFit(exp.5,design)
##eBayes()使用trend=TRUE对标准误差进行经验贝叶斯平滑，计算每个对比中每个基因的moderated t-statistic和log-odds。
fit=eBayes(fit) 
##topTable()给出一个最有可能在给定对比下差异表达的基因列表。
allDiff=topTable(fit,coef=2,adjust='fdr',number=Inf) 
##topTable函数的coef参数，coef=2是指design的第2列，即tumour，即把tumour与normal进行对比
write.table(allDiff,file = "allDiff.txt",sep = "\t",col.names = NA)
write.csv(allDiff,file = 'allDiff.csv')

###############################################################
###有配对信息的差异分析
pairinfo = factor(rep(1:18,2))
design=model.matrix(~ pairinfo+group_list)
fit=lmFit(exprSet,design)
fit=eBayes(fit) 
allDiff_pair=topTable(fit,adjust='BH',coef="T",number=Inf,p.value=0.05)
##adjust="BH"表示对p值的校正方法，包括了："none", "BH", "BY" and "holm"。


###########################################################################################
##画个火山图

library(ggplot2)                         

library(ggrepel)                               

data<-allDiff
data$significant="stable"
log2(1.5)
data$significant[data$logFC>=0.585 & data$P.Value <0.05]="up"
data$significant[data$logFC<= -0.585 & data$P.Value <0.05]="down"

ggplot(data,aes(logFC,-1*log10(P.Value)))+xlim(-2,2)+ylim(0,6)+
  geom_point(aes(color=significant),size=0.8)+theme_classic()+
  scale_color_manual(values = c("#2a9d8f","#EE7AE9","#f8961e"))+
  geom_hline(yintercept = 1.3,linetype=4,size=0.3)+
  geom_vline(xintercept = c(-0.585,0.585),linetype=4,size=0.3)+
  theme(title=element_text(size = 18),text = element_text(size=18))+
  labs(x="log2(foldchange)",y="-log10(p_value)")

select.FPKM <- (data$AveExpr > 5 )
table(select.FPKM)

select.log2FC <- abs(data$logFC) >1
table(select.log2FC)
select.qval <- (data$adj.P.Val< 0.05)
table(select.qval)

select.vec=(select.FPKM & select.log2FC & select.qval) 
table(select.vec)

degs.list=as.character(rownames(data))[select.vec]

label.deg=sample(degs.list,20)

p=ggplot(data,aes(logFC,-1*log10(P.Value)))+xlim(-2,2)+ylim(0,6)+
  geom_point(aes(color=significant),size=0.8)+theme_classic()+
  scale_color_manual(values = c("#2a9d8f","#EE7AE9","#f8961e"))+
  geom_hline(yintercept = 1.3,linetype=4,size=0.3)+
  geom_vline(xintercept = c(-0.5,0.5),linetype=4,size=0.3)+
  theme(title=element_text(size = 18),text = element_text(size=18))+
  labs(x="log2(foldchange)",y="-log10(p_value)")
data_selected <- data[label.deg,]
p + geom_label_repel(data=data_selected,
                     aes(label=rownames(data_selected)))
###########################
data_selected <- data["PCNX2",]
p + geom_label_repel(data=data_selected,
                     aes(label=rownames(data_selected)))


####################################################3
##画个热图
# 列名注释
annotation_col1 = data.frame(
  # Database =c(rep("GEO",36)),
  CellType =c(rep("control",18),rep("treatment",18)) 
)
rownames(annotation_col1)=colnames(exp.5)
exp.6=as.data.frame(exp.5)
exprSet.map=exp.6[label.deg,]
#install.packages("pheatmap")
pheatmap::pheatmap(exprSet.map, #热图的数据
                   cluster_rows = F,#行聚类
                   cluster_cols =F,#列聚类，可以看出样本之间的区分度
                   annotation_col =annotation_col1,
                   show_colnames=F,
                   scale = "row", #以行来标准化，这个功能很不错
                   color =colorRampPalette(c("blue", "white","red"))(100))
####聚类
# pheatmap::pheatmap(exprSet.map, #热图的数据
#                    cluster_rows =T,#行聚类
#                    cluster_cols =F,#列聚类，可以看出样本之间的区分度
#                    annotation_col =annotation_col1,
#                    show_colnames=F,
#                    scale = "row", #以行来标准化，这个功能很不错
#                    color =colorRampPalette(c("blue", "white","red"))(100))

######################################################
###箱线图
label.deg[1:10]
exp.x=exp.6["TTC32",]

z=t(exp.x)
z=as.data.frame(z)
z$type=c(rep("control",18),rep("treatment",18))

colnames(z)
library(ggpubr)
library(ggsignif)
library(ggplot2)


ggboxplot(z,x="type",y="TTC32",
          width = 0.6,fill="type",
          notch = T,palette = c("#00AFBB", "red","#E7B800"),
          add = "jitter",shape="type")

###加个p值
p=ggboxplot(z,x="type",y="TTC32",
            width = 0.6,fill="type",
            notch = T,palette = c("#00AFBB", "red","#E7B800"),
            add = "jitter",shape="type")

p + stat_compare_means(aes(group =type))

####################################################
###箱线图
ggplot(z,aes(type,TTC32,fill=type)) +
  geom_boxplot() + stat_compare_means(aes(group=type))+
  geom_jitter(aes(color = type),position = position_jitter(0.2)) +
  geom_signif(comparisons = list(c('control','treatment')) ,step_increase = 0.1,map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05),test =t.test,size=1,textsize = 6)


#######匹配一下
z$pairinfo=pairinfo=rep(1:18,2)
ggplot(z, aes(type,TTC32,fill=type)) +
  geom_boxplot() +
  geom_point(size=2, alpha=0.5) +
  geom_line(aes(group=pairinfo), colour="black", linetype="11") +
  xlab("") +
  ylab(paste("Expression of ","TTC32"))+
  theme_classic()+
  theme(legend.position = "none")

##########################################################################################
#####小提琴

ggplot(z,aes(type,TTC32,fill=type)) +
  geom_violin()+ geom_jitter(shape=16, position=position_jitter(0.2))

library(ggsignif)
compaired <- list(c("control", "treatment"))


library(ggsci)
colnames(z)
ggplot(z,aes(type,TTC32,fill=type)) +
  geom_violin()+
  geom_signif(comparisons =compaired ,step_increase = 0.1,map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05),test =t.test,size=2,textsize = 6)+
  geom_jitter(shape=16, position=position_jitter(0.2))

help(package="ggsci")
