###让报错变成英文 方便google
Sys.setenv(LANGUAGE = "en")
#禁止chr转成factor
options(stringsAsFactors = FALSE)
###清空环境
rm(list=ls())
###设置路径
setwd("C:\\shangke\\lession6")
##加载包  这个包是用来加快数据读取的
library(data.table)
library(dplyr)
library(tidyverse)
###读取胃癌的临床信息
stad.phe=fread("TCGA-STAD.GDC_phenotype.tsv.gz",header = T, sep = '\t',data.table = F)
###查看一下数据类型
class(stad.phe)
###读取表达谱信息  rnaseq的数据
stad.fkpm=fread("TCGA-STAD.htseq_fpkm.tsv.gz",header = T, sep = '\t',data.table = F)
class(stad.fkpm)
###查看列名 方便merger的合并列设置
colnames(stad.fkpm)
####读取探针信息 ，目的是为了将ensemble名字转为基因名
stad.pro=fread("gencode.v22.annotation.gene.probeMap",header = T, sep = '\t',data.table = F)
##查看一下数据
View(stad.pro)
colnames(stad.pro)
###我们只要前两列进行转换
stad.pro=stad.pro[,c(1,2)]
colnames(stad.pro)
###用merge函数将探针转化的信息和表达谱信息进行合并
stad.fkpm.pro=merge(stad.pro,stad.fkpm,by.y ="Ensembl_ID",by.x = "id" )
?merge
dim(stad.fkpm.pro)
distinct()

##使用aggregate函数主要目的是对gene列中的相同基因进行合并 
stad.fkpm.pro.1 <- aggregate(stad.fkpm.pro[,-c(1,2)], list(stad.fkpm.pro$gene), FUN=sum)
##*aggregate函数的功能比较强大，它首先将数据进行分组（按行），
# 然后对每一组数据进行函数统计，最后把结果组合成一个比较nice的表格返回。
# 根据数据对象不同它有三种用法，分别应用于数据框（data.frame）、公式（formula）和时间序列（ts）。
# 这里我们的输入和输出都为data.frame。
# aggregate(x, by, FUN, ..., simplify = TRUE)——by参数也可以包含多个类型的因子，得到的就是每个不同因子组合的统计结果。
# 这里我们的by参数为因子，是通过list函数将gene这一列转化为因子，转化后便不会有重复基因存在。
####################################################################
stad.fkpm.pro=distinct(stad.fkpm.pro,gene,.keep_all = T)

stad.fkpm.pro=distinct(stad.fkpm.pro,gene,.keep_all = T)
dim(stad.fkpm.pro)
stad.fkpm.pro <- column_to_rownames(stad.fkpm.pro,"gene")
###查看合并后的数据  主要看第一列和最后一列
View(stad.fkpm.pro)
###因为tcga数据中有癌和癌旁，所以我们先根据临床信息把癌和癌旁区分一下
View(stad.phe)
##通过查看临床信息，我们发现在sample_type.samples列中，Primary Tumor为癌，Solid Tissue Normal为癌旁
rownames(stad.phe)=stad.phe$submitter_id.samples
colnames(stad.phe)
table(stad.phe$sample_type.samples)
stad.phe.t=filter(stad.phe,sample_type.samples=="Primary Tumor")
stad.phe.n=filter(stad.phe,sample_type.samples=="Solid Tissue Normal")
intersect(c("a","b"),c("b","c"))
z1=intersect(rownames(stad.phe.t),colnames(stad.fkpm.pro))
z2=intersect(rownames(stad.phe.n),colnames(stad.fkpm.pro))
stad.t=stad.fkpm.pro[,z1]
stad.n=stad.fkpm.pro[,z2]
colnames(stad.n)=paste0("N",1:32)
paste("xy",c(1,5,8,11,3),sep = "-")
?paste
colnames(stad.t)=paste0("T",1:375)
stad.exp=merge(stad.n,stad.t,by.x = 0,by.y = 0)
colnames(stad.exp)
stad.exp <- column_to_rownames(stad.exp,"Row.names")
###读取gtex的表达矩阵，注意解压和不解压都是可以读取的
gtex.exp=fread("gtex_RSEM_gene_fpkm.gz",header = T, sep = '\t',data.table = F)
###
gtex.exp[1:5,1:5]
###读取gtex的临床样本注释信息
gtex.phe=fread("GTEX_phenotype.gz",header = T, sep = '\t',data.table = F)
##查看一下
View(gtex.phe)
###读取gtex的基因注释信息 也就是探针信息
gtex.pro=fread("probeMap_gencode.v23.annotation.gene.probemap",header = T, sep = '\t',data.table = F)
####我们比较一下v23和v22的差异  一个是60498  一个是60483
###现在我们先合并gtex的基因信息
###合并之前先看一下列名  找到共同的合并列
colnames(gtex.pro)
colnames(gtex.exp)
gtex.pro=gtex.pro[,c(1,2)]
###我们发现sample和id是共同的列
###gtex.exp=merge(gtex.exp,gtex.pro,by.x ="sample",by.y = "id" )

######################################################################
save(stad.exp,file = "stad.exp.rdata")
rm(stad.exp,stad.fkpm,stad.fkpm.pro,stad.n,stad.t,stad.phe)


colnames(gtex.phe)
rownames(gtex.phe)=gtex.phe$Sample
table(gtex.phe$primary_site)
colnames(gtex.phe)=c("Sample","body_site_detail (SMTSD)","primary_site","gender","patient","cohort")
colnames(gtex.phe)
table(gtex.phe$primary_site)
gtex.phe.s=filter(gtex.phe,primary_site=="Stomach")
x1=intersect(rownames(gtex.phe.s),colnames(gtex.exp))
View(gtex.phe.s)
colnames(gtex.exp)[1:10]
gtex.s=gtex.exp[,c("sample",x1)]

gtex.exp=merge(gtex.pro,gtex.s,by.x ="id",by.y ="sample")


gtex.s1=distinct(gtex.exp,gene,.keep_all = T)
gtex.s2 <- column_to_rownames(gtex.s1,"gene")
gtex.s2 =gtex.s2[,-1]


###我们发现一共有209个胃组织
###我们从官网可以看到gtex是按照log2(fpkm+0.001)处理的  stad是按照log2(fpkm+1)，所以我们在合并之前 先要把他们的处理方式变成一样。
gtex.s3=2^gtex.s2

gtex.s5=log2(gtex.s3-0.001+1)
colnames(gtex.s5)=paste0("G",1:174)


###取一波交集  目的是为了决定之后的gtex和stad合并是按照symbol还是enseble来合并
length(intersect(gtex.exp$id,stad.fkpm.pro$id))
length(intersect(gtex.exp$gene,rownames(stad.fkpm.pro)))
###我们可以看到 如果按照基因合并会有57793个交集  如果按照Ensembl却只有42566个，所以最后还是按照gene来合并
###现在要提取正常的胃组织的表达矩阵，我们要根据gtex的临床信息来匹配胃组织的sample

###现在数据的处理方式都相同了 就有可比性了 将gtex的胃组织的数据与tcga的数据进行合并
all.data=merge(gtex.s5,stad.exp,by.x = 0,by.y = 0)
all.data <- column_to_rownames(all.data,"Row.names")
fwrite(all.data,file = "stomach.cancer.gtex.tcga.txt",sep = "\t",row.names = T)
write.table(all.data,file = "stomach.cancer.gtex.tcga.txt",sep = "\t")
all.data1<- as.matrix(all.data)

#################################################################################################
library(limma)
nromalized.data=normalizeBetweenArrays(all.data)

?normalizeBetweenArrays
# gs1  5个基因  A B C D E

# 1S: A 2   B 3   C 4    D 5    E 6
# 2S：A 3  B 4    C 5    D 6    E 2
library(GSVA)
library(GSEABase)
##BiocManager::install("clusterProfiler")
library(clusterProfiler)
#install.packages("devtools")
#devtools::install_github("GSEA-MSigDB/GSEA_R")
library(GSEA)

kegggmt2 <- read.gmt("c2.cp.kegg.v7.4.symbols.gmt")

kegg_list = split(kegggmt2$gene, kegggmt2$term)

##BiocManager::install("GSVA")
library(GSVA)

?gsva
expr=as.matrix(expr)
kegg2 <- gsva(nromalized.data, kegg_list, kcdf="Gaussian",method = "gsva")


?gsva
es.dif <- gsva(all.data1, gmt, mx.diff=TRUE)
es.max <- gsva(all.data1, gmt, mx.diff=FALSE)
es.dif.nromalized <- gsva(nromalized.data, gmt, mx.diff=TRUE,kcdf="Gaussian",parallel.sz=8)
es.max.nromalized <- gsva(nromalized.data, gmt, mx.diff=FALSE)


? data.frame
annotation_col = data.frame(
  Tissuetype =c(rep("Stomach",174),rep("Solid Tissue Normal",32),rep("Tumor",375)),
  Database =c(rep("GTEX",174), rep("TCGA",407))
)

rownames(annotation_col)=colnames(nromalized.data)

pheatmap::pheatmap(kegg2[1:20,], #热图的数据
                   cluster_rows =T,#行聚类
                   cluster_cols =F,#列聚类，可以看出样本之间的区分度
                   annotation_col = annotation_col,
                   show_colnames=F,
                   scale = "row", #以行来标准化，这个功能很不错
                   color =colorRampPalette(c("#FF7744", "white","#AAAAAA","#0044BB"))(100))


###############################################################################################
##不同物种

##install.packages("msigdbr")
library(msigdbr)
msigdbr_species()

h.kegg= msigdbr(species ="Homo sapiens", category = "C2", subcategory = "KEGG")

h.kegg=h.kegg[,c(3,4)]
h.kegg_list = split(h.kegg$gene_symbol, h.kegg$gs_name)

###############################################################

##自定义基因集、另外一种gmt


gene.set=read.table("GENE.txt",
                    header =F,sep = '\t',quote = '')
kegg.name=read.table("GENE.SET.NAME.txt",
                     header =F,sep = '\t',quote = '')
gene.set1=as.matrix(gene.set)
gene.set2=t(gene.set1)
gmt=list()
for (i in 1:18) {
  y=as.character(gene.set2[,i]) 
  b=which(y=="")
  gmt[[i]]=y[-b]
}
names(gmt)=kegg.name$V1
gmt=gmt[-18]


###################################################

y=as.character(gene.set2[,1]) 
b=which(y=="")
gmt[[1]]=y[-b]

View(gmt)
getwd()





###############特定基因分析
exprSet.all.r=all.data[c("RORC", "RORB", "RORA"),]

exprSet.all.r=t(exprSet.all.r)
exprSet.all.r=as.data.frame(exprSet.all.r)

x=c(rep("GTEX",174),rep("N",32),rep("T",375))
exprSet.all.r$Type=x

exprSet.rorc=exprSet.all.r[,c(1,4)]
exprSet.rorc$Gene=rep("RORC")
colnames(exprSet.rorc)[1]="Relative Expression"
exprSet.rorb=exprSet.all.r[,c(2,4)]
exprSet.rorb$Gene=rep("RORB") 
colnames(exprSet.rorb)[1]="Relative Expression"
exprSet.rora=exprSet.all.r[,c(3,4)] 
exprSet.rora$Gene=rep("RORA") 
colnames(exprSet.rora)[1]="Relative Expression"
x.all=rbind(exprSet.rorc,exprSet.rorb,exprSet.rora)
colnames(x.all)
library(ggsignif)
library(ggpubr)
library(ggplot2)
p <- ggboxplot(x.all, x = "Gene", y = "Relative Expression",
               color = "Type", palette = "Type",
               add = "Type")
p + stat_compare_means(aes(group = Type))

table(x.all$Gene)
my_comparisons <- list(c("RORA","RORB"), c("RORA","RORC"),c("RORB", "RORC"))


p +geom_signif(comparisons = my_comparisons,
               step_increase = 0.2,map_signif_level = F,
               test = t.test,size=0.8,textsize =4)
?geom_signif
x.c.b=cbind(exprSet.rorc,exprSet.rorb)
GSEA
GSVA
GO/KEGG
colnames(x.c.b)=c("RORC","Type","Gene", "RORB","Type","Gene" )
x.c.b=x.c.b[,c(1,4)]
library(ggplot2)
library(ggpubr)
## Loading required package: magrittr
p1 <- ggplot(data = x.c.b, mapping = aes(x = RORC, y = RORB)) +
  geom_point(colour = "red", size = 2) +
  geom_smooth(method = lm, colour='blue', fill='gray') #添加拟合曲线
p1
p1 + stat_cor(method = "pearson", label.x =3, label.y = 1) #添加pearson相关系数


##############################################################################################
#####生存曲线
stad.fkpm.pro.1=column_to_rownames(stad.fkpm.pro.1,colnames(stad.fkpm.pro.1)[1])

library(limma)
stad.fkpm.pro.1=normalizeBetweenArrays(stad.fkpm.pro.1)
kegg_list=kegg_list[1:10]
kegg2 <- gsva(stad.fkpm.pro.1, kegg_list,method = "gsva")


sur=fread("TCGA-STAD.survival.tsv",data.table = F)

sur1 <- sur[sur$sample %in% colnames(stad.fkpm.pro.1),]
##将TCGA表达矩阵的病人样本与临床数据样本进行整合

sur.kegg=as.data.frame(kegg2)

sur.kegg.1=sur.kegg[2,]
sur.kegg.2=as.data.frame(t(sur.kegg.1))
colnames(sur.kegg.2)="Glycan"

sur.glycan <- merge(sur1, sur.kegg.2, by.x="sample", by.y=0)


##中位数分组
sur.glycan$m = ifelse(sur.glycan$Glycan>median(sur.glycan$Glycan),'high','low')


################## 加载生存分析的r包

library(survival)
#install.packages("survminer")
library(survminer)


################################################################################
colnames(sur.glycan)
ggsurvplot(survfit(Surv(OS.time,OS)~m, data=sur.glycan), conf.int=T, pval=TRUE)

ggsurvplot(survfit(Surv(OS.time,OS)~m, data=sur.glycan), conf.int=F, pval=T)

#Surv()函数创建生存数据对象（主要输入生存时间和状态逻辑值）

#survfit()函数对生存数据对象拟合生存函数，创建KM(Kaplan-Meier)生存曲线

#survdiff()用于不同组的统计检验
fit=survfit(Surv(OS.time,OS)~m, data=sur.glycan)
ggsurvplot(fit,   
           title = "Survival curve", # 添加标题
           font.main = c(16, "bold", "darkblue"), # 设置标题字体大小、格式和颜色
           font.x = c(14, "bold.italic", "red"), # 设置x轴字体大小、格式和颜色
           font.y = c(14, "bold.italic", "darkred"), # 设置y轴字体大小、格式和颜色
           font.tickslab = c(12, "plain", "darkgreen")) # 设置坐标轴刻度字体大小、格式和颜色

###title = "Survival curve", # 添加标题
############################################################################


###四分位生存曲线
quantile(sur.glycan$Glycan)
s2.exp$zz=ifelse(s2.exp$EZH2>9.53,'high', 
                 ifelse( s2.exp$EZH2>9.23& s2.exp$EZH2<9.53,'h.stable', 
                         ifelse( s2.exp$EZH2>8.88& s2.exp$EZH2<9.23,'l.stable','down') )
)
