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
###一定要设置自己的路径
# setwd("C:\\LH\\6.5")
Sys.setenv("VROOM_CONNECTION_SIZE"=99999999)

####################################################################################################
### 获取表达矩阵
?getGEO
##下载数据 getGPL = F 这个代表不下载平台注释文件，因为有时候网络不稳定。后面我们会在网页中下载，然后读取。
gset = getGEO(GEO='GSE12417', destdir=".",getGPL = T)
gpl = getGEO('GPL96', destdir='.')

anno = gpl@dataTable@table
######################保存数据 保存为rdata的数据结构
#save(gset,file = "123gse12417.rdata")
###打开工作路径的文件夹 我们会发现 123gse12417.rdata出现在文件夹下
##保存下载的数据 保存为rdata 当然也可以保存为rds
###########  这时候可以清空环境 保存的数据可以用load函数加载

#load("123gse12417.rdata")
############加载完数据后 gset对象又出现在环境中

###########也可以保存环境中的全局变量 自己看视频

##########将gset中的第二部分数据提取出来
e2=gset[[2]]


##当然也可以用环境中的白色括号提取
exp2=e2@assayData[["exprs"]]


## 获取分组信息 同理也有三种方法
phe=pData(e2)
pdata1=e2@phenoData@data
####  s4    s3:matrix 矩阵 data.frame() 数据框  character() 向量字符型 list 列表

#########################################################################################################
### 整理探针转化文件

library(data.table)
###读取测序平台信息，为了进行探针的转化  这里的探针就是exprSe的行名比如：10344614
###平台文件在geo中下载，下载之后要放在指定的路径下面
###用fread函数读取的目的是多线程快速读取
# anno=fread("GPL96.txt",sep = "\t",header = T,data.table = F)

View(anno)
colnames(anno)
gene=anno[,c(1,11)]
gene.1=anno[,c("ID","Gene Symbol")]
###exp是矩阵 merge合并必须是数据框的合并
####################################################################
#####数据框合并演示

a1=data.frame(
  x1=c(1, 2, 3),
  x2=c("b", "c", "d"),
  stringsAsFactors=FALSE
)


a2=data.frame(
  x1=c(10, 7, 12),
  x2=c("k1", "k2", "k3"),
  stringsAsFactors=FALSE
)
a3=cbind(a1,a2)
a4=rbind(a1,a2)

##############################  merge用法演示
d1 <- data.frame(
  x= c("a", "b", "c", "d"),
  y= c(1, 3, 5, 7)
)

d2 <- data.frame(
  x1=c(1, 7, 5),
  x2=c("m1", "m2", "m3"),
  stringsAsFactors=FALSE
)

d3=merge(x=d1,y=d2,by.x = "y",by.y = "x1")
?merge()
################################################ merge用法演示



exp1=as.data.frame(exp2)
exp.anno=merge(x=gene,y=exp1,by.x=1,by.y=0)
x=gene$`Gene Symbol`


## "DDR1 /// MIR4640"  这个东西叫字符串
### 字符串切割  字符串的处理公众号专门有一篇推文在讲
a1=strsplit(x,split = " /// ",fixed = T)
##把分割后的字符串的第一个元素提取出来，合并成为一个新的向量
gene.all = sapply(a1,function(x){x[1]})
##其实就是等同于下面

###公众号有一篇图文专门讲apply家族的用法 包括sapply lapply tapply，感兴趣的可以自己去看公众号的推文

### 这时候 我们可以发现anno有22283行，gene.all也是有22283个元素

### 这时候的a3创建是否一一对应？  原理：向量是有序的
###向量的特点看r语言基础 第二节课的第6页  向量的特性
## "DDR1 /// MIR4640"  这个东西叫字符串

###查看前十个元素
gene.all[1:10]


exp.anno$`Gene Symbol`=gene.all
###到此为止  我们讲注释的gene symbol文件与表达矩阵都已经准备好了
###########################################################################################
###进行探针的转化

exp1=exp.anno
colnames(exp1)[2]="gene.all"
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
View(exp3)
###去除第一列和第二列
exp4=exp3[,-c(1,2)]

######################
####生存曲线  先要提取生存的数据从临床信息中
s=phe$characteristics_ch1
s[1:5]

s1=strsplit(s,split = "; ",fixed = T)
s1[[1]][3]
##把分割后的字符串的第一个元素提取出来，合并成为一个新的向量
os.time = sapply(s1,function(x){x[3]})
status= sapply(s1,function(x){x[4]})

os.time[1:5]
os.time1=strsplit(os.time,split = " ",fixed = T)
os.time2= sapply(os.time1,function(x){x[3]})
os.time2=as.numeric(os.time2)
status[1:5]

status1=strsplit(status,split = ": ",fixed = T)
status2= sapply(status1,function(x){x[2]})
status2=as.numeric(status2)

s2=data.frame(os.time2,status2)

rownames(s2)=rownames(phe)
######################################################
###提取某个基因作为作为分组

exp4.e=exp4["EZH2",]

exp4.e <- data.frame(t(exp4.e))

###s2.exp <- merge(s2,exp4.e, by.x=0, by.y=0)
s2.exp <- merge(s2,exp4.e, by=0)


##############################################################################
###四分位生存曲线
quantile(s2.exp$EZH2)
s2.exp$zz=ifelse(s2.exp$EZH2>9.53,'high', 
                 ifelse( s2.exp$EZH2>9.23& s2.exp$EZH2<9.53,'h.stable', 
                         ifelse( s2.exp$EZH2>8.88& s2.exp$EZH2<9.23,'l.stable','down') )
)
############################################################################
##中位数分组
s2.exp$EZH2 = ifelse(s2.exp$EZH2>median(s2.exp$EZH2),'high','low')

################## 加载生存分析的r包

library(survival)
#install.packages("survminer")
library(survminer)

############################################
Surv(os.time2,status2)

x1=survfit(Surv(os.time2, status2)~EZH2, data=s2.exp)
summary(x1)

plot(x1)

plot(x1,xlab="Time(Days)",ylab="Survival probability",
     col=c("blue","red"),lty=2:3,lwd=2) 
# 添加图例
legend("topright",c("high","low"),
       col=c("blue","red"),lty=2:3,lwd=2,cex=1)



################################################################################
fit=survfit(Surv(os.time2, status2)~zz, data=s2.exp)

ggsurvplot(fit)
?ggsurvplot

ggsurvplot(survfit(Surv(os.time2, status2)~EZH2, data=s2.exp), conf.int=T, pval=TRUE)

ggsurvplot(survfit(Surv(os.time2, status2)~EZH2, data=s2.exp), conf.int=F, pval=T)

#Surv()函数创建生存数据对象（主要输入生存时间和状态逻辑值）

#survfit()函数对生存数据对象拟合生存函数，创建KM(Kaplan-Meier)生存曲线

#survdiff()用于不同组的统计检验

ggsurvplot(fit,   
           main = "Survival curve", # 添加标题
           font.main = c(16, "bold", "darkblue"), # 设置标题字体大小、格式和颜色
           font.x = c(14, "bold.italic", "red"), # 设置x轴字体大小、格式和颜色
           font.y = c(14, "bold.italic", "darkred"), # 设置y轴字体大小、格式和颜色
           font.tickslab = c(12, "plain", "darkgreen")) # 设置坐标轴刻度字体大小、格式和颜色
######################################################

ggsurvplot(fit, 
           surv.median.line = "hv", # 添加中位数生存时间线
           
           # Change legends: title & labels
           legend.title = "EZH2", # 设置图例标题
           legend.labs = c("high", "low"), # 指定图例分组标签
           
           # Add p-value and tervals
           pval = TRUE, # 设置添加P值
           pval.method = TRUE, #设置添加P值计算方法
           conf.int = TRUE, # 设置添加置信区间
           
           # Add risk table
           risk.table = TRUE, # 设置添加风险因子表
           tables.height = 0.2, # 设置风险表的高度
           tables.theme = theme_cleantable(), # 设置风险表的主题
           
           # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
           # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
           # palette = c("#E7B800", "#cc6666"), # 设置颜色画板
           ggtheme = theme_bw() # Change ggplot2 theme
)


################################
####十六进制颜色

fit <- survfit(Surv(os.time2, status2)~EZH2, data=s2.exp)
ggsurvplot(fit,
             pval = TRUE, conf.int = T,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#FF0000", "#5599FF")
)

################################
####十六进制颜色

fit <- survfit(Surv(os.time2, status2)~EZH2, data=s2.exp)
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#FF0000", "#5599FF")
)

####################################################

# 使用ggsurvplot_facet()函数绘制分面生存曲线
s2.exp$sex=c(rep("m",100),rep("f",63))
fit=survfit(Surv(os.time2, status2)~EZH2, data=s2.exp)
ggsurvplot_facet(fit, s2.exp, 
                 facet.by = "sex", # 设置分面变量
                 palette = "jco", # 设置颜色画板
                 pval = TRUE) # 添加pvalue值
#########################换个配色
ggsurvplot_facet(fit, s2.exp,  
                 facet.by = c("sex"),
                 palette = "npg", 
                 pval = TRUE,
                 surv.median.line = "hv",  # 增加中位生存时间
                 conf.int = TRUE) # 增加置信区间)

###################################
# 拟合多个分组变量
fit2 <- survfit( Surv(os.time2, status2) ~ sex + EZH2, data = s2.exp )
fit2

ggsurvplot(fit2)

#####绘制累计风险曲线
ggsurvplot(fit, data = s2.exp,
           conf.int = TRUE, # 增加置信区间
           fun = "cumhaz") # 绘制累计风险曲线

#############添加总患者生存曲线

ggsurvplot(fit, # 创建的拟合对象
           data = s2.exp,  # 指定变量数据来源
           conf.int = TRUE, # 显示置信区间
           pval = TRUE, # 添加P值
           surv.median.line = "hv",  # 添加中位生存时间线
           add.all = TRUE) # 添加总患者生存曲线



ggsurvplot(fit2,
           conf.int = TRUE,
           risk.table.col = "strata", # Change risk table color by groups
           ggtheme = theme_bw(), # Change ggplot2 theme
           surv.median.line = 'hv',
           fun = "event")
res.sum <- surv_summary(fit)
##想看fit每个时间点的结果用 surv_summary（）函数
ggsurvplot(res.sum,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF")
)
?surv_summary()
