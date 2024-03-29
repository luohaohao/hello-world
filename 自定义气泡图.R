Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())
# setwd("C:\\shangke\\lession7")
kegg=read.table("zi_ding_yi.txt",header = T,sep = "\t")


k = data.frame(kegg)
library(ggplot2)
library(dplyr)
before <- as.numeric(sub("/\\d+$",replacement = "", k$GeneRatio))
after <- as.numeric(sub("^\\d+/", "", k$GeneRatio))
x=k$GeneRatio
x1=as.character(x)
a=strsplit(x1,split="/",fixed=T)
be=sapply(a,function(x){x[1]})
be1=as.numeric(be)
after=sapply(a,function(x){x[2]})
after=as.numeric(after)
?strsplit
k$GeneRatio = be1/after
font.size =10

k %>% 
  ## 对进行p值排序
  arrange(p.adjust) %>% 
  ##指定富集的通路数目
  slice(1:10) %>% 
  ## 开始ggplot2 作图，其中fct_reorder调整因子level的顺序
  ggplot(aes(GeneRatio,forcats::fct_reorder(Description,Count)))+ 
  ## 画出点图
  geom_point(aes(color=p.adjust, size = Count)) +
  ## 调整颜色，guide_colorbar调整色图的方向
  scale_color_continuous(low="red", high="blue", guide=guide_colorbar(reverse=TRUE))+
  ## 调整泡泡的大小
  scale_size_continuous(range=c(3, 8))+
  ## 如果用ylab("")或出现左侧空白
  labs(y=NULL) +
  ## 如果没有这一句，上方会到顶
  ggtitle("")+
  ## 设定主题
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black",
                                   size = font.size, vjust =1 ),
        axis.text.y = element_text(colour = "black",
                                   size = font.size, hjust =1 ),
        axis.title = element_text(margin=margin(10, 5, 0, 0),
                                  color = "black",size = font.size),
        axis.title.y = element_text(angle=90))

