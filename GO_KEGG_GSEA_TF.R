Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())
library(dplyr)
##BiocManager::install("Rgraphviz")
library(msigdbr)
library(GSEABase)
library(ggupset)
library(topGO)
library(Rgraphviz)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(data.table)
library(tidyverse)
library(ggnewscale)
library(ggridges)

# setwd("C:\\shangke\\lession7")

counts <- read.csv("airway_scaledcounts.csv")
metadata <-  read.csv("airway_metadata.csv")

BiocManager::install("DESeq2")
library(DESeq2)
# DESeqDataSetFromMatrix的countData参数为data.frame,第一列为基因名
dds <- DESeqDataSetFromMatrix(countData=counts, 
                              colData=metadata, 
                              design=~dex, 
                              tidy=TRUE)
# colData=metadata第一列为样本名
dds <- DESeq(dds)
res1 <- results(dds,tidy=TRUE)

anno <- read.csv("annotables_grch38.csv")
anno=anno[,c(1,3)]
res2=merge(anno,res1,by=1)

res3=distinct(res2,symbol,.keep_all = T)
res3=na.omit(res3)

# log2(1.5)   至少大于1.2
###基础表达量筛选
###一般来说 count数要在5--10 
select.FPKM <- (res3$baseMean > 10)
table(select.FPKM)
select.log2FC <- abs(res3$log2FoldChange) > 1
table(select.log2FC)
select.qval <- (res3$pvalue < 0.01)
table(select.qval)

select.vec=select.FPKM & select.log2FC & select.qval 

table(select.vec)

degs.list=as.character(res3$symbol)[select.vec]


############################################################################
##david富集分析

# getwd()
# a=degs.list[1:40]

# b<-a[-which(is.na(a))]
# write.table(b,file="deg.gene_symbol.txt",
            row.names = F,quote = F,col.names = F)

###quote指定是否为字符型变量添加引号""






#################################################################################
keytypes(org.Hs.eg.db)
?enrichGO

erich.go.BP = enrichGO(gene =degs.list,
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)


barplot(erich.go.BP,showCategory = 8)
#横轴为基因个数，纵轴为富集到的GO Terms的描述信息
#颜色对应p.adjust值，红色p值小，蓝色p值大
#showCategory指定展示的GO Terms的个数，默认为10，即p.adjust最小的10个
dotplot(erich.go.BP,showCategory = 15,x = 'Count')
#横轴为GeneRatio，代表该GO term下富集到的基因个数占列表基因总数的比例
#纵轴为富集到的GO Terms的描述信息，showCategory指定展示的GO Terms的个数



erich.go.BP=erich.go.BP@result
# write.table(erich.go.BP,"erich.go.BP.deg.con.hf.txt",sep = "\t",col.names = NA)

dotplot(kegg, showCategory=20) #气泡图
barplot(kegg,showCategory=20,drop=T) #柱状图

enrichMap(erich.go.BP)

erich.go.CC = enrichGO(gene =degs.list,
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont = "CC",
                       pvalueCutoff = 0.5,
                       qvalueCutoff = 0.5)

dotplot(erich.go.CC,showCategory = 8)
barplot(erich.go.CC,showCategory = 8)
#有向无环图 GO DAG graph
#investigate how the significant GO terms are distributed over the GO graph. 
#The goplot function shows subgraph induced by most significant GO terms.
goplot(erich.go.CC,showCategory=5)
##GO分析画成树形图，可以更加帮助我们理解。
plotGOgraph(erich.go.CC) 


erich.go.MF = enrichGO(gene =degs.list,
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont = "MF",
                       pvalueCutoff = 0.5,
                       qvalueCutoff = 0.5)

dotplot(erich.go.MF,showCategory = 8)

barplot(erich.go.MF,showCategory = 8)
##GO terms关系网络图 Enrichment Map
##GO terms组织在一个有向无环图（directed acyclic graph）中，terms之间的边表示父子关系（parent-child relationship）。
erich.go.MF <- pairwise_termsim(erich.go.MF)
emapplot(erich.go.MF,showCategory=10)

#对于富集到的GO terms之间的基因重叠关系进行展示
#每个节点是一个富集到的GO term，默认画top30个富集到的GO terms
#节点大小对应该GO terms下富集到的基因个数，节点的颜色对应p.adjust的值，红色小蓝色大
#如果两个GO terms的差异基因存在重叠，说明这两个节点存在overlap关系，用线条连接起来

#GO term与差异基因关系网络图 Gene-Concept Network
cnetplot(erich.go.MF,showCategory=5)
#对于基因和富集的GO terms之间的对应关系进行展示
#图中灰色的点代表基因，黄色的点代表富集到的GO terms
#如果一个基因位于一个GO Terms下，则将该基因与GO连线
#黄色节点的大小对应富集到的基因个数，默认画top5富集到的GO terms

#圆形布局，给线条上色
cnetplot(erich.go.MF,showCategory=10,foldChange=degs.list,circular=TRUE,colorEdge=TRUE)

p1 <- cnetplot(erich.go.MF,showCategory=2,node_label="category")
p2 <- cnetplot(erich.go.MF,showCategory=2,node_label="gene") 
p3 <- cnetplot(erich.go.MF,showCategory=2,node_label="all") 
p4 <- cnetplot(erich.go.MF,showCategory=2,node_label="none") 
cowplot::plot_grid(p1, p2, p3, p4, ncol=2, labels=LETTERS[1:4])
###################################################################################
##所以透露全部展示

egoall <- enrichGO(degs.list, OrgDb=org.Hs.eg.db, ont='ALL',
                   pAdjustMethod='BH', pvalueCutoff=0.05, 
                   qvalueCutoff=0.2, keyType='SYMBOL')

head(egoall,1);dim(egoall)
dotplot(egoall,title='Top5 GO terms of each sub-class',
        showCategory=5,split='ONTOLOGY')+ 
  facet_grid(ONTOLOGY~.,scale="free")


#####################################################################################
##kegg富集分析
keytypes(org.Hs.eg.db)
# 转换ID
DEG.entrez_id = mapIds(x = org.Hs.eg.db,
                       keys =  degs.list,
                       keytype = "SYMBOL",
                       column = "ENTREZID")
?enrichKEGG
erich.kegg.res <- enrichKEGG(gene = DEG.entrez_id,
                             organism = "hsa",
                             keyType = "kegg")

dotplot(erich.kegg.res)

browseKEGG(erich.kegg.res, 'hsa04216')
zzh=erich.kegg.res@result
write.table(zzh,"erich.kegg.res.6.19.txt",sep = "\t",col.names = NA)

kegg=read.table("自定义.txt",header = T,sep = "\t")
View(kegg)

k = data.frame(kegg)
library(ggplot2)
library(dplyr)
before <- as.numeric(sub("/\\d+$", "", k$GeneRatio))
after <- as.numeric(sub("^\\d+/", "", k$GeneRatio))
k$GeneRatio = before /after
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

###################################################################################





###enrichKEGG使用在线数据速度实在是太慢了，所以可以先使用createKEGGdb生成本地KEGG.db包。 
install.packages("remotes")
remotes::install_github("YuLab-SMU/createKEGGdb")
library(createKEGGdb)
###下载kegg数据库
create_kegg_db("hsa")
install.packages("KEGG.db_1.0.tar.gz",repos=NULL,type="source")
library(KEGG.db)
ekegg <- enrichKEGG(DEG.entrez_id, organism='hsa',keyType="kegg",
                    pvalueCutoff=0.5,pAdjustMethod='BH',qvalueCutoff=0.5,
                    use_internal_data=T)

##将通路中的gene名转化
ekegg.res=ekegg@result
z1=strsplit(ekegg.res$geneID[1], "/",fixed=TRUE)
z2=z1[[1]]
z3 = mapIds(x = org.Hs.eg.db,
            keys = z2,
            keytype = "ENTREZID",
            column = "SYMBOL")


####################################################################################
#####GSEA

select.FPKM <- (res2$baseMean > 10)
table(select.FPKM)
rownames(res2)
gs.list=as.character(res2$symbol)[select.FPKM]
res3=distinct(res2,symbol,.keep_all = T)
rownames(res3)=res3$symbol
f1=res3[gs.list,]
exp.fc=f1[,c(2,4)]

###id转化
# mapID,bitr,soft
expm.id <- bitr(exp.fc$symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
head(expm.id)
exp.fc.id <- merge(exp.fc, expm.id,by.x="symbol", by.y="SYMBOL", all=F)
head(exp.fc.id)
exp.fc.id=na.omit(exp.fc.id)
exp.fc.id.sorted <- exp.fc.id[order(exp.fc.id$log2FoldChange, decreasing = T),]
head(exp.fc.id.sorted)
id.fc <- exp.fc.id.sorted$log2FoldChange
names(id.fc) <- exp.fc.id.sorted$ENTREZID
head(id.fc)
gsea <- gseKEGG(id.fc, organism = "hsa",pvalueCutoff = 0.5)
# gsgo <- gseGO(geneList = id.fc, ont = 'BP', OrgDb = org.Hs.eg.db, pvalueCutoff = 0.5)
####查看一下富集结果
View(gsea@result)
# NES一般认为abs大于1.5, pvalue小于0.01
?gseaplot2
# gseaplot2(gsea, geneSetID = 1,title = "Focal adhesion", pvalue_table = T,)
# gseaplot2(gsea,2:4,subplots = 1)
gseaplot2(gsea, geneSetID = c(2,3))

##峰峦图
ridgeplot(gsea)

##UpSet Plot
#着重于不同基因集间基因的重叠情况
upsetplot(erich.go.MF)
#对于GSEA结果将绘制不同类别的fold change分布
upsetplot(gsea) 

##Heatmap-like functional classification
heatplot(erich.go.MF)
heatplot(gsea,foldChange=id.fc)

###
p1 <- gseaplot(gsea,geneSetID=1,by="runningScore",title=gsea$Description[1])
p2 <- gseaplot(gsea,geneSetID=1,by="preranked",title=gsea$Description[1])
p3 <- gseaplot(gsea,geneSetID=4,title=gsea$Description[4])
cowplot::plot_grid(p1, p2, p3, ncol=2, labels=LETTERS[1:3])

dotplot(gsea)
####################################################################################
###msigdbr | 多个物种的通路数据集合
select.FPKM <- (res2$baseMean > 10)
table(select.FPKM)
rownames(res2)
degs.list=as.character(res2$symbol)[select.FPKM]
res3=distinct(res2,symbol,.keep_all = T)
rownames(res3)=res3$symbol
f1=res3[degs.list,]
exp.fc=f1[,c(2,4)]


deg = exp.fc$log2FoldChange
names(deg) = exp.fc$symbol
deg= sort(deg,decreasing = T)

h.kegg= msigdbr(species ="Homo sapiens", category = 'C2', subcategory = 'KEGG')
length(unique(h.kegg$gs_exact_source))
h.kegg= h.kegg[,c(3,4)]
gsea <- GSEA(deg,pvalueCutoff = 0.5, TERM2GENE =h.kegg)
gseaplot2(gsea, geneSetID = 6, title = gsea$Description[6])
gsea1=gsea@result

write.table(gsea1,"gsea.res.txt",sep = "\t",col.names = NA)

#####################################################################################
a=degs.list[1:40]

b<-a[!is.na(a)]
write.table(b,file="deg.gene_symbol.txt",
            row.names = F,quote = F,col.names = F)
getwd()
# 蛋白互做图
##https://string-db.org/cgi/network?pollingId=bUpew6ZPM4MT&sessionId=bq3ONpHVnJ16

library(ggraph)
help(package="ggraph")
library(igraph)
nodes<-read.csv("deg.gene_symbol.txt",header=F)
links<-read.table("string_interactions_short.tsv",header=F,sep="\t")
nodes$Name<-nodes$V1
nodes
links
net<-graph_from_data_frame(d=links,vertices=nodes,directed = T)
# plot(net)
# ggraph(net,layout = "linear",circular=T)+
#   geom_edge_link(color="blue")+
#   geom_node_text(aes(label=Name))+
#   geom_node_point(size=10,color="red",alpha=0.5)+
#   theme_void()

?ggraph

## layout = circlepack
# ggraph(net,layout = "treemap")+
#   geom_edge_link(color="blue")+
#   geom_node_text(aes(label=Name))+
#   geom_node_point(size=10,color="red",alpha=0.5)+
#   theme_void()

# 几个layout比较下来还是layout="kk"好看一点。
# ggraph(net,layout = "kk")+
#   geom_edge_link(color="blue")+
#   geom_node_text(aes(label=Name))+
#   geom_node_point(size=10,color="red",alpha=0.5)+
#   theme_void()

# 接下来试着添加基因名字并按照上调和下调添加颜色
library(ggraph)
help(package="ggraph")
library(igraph)
nodes
links
nodes$Name<-nodes$V1
nodes$Group<-c(rep("Up",5),rep("Down",5))
net<-graph_from_data_frame(d=links,vertices=nodes,directed = T)
ggraph(net,layout = "kk")+
  geom_edge_link(color="yellow")+
  geom_node_point(size=10,aes(color=Group),alpha=0.5)+
  geom_node_text(aes(label=Name))+
  theme_void()

################################################################################
# 转录因子富集
library(RcisTarget)
geneLists <- list(hypoxia= degs.list ) 
geneLists
data(motifAnnotations_hgnc)
motifAnnotations_hgnc

#### https://resources.aertslab.org/cistarget/
featherURL="https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-500bp-upstream-7species.mc9nr.feather"
?download.file()

download.file(featherURL, destfile="abcd")
motifRankings <- importRankings("H:/rcistarget/hg19-tss-centered-10kb-7species.mc9nr.feather")
# Motif enrichment analysis:
motifEnrichmentTable_wGenes <- cisTarget(geneLists, motifRankings,
                                         motifAnnot=motifAnnotations_hgnc)
motifs_AUC <- calcAUC(geneLists, motifRankings)
motifEnrichmentTable <- addMotifAnnotation(motifs_AUC, 
                                           motifAnnot=motifAnnotations_hgnc)
motifEnrichmentTable_wGenes <- addSignificantGenes(motifEnrichmentTable, 
                                                   geneSets=geneLists,
                                                   rankings=motifRankings, 
                                                   nCores=1,
                                                   method="aprox")
write.table(motifEnrichmentTable_wGenes,"10-0-deg.537gene.auc.txt",sep = "\t")
