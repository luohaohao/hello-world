
rm(list = ls())
options(stringsAsFactors = F)

library(clusterProfiler)
library(org.Hs.eg.db)

# 读取3个软件的差异分析结果
load("data/Step01-airwayData.Rdata")
load(file = "data/Step03-limma_voom_nrDEG.Rdata")
load(file = "data/Step03-DESeq2_nrDEG.Rdata")
load(file = "data/Step03-edgeR_nrDEG.Rdata")
ls()

# 提取所有差异表达的基因名
limma_sigGene <- DEG_limma_voom[DEG_limma_voom$regulated!="normal",1]
edgeR_sigGene <- DEG_edgeR[DEG_edgeR$regulated!="normal",1]
DESeq2_sigGene <- DEG_DESeq2[DEG_DESeq2$regulated!="normal",1]


# 根据需要更改DEG的值
DEG <- limma_sigGene
gene_all <- rownames(filter_count)


#### 第一步，从org.Hs.eg.db提取ENSG的ID 和GI号对应关系
keytypes(org.Hs.eg.db)

# bitr in clusterProfiler     biological ID translator using OrgDb
allID <- bitr(gene_all, fromType = "ENSEMBL", toType = c( "ENTREZID" ), OrgDb = org.Hs.eg.db )
degID <- bitr(DEG, fromType = "ENSEMBL", toType = c( "ENTREZID" ), OrgDb = org.Hs.eg.db )
head(degID)


# KEGG analysis----
enrich <- enrichKEGG(gene =degID[,2],organism='hsa',universe=allID[,2],pvalueCutoff=1,qvalueCutoff=1)

GeneRatio <- GeneRatio <- as.numeric(apply(str_split(enrich$GeneRatio,'/',simplify = T),1,function(x) as.numeric(x[1])/as.numeric(x[2])))
BgRatio <- as.numeric(lapply(strsplit(enrich$BgRatio,split="/"),function(x) as.numeric(x[1])/as.numeric(x[2])  ))
enrich_factor <- GeneRatio/BgRatio
out <- data.frame(enrich$ID,enrich$Description,enrich$GeneRatio,enrich$BgRatio,round(enrich_factor,2),enrich$pvalue,enrich$qvalue,enrich$geneID)
colnames(out) <- c("ID","Description","GeneRatio","BgRatio","enrich_factor","pvalue","qvalue","geneID")
write.table(out,"deg_analysis/trut_VS_untrt_enrich_KEGG.xls",row.names = F,sep="\t",quote = F)

out_sig0.05 <- out[out$qvalue<0.05,]

# barplot
bar <- barplot(enrich,showCategory=10,title="KEGG Pathway",colorBy="p.adjust")
bar

# 保存
pdf(file = "deg_analysis/kegg_bar_plot.pdf",width = 8,height = 6)
print(bar)
dev.off()

# dotplot
dot <- dotplot(enrich,x="geneRatio",showCategory=10,font.size=12,title="KEGG Pathway")
dot

# 保存
pdf(file = "deg_analysis/kegg_dot_plot.pdf",width = 8,height = 6)
print(dot)
dev.off()


# GO 
enrich <- enrichGO(gene =degID[,2],OrgDb='org.Hs.eg.db',ont="ALL",universe=allID[,2],pvalueCutoff=1,qvalueCutoff=1)

GeneRatio <- as.numeric(lapply(strsplit(enrich$GeneRatio,split="/"),function(x) as.numeric(x[1])/as.numeric(x[2])))
BgRatio <- as.numeric(lapply(strsplit(enrich$BgRatio,split="/"),function(x) as.numeric(x[1])/as.numeric(x[2])))
enrich_factor <- GeneRatio/BgRatio
out <- data.frame(enrich$ID,enrich$Description,enrich$GeneRatio,enrich$BgRatio,round(enrich_factor,2),enrich$pvalue,enrich$qvalue,enrich$geneID)
colnames(out) <- c("ID","Description","GeneRatio","BgRatio","enrich_factor","pvalue","qvalue","geneID")
write.table(out,"deg_analysis/trut_VS_untrt_enrich_GO.xls",row.names = F,sep="\t",quote = F)

# barplot
bar <- barplot(enrich,showCategory=10,title="Go",colorBy="p.adjust")
bar

# 保存
pdf(file = "deg_analysis/BP_bar_plot.pdf",width = 6,height = 6)
print(bar)
dev.off()

# dotplot
dot <- dotplot(enrich,x="geneRatio",showCategory=10,font.size=12,title="Biological Pathway")
dot

# 保存
pdf(file = "deg_analysis/BP_dot_plot.pdf",width = 6,height = 6)
print(dot)
dev.off()



#获取某通路基因名
if(!requireNamespace("BiocManager", quietly = TRUE)){install.packages("BiocManager")}
#BiocManager::install("KEGGREST", version = "3.10") 
library(KEGGREST)
listDatabases() 
gs<-keggGet('hsa04728')
gs
#获取通路中gene信息 
gs[[1]]$GENE 
#查找所有基因 
genes<-unlist(lapply(gs[[1]]$GENE,function(x) strsplit(x,';'))) 
genelist <- genes[1:length(genes)%%3 ==2] 
genelist <- data.frame(genelist)  
#把结果写入表格中 
write.table(genelist, "hsa04728.csv",            
            row.names=FALSE,col.names=TRUE,sep=",") 


