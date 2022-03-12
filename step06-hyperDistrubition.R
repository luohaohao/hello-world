# 2017.8.9 BY Juan Zhang

rm(list=ls())
options(stringsAsFactors = F)


## 加载差异分析结果
lnames <- load("../Analysis/deg_analysis/Step03-limma_voom_nrDEG.Rdata")
lnames

# 提取差异表达的gene
DEG <- DEG_limma_voom$Geneid[which(DEG_limma_voom$regulated!="normal")]

# 将symbol转换成gene id
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)

library(clusterProfiler)
id2ENTREZ <- bitr(DEG, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db )
DEGs <- unique(id2ENTREZ$ENTREZID)

# 分析时所用表达谱所有基因
id2ENTREZ <- bitr(DEG_limma_voom$Geneid, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db )
ref_gene <- unique(id2ENTREZ$ENTREZID)



## 加载通路数据，
# kegg_gid为KEGG数据库里面所有基因
# non_dis_pathway，通路数据，即gene set，每一个通路由一组功能相关的基因组成一个基因集合，行使特定生物学功能。
lnames <- load("../Analysis/data/kegg_pathway240_2017.6.Rdata") 
lnames
class(kegg_gid)
kegg_gid




## 实现超几何检验原理
diff_gene <- intersect(kegg_gid,DEGs) # 与数据库取交集后的差异表达基因
backg_gene <- intersect(kegg_gid,ref_gene) # 数据库中的基因与表达谱所有基因取交集作为背景

K <- length(diff_gene) # 差异表达基因数L2
N <- length(backg_gene) # 背景基因数L


p_value <- NULL
path_hsa <- NULL
path_name <- NULL

for(i in 1:length(non_dis_pathway)){  
   path <- non_dis_pathway[[i]]
   path_hsa[i] <- path$pathwayId  # 得到此通路的KEGG ID
   path_name[i] <- path$pathwayName  # 得到通路的名字
   M <- length(intersect(path$genesInPathway,backg_gene)) # 一条通路中的基因数L1
   X <- length(intersect(diff_gene,path$genesInPathway)) # 这条通路中感兴趣的基因数K
   if(X==0){
     p_value[i] <- 1
    }else{
         p_value[i] <- (1-phyper(X-1,M,N-M,K))
        }
}

fdr <- p.adjust(p_value,method="BH",length(p_value))
res <- data.frame(path_hsa,path_name,p_value,fdr)
res <- res[order(res$p_value),]
write.table(res,file="../Analysis/hyperDistribution_DEG_KEGG.xls",sep="\t",row.name=F,quote=F)
