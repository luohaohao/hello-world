
rm(list = ls())

## Installing R packages
bioPackages <-c( 
  "corrplot","ggrepel", #绘制相关性图形
  "stringr", #处理字符串的包
  "readr","tximport","dplyr", #处理salmon表达量的扩展包
  "FactoMineR","factoextra", #PCA分析软件
  "limma","edgeR","DESeq2", #差异分析的三个软件包
  "clusterProfiler", "org.Hs.eg.db", #安装进行GO和Kegg分析的扩展包
  "GSEABase","GSVA", #安装进行GSEA分析的扩展包
  "airway" # 包含数据集的bioconductor软件包
  )


## If you are in China, run the command below
local({
  r <- getOption( "repos" );# set CRAN mirror for users in China
  r[ "CRAN" ] <- "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"; # CRAN的镜像地址
  options( repos = r )
  
  BioC <- getOption( "BioC_mirror" ); # set bioconductor mirror for users in China
  BioC[ "BioC_mirror" ] <- "https://mirrors.ustc.edu.cn/bioc/"; # bioconductor的镜像地址
  options( BioC_mirror = BioC )
})

# 检查是否设定完毕
getOption("BioC_mirror")
getOption("CRAN")

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# 安装devtools管理github上的软件包
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")


## Installing missing packages
lapply( bioPackages, 
        function( bioPackage ){
          if(!bioPackage %in% rownames(installed.packages())){
              CRANpackages <- available.packages()

              if(bioPackage %in% rownames(CRANpackages)){
                install.packages( bioPackage)
              }else{
                  BiocManager::install(bioPackage,suppressUpdates=F,ask=F)
              }
          }
        })


## 验证R扩展包是否安装成功
library(limma)
library(edgeR)
library(DESeq2)
library(FactoMineR)
library(factoextra)
library(clusterProfiler)
library(org.Hs.eg.db)

# 不显示加载信息
suppressMessages(library(limma))





