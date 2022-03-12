
rm(list = ls())
options(stringsAsFactors = F)

read <- read.table("../Analysis/data/SRR1039510.seq",header = F,sep = "\t",stringsAsFactors = )

library(stringr)
gc <- NULL
for(i in 1:nrow(read)){
  gc[i] <- round(sum(str_count(read$V1[i], c("G", "C")))/nchar(read$V1[i]),4)*100
  print(i)
}

plot(density(gc,bw = 1))

read$gc <- gc
write.table(read,"data/rawdata/SRR6236728_1.fq.addGC.xls",col.names = F,row.names = F,sep = "\t")

temp <- "TGGGAGGCTGAGGCAGGAGAATCACTTAAACCTGGGAGGCAGAGGTTACAGTGAGCCGAGATT"
round(sum(str_count(temp, c("G", "C")))/nchar(temp),4)*100


id <- NULL
name <- NULL

for(i in 1:nrow(gff_gene)){
  temp <- unlist(strsplit(gff_gene[i,9],";"))
  id[i] <- substr(temp[1],9,nchar(temp[1]))
  name[i] <- substr(temp[2],6,nchar(temp[2]))
  print(i)
}

id2name <- cbind(id,name)











