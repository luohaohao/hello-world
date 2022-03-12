

group_list=c(rep('normal',3),rep('tumour_1',3),rep('tumour_2',3),rep('tumour_3',3))
print(group_list)

group_list1 <-factor(group_list)
design=model.matrix(~group_list)
colnames(design) <- levels(group_list1)
rownames(design) <- paste0("sample",1:12)
View(design)
##这跟design的设置有关，design=model.matrix(~group_list)会默认第一列为“截距”
#而后面的分组都会跟这个“截距”进行对比，这个“截距”会被默认为分母
##而design的列名的顺序是取决于group_list的level顺序
#######如果想改变截距怎么办呢？
group_list2 <- factor(group_list,levels = c('tumour_1','normal','tumour_2','tumour_3'))
design1=model.matrix(~group_list2)
colnames(design1) <- levels(group_list2)
rownames(design1) <- paste0("sample",1:12)
View(design)
##这时的截距就变成了tumor_1
allDiff=topTable(fit2,adjust='fdr',coef="tumour_1",number=Inf)
allDiff=topTable(fit2,adjust='fdr',coef=2,number=Inf)  
###这里coef=2也可以写成coef="tumour_1"
#topTable函数的coef参数，coef=2是指design的第2列，即tumour_1，即把tumour_1与normal进行对比
##coef=2是指design1的第3列，即把tumour_2和tumour_1进行对比

###########################################################################################

group_list3 <-factor(group_list)
design3 <- model.matrix(~0+group_list3)
colnames(design3) <- levels(group_list3)
rownames(design3) <- paste0("sample",1:12)
View(design3)
library(edgeR)
contrast.matrix <- makeContrasts(contrasts=c('tumour_2-normal',
                                             'tumour_1-normal',
                                             'tumour_2-tumour_1'),levels = design3)

fit=lmFit(r,design)
fit1 <- contrasts.fit(fit, contrast.matrix)
fit2=eBayes(fit1) 
allDiffc=topTable(fit2,adjust='fdr',coef=2,number=Inf,p.value=0.05)
#coef=2，即指contrast.matrix的第2列
#所以这个方法是自由度最大的一种