#### GSEA ####
### GSEA——将差异基因按照logFC的大小排序，看差异分数的正负（正值为目标基因集为上调，负值为目标基因集下调）
## 优势：1.纳入所有基因，没有人为设定阈值
##       2.结果区分正反，有利于指导科学实验
##       3.基因集种类齐全，可自行订制——基因集网站：Board;Enrichr
## 输入的数据——向量
## 向量内容：从大到小排列的logFC
## 向量名称：基因的ENTREZID
library(clusterProfiler)
load(file = "output/allDiff.Rdata") #只需要将所有的差异基因名称按照大小排序就可以，不用关心上调还是下调
#获得基因列表
gene <- rownames(allDiff)
## 转换
gene = bitr(gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")#鼠：org.Mm.eg.db

## 去重
gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)
gene_df <- data.frame(logFC=allDiff$logFC,
                      SYMBOL = rownames(allDiff))
gene_df <- merge(gene_df,gene,by="SYMBOL")

## geneList 三部曲
## 1.获取基因logFC
geneList <- gene_df$logFC
## 2.命名
names(geneList) = gene_df$SYMBOL
## 3.排序很重要
geneList = sort(geneList, decreasing = TRUE)

head(geneList)

## 上方为数据处理步骤，以下的内容就是流程性步骤
################################################################
#################################################################
### GSEA 变化在于gene set
### 选择不同的gene set也就对应不同的分析结果——构建自己的基因集
### 选择自己实验室可以做的表型（铁死亡、双硫死亡、GPCR）——看处理后的基因表型变化
### 缩小要做的基因范围——1.通过看某个基因在多个通路中均有表达变化；2.通过看MAplot的表达量；3.通过转录因子下游基因的集体变化来判断
#################################################################
### 1.hallmarks gene set——此基因集是与肿瘤最相关的50个通路
# 不同基因集对应的不同基因特点
# https://mp.weixin.qq.com/s/YVNOeYEIvh7adXb4nSOO5A
## 读入hallmarks gene set，从哪来？
## GSEA-download-Molecular Signatures Database
hallmarks <- read.gmt("resource/h.all.v2022.1.Hs.symbols.gmt")
table(hallmarks$term) #看每个基因集中对应多少基因
length(table(hallmarks$term)) #看一共有多少个基因集

### 主程序GSEA
y <- GSEA(geneList,TERM2GENE =hallmarks)

yd <- as.data.frame(y)
### 看整体分布
library(ggplot2)
library(aPEAR)
enrichmentNetwork(y@result, repelLabels = TRUE, drawEllipses = TRUE)
dotplot(y,showCategory=30,
        split=".sign", #按照正负来区分激活/抑制
        font.size = 8,
        label_format = 60)+facet_grid(~.sign)

#但是这种图片需要将每个通路的变化都解释，所以选择有需要的通路单独展示并解释。

#df <- ggplot2::fortify(y, showCategory = 30, split=".sign")
### https://github.com/YuLab-SMU/enrichplot/blob/master/R/method-fortify.R

### 可以修改标签长度
library(stringi)
library(ggplot2)
dotplot(y,showCategory=12,split=".sign")+
  facet_grid(~.sign)+
  scale_y_discrete(labels=function(x) stri_sub(x,10))

### 选择需要呈现的来作图
library(enrichplot)
gseaplot2(y,"HALLMARK_E2F_TARGETS",color = "red",pvalue_table = T)
gseaplot2(y,10,color = "red",pvalue_table = T)
gseaplot2(y, geneSetID = 1:3) #多个通路合并在一起作图

### cutting edge作图
### 看不同通路中的共有的基因
if(!requireNamespace("ggnewscale",quietly = TRUE)){
  options("repos"=c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
  install.packages("ggnewscale",update = F,ask = F)
}
cnetplot(y,showCategory = 4,foldChange = geneList,colorEdge = T)

#################################################################
### 2.kegg 通路
## 读入kegg gene set
kegg <- read.gmt("resource/c2.cp.kegg.v2022.1.Hs.symbols.gmt")
y <- GSEA(geneList,TERM2GENE =kegg)
yd <- as.data.frame(y)
### 看整体分布
dotplot(y,showCategory=12,
        font.size = 9,
        split=".sign",
        label_format = 60)+
  facet_grid(~.sign)
## KEGG的优势：可以看整个通路的基因图，表示出的是上调还是下调

#################################################################
### 3.转录因子
## 读入转录因子,高能，需要专业知识消化
## Chip-seq（染色质免疫共沉淀）-可用于结合转录因子
ENCODE_TF <- read.gmt("resource/ENCODE_TF_ChIP-seq_2015.txt")
tfbs <- GSEA(geneList,TERM2GENE = ENCODE_TF)
tfbsd <- as.data.frame(tfbs)
### 看整体分布
dotplot(tfbs,showCategory=30,split=".sign")+facet_grid(~.sign)

gseaplot2(tfbs,"E2F4 MEL cell line mm9",color = "red",pvalue_table = T)
### 这些gene set从哪来的
### 1.board 官网
### 2.https://amp.pharm.mssm.edu/Enrichr/#stats
### 该部分所有数据集已经下载，在附赠教程中
### 3.自己自作或者收集

#################################################################
## GZ06_必备技能GSEA，富集分析的神器，量化一切通路。 (课程赠送，微店领取) 
## https://weidian.com/item.html?itemID=3749911422
## 其他物种的GSEA，比如老鼠该怎么做呢
### 使用msigdbr包
### https://cran.r-project.org/web/packages/msigdbr/vignettes/msigdbr-intro.html
## 下游分析特别好的教程！！！！
## https://hbctraining.github.io/DGE_workshop_salmon/lessons/functional_analysis_2019.html
#####################################################
### PPi 网络构建，参考附赠视频
### 以一款大杀器，结束GEO的教程
### https://amp.pharm.mssm.edu/Enrichr/

## GEO教程长期更新的链接是这个:
## https://codingsoeasy.com/archives/geo
