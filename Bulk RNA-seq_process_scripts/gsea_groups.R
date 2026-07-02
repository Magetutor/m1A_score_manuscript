library(clusterProfiler)
library(org.Mm.eg.db)
library(GseaVis)
#获得基因列表
gene <- rownames(allDiff_7d)
## 转换
gene = bitr(gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")#鼠：org.Mm.eg.db

## 去重
gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)
gene_df <- data.frame(logFC=allDiff_7d$logFC,
                      SYMBOL = rownames(allDiff_7d))
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
hallmarks <- read.gmt("data/m5.go.bp.v2023.2.Mm.symbols.gmt")
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
gseaNb(object= y,geneSetID="GOBP_MACROPHAGE_ACTIVATION",
       subPlot=3,
       addPval=T,
       pvalX=0.95,
       pvalY=0.8)

### cutting edge作图
### 看不同通路中的共有的基因
if(!requireNamespace("ggnewscale",quietly = TRUE)){
  options("repos"=c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
  install.packages("ggnewscale",update = F,ask = F)
}
cnetplot(y,showCategory = 4,foldChange = geneList,colorEdge = T)

