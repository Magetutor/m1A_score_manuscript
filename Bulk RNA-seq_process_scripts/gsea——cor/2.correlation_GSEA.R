# 基于单基因批量相关性分析的GSEA #
## 在相关性分析时做过一对多基因的相关性，
## 但是这个方法有个小缺陷，并不知道最后富集的通路是正向影响还是反向影响。
## 也就是无法判断方向。判断方向的工具也不是没有，GSEA就是一个。
## 所以，我想能不能把批量相关性分析和GSEA结合起来。
## GSEA需要的gene set是现成的没有问题，但是genelist没有，
## 这里我们可以把所有基因跟单个基因的相关性系数当做LogFC，有正有负，
## 就解决了geneList的问题。
## 输入的为基因表达矩阵——vst/TPM
## 行名为基因名，列名为样本名
load(file = "BRCA_mRNA_exprSet.Rdata")
exprSet <- spinal_Exp_vst
test <- exprSet[1:10,1:10]

#这个函数只要输入一个基因，他就会批量计算这个基因跟其他编码基因的相关性，返回相关性系数和p值。
batch_cor <- function(gene){
  y <- as.numeric(exprSet[gene,])
  rownames <- rownames(exprSet)
  do.call(rbind,future_lapply(rownames, function(x){
    dd  <- cor.test(as.numeric(exprSet[x,]),y,method="spearman")
    data.frame(gene=gene,genelist=x,cor=dd$estimate,p.value=dd$p.value )
  }))
}
# 以PCDC1这个基因为例
library(future.apply)
plan(multiprocess)
system.time(dd <- batch_cor("Trmt10c"))
# 制作genelist
gene <- dd$genelist
## 转换
library(clusterProfiler)
gene = bitr(gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db") #鼠：org.Mm.eg.db
## 去重
gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)

gene_df <- data.frame(logFC=dd$cor,
                      SYMBOL = dd$genelist)
gene_df <- merge(gene_df,gene,by="SYMBOL")

## geneList 三部曲
## 1.获取基因logFC
geneList <- gene_df$logFC
## 2.命名
names(geneList) = gene_df$SYMBOL
## 3.排序很重要
geneList = sort(geneList, decreasing = TRUE)

library(clusterProfiler)
## 读入hallmarks gene set，从哪来？
hallmarks <- read.gmt("data/m5.go.v2023.2.Mm.symbols (1).gmt")
# 需要网络
y <- GSEA(geneList,TERM2GENE =hallmarks)

### 看整体分布
library(ggplot2)
dotplot(y,showCategory=12,split=".sign")+facet_grid(~.sign)
### 特定通路作图
yd <- data.frame(y)
library(enrichplot)
gseaplot2(y,"HALLMARK_INTERFERON_ALPHA_RESPONSE",color = "red",pvalue_table = T)
