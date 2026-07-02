library(tibble)
library(tidyr)
library(clusterProfiler)
load("E:/R stuation/py/public_mice/output/2d/diffgene.Rdata")
load("E:/R stuation/py/public_mice/output/3d_35d/diffgene_3d_limma.Rdata")
load("E:/R stuation/py/public_mice/output/3d_35d/diffgene_35d_limma.Rdata")
load("E:/R stuation/py/public_mice/output/7d_2/diffgene_7d_limma.Rdata")
## 2d
diffgene_2d <- diffgene
diffgene_2d$group <- ifelse(diffgene_2d$logFC > 0,"Up","Down")
diffgene_2d <- diffgene_2d %>%
  rownames_to_column("gene_symbol")
mydata_2d <- diffgene_2d[,c("gene_symbol","group")]
group_2d <- data.frame("gene_symbol" = diffgene_2d$gene_symbol,"group" = rep("All",nrow(mydata_2d)))
mydata_2d <- rbind(mydata_2d,group_2d)
mydata_2d <- cbind("cluster" = rep("2d",nrow(mydata_2d)),mydata_2d)                      

## 3d
diffgene_3d$group <- ifelse(diffgene_3d$logFC > 0,"Up","Down")
diffgene_3d <- diffgene_3d %>%
  rownames_to_column("gene_symbol")
mydata_3d <- diffgene_3d[,c("gene_symbol","group")]
group_3d <- data.frame("gene_symbol" = diffgene_3d$gene_symbol,"group" = rep("All",nrow(mydata_3d)))
mydata_3d <- rbind(mydata_3d,group_3d)
mydata_3d <- cbind("cluster" = rep("3d",nrow(mydata_3d)),mydata_3d)

## 7d
diffgene_7d$group <- ifelse(diffgene_7d$logFC > 0,"Up","Down")
diffgene_7d <- diffgene_7d %>%
  rownames_to_column("gene_symbol")
mydata_7d <- diffgene_7d[,c("gene_symbol","group")]
group_7d <- data.frame("gene_symbol" = diffgene_7d$gene_symbol,"group" = rep("All",nrow(mydata_7d)))
mydata_7d <- rbind(mydata_7d,group_7d)
mydata_7d <- cbind("cluster" = rep("7d",nrow(mydata_7d)),mydata_7d)

## 35d
diffgene_35d$group <- ifelse(diffgene_35d$logFC > 0,"Up","Down")
diffgene_35d <- diffgene_35d %>%
  rownames_to_column("gene_symbol")
mydata_35d <- diffgene_35d[,c("gene_symbol","group")]
group_35d <- data.frame("gene_symbol" = diffgene_35d$gene_symbol,"group" = rep("All",nrow(mydata_35d)))
mydata_35d <- rbind(mydata_35d,group_35d)
mydata_35d <- cbind("cluster" = rep("35d",nrow(mydata_35d)),mydata_35d)

## 合并
mydata <- rbind(mydata_2d,mydata_3d,mydata_7d,mydata_35d)

## 加载数据集
m_df <- read.gmt("data/m5.go.bp.v2023.2.Mm.symbols.gmt")
xx <- compareCluster(gene_symbol~group+cluster, 
                     data=mydata, 
                     fun="enricher",
                     TERM2GENE = m_df,
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.05)

## 将通路里面的Description变为首字母大写，去掉GO；添加富集分数
dd <- xx@compareClusterResult
for (i in 1:length(dd$Description)) {
  print(i)
  str = unlist(strsplit(dd$Description[i],split = "_"))[-1]
  str = paste(stringr::str_to_title(str),collapse = " ")
  dd$Description[i] = str
}
head(dd$BgRatio)
as.numeric(sub("/\\d+", "", dd$BgRatio))
dd$richFactor =dd$Count / as.numeric(sub("/\\d+", "", dd$BgRatio))
xx@compareClusterResult = dd
xx@compareClusterResult$cluster <- factor(xx@compareClusterResult$cluster,
                                      levels = c("2d","3d","7d","35d")) #排序
xx@compareClusterResult$group <- factor(xx@compareClusterResult$group,
                                        levels = c("Up","Down","All")) #排序
## 画图
library(ggplot2)
dd <- xx@compareClusterResult
go_n <- go_f[c(grep("mitochondrial",go_f$Description),
          grep("ATP",go_f$Description),
          grep("energy",go_f$Description),
          grep("oxidative phosphorylation",go_f$Description),
          grep("aerobic respiration",go_f$Description),
          grep("respiratory electron transport chain",go_f$Description)),]

ggplot(go_n,aes(x=richFactor, y=Description),
       label_format = 1) + #多个分组时需要选取Cluster
  geom_point(aes(size = Count,color = -log10(p.adjust))) + # 气泡大小及颜色设置
  #facet_grid(~cluster) +
  labs(x = "Rich Factor",
       y = "Description",
       title = "GO:CC Enrichment Dotplot", # 设置坐标轴标题及图标题
       size = "Count",
       fill = "-Log10(Adjusted p value)") +
  theme_bw() + 
  scale_color_distiller(palette = "YlOrBr",direction = 1) + 
  scale_size_continuous(range=c(5,8))+
  theme (text = element_text (size = 15,color = "black"),
         axis.text = element_text(size = 13,colour = "black"),
         axis.title = element_text(size = 15))


aa <- dd


go_n$Description <- factor(go_n$Description,
                         levels = c(go_n$Description[order(go_n$richFactor)]))

ggsave(file = "output/GO_BP/GO_BP.pdf")



