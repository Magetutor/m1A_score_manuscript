#### 富集分析 ####
## 1.基因集网站：Board;Enrichr
## 分析方法ORA/GSEA
### ORA(GO/KEGG)——提取出差异基因来看其分布于哪些通路中
## 优势：提取出上调下调基因分面来看不同
## 局限：人为规定了差异基因的选择范围，会忽略一部分的信息

### GSEA——将差异基因按照logFC的大小排序，看差异分数的正负（正值为目标基因集为上调，负值为目标基因集下调）
## 优势：1.纳入所有基因，没有人为设定阈值
##       2.结果区分正反，有利于指导科学实验
##       3.基因集种类齐全，可自行订制——基因集网站：Board;Enrichr

## 输入的数据
## 差异基因：差异基因的名称
##### 1.KEGG 分析 #####


## 1.1 普通分析——不区分上下调

# 1.1.1 加载差异基因数据或根据allDiff进行筛选
library(clusterProfiler)
load(file = "output/diffgene.Rdata")
#diffgene <- allDiff %>% 
#  filter(gene_symbol !="") %>% 
#  filter(adj.P.Val < 0.05) %>% 
#  filter(abs(logFC) > 1) %>%
#  filter(cluster == "7d") %>%
#  arrange(desc(logFC)) #从大到小排序，后面提取上调/下调基因时也是这个顺序
### 这个分析需要什么数据？——差异基因的名称
### 将所需要的差异基因提取出来，那么分析中的想要通路富集比例会升高
### 获得基因列表
gene <- unlist(cor_data_sig %>%
  filter(cor > 0.7) %>%
  dplyr::select(gene2))
#基因名称转换，返回的是数据框-需要的是基因ENRIZID名称
gene = bitr(gene, 
            fromType="SYMBOL", 
            toType="ENTREZID", 
            OrgDb="org.Mm.eg.db") #鼠：org.Mm.eg.db
head(gene)
which(is.na(gene$ENTREZID)) #查找缺失值的位置
# 如果你想把KEGG 本地化——不建议，一直在更新
# organism = 'hsa'
# http://rest.kegg.jp/list/organism
# https://mp.weixin.qq.com/s/PwrdQAkG3pTlwMB6Mj8wXQ
# library(createKEGGdb)
# species <-c("mmu","hsa")
# create_kegg_db(species)
# library(KEGG.db)
# EGG <- enrichKEGG(gene = gene$ENTREZID,
#                   organism = 'hsa',
#                   pvalueCutoff = 0.05,
#                   use_internal_data =T)

## 在线分析
EGG <- enrichKEGG(gene = gene$ENTREZID,
                  organism = 'mmu', #"mmu"
                  pvalueCutoff = 0.05)


## 画图 ?可以显示函数的相关调节参数
library(ggplot2)
barplot(EGG)
dotplot(
  EGG,
  x = "GeneRatio",#"Count"
  color = "p.adjust",
  showCategory = 20, #展示出的通路结果—— test <- as.data.frame(EGG) 多少行就是多少通路
  font.size = 12,
  title = "",
  orderBy = "x",
  label_format = 30,
) +   #颜色变化
  scale_colour_gradient(  
    low = "#23A8DF",
    high = "#F18F17"
  )

## 分析出的差异基因富集通路为处理后的样本影响——可以将无关的通路删除，不展示出来
### KEGG的富集分析比较特殊，他的背后是个网站 
KEGG_df <- as.data.frame(EGG)
symboldata <- setReadable(EGG, OrgDb="org.Hs.eg.db", keyType = "ENTREZID") #基因ID转化为gene symbol
symboldata  <- as.data.frame(symboldata)

#输入标号可以看到在网站中总结好的某一通路中的差异基因,但是忽略了变化倍数与基因的上下调
browseKEGG(EGG, 'hsa04110') 
save(EGG,file = "output/EGG.Rdata")





##### 2.GO 分析 #####
library(clusterProfiler)
load(file = "data/diffgene.Rdata")
### 这个分析需要什么数据？
### 获得基因列表
gene <- rownames(diffgene)
#基因名称转换，返回的是数据框
gene = bitr(gene, 
            fromType="SYMBOL", 
            toType="ENTREZID", 
            OrgDb="org.Hs.eg.db")
head(gene)
#GO分析的细胞组分 CC
ego_CC <- enrichGO(gene = gene,
                   OrgDb= "org.Mm.eg.db",
                   keyType= "SYMBOL",
                   ont = "CC",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = TRUE)

#**作图**
#条带图
barplot(ego_CC)
# 点图
dotplot(ego_CC)
# GO 作图
goplot(ego_CC)
### 重点是如何解释图：例如富集到的位置在染色质上，那么可以推断其主要发生的位置在细胞核内。

#GO分析的生物过程BP
ego_BP <- enrichGO(gene = gene,
                   OrgDb= "org.Mm.eg.db",
                   keyType= "SYMBOL",
                   ont = "BP",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = TRUE)

#**作图**
#条带图
barplot(ego_BP)
# 点图
dotplot(ego_BP)
# GO 作图
goplot(ego_BP)

#GO分析分子功能MF：

ego_MF <- enrichGO(gene = gene$ENTREZID,
                   OrgDb= "org.Hs.eg.db",
                   ont = "MF",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.01,
                   readable = TRUE)

#**作图**
#条带图
barplot(ego_MF)
# 点图
dotplot(ego_MF)
# GO 作图
goplot(ego_MF)

############################################
### 尝试新的GO聚类的方法
### 三个一起做，需要时间，并且如果其中一个没富集出来的话就没有图
### 其中一个方面没富集出来的话，可以分开来做
go <- enrichGO(gene = gene$ENTREZID, OrgDb = "org.Hs.eg.db", ont="all")
#save(go,file = "output/go.Rdata")
library(ggplot2)
p <- barplot(go, label_format = 60,split="ONTOLOGY") +
  facet_grid(ONTOLOGY~., scale="free")
p
