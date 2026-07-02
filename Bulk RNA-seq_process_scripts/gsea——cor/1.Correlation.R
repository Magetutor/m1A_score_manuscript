#### 6.一对多基因 ####
### 4.单基因批量相关性分析-附赠课程-guilt by association
### 1.单个基因与其他所有基因检测相关性，推测出与其最相关联的基因
### 2.目标基因与其他基因相关性分析，对与其相关的基因进行功能富集从而间接反应该基因的功能
### 注意:此处技能在后续课程，课上只作为演示！
### 相关性分析单次操作
gene1 = as.numeric(exprSet["Trmt10c",])
gene2 = as.numeric(exprSet["Gbe1",])
aa = cor.test(gene1,gene2,method="spearman")
### 提取p值和相关性系数
aa$p.value
aa$estimate

### 能单次操作就能批量操作
##1.设定容器
correlation <- data.frame()
##2.准备数据
data <- exprSet
##3.获取基因列表
genelist <- colnames(data)
##4.指定基因
gene <- "Trmt10c"
genedata <- as.numeric(data[,gene])
##5.开始for循环
for(i in 1:length(genelist)){
  ## 1.指示
  print(i)
  ## 2.计算
  dd = cor.test(genedata,as.numeric(data[,i]),method="spearman")
  ## 3.填充
  correlation[i,1] = gene
  correlation[i,2] = genelist[i]
  correlation[i,3] = dd$estimate
  correlation[i,4] = dd$p.value
}

colnames(correlation) <- c("gene1","gene2","cor","p.value")
## 6.p值矫正
### 校正后的P值会扩大，其更可靠
correlation$padjust = p.adjust(correlation$p.value,method = "BH")

## 7.筛选p值小于0.05，按照相关性系数绝对值选前500个的基因， 数量可以自己定
library(dplyr)
library(tidyr)
cor_data_sig <- correlation %>% 
  filter(padjust < 0.05 & abs(cor) >= 0.7) 
## 横向柱状图展示
#添加上下调分组标签：
dt <- cor_data_sig
dt$group <- case_when(dt$cor > 0 ~ 'Positive correlation',
                      dt$cor < 0 ~ 'Negative correlation')
dt$gene2 <- factor(dt$gene2,levels = rev(c("Uqcrb","mt-Nd4","Zfp655","Gbe1",
                                       "Dld","Cfl2","Crls1","Trrap","Nisch")))
p <- ggplot(dt,
            aes(x =cor, y = gene2, fill = group)) + #数据映射
  geom_col() + #绘制添加条形图
  theme_bw()
mytheme <- theme(
  legend.position = 'none',
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  axis.line.x = element_line(color = 'black',size = 1.1),
  axis.text = element_text(size = 12)
)
p1 <- p + mytheme+
  scale_x_break(c(0.1,0.7),scales =1,space = 0.1)+
  scale_x_break(c(-0.1,-0.7),scales =1,space = 0.1)
p1

#先根据上下调标签拆分数据框：
up <- dt[which(dt$cor > 0),]
down <- dt[which(dt$cor < 0),]
#添加上调pathway标签：
p2 <- p1 +
  geom_text(data = up,
            aes(x = -0.01, y = gene2, label = gene2),
            size = 5,
            hjust = 1) #标签右对齐
p2
#添加下调pathway标签：
p3 <- p2 +
  geom_text(data = down,
            aes(x = 0.01, y = gene2, label = gene2),
            size = 5,
            hjust = 0) #标签左对齐
p3
#继续调整细节：
p4 <- p3 +
  labs(x = 'Correlation coefficient', y = ' ', title = 'Correlation analysis') + #修改x/y轴标签、标题添加
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14)) #主标题居中、字号调整
p4
p5 <- p4 +
  scale_fill_manual(values = c("#4978bc", "#ea604c"))
p5
p6 <- p5 +
  geom_text(x = -0.705, y = 6.5, label = "Positive correlation", size = 6, color = "#ea604c") +
  geom_text(x = 0.72, y = 1.5, label = "Negative correlation", size = 6, color = "#4978bc")
p6
library(export)
## 导成PPT可编辑的格式
graph2ppt(file="output/发表级别图-初版/ourdata/barplot.pptx")




## 8.随机选取正的和负的分别作图验证
## 正相关的选取IL2RG
library(ggstatsplot)
ggscatterstats(data = exprSet, 
               y = PDCD1, 
               x = IL2RG,
               centrality.para = "mean",                              
               margins = "both",                                         
               xfill = "#CC79A7", 
               yfill = "#009E73", 
               marginal.type = "histogram",
               title = "Relationship between PDCD1 and IL2RG")
## 负相关的选取MARK1
ggscatterstats(data = exprSet, 
               y = PDCD1, 
               x = MARK1,
               centrality.para = "mean",                              
               margins = "both",                                         
               xfill = "#CC79A7",
               yfill = "#009E73", 
               marginal.type = "histogram",
               title = "Relationship between PDCD1 and IL2RG")
## 我们还可以用cowplot拼图
p1 <- ggscatterstats(data = exprSet, 
                     y = PDCD1, 
                     x = IL2RG,
                     centrality.para = "mean",                              
                     margins = "both",                                         
                     xfill = "#CC79A7", 
                     yfill = "#009E73", 
                     marginal.type = "histogram",
                     title = "Relationship between PDCD1 and IL2RG")

p2 <- ggscatterstats(data = exprSet, 
                     y = PDCD1, 
                     x = MARK1,
                     centrality.para = "mean",                              
                     margins = "both",                                         
                     xfill = "#CC79A7", 
                     yfill = "#009E73", 
                     marginal.type = "histogram",
                     title = "Relationship between PDCD1 and IL2RG")
plot_grid(p1,p2,nrow = 1,labels = LETTERS[1:2])

## 5.下面进行聚类分析
### 既然确定了相关性是正确的，那么我们用我们筛选的基因进行富集分析就可以反推这个基因的功能
library(clusterProfiler)
#获得基因列表
library(stringr)
gene <- str_trim(cor_data_sig$gene2,'both')
#基因名称转换，返回的是数据框
gene = bitr(gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
go <- enrichGO(gene = gene, keyType = "SYMBOL",OrgDb = "org.Mm.eg.db", ont="all")
go_s <- simplify(
  go,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang",
  semData = NULL
)
# 条形图
go_f <- as.data.frame(ego_CC)
go_f$richFactor =go_f$Count / as.numeric(sub("/\\d+", "", go_f$BgRatio))
go_f$log.adj.p <- -log10(go_f$p.adjust)
go_n <- go_f %>%
  top_n(n=10,wt=log.adj.p)
ggbarplot(go_n, x="Description", y="log.adj.p", fill = "ONTOLOGY", color = "white", 
          orientation = "horiz",   #横向显示
          palette = "nejm",        #配色方案，常用还有npg，aaas，jama，jco
          legend = "right",        #图例位置
          sort.val = "asc",        #倒序，顺序改为desc
          sort.by.groups=TRUE)+    #按组排序
  scale_y_continuous(expand=c(0, 0)) + scale_x_discrete(expand=c(0,0))+
  ggtitle("Negative correlation") +
  labs(x="Description",y="-Log10(Adjusted p value)",fill=" ")+
  theme(plot.title = element_text(size = 20, face = "bold",hjust=-0.2),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16))
ggsave(file = "Negative_correlation_0.7.pdf",width=15, height=8)


#  气泡图
dotplot(go, split="ONTOLOGY")+ facet_grid(ONTOLOGY~., scale="free")

#但是这个方法有个小缺陷，并不知道最后富集的通路是正向影响还是反向影响，
#也就是无法判断方向。判断方向的工具也不是没有，GSEA就是一个。
#所以，我想能不能把批量相关性分析和GSEA结合起来。
# GSEA
#https://mp.weixin.qq.com/s/sZJPW8OWaLNBiXXrs7UYFw







