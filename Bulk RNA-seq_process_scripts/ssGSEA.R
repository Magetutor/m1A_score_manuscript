setwd("D:/bilibiliR/11、ssGSEA免疫浸润")
rm(list=ls(all=TRUE))

# 加载所需要的R包 ----------------------------------------------------------------
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GSVA")

library(GSVA)#算法
library(tidyverse)#数据处理
library(ggpubr)#绘图
library(ggplot2)#绘图
library(pheatmap)#绘制热图
library(readxl)
library(export)
#ssGSEA是一种算法，根据我们输入的参考基因集进行相应的分析。

# 准备数据:转录组+免疫细胞类型列表+分组数据 -------------------------------------
rownames(exprSet) <- exprSet[,1]
DEG_expr <- exprSet
group <- metadata_merge
group <- group[c(1:4,8:10,5:7,11:13),]
identical(colnames(DEG_expr),group$sample)
group$group <- factor(group$group,levels = c("sham","7d"))
table(markergenes_mice$Cell_type)
#数据处理
#盛夏的果实函数科普：
#lapply 函数是 R 语言的一个内置函数，用于对列表或向量中的每一个元素执行函数。
geneset <- split(markergenes_mice,markergenes_mice$Cell_type)

im_geneset <- lapply(geneset, function(x){
  gene = x$gene
  unique(gene)
})
lapply(im_geneset[1:3], head)

DEG_expr <- as.matrix(DEG_expr) 
# 开始进行ssGSEA分析 ------------------------------------------------------------
result <- gsva(DEG_expr,im_geneset,method = "ssgsea",kcdf = "Gaussian") #“kdcf”参数默认为 "Gaussian"，适用于经过对数转换的芯片数据和测序数据的log-CPMs, log-RPKMs 或者 log-TPMs，"Poisson"适用于测序数据的count格式。

#Min-Max标准化：
#将每个样本中不同的免疫细胞富集分数标准化到0-1间；
result1 <- result
for (i in colnames(result)) {
  result1[,i] <- (result[,i] -min(result[,i]))/(max(result[,i] )-min(result[,i] ))
}


# 结果可视化 -------------------------------------------------------------------
# 相关R包载入：
library(stringr) #数据清洗
library(tidyverse) #数据清洗
library(pheatmap) #热图
library(ggplot2) #箱线图等
library(ggpubr) #箱线图等
library(paletteer) #配色
library(corrplot) #相关性分析
library(rstatix) #添加显著性标记
# 热图
annotation <- data.frame("group" = group$group)
rownames(annotation) <- rownames(group)
p <- pheatmap(result,
         show_colnames = F,
         cluster_cols = F,
         scale = "row", #以行来标准化，这个功能很不错
         annotation_col = annotation,
         color = colorRampPalette(c("#fdebac","white" ,"#582e8c"))(200),#调色
         cellwidth = 25, # 格子宽度
         cellheight = 20,
         fontsize = 12)
graph2pdf(p,file="output/heatmap.pdf")

# 箱式图 -------------------------------------------------------------------
# 比较不同类型免疫细胞的表达情况,剔除正常组，只展示脊髓损伤组,从大到小排列
# 数据处理 ——损伤组/正常组
result2 <- result1[,c(1:12)]
#宽数据转换为长数据(ggplot2绘图格式)：
dt <- result2 %>% t() %>% as.data.frame() %>%
  rownames_to_column("sample") %>%
  gather(key = cell_type,
         value = value, -sample)
head(dt)
dt <- dt %>%
  filter(!(cell_type == "Pyroptosis"))
#将value根据样本转换为百分比形式(新增一列)：
dtt <- dt %>%
  group_by(sample) %>%
  mutate(proportion = round(value/sum(value),3))
head(dtt)

#重新指定箱线图排序（按相对丰度中位数从大到小）：
dtt_arrange <- dtt %>%
  group_by(cell_type) %>%
  summarise(de = median(proportion)) %>%
  arrange(desc(de)) %>%
  pull(cell_type)

dtt$cell_type <- factor(dtt$cell_type,levels = unique(dtt_arrange))

#重新绘图(代码相同)：
#自定义主题：
mytheme <- theme(axis.title = element_text(size = 12),
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 plot.title = element_text(size = 13,
                                           hjust = 0.5,
                                           face = "bold"),
                 legend.text = element_text(size = 10),
                 legend.position = "bottom")

#配色挑选
col1 <-  colorRampPalette(c("darkred","white","darkblue"))
p1 <- ggplot(dtt,
             aes(x = cell_type,y = proportion,fill = cell_type)) +
  geom_boxplot(color = "black",alpha = 0.6,outlier.shape = 21,outlier.size = 1.2) +
  scale_fill_manual(values = col1(14)) +
  labs(x = "cell type", y = "proportion") +
  theme_bw() + mytheme
p1
graph2pdf(p1,file="output/boxplot_propotion_group.pdf")

### 分组箱线图+P值
## 比较不同分组（肿瘤样本和正常样本）内不同免疫细胞的表达情况，同时组间进行两两差异比较。
dt <- result1 %>% t() %>% as.data.frame() %>%
  rownames_to_column("sample") %>%
  gather(key = cell_type,
         value = value, -sample)
head(dt)
#将value根据样本转换为百分比形式(新增一列)：
dtt <- dt %>%
  group_by(sample) %>%
  mutate(proportion = round(value/sum(value),3))
head(dtt)
dtt <- dtt %>%
  filter(!(cell_type == "CD8 T cells")) %>%
  filter(!(cell_type == "T cells")) %>%
  filter(!(cell_type == "Mast cells")) %>%
  filter(!(cell_type == "Neutrophils")) %>%
  filter(!(cell_type == "Granulocytes")) %>%
  filter(!(cell_type == "B derived")) %>%
  filter(!(cell_type == "Eosinophils")) %>%
  filter(!(cell_type == "Neutrophil")) %>%
  filter(!(cell_type == "Microglia"))
#重新指定箱线图排序（按相对丰度中位数从大到小）：
dtt_arrange <- dtt %>%
  group_by(cell_type) %>%
  summarise(de = median(proportion)) %>%
  arrange(desc(de)) %>%
  pull(cell_type)

dtt$cell_type <- factor(dtt$cell_type,levels = unique(dtt_arrange))
# 添加组别
dtt$group <- group$group[match(dtt$sample,group$sample)]
dtt$group <- factor(dtt$group,levels = c("sham","7d"))
# 计算P值
#使用t test或wilcox test进行两两比较(T检验为例)：
t <- t_test(group_by(dtt, cell_type), proportion ~ group)
tj <- adjust_pvalue(t, method = 'fdr') #p值矫正；
tj
#根据p.adj添加显著性标记符号；
tj <- add_significance(tj, 'p.adj')
tj
#在图表中添加 p 值或者显著性标记；
lab <- add_xy_position(tj, x = 'cell_type', dodge = 0.65)
#ggpubr绘图：
p3 <- ggboxplot(dtt, 
                x = "cell_type", 
                y = "proportion",
                fill = "group",
                alpha = 0.8,
                color = "black") +
  scale_fill_manual(values = c("navy","firebrick3")) +
  labs(x = "", y = "proportion") +
  theme_bw() + 
  mytheme + 
  theme(axis.text.x = element_text(angle = 45,size = 11),
        axis.title.y = element_text(size = 11),
        axis.title  = element_text(size = 12)) +
  stat_pvalue_manual(lab, label = 'p.adj.signif', label.size=4, bracket.size=0.5, tip.length = 0.02)
p3
graph2pdf(file="output/boxplot_propotion_group.pdf")
# 相关性分析
#计算相关性系数：
resmcor <- cor(t(result1), method = "pearson")
View(resmcor)
corrplot(resmcor,
         method = "square",
         order = "hclust",
         tl.cex = 0.6,
         tl.col = "black")
## 添加显著性
resmorp <- cor.mtest(resmcor, method = "pearson",conf.level = 0.95) #使用cor.mtest做显著性检验;

#提取p值矩阵；
p.mat <- resmorp$p
View(p.mat)
#相关性热图中展示显著性标记：
col2 <-  colorRampPalette(c("darkblue", "white", "darkred"))
corrplot(resmcor,
         method = "color",
         order = "hclust",
         tl.cex = 1.1,
         tl.col = "black",
         tl.srt = 45,
         col = col2(1000),
         p.mat = resmorp$p, sig.level = c(.001, .01, .05),outline="white",
         insig = "label_sig",pch.cex = 0.8, pch.col = "white")
graph2pdf(file="output/corrplot.pdf")

## 小胶质细胞与焦亡基因集相关性
result_t <- t(result)
result_t <- data.frame(result_t)
ggscatter(result_t, x = "Pyroptosis", y = "Monocyte",
                size = 1.5,
                add = "reg.line",  # 添加回归线
                add.params = list(color = "#0AA1FF", fill = "#a5dff9",size = 1),  # 自定义回归线的颜色
                conf.int = TRUE) +  # 添加置信区间
stat_cor(method = "pearson") +
xlab("Pyroptosis") +
ylab("Monocyte") +
ggtitle("The correlation between Pyroptosis and Monocyte") +
theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text = element_text(size = 12),
        axis.title  = element_text(size = 13))
graph2pdf(file="output/corrplot_2.pdf")

