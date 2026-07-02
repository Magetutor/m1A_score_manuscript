setwd('/home/bio/Projects/NC2022/')

rm(list = ls())
# remove.packages("UCell")
# BiocManager::install("UCell")
library(Seurat)
library(tidyverse)
library(readr)
library(MAST)
library(UCell)
library(ggpubr)
library(ggthemes)
library(ggsci)
library(tidydr)
library(grid)

# 导入数据
scobj_nc <- qs::qread("/home/bio/Projects/NC2022/data/scobj_m1a_240325.qs")
m1Aset <- read.csv("./m1A_genesets.csv")

scobj@meta.data <- scobj@meta.data %>%
  mutate(celltype = case_when(
    coarse_clusters == "Oligodendrocytes" ~ 'Oligodendrocytes',
    coarse_clusters == "Schwann" ~ 'Schwann',
    coarse_clusters == "Oligodendrocyte Progenitors/Precursors" ~ 'OPC',
    coarse_clusters == "Neurons" ~ 'Neurons',
    coarse_clusters == "Microglia/Hematopoietic" ~ 'Microglia',
    coarse_clusters == "Endothelial" ~ 'Endothelial',
    coarse_clusters == "Leptomeninges" ~ 'Leptomeninges',
    coarse_clusters == "Ependymal" ~ 'Ependymal',
    coarse_clusters == "Astrocytes" ~ 'Astrocytes',
    coarse_clusters == "Pericytes" ~ 'Pericytes'
  ))
scobj@meta.data$celltype <- factor(scobj@meta.data$celltype,
                                   levels = c('Oligodendrocytes','Schwann','OPC','Neurons','Microglia',
                                              'Endothelial','Leptomeninges','Ependymal','Astrocytes','Pericytes'))
Idents(scobj) <- 'celltype'

qs::qsave(scobj,file = './data/scobj_m1a_240325.qs')

Cell.col <-  c('#7194b8', '#f9a657', '#94c6c3', '#e77a7b',
               '#78b572', '#f4d76f','#c294b6', '#ffb3ba', '#b2937f','#5ec6e7')

# 先看看有哪些条目
colnames(scobj@meta.data)

# 看看有多少种类细胞
table(scobj@meta.data$celltype) # 10种类细胞

table(scobj@meta.data$condition) 

##### Figure1A dimplot 改变颜色 ########
# scobj@meta.data$condition <- factor(scobj@meta.data$time,
#                                    levels = c('uninjured','4hr','1d','3d','7d','14d','38d'))
Figure1A <- DimPlot(
  scobj,
  reduction = "umap",
  group.by = "celltype",
  label = T, # 是否标注细胞类型
  label.size = 4,
  cols= Cell.col, #自定义颜色
  pt.size = 0.5, # 原点的大小
  repel = T) + #标注有点挤，repel=T可以让排列更加合理)  
  tidydr::theme_dr(xlength = 0.1, # 添加箭头
                   ylength = 0.1,
                   arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(#aspect.ratio = 1,
    panel.grid = element_blank())+
  NoLegend() + 
  labs(title = "Cell Type")

# 保存为pdf
pdf('./figures/Figures_1/Figure1A_umap_all.pdf',width = 6,height = 6)
print(Figure1A)
dev.off()

##### Figure1B m1A调节因子的表达情况-细胞类型 ########

expr <-AverageExpression(scobj, group.by = "celltype",
                         assays = "RNA")[["RNA"]]#orig.ident  lable  cell_type

features <- c("Alkbh1", "Alkbh3", "Fto", "Trmt10c", "Trmt6", "Trmt61a", 
"Ythdc1", "Ythdf1", "Ythdf2", "Ythdf3")
features <-c('mt-Atp6','mt-Nd4','mt-Co1','mt-Cytb')
complex <- c('Ndufb8','Sdhb')

VlnPlot(scobj,
        features = features,
        group.by = "time",
        pt.size = 0,
       # split.by = "time",
        sort = F) +
  #NoLegend() +
  labs(title = "m1A Score") +
  #scale_fill_manual(values = c("#5797bc","#e44349","#84bb7a"))+
  theme(aspect.ratio = 0.3,
        axis.title.x = element_blank())


allgene <- data.frame(rownames(expr))
features[!(features %in% rownames(expr))]


p1 <- pheatmap::pheatmap(
  expr[features,],
  scale = "row",
  border_color = "white",
  #cellwidth = 10,
  #cellheight = 10,
  fontsize = 10,
  cluster_rows = F,
  cluster_cols = F,
  number_color = "white",
  color = colorRampPalette(c("#69b9e8", "white", "#ef6873"))(50)
)
ggsave('./figures/Figures_1/Figure1B_celltype_m1A_exp_heat.pdf',
       p1,
       height=5,
       width=6)

library(scRNAtoolVis)
Idents(scobj) <- 'celltype'
jjDotPlot(object = scobj,
          gene = features,
          xtree = F,
          ytree = F,
          rescale = T,
          rescale.min = 0,
          rescale.max = 1
          # point.geom = F,
          # tile.geom = T
) + coord_flip()
ggsave('./figures/Figures_1/Figure1B_m1A_exp_dot_celltype.pdf',height = 8,width = 10)


##### Figure1C m1A调节因子的表达情况-组别 ########
Idents(scobj) <- 'condition'
DefaultAssay(scobj)
jjDotPlot(object = scobj,
          gene = features,
          xtree = F,
          ytree = F,
          rescale = T,
          rescale.min = 0,
          rescale.max = 1
          # point.geom = F,
          # tile.geom = T
) + coord_flip()

ggsave('./figures/Figures_1/Figure1C_m1A_exp_dot_time.pdf',height =5,width = 7)

#### Figure1D m1A调节因子的表达相关性 ####

m1a_expr <-
  AverageExpression(
    scobj,
    assays = "RNA",
    slot = "data",
    group.by = "celltype",
    features = features
  )[["RNA"]]

#bk <- c(seq(0.8, 0.8999, by = 0.01), seq(0.9, 1, by = 0.01))

cor(m1a_expr) %>%
  pheatmap::pheatmap(
    # cellwidth = 8,
    # cellheight = 8,
    border_color = "white",
    treeheight_row = 0.5,
    treeheight_col = 0.5,
    fontsize = 8,
    # annotation_col = annotation_col,
    color = c(
      colorRampPalette(colors = c("#227e3b", "#f5f4f6"))(20),
      colorRampPalette(colors = c("#f5f4f6", "#7e3b8d"))(20)
    ),
    #legend_breaks = seq(0.8, 1, 0.05),
    #breaks = bk,
    main = "m1A",
    filename = file.path("./figures/Figures_1/Figure1D_m1a_heatmap.pdf")
  )

#### Figure1E m1A调节因子的表达featureplot ####
p <- FeaturePlot(
  scobj,
  features = features,
  order = T,
  #raster = T,
  pt.size = 0.5,
  ncol = 5,
  reduction = "umap",cols = c('lightgrey', 'blue','seagreen2'), #color <- c('lightgrey', 'blue','seagreen2') c("lightgrey", "red")
  # cols = c("lightgrey", "firebrick3"),
  #slot = "scale.data"
) 
print(p)
ggsave('./figures/Figures_1/Figure1E_m1A_feature1.pdf',height =10,width = 25)



###### Figure1F 计算m1A评分##########
signatures <- list()
signatures$m1a <- m1Aset$gene
features <- c("Alkbh1", "Alkbh3", "Fto", "Trmt10c", "Trmt6", "Trmt61a", 
              "Ythdc1", "Ythdf1", "Ythdf2", "Ythdf3")
DefaultAssay(scobj)

scobj <- AddModuleScore_UCell(scobj, features = signatures, name = "_ucell")

colnames(scobj@meta.data)

Idents(scobj) <- 'celltype'

DimPlot(scobj, reduction = "umap",split.by = "condition",group.by = 'celltype',
        cols= Cell.col)+
  guides(color = guide_legend(ncol = 10,override.aes = list(size = 2))) +
  theme(legend.position = "top", legend.direction = "horizontal")

ggsave('/home/bio/Projects/NC2022/figures/Figures_1/Figure1_times_split.pdf',height =5,width = 16)


FeaturePlot(scobj, features = 'm1a_ucell',order = T,split.by = "condition",ncol = 5)+
  theme(legend.position = "right")
ggsave('/home/bio/Projects/NC2022/figures/Figures_1/Figure1_m1Ascore_split.pdf',height =4,width = 18)

# 组间以及细胞间画图
plot.data <- scobj@meta.data[,c('condition',"coarse_clusters",'m1a_ucell')]

comparisons <- list(c("B_1dpi","A_Uninj"),c("C_1wpi", "A_Uninj"), 
                    c("D_3wpi", "A_Uninj"), c("E_6wpi","A_Uninj"))

ggboxplot(plot.data, x = "condition", y = 'm1a_ucell', fill = "condition", 
          palette = c('#4682b4','#33a02c','#ffa500','#ff4500',
                      "#b16268",  "#288c66", 
                      "#264939",  "#6b76ae"),
          error.plot = "errorbar") + 
  stat_compare_means(comparisons = comparisons) + # label = 'p.signif'
  theme_base()+
  ylab('m1A Score')+
  xlab('')

ggsave('./figures/Figures_1/Figure1F_m1A_score_time.pdf',
       ggplot2::last_plot(),
       width = 6,
       height = 5)

## 每个细胞类型中的图

ggplot(plot.data, aes_string(x = 'condition', # aes_string这个函数在批量中很有用
                              y = 'm1a_ucell', 
                              fill = 'condition')) +
  geom_boxplot() + #width = 0.9默认0.9调剂柱子宽度
  # geom_jitter(size = 0.4, width = 0.5)+
  scale_fill_manual(values =  c('#4682b4','#33a02c','#ffa500','#ff4500',
                                "#b16268",  "#288c66", 
                                "#264939",  "#6b76ae")) +
  theme_classic2() + 
  facet_wrap(~coarse_clusters,scales='free_x',ncol = 5)+
  stat_compare_means(comparisons = comparisons,
                     label =  "p.signif",
                     method="wilcox.test")+
  scale_y_continuous(limits = c(0,0.5),breaks = seq(0,0.5,0.1))+
  ylab('m1A Score')+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x =  element_blank(),
        strip.background = element_rect(
          color = "white", fill = "white"),# 分面标题灰色的修改
        strip.text.x = element_text(size = 10), #分面标题字体大小,
        #legend.position = 'top', aspect.ratio = 0.4
        axis.text.y = element_text(size = 10, color = 'black'),
        axis.title.y=element_text(size=10),legend.position = 'bottom')

ggsave('./figures/Figures_1/Figure1F_m1A_score_celltype.pdf',
       ggplot2::last_plot(),
       width = 10,
       height = 6)

###### Figure1G 每个基因的表达比例情况##########
features[!(features %in% rownames(scobj@assays$RNA@counts))]
prop <- data.frame()
for (feature in features) {
  if (feature %in% rownames(scobj@assays$RNA@data)) {
    pn <-
      as.data.frame(t(ifelse(
        as.matrix(scobj[feature,]@assays$RNA@data) > 0,
        "postive",
        "negative"
      )))
    scobj <-
      AddMetaData(scobj,
                  metadata = pn,
                  col.name = paste0("pn_", feature))
    df <-
      as.data.frame(prop.table(
        table(scobj$coarse_clusters,
              scobj@meta.data[, paste0("pn_", feature)]),
        margin = 1
      ))
    df$feature <- feature
    prop <- rbind(prop, df)
  }
}
head(prop)
colnames(prop) <- c("celltype", "group", "freq", "feature")


ggplot(prop[prop$group == "postive",],
       aes(x = freq * 100, y = feature, fill = celltype)) +
  geom_bar(stat = "identity") +
  facet_grid(. ~ celltype) +
  theme_bw() +
  scale_fill_manual(values = Cell.col)+
  theme(
    axis.text.y = element_text(
      face = "italic",
      colour = "black",
      size = 10,
      family = "Arial"
    ),
    axis.title.y = element_blank(),
    axis.text.x = element_text(
      colour = "black",
      size = 10,
      family = "Arial"
    ),
    strip.background = element_blank(),
    strip.text = element_text(size = 10, family = "Arial"),
    axis.title.x = element_text(size = 10, family = "Arial"),
    panel.grid = element_blank()
  ) +
  labs(x = "Cell proportion (%)") +
  guides(fill = "none")

ggsave('./figures/Figure_1/Figure1G_m1A_postive_gene_bar.pdf',
       ggplot2::last_plot(),#最后一张图
       height=6,
       width=18)

install.packages("showtext")
library(showtext)
font_add('Arial', '~/Public/arial.ttf')#使用自己的路径
showtext_auto()


VlnPlot(scobj,
        features = "m1a_ucell",
        # group.by = "group",
        pt.size = 0,
        # split.by = "Disease",
        sort = F) +
  #NoLegend() +
  labs(title = "m1A Score") +
  #scale_fill_manual(values = c("#5797bc","#e44349","#84bb7a"))+
  theme(aspect.ratio = 0.3,
        axis.title.x = element_blank())

####y提取一下神经元看看有多少##########
neuro <- subset(scobj, coarse_clusters == 'Neurons')
table(neuro@meta.data$condition)

# 将m1A打分的数据进行保存
qs::qsave(scobj,file = './data/GSE172167_SCI_all_nuclei_m1ascore.qs')

### 推荐阅读
### https://carmonalab.github.io/UCell_demo/UCell_Seurat_vignette.html
