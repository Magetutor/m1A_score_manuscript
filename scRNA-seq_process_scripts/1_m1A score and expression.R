setwd('/home/bio/Projects/JEM_SCI2021_GSE162610/biyelunwen/')

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
library(tidydr)
library(grid)

# 导入数据
scobj <- readRDS("../completedata/sci.rds")
m1Aset <- read.csv("./m1A_genesets.csv")

# 定义颜色
Cell.col <-c("#9479ad", #Astrocyte
             "#cc4a50", #Neuron
             "#f8766d", # NSC
             "#83ccd2", # Oligodendrocyte
             "#47b9c9", #OPC
             "#e77e2c", # Bcell
             "#9b9c9b", # Granulocyte
             "#7fb687", # Macrophage
             "#1d933f", # Microglia
             "#66c2a5", # Monocyte
             "#caa980", # Neutrophils
             "#00b6eb", # NK
             "#986156", # Endothelial
             "#1271b4", # Microglia
             "#bf783e"# Ependymal
)

##### Figure1A dimplot 改变颜色 ########
Figure1A <- DimPlot(
  scobj,
  reduction = "umap",
  group.by = "celltype",
  label = T, # 是否标注细胞类型
  label.size = 4,
  cols= Cell.col, #自定义颜色
  pt.size = 0.1, # 原点的大小
  repel = T) + #标注有点挤，repel=T可以让排列更加合理)  
  tidydr::theme_dr(xlength = 0.1, # 添加箭头
                   ylength = 0.1,
                   arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(#aspect.ratio = 1,
    panel.grid = element_blank())+
  NoLegend() + 
  labs(title = "Cell Type")

# 保存为pdf
pdf('./figures/Figure_1/Figure1A_umap_all.pdf',width = 6,height = 6)
print(Figure1A)
dev.off()

##### Figure1B m1A调节因子的表达情况-细胞类型 ########

expr <-AverageExpression(scobj, group.by = "celltype",
                         assays = "RNA")[["RNA"]]#orig.ident  lable  cell_type

features <- c("Alkbh1", "Alkbh3", "Fto", "Trmt10c", "Trmt6", "Trmt61a", 
"Ythdc1", "Ythdf1", "Ythdf2", "Ythdf3")
features[!(features %in% rownames(expr))]

features1 <- c('Cpsf6','Cpsf7','Nudt21','Cstf2')

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
ggsave('./figures/Figure_1/Figure 1B_celltype_m1A_exp_heat.pdf',
       p1,
       height=4,
       width=7)

library(scRNAtoolVis)
Idents(scobj) <- 'celltype'
Idents(scobj) <- 'time'
Idents(sub_seurat_obj1) <- 'condition'
DotPlot(object = sub_seurat_obj1,features = features1) + coord_flip()

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
ggsave('./figures/Figure_1/Figure1B_m1A_exp_dot_celltype.pdf',height = 8,width = 10)

VlnPlot(scobj, features = 'Rcn1', 
       split.by = 'time',
       group.by = "celltype",pt.size = 0 ) + NoLegend()

##### Figure1C m1A调节因子的表达情况-组别 ########
Idents(scobj) <- 'time'
DefaultAssay()
jjDotPlot(object = scobj,
          gene = 'Rcn1',
          xtree = F,
          ytree = F,
          rescale = T,
          rescale.min = 0,
          rescale.max = 1
          # point.geom = F,
          # tile.geom = T
) + coord_flip()
ggsave('./figures/Figure_1/Figure1C_m1A_exp_dot_time.pdf',height =5,width = 7)

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
    filename = file.path("./figures/Figure_1/Figure1D_m1a_heatmap.pdf")
  )

#### Figure1E m1A调节因子的表达featureplot ####
p <- FeaturePlot(
  scobj,
  features = 'Rcn1',
  order = T,
  split.by = 'time',
  #raster = T,
  pt.size = 0.5,
  ncol = 5,
  reduction = "umap",cols = c("lightgrey", "red"),
  # cols = c("lightgrey", "firebrick3"),
  #slot = "scale.data"
) 
print(p)
ggsave('./figures/Figure_1/Figure1E_m1A_feature.pdf',height =7,width = 15)


###### Figure1F 计算m1A评分##########
signatures <- list()
signatures$m1a <- m1Aset$gene
features <- c("Alkbh1", "Alkbh3", "Fto", "Trmt10c", "Trmt6", "Trmt61a", 
              "Ythdc1", "Ythdf1", "Ythdf2", "Ythdf3")
DefaultAssay(scobj)

scobj <- AddModuleScore_UCell(scobj, features = signatures, name = "_ucell")

colnames(scobj@meta.data)

DimPlot(scobj, reduction = "umap",split.by = "time")

FeaturePlot(scobj, features = 'm1a_ucell',order = T,split.by = "time")

# 组间以及细胞间画图
plot.data <- scobj@meta.data[,c('time',"celltype",'m1a_ucell')]

comparisons <- list(c("1dpi", "Uninjured"), c("3dpi", "Uninjured"), c("7dpi","Uninjured"))

ggboxplot(plot.data, x = "time", y = 'm1a_ucell', fill = "time", 
          palette = c('#4682b4','#33a02c','#ffa500','#ff4500'),
          error.plot = "errorbar") + 
  stat_compare_means(comparisons = comparisons) + # label = 'p.signif'
  theme_base()+
  ylab('m1A Score')+
  xlab('')

ggsave('./figures/Figure_1/Figure1F_m1A_score_time.pdf',
       ggplot2::last_plot(),
       width = 5,
       height = 4)

## 每个细胞类型中的图

ggplot(plot.data, aes_string(x = 'time', # aes_string这个函数在批量中很有用
                              y = 'm1a_ucell', 
                              fill = 'time')) +
  geom_boxplot() + #width = 0.9默认0.9调剂柱子宽度
  # geom_jitter(size = 0.4, width = 0.5)+
  scale_fill_manual(values =  c('#4682b4','#33a02c','#ffa500','#ff4500')) +
  theme_classic2() + 
  facet_wrap(~celltype,scales='free_x',ncol = 8)+
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

ggsave('./figures/Figure_1/Figure1F_m1A_score_celltype.pdf',
       ggplot2::last_plot(),
       width = 13,
       height = 6)

###### Figure1G 每个基因的表达比例情况##########
features[!(features %in% rownames(scobj@assays$SCT@counts))]
prop <- data.frame()
for (feature in features) {
  if (feature %in% rownames(scobj@assays$SCT@data)) {
    pn <-
      as.data.frame(t(ifelse(
        as.matrix(scobj[feature,]@assays$SCT@data) > 0,
        "postive",
        "negative"
      )))
    scobj <-
      AddMetaData(scobj,
                  metadata = pn,
                  col.name = paste0("pn_", feature))
    df <-
      as.data.frame(prop.table(
        table(scobj$celltype,
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



### 推荐阅读
### https://carmonalab.github.io/UCell_demo/UCell_Seurat_vignette.html
