################### After SCENIC #######################
## Aims:
## 从SCENIC outputs到生物学假设
## Input data: PBMC (ctrl vs IFNB simulated)
## => 1. 受到IFNB影响最大的PBMC细胞类型是什么?        | cellular level
##    2. 哪些TF驱动了IFNB simulated PBMC的转录组变化? | Molecular level
##    3. 哪些TF驱动了哪些细胞类型的什么样的变化?
##       (下游的基因, related to某些生物学功能)       | Functional level

library(Seurat)
library(tidyverse)
library(patchwork)
library(AUCell)

# 定义颜色
group <- c('#4682b4','#33a02c','#ffa500','#ff4500')
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
setwd(here::here())
source("./R_pyscenic/compute_module_score.R")

## 导入Seurat对象
myeloid <- readRDS("/home/bio/Projects/JEM_SCI2021_GSE162610/completedata/myeloid.rds")
seu <- myeloid

# seu <- subset(seu, celltype %in% c("Mono/Mk Doublets", "Eryth"), invert = TRUE)

DimPlot(seu, group.by = "celltype", reduction = "umap", label = T)

celltype.levels <- c("CD14 Mono", "CD16 Mono", "DC", "pDC", "B cell", "B Activated",
                     "NK", "CD8 T", "CD4 Memory T", "T activated", "CD4 Naive T", "Mk")

seu$celltype <- factor(seu$celltype, levels = celltype.levels)


## 导入regulon (gene list)
regulons <- clusterProfiler::read.gmt("output/myeloid_pyscenic/02-meyloid.regulons.gmt")
## data.frame -> list, list中的每个元素为一个gene set
rg.names <- unique(regulons$term)
regulon.list <- lapply(rg.names, function(rg) {
  subset(regulons, term == rg)$gene
})
names(regulon.list) <- sub("[0-9]+g", "\\+", rg.names)
summary(sapply(regulon.list, length))
print(regulon.list[1])
saveRDS(regulon.list, "output/myeloid_pyscenic/03-1.meyloid.regulons.rds")

## 用AUCell计算RAS matrix
## RAS = regulon activity score
seu <- ComputeModuleScore(seu, gene.sets = regulon.list, min.size = 10, cores = 1)
seu
DefaultAssay(seu) <- "AUCell"

p1 <- FeaturePlot(seu, features = "Etv6(+)", split.by = "time")
DefaultAssay(seu) <- "RNA"

p2 <- FeaturePlot(seu, features = "Etv6", slot = 'data',split.by = "time")
(p1 / p2) & scale_color_viridis_c()

VlnPlot(seu, group.by = "celltype", features = "Etv6(+)", pt.size = 0,
        split.by = "time", split.plot = TRUE, cols = group) + ylab("TF activity") +
  geom_hline(aes(yintercept = 0.07), colour = 'black', linetype = "dashed")

ggsave('./figures/Figure_2/Figure2A_Etv6_TFactivity.pdf',width = 8,height =4 )

VlnPlot(seu, group.by = "celltype", features = "Etv6", pt.size = 0,
        split.by = "time", split.plot = TRUE, cols = group)
ggsave('./figures/Figure_2/Figure2B_Etv6_TFexpress.pdf',width = 8,height =4 )


## 用RAS matrix计算UMAP
seu <- RunUMAP(object = seu,
               features = rownames(seu),
               metric = "correlation", # 注意这里用correlation效果最好
               reduction.name = "umapRAS",
               reduction.key = "umapRAS_")

 ## 可视化：UMAP on harmony
p1 <- DimPlot(seu, reduction = "umap", group.by = "celltype") + ggsci::scale_color_d3("category20") + NoLegend()
p2 <- DimPlot(seu, reduction = "umap", group.by = "time") + NoLegend()

## 可视化：UMAP on RAS
p3 <- DimPlot(seu, reduction = "umapRAS", group.by = "celltype") + ggsci::scale_color_d3("category20")
p4 <- DimPlot(seu, reduction = "umapRAS", group.by = "time")

pdf('./figures/Figure_2/Figure2C_4umap.pdf',width = 10,height =8 )
(p1 + p3) / (p2 + p4)
dev.off()

## 推测：INFB对髓系细胞的影响更大

## 换一种方式：PCA'
DefaultAssay(seu) <- "AUCell"
seu <- ScaleData(seu)
seu <- RunPCA(object = seu,
              features = rownames(seu),
              reduction.name = "pcaRAS",
              reduction.key = "pcaRAS_")

## 可视化：PCA on RAS
p3 <- DimPlot(seu, reduction = "pcaRAS", group.by = "celltype") + ggsci::scale_color_d3("category20")
p4 <- DimPlot(seu, reduction = "pcaRAS", group.by = "time")
p3 + p4

pdf('./figures/Figure_2/Figure2D_4pca.pdf',width = 10,height =4)
(p3 + p4)
dev.off()
## PC1 encoding the regulons related to cell type
## PC2 encoding the regulons affected by INFB treatment
## The INFB induced transcriptome shift is orthogonal to the cell identity transcriptional programs.

VlnPlot(seu, group.by = "celltype", features = "pcaRAS_1", pt.size = 0,
        split.by = "time", split.plot = TRUE, cols = c("blue", "red"))

VlnPlot(seu, group.by = "celltype", features = "pcaRAS_2", pt.size = 0,
        split.by = "time", split.plot = TRUE, cols = c("blue", "red"))

qs::qsave(seu, "output/03-2.myeloid.seurat.aucell.qs")

###### Figure2 计算meyloid m1A评分##########
signatures <- list()
signatures$m1a <- features

DefaultAssay(seu)

seu <- AddModuleScore_UCell(seu, features = signatures, name = "_ucell")

Idents(seu)

colnames(seu@meta.data)
Idents(seu) <- 'celltype'

DimPlot(seu, reduction = "umap",split.by = "time")+ ggsci::scale_color_d3("category20")
ggsave('./figures/Figure_2/Figure2E_umap_times.pdf',width = 12,height =4)


FeaturePlot(seu, features = 'm1a_ucell',order = T,split.by = "time")
ggsave('./figures/Figure_2/Figure2E_umap_m1aucell.pdf',width = 12,height =3)

# 组间以及细胞间画图
plot.data <- seu@meta.data[,c('time',"celltype",'m1a_ucell')]

comparisons <- list(c("1dpi", "Uninjured"), c("3dpi", "Uninjured"), c("7dpi","Uninjured"))

ggboxplot(plot.data, x = "time", y = 'm1a_ucell', fill = "time", 
          palette = c('#4682b4','#33a02c','#ffa500','#ff4500'),
          error.plot = "errorbar") + 
  stat_compare_means(comparisons = comparisons) + # label = 'p.signif'
  theme_base()+
  ylab('m1A Score')+
  xlab('')

ggsave('./figures/Figure_2/Figure2F_m1A_score_time.pdf',
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
  scale_y_continuous(limits = c(0,0.45),breaks = seq(0,0.5,0.1))+
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

ggsave('./figures/Figure_2/Figure2F_m1A_score_celltype.pdf',
       ggplot2::last_plot(),
       width = 13,
       height = 6)

qs::qsave(seu, "output/03-2.myeloid.seurat.aucell.qs")

######计算m1A评分和转录因子的相关性########
