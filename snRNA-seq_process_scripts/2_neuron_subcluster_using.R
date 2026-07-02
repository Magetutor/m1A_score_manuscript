
# 设置工作目录
setwd('/home/bio/Projects/NC2022/neuro_sub/')
# 清空变量 加载包
# rm(list = ls())
library(Seurat)
# packageVersion('sub_seurat_obj1')
library(DESeq2)
library(tidyverse)
library(ggthemes)
# 加载数据
sub_seurat_obj1 <- qs::qread("../data/scobj_m1a_240325.qs")
m1Aset <- read.csv("../m1A_genesets.csv")


sub_seurat_obj1 <- subset(sub_seurat_obj1,celltype %in% c("Neurons"))

# sub_seurat_obj1_obj$Group[sub_seurat_obj1_obj$Group %in% c("MOCK")] <- "PBS"
# sub_seurat_obj1_obj$Group[sub_seurat_obj1_obj$Group %in% c("V_FS")] <- "FSS13025"
# sub_seurat_obj1_obj$Group[sub_seurat_obj1_obj$Group %in% c("V_59")] <- "PRVABC59"

######直接分群

DefaultAssay(sub_seurat_obj1) <- 'integrated'

sub_seurat_obj1 <- ScaleData(sub_seurat_obj1, verbose = FALSE)

sub_seurat_obj1 <- RunPCA(sub_seurat_obj1, npcs = 50, verbose = FALSE)
ElbowPlot(sub_seurat_obj1,ndims = 30)

sub_seurat_obj1 <- RunUMAP(sub_seurat_obj1, reduction = "pca", dims = 1:10)

#sub_seurat_obj1 <- RunTSNE(sub_seurat_obj1, reduction = "pca", dims = 1:10)

sub_seurat_obj1 <- FindNeighbors(sub_seurat_obj1, reduction = "pca", dims = 1:10)

sub_seurat_obj1 <- FindClusters(sub_seurat_obj1, resolution = 0.05)

DimPlot(sub_seurat_obj1, reduction = "umap",pt.size = 0.5,
        # split.by = 'Group',
        label = T,label.box = T,label.size = 4,repel = T) + 
  theme_bw() +
  labs( x= "UMAP_1",y= "UMAP_2",title = "") +
  theme(panel.grid = element_blank()
        #,aspect.ratio = 1
  )+NoLegend()

# 看一下每个群的marker
DefaultAssay(sub_seurat_obj1) <- 'integrated'

markers <- FindAllMarkers(
  sub_seurat_obj1,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

top5.markers <- markers %>%  group_by(cluster) %>% slice_max(n = 5, order_by = avg_log2FC)

top10.markers <- markers %>%  group_by(cluster) %>% slice_max(n = 10, order_by = avg_log2FC)

top20.markers <- markers %>%  group_by(cluster) %>% slice_max(n = 20, order_by = avg_log2FC)

top50.markers <- markers %>%  group_by(cluster) %>% slice_max(n = 50, order_by = avg_log2FC)

top100.markers <- markers %>%group_by(cluster) %>% slice_max(n = 100, order_by = avg_log2FC)

source("../selfcode/R/data.R")
source("../selfcode/R/objects.R")
source("../selfcode/R/rip-seq.R")
source("../selfcode/R/rna-seq.R")
source("../selfcode/R/scrna-seq.R")
source("../selfcode/R/utilities.R")
source("../selfcode/R/visualization.R")

DefaultAssay(sub_seurat_obj1) <- 'RNA'


top5.markers$gene
DotPlot(sub_seurat_obj1, features = unique(top5.markers$gene), col.min = 0) + 
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))

# 确定每个亚群是什么
features <- c('Snhg11','Rbfox1', 'Slc17a6',# 
              'Gad1','Gad2','Sox5','Lyz2', # 
              'Stmn1','Pclaf','Hmgb2','Ube2c','Tubb5', # 
              'Fabp5','Spp1','Gpnmb','Ccl9'# 
)

features <- c('Slit2','Esrrg', 'Cdh18', 'Cntn5',# Slit2_IN 0
              'Npy','Gpc5','Chrm3','Nxph1', # Npy_IN 1
              'Cck','Nts','Adarb2','Erbb4', # Cck_EN 2
              'Sox5','Tac1','Fam19a1','Lmo3',# Sox5_EN 3
              'Pde11a','Cpne8','Prex2','Fbn2', # Pde11a_EN 4
              'Tac2','Nmu','Sst','Maml3', # Tac2_EN 5 
              'Gal','Pnoc','Kcnt2','Cartpt' # Gal_IN 6
              #'Pkd1l2','Rnf220', 'Meis2','Myo3b'#  7
              )

DotPlot(sub_seurat_obj1, features = unique(features), col.min = 0) + 
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))

Idents(sub_seurat_obj1) <- 'seurat_clusters'
DotPlot(sub_seurat_obj1, dot.scale = 5, col.min = 0, # col.max = 1., # scale = F,
        features = unique(features)) +
  scale_color_distiller(palette = "RdBu", direction = -1) + # YlOrRd
  coord_fixed() +
  # theme_cat() +
  theme(
    axis.title = element_blank(),
    legend.margin = margin(b = -8),
    legend.position = "top",
    axis.text.x = element_text(face = "italic")
  ) +
  guides(
    x = guide_axis(angle = 90),
    color = guide_colorbar(
      title = "Expression",
      frame.colour = "black",
      ticks.colour = "black"
    )
  ) +
  scale_y_discrete(limits = rev)

ggsave('./figure/Figure_2/A_marker_dot.pdf',width = 10,height =4)

# 命名 cell_type
new.cluster.ids <- c(
  "c0",# 00
  "c1",#01
  "c2",#02
  "c3",#03
  "c4",#04
  "c5",#05
  "c6",#06
  "c7")

names(new.cluster.ids) <- levels(sub_seurat_obj1)
sub_seurat_obj1 <- RenameIdents(sub_seurat_obj1,
                       new.cluster.ids)
sub_seurat_obj1$cell_type <- Idents(sub_seurat_obj1)

cell_type <- unique(sub_seurat_obj1@meta.data$cell_type)

# 使用top100marker功能富集来看看
# 给markers添加c1标签
markers1 <- top100.markers
markers1$cluster <- as.character(markers1$cluster)
markers1 <- markers1 %>% mutate(cell_type = case_when(
  cluster == "0" ~ "c0",
  cluster == "1" ~ "c1",
  cluster == "2" ~ "c2",
  cluster == "3" ~ "c3",
  cluster == "4" ~ "c4",
  cluster == "5" ~ "c5",
  cluster == "6" ~ "c6",
  cluster == "7" ~ "c7"
))

library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)
for(i in cell_type){
  path <- paste0('/home/bio/Projects/NC2022/neuro_sub/figure/Figure_2/', i, '/')
  dir.create(path,recursive = T)
  setwd(path)
  outdir <- "./files/"
  dir.create(outdir,recursive = T)
  outdir2 <- "./plots/"
  dir.create(outdir2,recursive = T)
  
  dataf <- markers1[markers1$cell_type == i,]
  
  
  #### GO功能富集 ----
  
  #### 1.GO-UP ----
  bp <-
    enrichGO(
      dataf$gene,
      OrgDb = org.Mm.eg.db,
      keyType = 'SYMBOL',
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2)
  
  term <- bp@result
  #保存结果
  write.table(term,
              file = file.path(outdir , paste0(i,'_GOBP.csv')),
              quote = F,
              sep = ",",
              row.names = F
  )
  #绘制TOP20
  df <- term[1:30,]
  df$labelx=rep(0,nrow(df))
  df$labely=seq(nrow(df),1)
  ggplot(data = df,
         aes(x = -log10(pvalue),
             y = reorder(Description,-log10(pvalue))))  +
    geom_bar(stat="identity",
             alpha=1,
             fill= "#ee6470",
             width = 0.8) +
    geom_text(aes(x=labelx,
                  y=labely,
                  label = df$Description),
              size=3.5,
              hjust =0)+
    theme_classic()+
    theme(axis.text.y = element_blank(),
          axis.line.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.x = element_line(colour = 'black', linewidth = 1),
          axis.text.x = element_text(colour = 'black', size = 10),
          axis.ticks.x = element_line(colour = 'black', linewidth = 1),
          axis.title.x = element_text(colour = 'black', size = 12))+
    xlab("-log10(pvalue)")+
    ggtitle(i)+
    scale_x_continuous(expand = c(0,0))
  #保存
  ggsave(file.path(outdir2, paste0( i, '_GO_Barplot.pdf')),
         ggplot2::last_plot(),#最后一张图
         height=6,
         width=6)
  
  
  #### KEGG功能富集 ----
  
  top.genes <-
    dataf$gene
  
  
  convert <- bitr(
    top.genes,
    fromType = "SYMBOL",
    toType = c("ENTREZID"),
    OrgDb = org.Mm.eg.db
  )
  top.genes <- convert$ENTREZID
  
  ego <- enrichKEGG(
    gene = top.genes,
    keyType = "kegg",
    organism  = 'mmu',
    # human: hsa, mouse: mmu
    pvalueCutoff  = 0.05,
    pAdjustMethod  = "BH",
    qvalueCutoff  = 0.2,
  )
  
  term <- ego@result
  term$celltype <- i
  
  write.table(term,
              file = file.path(outdir , paste0(i,'_KEGG.csv')),
              quote = F,
              sep = ",",
              row.names = F
  )
  #绘制TOP20
  df <- term[1:30,]
  df$labelx=rep(0,nrow(df))
  df$labely=seq(nrow(df),1)
  ggplot(data = df,
         aes(x = -log10(pvalue),
             y = reorder(Description,-log10(pvalue))))  +
    geom_bar(stat="identity",
             alpha=1,
             fill= "#ee6470",
             width = 0.8) +
    geom_text(aes(x=labelx,
                  y=labely,
                  label = df$Description),
              size=3.5,
              hjust =0)+
    theme_classic()+
    theme(axis.text.y = element_blank(),
          axis.line.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.x = element_line(colour = 'black', linewidth = 1),
          axis.text.x = element_text(colour = 'black', size = 10),
          axis.ticks.x = element_line(colour = 'black', linewidth = 1),
          axis.title.x = element_text(colour = 'black', size = 12))+
    xlab("-log10(pvalue)")+
    ggtitle(i)+
    scale_x_continuous(expand = c(0,0))
  #保存
  ggsave(file.path(outdir2, paste0( i, '_UP_Barplot_KEGG.pdf')),
         ggplot2::last_plot(),#最后一张图
         height=6,
         width=6)
  
}
setwd('/home/bio/Projects/NC2022/neuro_sub/')
save(markers,file = './data/each_cluster_markers.Rda')

## 命名为生物学意义亚群

Idents(sub_seurat_obj1) <- 'seurat_clusters'

features <- c('Slit2','Esrrg', 'Cdh18', 'Cntn5',# Slit2_IN 0
              'Npy','Gpc5','Chrm3','Nxph1', # Npy_IN 1
              'Cck','Nts','Adarb2','Erbb4', # Cck_EN 2
              'Sox5','Tac1','Fam19a1','Lmo3',# Sox5_EN 3
              'Pde11a','Cpne8','Prex2','Fbn2', # Pde11a_EN 4
              'Tac2','Nmu','Sst','Maml3', # Tac2_EN 5 
              'Gal','Pnoc','Kcnt2','Cartpt', # Gal_IN 6
              'Pkd1l2','Rnf220', 'Meis2','Myo3b'#  7
)

new.cluster.ids <- c(
  "Slit2_IN",# 00
  "Npy_IN",#01
  "Cck_EN",#02
  "Sox5_EN",#03
  'Pde11a_EN',#04
  "Tac2_EN",#05
  "Gal_IN", #06
  "c7"
  )

names(new.cluster.ids) <- levels(sub_seurat_obj1)
sub_seurat_obj1 <- RenameIdents(sub_seurat_obj1,
                       new.cluster.ids)
sub_seurat_obj1$celltype <- Idents(sub_seurat_obj1)

DimPlot(sub_seurat_obj1 , reduction = "umap",raster=FALSE,label = T,pt.size = 0.4)

DimPlot(sub_seurat_obj1, reduction = "umap",raster=FALSE,
              label = T,label.box = T,label.size = 2,repel = T) + theme_bw() +
  labs( x= "UMAP_1",y= "UMAP_2",title = "") +
  theme(panel.grid = element_blank(),aspect.ratio = 1)+NoLegend()




DimPlot(sub_seurat_obj1, reduction = "umap",raster=FALSE,split.by = "condition",
        label = T,repel = T) + theme_bw() +
  labs( x= "UMAP_1",y= "UMAP_2",title = "") +
  theme(panel.grid = element_blank(),aspect.ratio = 1)

# 因为c7太少，直接去掉c7
Idents(sub_seurat_obj1) <- 'cell_type'

sub_seurat_obj1 <- subset(x = sub_seurat_obj1, idents = c('c7'), invert = TRUE)

sub_seurat_obj1@meta.data$celltype <- factor(sub_seurat_obj1@meta.data$celltype,
                                             levels = c("Slit2_IN", 
                                                        "Npy_IN", "Cck_EN", "Sox5_EN", 
                                                        "Pde11a_EN", "Tac2_EN", "Gal_IN"))
#save(sub_seurat_obj1,"./data/neuron_subcluster_20250622.Rda")
qs::qsave(sub_seurat_obj1,"./data/neuron_subcluster_20250622.qs")
#### 细胞比例计算 ----
# 查看整体每个细胞类型的数量
table(Idents(sub_seurat_obj1))
# 计算整体每个细胞类型所占百分比
prop.table(table(Idents(sub_seurat_obj1))) * 100

# 查看各个分组每个细胞类型数量
table(Idents(sub_seurat_obj1), sub_seurat_obj1$condition)
table(sub_seurat_obj1$condition)
# 计算各个分组每个细胞类型所占百分比
prop.table(table(Idents(sub_seurat_obj1), sub_seurat_obj1$condition), margin = 2) * 100
#### 可视化 
meta <- sub_seurat_obj1@meta.data

groupcol <- c("#E64B35CC", "#4DBBD5CC", "#00A087CC", "#3C5488FF",
               "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "#DC0000FF" ,"#7e6148")

p <- ggplot(data = meta, aes(x = celltype, fill = condition)) +
  geom_bar(position = "fill",
           width = 0.7,
           size = 0.3,
           color = "black") +
  labs(title = "Group", x = "", y = "Fraction of Cell") +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 8),
    legend.title = element_blank(),
    legend.background = element_blank(),
    axis.title.x = element_text(size = 8),
    axis.line.x = element_line(size = 0.5),
    axis.ticks.x = element_line(size = 0.5, colour = "black"),
    axis.ticks.length = unit(.5, "lines"),
    axis.ticks.y = element_blank(),
    axis.text = element_text(size = 8, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 8),
    strip.background = element_rect(color = "white", fill = "white"),
    strip.text.y = element_text(size = 8)
  ) +
  # coord_flip() +
  scale_y_continuous(expand = c(0, 0)) +
  # scale_fill_manual(values = c("#ff6347", "#4f94cd","red","blue","green")) +
  scale_fill_manual(values = groupcol) +
  # geom_hline(yintercept = 0.25,
  #            linetype = "dashed",
  #            size = 0.5) +
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        axis.text = element_text(colour = 'black',size = 10),
        legend.title = element_blank()
        #axis.text.y = element_text(face = 'italic')
        )


ggsave(
  file.path( "./figure/Figure_2/B_percentage-condition.pdf"),
  plot = p,
  height = 4.5,
  width = 6
)



# 填充细胞类型，
# cellcol <- c('#80120b', '#252662', '#9a9a9a', '#fbd9b4', 
#          '#e70c4c', '#5ec6e7','#f18433', '#fce218', '#b3d8b6', '#78398f')

cellcol <- c('#7194b8', '#f9a657','#e77a7b', '#94c6c3', 
         '#78b572', '#f4d76f','#c294b6', '#ffb3ba', '#b2937f','#5ec6e7')

# cellcol = c("#9ACBDE","#F47892","#64AE59","#6DCCDD",
#                  "#ad6593","#F6A395","#C4A5DE","#B9C984",
#                  "#F8AD77","#FF6347","#EDC6DD")

p2 <- ggplot(data = meta, 
             aes(x = condition,
                 fill = celltype)) +
  geom_bar(position = "fill",
           width = 0.7,
           size = 0.3,
           color = "black") +
  labs(title = "", x = "", y = "Fraction of Cell") +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 8),
    legend.title = element_blank(),
    legend.background = element_blank(),
    axis.title.x = element_text(size = 8),
    axis.line.x = element_line(size = 0.5),
    axis.ticks.x = element_line(size = 0.5, colour = "black"),
    axis.ticks.length = unit(.5, "lines"),
    axis.ticks.y = element_blank(),
    axis.text = element_text(size = 8, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 8),
    strip.background = element_rect(color = "white", fill = "white"),
    strip.text.y = element_text(size = 8)
  ) +
  # coord_flip() +
  scale_y_continuous(expand = c(0, 0)) +
  # scale_fill_manual(values = c("#ff6347", "#4f94cd","red","blue","green")) +
  scale_fill_manual(values = cellcol) +
  # geom_hline(yintercept = 0.25,
  #            linetype = "dashed",
  #            size = 0.5) +
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))


p2


ggsave(
  file.path( "./figure/Figure_2/C_percentage-celltype.pdf"),
  plot = p2,
  height = 4.5,
  width = 6)

#### A图 画umap图 ###########
Idents(sub_seurat_obj1) <- 'celltype'

DimPlot(sub_seurat_obj1, label = T, pt.size = 1)+
  NoLegend()+labs(x = "UMAP1", y = "UMAP2",title = "") +
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) + scale_color_manual(values = cellcol)

ggsave( "./figure/Figure_2/A1_umap.pdf", height = 6,width = 6)

#### D图 细胞比例计算 ###########
source('../selfcode/Singlecellratio_plotstat.R')
my_comparisons <- list( c("B_1dpi", "A_Uninj"),
                        c("C_1wpi", "A_Uninj"),
                        c("D_3wpi", "A_Uninj"),
                        c("E_6wpi", "A_Uninj"))


Singlecellratio_plotstat(sub_seurat_obj1, group_by = "condition",
                         meta.include = c("condition","orig.ident"),
                         comparisons = my_comparisons, color_by = 'condition',
                         group_by.point = "orig.ident",label.x = 1, pt.size = 1,
                         label = 'p.signif', ncol =3) + scale_color_manual(values = groupcol)+
  theme(legend.position = "top", legend.direction = "horizontal",
        legend.text = element_text(size = 10),
        legend.title=element_blank() # 去掉legend title
  )

ggsave( "./figure/Figure_2/D_cellpercentage_2.pdf", height = 8,width = 8)


my_comparisons1 <- list(c("Npy_IN", "Slit2_IN"),
                        c("Cck_EN", "Slit2_IN"),
                        c("Sox5_EN", "Slit2_IN"),
                        c("Pde11a_EN", "Slit2_IN"),
                        c("Tac2_EN", "Slit2_IN"),
                        c("Gal_IN", "Slit2_IN"))

Idents(sub_seurat_obj1) <- 'condition'
Singlecellratio_plotstat1(sub_seurat_obj1, group_by = "celltype",
                         meta.include = c("celltype","orig.ident"),
                         comparisons = my_comparisons1, color_by = 'celltype',
                         group_by.point = "orig.ident",label.x = 1, pt.size = 1,
                         label = 'p.signif', ncol =4) + scale_color_manual(values = cellcol)+
  theme(legend.position = "top", legend.direction = "horizontal",
        legend.text = element_text(size = 10),
        legend.title=element_blank() # 去掉legend title
  )

ggsave( "./figure/Figure_2/D_cellpercentage_3.pdf", height = 8,width = 8)


#### E图 神经元亚群m1A打分比较 ###########
features <- c("Alkbh1", "Alkbh3", "Fto", "Trmt10c", "Trmt6", "Trmt61a", 
              "Ythdc1", "Ythdf1", "Ythdf2", "Ythdf3")
signatures <- list()
signatures$m1a <- features


DefaultAssay(sub_seurat_obj1) <- 'RNA'

sub_seurat_obj1 <- AddModuleScore_UCell(sub_seurat_obj1, features = signatures, name = "_ucell")

colnames(sub_seurat_obj1@meta.data)

Idents(sub_seurat_obj1) <- 'celltype'

DimPlot(sub_seurat_obj1, reduction = "umap",split.by = "condition",cols= cellcol)+
  guides(color = guide_legend(ncol = 10,override.aes = list(size = 2))) +
  theme(legend.position = "top", legend.direction = "horizontal")

ggsave('./figure/Figure_2/A1_umap_split.pdf',height =5,width = 16)


FeaturePlot(sub_seurat_obj1, features = 'm1a_ucell', cols = c("gray", "red"),
            order = T,min.cutoff = 0.01,
            split.by = "condition",ncol = 5) + theme(legend.position = "right")

ggsave('./figure/Figure_2/E1_m1ascore_split.pdf',height =4,width = 20)

# 组间以及细胞间画图
plot.data <- sub_seurat_obj1@meta.data[,c('condition',"celltype",'m1a_ucell')]

comparisons <- list(c("B_1dpi","A_Uninj"),c("C_1wpi", "A_Uninj"), 
                    c("D_3wpi", "A_Uninj"), c("E_6wpi","A_Uninj"))

ggboxplot(plot.data, x = "condition", y = 'm1a_ucell', fill = "condition", 
          palette = groupcol,
          error.plot = "errorbar") + 
  stat_compare_means(comparisons = comparisons) + # label = 'p.signif'
  theme_base()+
  ylab('m1A Score')+
  xlab('')

ggsave('./figure/Figure_2/E2_m1ascore_group.pdf',
       ggplot2::last_plot(),
       width = 6,
       height = 5)

## 每个细胞类型中的图

ggplot(plot.data, aes_string(x = 'condition', # aes_string这个函数在批量中很有用
                             y = 'm1a_ucell', 
                             fill = 'condition')) +
  geom_boxplot() + #width = 0.9默认0.9调剂柱子宽度
  # geom_jitter(size = 0.4, width = 0.5)+
  scale_fill_manual(values =  groupcol) +
  theme_classic2() + 
  facet_wrap(~celltype,scales='free_x',ncol = 4)+
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

ggsave('./figure/Figure_2/F_m1A_score_celltype.pdf',
       ggplot2::last_plot(),
       width = 10,
       height = 6)

save(sub_seurat_obj1,"./data/neuron_subcluster_m1ascore.Rda")
save(sub_seurat_obj,"./data/neuron_SCTsubcluster_uncomplete.Rda")
qs::qsave(sub_seurat_obj1,"./data/neuron_subcluster_m1ascore.qs")

### 批量找各个组与对照组相比的差异基因
# 批量找各个细胞在V59中的差异基因
dput(unique(sub_seurat_obj1@meta.data$celltype))
dput(unique(sub_seurat_obj1@meta.data$condition))

Cells <- c("Slit2_IN",  "Npy_IN", "Cck_EN", "Sox5_EN", "Pde11a_EN", "Tac2_EN", "Gal_IN")
group <- c("B_1dpi", "C_1wpi", "D_3wpi", "E_6wpi")


for (j in group) {
  path <- paste0('/home/bio/Projects/NC2022/neuro_sub/figure/Figure_3/', j, '/')
  dir.create(path,recursive = T)
  setwd(path)
  outdir <- "./files/"
  dir.create(outdir,recursive = T)
  outdir2 <- "./plots/"
  dir.create(outdir2,recursive = T)
  
  diff <- list()
  for (i in 1:length(Cells)) {
    cells <- subset(sub_seurat_obj1,
                    celltype == Cells[i])
    Idents(cells) <- cells@meta.data$condition
    cells.diff <- FindMarkers(cells,
                              ident.1 = j,
                              ident.2 = "A_Uninj",
                              test.use = "MAST",
                              min.pct = 0.25, # 过滤掉那些在50%以下细胞中检测到的基因
                              logfc.threshold = 0.10 # 过滤掉那些在不同组之间平均表达的差异倍数低于2的基因
    )
    diff[[i]] <- cells.diff
  }
  names(diff) <- Cells
  
# 分成上下调基因
  # i = 'Slit2_IN'
  updiff <- list()
  for (i in 1:length(Cells)) {
    up <- diff[[i]] %>% 
      filter(avg_log2FC > 0)
    updiff[[i]] <- up
  }
  names(updiff) <- Cells
  
  dndiff <- list()
  for (i in 1:length(Cells)) {
    dn <- diff[[i]] %>% 
      filter(avg_log2FC < 0)
    dndiff[[i]] <- dn
  }
  names(dndiff) <- Cells
  
  save(diff, file = './files/diff.Rda')
  save(updiff, file = './files/updiff.Rda')
  save(dndiff, file = './files/updiff.Rda')
  
  
  # 每个细分亚群的功能富集结果
  # ######富集分析#######
  # 基因ID转换
  library(tidyverse)
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
  library(org.Mm.eg.db)
  
  ####### 上调基因
  # 批量转换
  scRNA.diff.up <- list()
  for (i in 1:length(updiff)) {
    ENTREZID <- bitr(rownames(updiff[[i]]), 
                     fromType = "SYMBOL", 
                     toType = "ENTREZID", 
                     OrgDb = "org.Mm.eg.db",
                     drop = T)
    ENTREZID <- ENTREZID$ENTREZID
    scRNA.diff.up[[i]] <- ENTREZID
  }
  names(scRNA.diff.up) <- names(updiff)
  
  # load('scRNA.diff_V59_Endo.Rdata')
  
  # 多个cluster比较
  # KEGG分析
  xx <- compareCluster(scRNA.diff.up, 
                       fun="enrichKEGG",
                       organism="mmu",
                       pvalueCutoff=0.1,
                       qvalueCutoff=0.1)
  xx <- pairwise_termsim(xx) # 获取相似性矩阵
  # dotplot(xx, showCategory = 5) #,label_format = 60 可以修改通路一行还是两行
  
  xx <- setReadable(xx,OrgDb = "org.Mm.eg.db", keyType="ENTREZID")
  
  save(xx, file = "./files/up_compareCluster.kegg.Rda")
  
  dotplot(xx, showCategory = 5,label_format = 60) +
    theme_clean() +
    scale_color_gradientn(colors =  rev(c('#318fc4', '#ca2b2b')))
  
  ggsave("./plots/up_Diffgene_kegg.pdf",
         ggplot2::last_plot(),
         width = 10,
         height = 6)
  
  # GO BP分析
  xx_BP <- compareCluster(scRNA.diff.up, 
                          fun = "enrichGO",
                          OrgDb = "org.Mm.eg.db",
                          ont = "BP",
                          pvalueCutoff=0.01,
                          qvalueCutoff=0.01)
  
  xx_BP <- setReadable(xx_BP,OrgDb ="org.Mm.eg.db" )
  
  save(xx_BP, file = "./files/up_compareCluster.bp.Rda")
  
  dotplot(xx_BP, showCategory = 5,label_format = 60) +
    theme_clean() +
    scale_color_gradientn(colors =  rev(c('#318fc4', '#ca2b2b')))
  
  ggsave("./plots/up_Diffgene_GO.pdf",
         ggplot2::last_plot(),
         width = 8,
         height = 6)
  
  
  ####### 上调基因
  # 批量转换
  scRNA.diff.dn <- list()
  for (i in 1:length(dndiff)) {
    ENTREZID <- bitr(rownames(dndiff[[i]]), 
                     fromType = "SYMBOL", 
                     toType = "ENTREZID", 
                     OrgDb = "org.Mm.eg.db",
                     drop = T)
    ENTREZID <- ENTREZID$ENTREZID
    scRNA.diff.dn[[i]] <- ENTREZID
  }
  names(scRNA.diff.dn) <- names(dndiff)
  
  # load('scRNA.diff_V59_Endo.Rdata')
  
  # 多个cluster比较
  # KEGG分析
  xx <- compareCluster(scRNA.diff.dn, 
                       fun="enrichKEGG",
                       organism="mmu",
                       pvalueCutoff=0.1,
                       qvalueCutoff=0.1)
  xx <- pairwise_termsim(xx) # 获取相似性矩阵
  # dotplot(xx, showCategory = 5) #,label_format = 60 可以修改通路一行还是两行
  
  xx <- setReadable(xx,OrgDb = "org.Mm.eg.db", keyType="ENTREZID")
  
  save(xx, file = "./files/dn_compareCluster.kegg.Rda")
  
  dotplot(xx, showCategory = 5,label_format = 60) +
    theme_clean() +
    scale_color_gradientn(colors =  rev(c('#318fc4', '#ca2b2b')))
  
  ggsave("./plots/dn_Diffgene_kegg.pdf",
         ggplot2::last_plot(),
         width = 10,
         height = 6)
  
  # GO BP分析
  xx_BP <- compareCluster(scRNA.diff.dn, 
                          fun = "enrichGO",
                          OrgDb = "org.Mm.eg.db",
                          ont = "BP",
                          pvalueCutoff=0.01,
                          qvalueCutoff=0.01)
  
  xx_BP <- setReadable(xx_BP,OrgDb ="org.Mm.eg.db" )
  
  save(xx_BP, file = "./files/dn_compareCluster.bp.Rda")
  
  dotplot(xx_BP, showCategory = 5,label_format = 60) +
    theme_clean() +
    scale_color_gradientn(colors =  rev(c('#318fc4', '#ca2b2b')))
  
  ggsave("./plots/dn_Diffgene_GO.pdf",
         ggplot2::last_plot(),
         width = 8,
         height = 6)
}

setwd('/home/bio/Projects/NC2022/neuro_sub')
