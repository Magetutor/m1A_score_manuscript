# rm(list = ls())
library(Seurat)
library(DESeq2)
library(tidyverse)
library(ggthemes)
library(monocle)
library(patchwork)
library(ggpubr)
library(ggsci)
library(ggthemes)

type_colors <- c("#E64B35CC", "#4DBBD5CC", "#00A087CC", "#3C5488FF", "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "#DC0000FF" ,"#7e6148")


group_colors <-  c('#4682b4','#33a02c','#ffa500','#ff4500')


# 加载数据

seu <- qs::qread("output/03-2.myeloid.seurat.aucell.qs")

dim(seu)

table(seu@meta.data$celltype)

# 单独提取出microglia
seurat.obj <- subset(seu,myeloid_subcluster %in% c('Homeostatic Microglia','Inflammatory Microglia','Dividing Microglia',
                                         'Migrating Microglia','Interferon Myeloid'))
table(seurat.obj@meta.data$myeloid_subcluster)

Idents(seurat.obj) <- "myeloid_subcluster"

# seurat.obj <- subset(seurat.obj,idents = 'c5', invert = TRUE)

seurat.obj@meta.data$cell_type2 <- seurat.obj@meta.data$myeloid_subcluster

DefaultAssay(seurat.obj) <- 'RNA'

#寻找高变基因
markers <- FindAllMarkers(seurat.obj, only.pos = T)
table(markers$cluster)

saveRDS(markers, "./output/monocle/micro_markers.rds")
markers <-readRDS("./biyelunwen/output/monocle/micro_markers.rds")
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top20 <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

#### 2.对象转化 ----
data <- as(as.matrix(seurat.obj@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = seurat.obj@meta.data)

fData <-
  data.frame(gene_short_name = row.names(data),
             row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
data <- data[,rownames(pd)]
cds <- newCellDataSet(
  data,
  phenoData = pd,
  featureData = fd,
  lowerDetectionLimit = 0.5,
  expressionFamily = negbinomial.size())

#### 3.质控 ----
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <-
  detectGenes(cds, min_expr = 0.1) #这一操作会在fData(cds)中添加一列num_cells_expressed
head(fData(cds))
expressed_genes <- row.names(subset(fData(cds),
                                    num_cells_expressed >= 10)) #过滤掉在小于10个细胞中表达的基因。
saveRDS(cds, "./output/monocle/Micro_cds_ori.rds")
cds0 <- cds
n_markers <- as.numeric(c(50,60,70,80,85)) 

for(i in n_markers){ 
  #### 4.设置路径 ----
  path <- paste0("/home/bio/Projects/JEM_SCI2021_GSE162610/biyelunwen/output/monocle/", i, '/')
  dir.create(path,recursive = T)
  setwd(path)
  outdir <- './files/'
  dir.create(outdir,recursive = T)
  outdir2 <- './plots/'
  dir.create(outdir2,recursive = T)
  
  #### 5.差异基因选择 ----
  # 可输入seurat筛选出的高变基因：
  # expressed_genes <- VariableFeatures(seurat.obj)
  # diff <-
  #   differentialGeneTest(cds[expressed_genes,],
  #                        fullModelFormulaStr = "~group",
  #                        cores = 10)
  # # ~后面是表示对谁做差异分析的变量，理论上可以为p_data的任意列名
  # head(diff)
  
  ## 差异表达基因作为轨迹构建的基因,差异基因的选择标准是qval<0.01,decreasing=F表示按数值增加排序
  top_markers <- markers %>% 
    filter(p_val <= 0.05) %>% 
    arrange("avg_log2FC") %>% 
    group_by(cluster) %>% 
    top_n(i, wt = avg_log2FC) %>% 
    pull(gene)
  # markers <- subset(diff, qval < 0.05)
  # markers <- markers[order(markers$qval, decreasing = F),]
  # head(markers)
  # dim(markers)
  # ## 差异基因的结果文件保存
  saveRDS(top_markers, file.path(outdir, "top_markers.rds"))
  
  #### 轨迹构建基因可视化 ----
  ordergene <- top_markers
  cds <- setOrderingFilter(cds0, ordergene)
  # #这一步是很重要的，在我们得到想要的基因列表后，我们需要使用setOrderingFilter将它嵌入cds对象，后续的一系列操作都要依赖于这个list。
  # #setOrderingFilter之后，这些基因被储存在cds@featureData@data[["use_for_ordering"]]，可以通过table(cds@featureData@data[["use_for_ordering"]])查看
  # plot_ordering_genes(cds)
  # #出的图黑色的点表示用来构建轨迹的差异基因
  # # 灰色表示背景基因。
  # # 红色的线是根据第2步计算的基因表达大小和
  # # 离散度分布的趋势(可以看到，
  # # 找到的基因属于离散度比较高的基因)
  # ggsave(file.path(outdir2, 'plot_ordering_genes.pdf'),
  #        ggplot2::last_plot(),
  #        width = 4,height = 3)
  
  # 降维
  cds <- reduceDimension(cds, 
                         max_components = 2, 
                         #num_dim = 6, 
                         #sigma=0.1,
                         reduction_method = 'DDRTree', 
                         #residualModelFormulaStr = "~orig.ident", #去除样本影响
                         verbose = F)
  # 拟时间轴轨迹构建和在拟时间内排列细胞
  cds <- orderCells(cds)
  #⚠️使用root_state参数可以设置拟时间轴的根，如下面的拟时间着色图中可以看出，左边是根。根据state图可以看出，根是State1，若要想把另一端设为根，可以按如下操作
  #cds <- orderCells(cds, root_state = 2) #把State5设成拟时间轴的起始点
  saveRDS(cds, file.path(outdir, "cds.rds"))
  
  #### 绘图 ----
  # 以时间上色
  p1 <- plot_cell_trajectory(cds,
                             color_by = "Pseudotime",
                             cell_size = 0.1,
                             size = 1,
                             show_backbone = TRUE)
  # 以状态上色
  p2 <- plot_cell_trajectory(cds,
                             color_by = "State",
                             cell_size = 0.1,
                             size = 1,
                             show_backbone = TRUE)+
    guides(color=guide_legend( override.aes = list(size=4,alpha=1)))
  # 以分组上色
  p3 <- plot_cell_trajectory(cds,
                             color_by = "time",
                             cell_size = 0.1,
                             size = 1,
                             show_backbone = TRUE)+
    scale_color_manual(values = group_colors)+
    guides(color=guide_legend(title = "time",nrow = 3, override.aes = list(size=4,alpha=1)))
  
  p4 <- plot_cell_trajectory(cds) + 
    facet_wrap(~time, nrow = 1)+
    guides(color=guide_legend( override.aes = list(size=4,alpha=1)))
  
  # 以细胞类型上色
  p5 <- plot_cell_trajectory(cds,
                             cell_size = 0.1,
                             show_branch_points = F,
                             color_by = "myeloid_subcluster",
                             size = 1,
                             show_backbone = TRUE)+
    scale_color_manual(values = type_colors)+
    guides(color=guide_legend(title = "cell_type2", nrow = 2, override.aes = list(size=4,alpha=1)))
  
  p6 <- plot_cell_trajectory(cds,cell_size = 0.1,) +
    facet_wrap(~cell_type2, nrow = 1)+
    guides(color=guide_legend(override.aes = list(size=4,alpha=1)))
  
  p7 <- plot_cell_trajectory(cds,
                             color_by = "m1a_ucell",
                             cell_size = 0.1,
                             size = 1,
                             show_backbone = TRUE)+ viridis::scale_color_viridis()
  
  
  p <- p1 + p2 + p3 + p5
  plot_layout(nrow = 3,byrow = T)
  p
  ggsave(file.path(outdir2, 'plot_cell_trajectory1.pdf'),
         p, width = 8,height = 8)
  p <- p4/p6
  ggsave(file.path(outdir2, 'plot_cell_trajectory2.pdf'),
         p, width = 10,height =7)
  
  ggsave(file.path(outdir2, 'Pseudotime.pdf'),
         p1, width = 4,height =4 )
  
  ggsave(file.path(outdir2, 'cell_type.pdf'),
         p5, width = 4,height =4 )
  
  ggsave(file.path(outdir2, 'Pseudo_m1A.pdf'),
         p7, width = 4,height =4 )
  
  
  # 密度图展示分群效果
  library(ggpubr)
  df <- pData(cds) # pData(cds) cds@phenoData@data的内容
  p1 <- ggplot(df,aes(Pseudotime, fill=cell_type2)) +
    geom_density(bw=0.5,size=1,alpha =0.5)+theme_classic2()
  p2 <- ggplot(df,aes(Pseudotime,fill=time)) +
    geom_density(bw=0.5,size=1,alpha =0.5)+theme_classic2()
  p <- p1+p2
  ggsave(file.path(outdir2, 'cell_typeandGroup.pdf'),
         p, width = 12,height =4 )
}

### 选择n_markers = 60进行重新画图

# 以时间上色
cds <- readRDS('../60/files/cds.rds')

p1 <- plot_cell_trajectory(cds,
                           color_by = "Pseudotime",
                           cell_size = 0.1,
                           size = 1,
                           show_backbone = TRUE)
# 以状态上色
p2 <- plot_cell_trajectory(cds,
                           color_by = "State",
                           cell_size = 0.1,
                           size = 1,
                           show_backbone = TRUE)+
  guides(color=guide_legend( override.aes = list(size=4,alpha=1)))
# 以分组上色
p3 <- plot_cell_trajectory(cds,
                           color_by = "Group",
                           cell_size = 0.1,
                           size = 1,
                           show_backbone = TRUE)+
  scale_color_manual(values = group_colors)+
  guides(color=guide_legend(title = "Group",nrow = 3, override.aes = list(size=4,alpha=1)))

p4 <- plot_cell_trajectory(cds) + 
  facet_wrap(~time, nrow = 1)+
  guides(color=guide_legend( override.aes = list(size=4,alpha=1)))
# 以细胞类型上色
p5 <- plot_cell_trajectory(cds,
                           cell_size = 0.1,
                           show_branch_points = F,
                           color_by = "myeloid_subcluster",
                           size = 1,
                           show_backbone = TRUE)+
  scale_color_manual(values = type_colors)+
  guides(color=guide_legend(title = "myeloid_subcluster", nrow = 2, override.aes = list(size=4,alpha=1)))

p6 <- plot_cell_trajectory(cds,cell_size = 0.1,) +
  facet_wrap(~myeloid_subcluster, nrow = 1)+
  guides(color=guide_legend(override.aes = list(size=4,alpha=1)))


p <- p1 + p2 + p3 + p5
plot_layout(nrow = 3,byrow = T)
p
ggsave(file.path(outdir2, 'plot_cell_trajectory1.pdf'),
       p, width = 8,height = 10)
p <- p4/p6
ggsave(file.path(outdir2, 'plot_cell_trajectory2.pdf'),
       p, width = 7,height =6)
ggsave(file.path(outdir2, 'Pseudotime.pdf'),
       p1, width = 3.5,height =4 )
ggsave(file.path(outdir2, 'cell_type.pdf'),
       p5, width = 3.5,height =4 )

# 密度图展示分群效果
library(ggpubr)
df <- pData(cds) # pData(cds) cds@phenoData@data的内容
p1 <- ggplot(df,aes(Pseudotime, fill=cell_type2)) +
  geom_density(bw=0.5,size=1,alpha =0.5)+theme_classic2()
p2 <- ggplot(df,aes(Pseudotime,fill=Group)) +
  geom_density(bw=0.5,size=1,alpha =0.5)+theme_classic2()
p <- p1+p2
ggsave(file.path(outdir2, 'cell_typeandGroup.pdf'),
       p, width = 12,height =4 )
