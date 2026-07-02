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
m1Aset <- read.csv("../m1A_genesets.csv")
#load("./data/neuron_subcluster_20250622.Rda")
sub_seurat_obj1 <- qs::qread("./data/neuron_subcluster_20250622.qs")
# 使用的数据就是神经元的亚分群数据 sub_seurat_obj1
# 使用SuperCell包进行
# if (!requireNamespace("remotes")) install.packages("remotes")
# remotes::install_github("GfellerLab/SuperCell")

library(SuperCell)

# 定义参数
MC_tool = "SuperCell"
proj_name <- "subneuron"
gamma = 20 # the requested graining level.
k_knn = 30 # the number of neighbors considered to build the knn network.
nb_var_genes = 2000 # number of the top variable genes to use for dimensionality reduction 
nb_pc = 30 # the number of principal components to use.  

# Metacells identification
MC <- SuperCell::SCimplify(Seurat::GetAssayData(sub_seurat_obj1, slot = "data"),  # single-cell log-normalized gene expression data
                           k.knn = k_knn,
                           gamma = gamma,
                           # n.var.genes = nb_var_genes,  
                           n.pc = nb_pc
                           #genes.use = Seurat::VariableFeatures(sub_seurat_obj1)
)

MC.GE <- supercell_GE(Seurat::GetAssayData(sub_seurat_obj1, slot = "counts"),
                      MC$membership,
                      mode =  "sum"
)
dim(MC.GE) 

# 注释细胞
print(annotation_label)
#> [1] "celltype_simplified"
MC$annotation <- supercell_assign(clusters = sub_seurat_obj1@meta.data[, "celltype"], # single-cell annotation
                                  supercell_membership = MC$membership, # single-cell assignment to metacells
                                  method = "absolute"
)

MC$condition <- supercell_assign(clusters = sub_seurat_obj1@meta.data[, "condition"], # single-cell annotation
                                  supercell_membership = MC$membership, # single-cell assignment to metacells
                                  method = "absolute"
)

head(MC$annotation)
head(MC$condition)
# 可以先画个图看看
# plot network of metacells
supercell_plot(
  MC$graph.supercells, 
  group = MC$annotation, 
  lay.method = 'drl',
  seed = 1, 
  alpha = -pi/2,
  main  = "Metacells colored by sc assignment"
)

pdf("./figure/Figure_4_supercell/A_spercell_origin_celltype.pdf",height = 6,width = 10)
mysupercell(
  MC$graph.supercells, 
  group = MC$annotation, 
  lay.method = 'drl',
  seed = 1, 
  alpha = -pi/2,
  main  = "Metacells colored by sc_celltype assignment"
)
dev.off()

pdf("./figure/Figure_4_supercell/A_spercell_origin_group.pdf",height = 6,width = 10)
mysupercell(
  MC$graph.supercells, 
  group = MC$condition, 
  lay.method = 'drl',
  seed = 1, 
  alpha = -pi/2,
  main  = "Metacells colored by sc_group assignment"
)
dev.off()


# 查看每个元细胞含有的细胞系类别的纯度
purity <- supercell_purity(clusters = sub_seurat_obj1@meta.data[, "celltype"], 
                           supercell_membership =  MC$membership, method = 'entropy')

pdf("./figure/Figure_4_supercell/B_metacell_purity.pdf",height = 5,width = 6)
hist(purity, main = "Purity of metacells in terms of composition")
dev.off()

# MC$purity <- purity

# 直接看基因表达
genes.to.plot <- c("Alkbh1","Alkbh3","Fto","Trmt10c","Trmt6")

# 导出成seurat对象
colnames(MC.GE) <- as.character(1:ncol(MC.GE))
MC.seurat <- CreateSeuratObject(counts = MC.GE, 
                                meta.data = data.frame(size = as.vector(table(MC$membership)))
)
# 添加细胞类型
MC.seurat[["celltype_simplified"]] <- MC$annotation
MC.seurat[["celltype"]] <- MC$annotation
# 添加分组信息
MC.seurat[["condition"]] <- MC$condition

# save single-cell membership to metacells in the MC.seurat object
MC.seurat@misc$cell_membership <- data.frame(row.names = names(MC$membership), membership = MC$membership)
MC.seurat@misc$var_features <- MC$genes.use 

# Save the PCA components and genes used in SCimplify  
PCA.res <- irlba::irlba(scale(Matrix::t(sub_seurat_obj1@assays$RNA@data[MC$genes.use, ])), nv = nb_pc)
pca.x <- PCA.res$u %*% diag(PCA.res$d)
rownames(pca.x) <- colnames(sub_seurat_obj1@assays$RNA@data)
MC.seurat@misc$sc.pca <- CreateDimReducObject(
  embeddings = pca.x,
  loadings = PCA.res$v,
  key = "PC_",
  assay = "RNA"
)
if(packageVersion("Seurat") >= 5) {
  MC.seurat[["RNA"]] <- as(object = MC.seurat[["RNA"]], Class = "Assay")
}
print(paste0("Saving metacell object for the ", proj_name, " dataset using ", MC_tool))
#> [1] "Saving metacell object for the bmcite dataset using SuperCell"

# 保存数据
save(MC.seurat, file = './figure/Figure_4_supercell/data/subneuron_MC.seurat.Rda')

#### 上面构建了metacell，接下来就是继续处理，然后按照seurat流程走
# 先来看看 这个时候是count，所以画图都是平行的
FeatureScatter(object = MC.seurat, feature1 = 'Il6', feature2 = 'Ythdc1')

MC_tool = "SuperCell"
proj_name = "subneuron"
annotation_column = "celltype"

celltypes <- c("Slit2_IN", "Npy_IN","Gal_IN", "Cck_EN", "Sox5_EN", "Pde11a_EN", "Tac2_EN")

celltype_colors <- c("#3477a9", "#96c3d8","#1E88E5", "#d62e2d", "#f47d2f",
                     "#F06292", "#4a9d47")
# c("#1E88E5", "#FFC107", "#004D40", "#9E9D24",
#   "#F06292", "#546E7A", "#D4E157", "#76FF03", 
#   "#26A69A", "#AB47BC", "#D81B60", "#42A5F5",
#   "#2E7D32", "#FFA726", "#5E35B1", "#EF5350","#6D4C41")
names(celltype_colors) <- celltypes

# 定义好细胞ident
MC.seurat@meta.data$celltype <- MC.seurat@meta.data$celltype_simplified

MC.seurat@meta.data$celltype <- factor(MC.seurat@meta.data$celltype,
                                       levels = celltypes )
Idents(MC.seurat) <- 'celltype'

##### Seurat降维流程#####
MC.seurat <- NormalizeData(MC.seurat)
MC.seurat <- FindVariableFeatures(MC.seurat, selection.method = "vst", nfeatures = 2000)
MC.seurat <- ScaleData(MC.seurat)
#> Centering and scaling data matrix
MC.seurat <- RunPCA(MC.seurat, verbose = F)
MC.seurat <- RunUMAP(MC.seurat, dims = 1:30, verbose = F, min.dist = 1)
#> Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
#> To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
#> This message will be shown once per session

data <- cbind(Embeddings(MC.seurat, reduction = "umap"),
              data.frame(size = MC.seurat$size,
                         cell_type = MC.seurat@meta.data[, annotation_column]))
colnames(data)[1:2] <- c("umap_1", "umap_2")

data$cell_type <- factor(data$cell_type,levels = celltypes)
p_annot <- ggplot(data, aes(x= umap_1, y=umap_2, color = cell_type)) + geom_point(aes(size=size)) +
  ggplot2::scale_size_continuous(range = c(0.5,  0.5*max(log((data$size))))) +
  ggplot2::scale_color_manual(values = celltype_colors) +
  theme_classic() + guides(color=guide_legend(ncol=1))+
  theme(axis.title = element_text(size = 12, colour = "black"),
        axis.text = element_text(size = 10, colour = "black"))+
  guides(color=guide_legend(override.aes = list(size=3,alpha=1))) # 注意fill和color的变化
p_annot

ggsave(
  file.path( "./figure/Figure_4_supercell/C_Supercell_Umap.pdf"),
  plot = p_annot,
  height = 5,
  width = 6.5)

##### Seurat聚类流程#####
MC.seurat <- FindNeighbors(MC.seurat, reduction = "pca", dims = 1:30)

MC.seurat <- FindClusters(MC.seurat, resolution = 0.1)

data <- cbind(Embeddings(MC.seurat, reduction = "umap"),
              data.frame(size = MC.seurat$size,
                         cluster = MC.seurat$seurat_clusters))

colnames(data)[1:2] <- c("umap_1", "umap_2")
p_cluster <- ggplot(data, aes(x= umap_1, y=umap_2, color = cluster)) + geom_point(aes(size=size)) +
  ggplot2::scale_size_continuous(range = c(0.5, 0.5*max(log1p((data$size))))) +
  theme_classic() + guides(color=guide_legend(ncol=1))+
  theme(axis.title = element_text(size = 12, colour = "black"),
        axis.text = element_text(size = 10, colour = "black"))+
  guides(color=guide_legend(override.aes = list(size=3,alpha=1))) # 注意fill和color的变化
p_cluster

##### 差异分析 ####
# Set idents to metacell clusters

Idents(MC.seurat) <- "seurat_clusters"
cells_markers <- FindMarkers(MC.seurat, ident.1 = "0", only.pos = TRUE,
                             logfc.threshold = 0.25, min.pct = 0.1, 
                             test.use = 'wilcox', pseudocount.use = 1)
#genes.to.plot <- c("Alkbh1","Alkbh3","Fto","Trmt10c","Trmt6")
test_marker <- c("Slit2","Esrrg","Cdh18",'Snrpn')
cells_markers[test_marker, ]

VlnPlot(MC.seurat, test_marker, ncol = 3, pt.size = 0.0)
p_cluster + p_annot

##### Visualize gene-gene correlation ####
# 单细胞水平相关性
# cells_markers <- cells_markers[order(cells_markers$avg_log2FC, decreasing = T),]
gene_x <- test_marker[1:3] 
gene_y <- test_marker[4]

alpha <- 0.7

p.sc <- SuperCell::supercell_GeneGenePlot(
  GetAssayData(sub_seurat_obj1, slot = "data"),
  gene_x = gene_x,
  gene_y = gene_y,
  clusters = sub_seurat_obj1@meta.data[, annotation_column],
  sort.by.corr = F,
  alpha = alpha,
  color.use = celltype_colors
)
p.sc$p 

# metacell水平相关性

p.MC <- SuperCell::supercell_GeneGenePlot(GetAssayData(MC.seurat, slot = "data"),
                                          gene_x = gene_x,
                                          gene_y = gene_y,
                                          clusters = MC.seurat@meta.data[, annotation_column],
                                          sort.by.corr = F, supercell_size = MC.seurat$size,
                                          alpha = alpha,
                                          color.use = celltype_colors)
p.MC$p 

# 只看一个细胞类型
Slit2_IN <- subset(MC.seurat,celltype =='Slit2_IN')

p.MC <- SuperCell::supercell_GeneGenePlot(GetAssayData(Slit2_IN, slot = "data"),
                                          gene_x = 'Trmt10c',
                                          gene_y = 'Fto',
                                         # clusters = MC.seurat@meta.data[, annotation_column],
                                          # sort.by.corr = F, supercell_size = MC.seurat$size,
                                          alpha = alpha,
                                          color.use = celltype_colors)
p.MC$p 

# 保存数据
save(MC.seurat, file = './figure/Figure_4_supercell/data/subneuron_MC.seurat.Rda')
