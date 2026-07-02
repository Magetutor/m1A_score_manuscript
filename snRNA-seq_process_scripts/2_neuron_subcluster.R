# 设置工作目录
setwd('/home/bio/Projects/NC2022/neuro_sub/')
# 清空变量 加载包
# rm(list = ls())
library(Seurat)
packageVersion('Seurat')
library(DESeq2)
library(tidyverse)
library(ggthemes)
# 加载数据
scobj <- qs::qread("./data/scobj_m1a_240325.qs")
m1Aset <- read.csv("./m1A_genesets.csv")


table(scobj@meta.data$condition)
table(scobj@meta.data$celltype)


sub_seurat_obj <- subset(scobj,celltype %in% c("Neurons"))

# 挑选T细胞进⾏进⼀步分群
sub_seurat_obj <- subset(scobj,celltype %in% c("Neurons"))
sub_seurat_obj
sub_seurat_obj <- SplitObject(sub_seurat_obj,
                         split.by = "orig.ident")
sub_seurat_obj
# 第一次做过线粒体、核糖体、细胞周期质控了所以直接从SCT开始
# SCTransform
# 速度较慢，样本越多越慢
# 循环运⾏
for(i in 1:length(sub_seurat_obj)){
  sub_seurat_obj[[i]] <- SCTransform(
    sub_seurat_obj[[i]],
    variable.features.n = 3000,
    #vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"), # 去除线粒体和细胞周期的影响
    verbose = FALSE)
}

# 数据整合
# load("Endo.data.SCT.Rda")
features <- SelectIntegrationFeatures(object.list = sub_seurat_obj,
                                      nfeatures = 3000)
# 准备整合
sub_seurat_obj <- PrepSCTIntegration(object.list = sub_seurat_obj,
                                anchor.features = features)
names(sub_seurat_obj)

# 寻找锚定基因(速度很慢)
AnchorSet <- FindIntegrationAnchors(object.list = sub_seurat_obj,
                                    reference = 1,
                                    normalization.method = "SCT",
                                    anchor.features = features
)
# save(AnchorSet, file = "AnchorSet.Rda")
# 整合数据
# load("AnchorSet.Rda")
sub_seurat_obj <- IntegrateData(anchorset = AnchorSet,
                           normalization.method = "SCT")


######直接分群


DefaultAssay(sub_seurat_obj) <- 'integrated'

sub_seurat_obj <- ScaleData(sub_seurat_obj, verbose = FALSE)

sub_seurat_obj <- RunPCA(sub_seurat_obj, npcs = 50, verbose = FALSE)
ElbowPlot(sub_seurat_obj,ndims = 50)
sub_seurat_obj <- RunUMAP(sub_seurat_obj, reduction = "pca", dims = 1:20)
sub_seurat_obj <- RunTSNE(sub_seurat_obj, reduction = "pca", dims = 1:30)

sub_seurat_obj <- FindNeighbors(sub_seurat_obj, k.param = 20, dims = 1:20)

sub_seurat_obj <- FindClusters(sub_seurat_obj, 
                               reduction="umap", resolution = 0.1)

sub_seurat_obj <- sub_seurat_obj %>% FindClusters(resolution =seq(from = 0.1,to = 1.2, by = 0.1))



library(clustree)
pdf('clustree.pdf',width = 14,height = 8)
clustree(sub_seurat_obj) + theme(legend.position = "bottom") +
  scale_color_manual(values = c(pal_npg()(10),pal_d3()(8))) +
  scale_edge_color_continuous(low = "grey80", high = "red")
dev.off()


DimPlot(sub_seurat_obj, reduction = "umap",pt.size = 0.5,
        # split.by = 'condition',
        label = T,label.box = T,label.size = 4,repel = T) + 
  theme_bw() +
  labs( x= "UMAP_1",y= "UMAP_2",title = "") +
  theme(panel.grid = element_blank()
        #,aspect.ratio = 1
  )+NoLegend()


DefaultAssay(sub_seurat_obj) <- 'integrated'
markers <- FindAllMarkers(
  sub_seurat_obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)


top5.markers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)



top.markers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)


top.markers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC)


top.markers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 50, order_by = avg_log2FC)


top100.markers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = avg_log2FC)

write.table(top.markers,
            './results/seurat/NK/dim30r0点4/top.markers100.csv',
            quote = F,
            sep = ",",
            row.names = F)

source("/home/lixy/R/x86_64-pc-linux-gnu-library/4.3/hephaestus/R/data.R")
source("/home/lixy/R/x86_64-pc-linux-gnu-library/4.3/hephaestus/R/objects.R")
source("/home/lixy/R/x86_64-pc-linux-gnu-library/4.3/hephaestus/R/rip-seq.R")
source("/home/lixy/R/x86_64-pc-linux-gnu-library/4.3/hephaestus/R/rna-seq.R")
source("/home/lixy/R/x86_64-pc-linux-gnu-library/4.3/hephaestus/R/scrna-seq.R")
source("/home/lixy/R/x86_64-pc-linux-gnu-library/4.3/hephaestus/R/utilities.R")
source("/home/lixy/R/x86_64-pc-linux-gnu-library/4.3/hephaestus/R/visualization.R")

DefaultAssay(sub_seurat_obj) <- 'RNA'


features <- c("P2ry12",'Igf1',"Msr1","Cdk1")

features <- c('Jun','P2ry12','Siglech','Selplg','Btg2','Fos','Tmem119','Egr1','Klf2','Hspa1a', # Homeostatic_micro
              'Il1b','Apoe','Igf1','Apoc4','Apoc1','Cxcl2','Lyz2', # inflammatory_micro
              'Stmn1','Tuba1b','Pclaf','H2az1','Hmgb2','Birc5','Ube2c','Tubb5','Ran','Cks1b', #Dividing_micro
              'Lgals3','Npy','Vim','Plin2','Adam8','Lgals1','Fabp5','Spp1','Gpnmb','Ccl9',# Migrating_micro
              'Isg15','Ifitm3','Ifit3','Irf7','Rsad2','Ifi204','Rtp4','Bst2','Ccl5','Cxcl10' #Interferon_mmicro
)
# 筛选一下
features <- c('P2ry12','Btg2','Fos','Tmem119', # Homeostatic_micro
              'Igf1','Apoc1','Cxcl2','Lyz2', # inflammatory_micro
              'Stmn1','Pclaf','Hmgb2','Ube2c','Tubb5', #Dividing_micro
              'Fabp5','Spp1','Gpnmb','Ccl9'# Migrating_micro
)

feture22 <- c('P2ry12','Btg2','Fos','Tmem119', # Homeostatic_micro
              'Stmn1','Pclaf','Hmgb2','Ube2c','Tubb5', #Dividing_micro
              'Ear2','Trem1','Thbs4','Hp','Trem3',
              'Ccnb2','Hcar2','Tnfsf8',
              'Spib','Mzb1','Cd7','Rell1',
              'Celf5','Cep290','Crip2','Arap2'
)
top5.markers$gene
DotPlot(sub_seurat_obj, features = unique(features), col.min = 0) + 
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))


DotPlot(sub_seurat_obj, dot.scale = 4, col.min = 0, # col.max = 1., # scale = F,
        features = unique(top5.markers$gene)) +
  scale_color_distiller(palette = "RdBu", direction = -1) + # YlOrRd
  coord_fixed() +
  theme_cat() +
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

# 分配细胞名字
merged <- sub_seurat_obj

new.cluster.ids <- c(
  "c0",# 00
  "c1",#01
  "c2",#02
  "c3",#03
  "c4",#04
  "c5")#05

names(new.cluster.ids) <- levels(merged)
merged <- RenameIdents(merged,
                       new.cluster.ids)
merged$cell_type <- Idents(merged)

celltype <- unique(merged@meta.data$cell_type)

# 给markers添加c1标签
markers1 <- markers
markers1$cluster <- as.character(markers1$cluster)
markers1 <- markers1 %>% mutate(cell_type = case_when(
  cluster == "0" ~ "c0",
  cluster == "1" ~ "c1",
  cluster == "2" ~ "c2",
  cluster == "3" ~ "c3",
  cluster == "4" ~ "c4"
))
paste(markers1$gene[markers1$cell_type=='c0'],collapse = ',')

paste(my_vector, collapse = ";")

write.csv(markers1,'markers1.csv',row.names = F)

# 使用top100marker功能富集来看看
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)
for(i in celltype){
  path <- paste0('/home/zhangc/analysis/zika_new/Figure_code/20240101_microglia_resub', i, '/')
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
              file = file.path(outdir , paste0(i,'_UP_GOBP.csv')),
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
setwd('/home/zhangc/analysis/zika_new/Figure_code/20240101_microglia_resub/')
write.csv(markers1,'microglia_sub_marker.csv',row.names = F)

# 分配细胞名字
merged <- sub_seurat_obj
Idents(merged) <- 'seurat_clusters'

new.cluster.ids <- c(
  "Homeo",# 
  "Dividing",#01
  "Adhesion",#02
  "Genesis",#03
  "Regulatory",#04
  "c5")#05

names(new.cluster.ids) <- levels(merged)
merged <- RenameIdents(merged,
                       new.cluster.ids)
merged$cell_type2 <- Idents(merged)

DimPlot(merged , reduction = "umap",raster=FALSE,label = T,pt.size = 0.4)

p1 <- DimPlot(merged, reduction = "umap",raster=FALSE,
              label = T,label.box = T,label.size = 2,repel = T) + theme_bw() +
  labs( x= "UMAP_1",y= "UMAP_2",title = "") +
  theme(panel.grid = element_blank(),aspect.ratio = 1)+NoLegend()
p2 <- DimPlot(merged, reduction = "tsne",#split.by = "Group",#Group.by = "orig.ident",
              label = T,label.box = T,raster=FALSE,label.size = 2,repel = T) + theme_bw() +
  labs( x= "TSNE_1",y= "TSNE_2",title = "") +
  theme(panel.grid = element_blank(),aspect.ratio = 1)+NoLegend()
p1|p2

DimPlot(merged, reduction = "umap",raster=FALSE,group.by = "Group",
        label = T,label.box = T,label.size = 2,repel = T) + theme_bw() +
  labs( x= "UMAP_1",y= "UMAP_2",title = "") +
  theme(panel.grid = element_blank(),aspect.ratio = 1)+NoLegend()

DimPlot(merged, reduction = "umap",raster=FALSE,split.by = "Group",
        label = T,repel = T) + theme_bw() +
  labs( x= "UMAP_1",y= "UMAP_2",title = "") +
  theme(panel.grid = element_blank(),aspect.ratio = 1)


saveRDS(merged,"./20240108_microglia.Rds")

####============== Figure S1 细胞类型用疾病填充  ===================
seurat <- merged

#### 各疾病细胞比例计算 ----
# 查看整体每个细胞类型的数量
table(Idents(seurat))
# 计算整体每个细胞类型所占百分比
prop.table(table(Idents(seurat))) * 100

# 查看各个分组每个细胞类型数量
table(Idents(seurat), seurat$Group)
table(seurat$Group)
# 计算各个分组每个细胞类型所占百分比
prop.table(table(Idents(seurat), seurat$Group), margin = 2) * 100
#### 可视化 
meta <- seurat@meta.data

allcolour <- c("#E64B35CC", "#4DBBD5CC", "#00A087CC", "#3C5488FF",
               "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "#DC0000FF" ,"#7e6148")

meta2 <- subset(meta, meta$Group == c("MOCK","V_FS","V_59"))
p <- ggplot(data = meta, aes(x = cell_type2, fill = Group)) +
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
  scale_fill_manual(values = allcolour) +
  # geom_hline(yintercept = 0.25,
  #            linetype = "dashed",
  #            size = 0.5) +
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))


ggsave(
  file.path( "./percentage-Group.pdf"),
  plot = p,
  height = 4.5,
  width = 6
)

####============== Figure S1 疾病用细胞类型填充  ===================

type_colors2 = c("#9ACBDE","#F47892","#64AE59","#6DCCDD",
                 "#ad6593","#F6A395","#C4A5DE","#B9C984",
                 "#F8AD77","#FF6347","#EDC6DD")


type_colors2=c("#9ACBDE","#F47892","#64AE59","#ad6593","#6DCCDD",
               "#F6A395","#C4A5DE","#B9C984",
               "#F8AD77","#FF6347","#EDC6DD")

x = factor(meta$Group,levels=c("MOCK","V_FS","V_59"))



p2 <- ggplot(data = meta, 
             aes(x = factor(Group, levels=c("MOCK","V_FS","V_59")),
                 fill = cell_type2)) +
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
  scale_fill_manual(values = type_colors2) +
  # geom_hline(yintercept = 0.25,
  #            linetype = "dashed",
  #            size = 0.5) +
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))


p2


ggsave(
  file.path( "./ratio-celltype.pdf"),
  plot = p2,
  height = 6,
  width = 4
)

####cluster相似性----
#### 导入包 ----
#file.choose()
####样本相似性----
library(Seurat)
library(dplyr)
library(ggplot2)
library("corrplot")
av_expr <-
  AverageExpression(
    merged,
    assays = "RNA",
    group.by = "cell_type"
  )[["RNA"]] %>%
  as.data.frame()

corr<-cor(av_expr,method="spearman")
#pearson相关系数：适用于连续型变量，且变量服从正态分布的情况，为参数性的相关系数。
#spearman等级相关系数：适用于连续型及分类型变量，为非参数性的相关系数。
#kendall秩相关系数：适用于定序变量或不满足正态分布假设的等间隔数据。

library(pheatmap)
p <- pheatmap(corr, shape = "circle",cellwidth = 15,cellheight =15,
              # color = colorRampPalette(c("navy", "white", "firebrick3"))(50)
              color = colorRampPalette(c("#B9d6e9","#2e86c1"))(50),
              # cluster_rows = F,cluster_cols = F
)
p

ggsave(file.path(outdir,'corr.pdf'),
       p,
       height=5,
       width=5)



corrplot(corr,type="upper", #保留右上部分图形
         addCoef.col="black",#添加相关系数，颜色为黑色
         diag=F)#去掉自身相关

col <- colorRampPalette(c("#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
                          "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080"))
corrplot(corr,
         method="pie",#method="color",#调整为正方形
         col=col(500000), #颜色调整
         type="upper", #保留右上部分
         order="hclust", #层次聚类
         addCoef.col = "black", #添加相关系数
         tl.col="black", tl.srt=45 #修改字体
         #diag=FALSE #去除自身相关
)

library(pheatmap)
p <- heatmap(corr)

ggsave(p,'./data/intermediate/cor_cluster.pdf',width = 5,height =5 )

ggsave(file.path(outdir,paste0(case, '样本相似性.pdf')),
       p,#最后一张图
       height=5,
       width=5)



feture22 <- c('P2ry12','Btg2','Fos','Tmem119', # Homeostatic_micro
              'Stmn1','Pclaf','Hmgb2','Ube2c','Tubb5', #Dividing_micro
              'Ear2','Trem1','Thbs4','Hp','Trem3',
              'Ccnb2','Hcar2','Tnfsf8',
              'Spib','Mzb1','Cd7','Rell1',
              'Celf5','Cep290','Crip2','Arap2'
)

DotPlot(seurat, #dot.scale = 4, 
        #col.min = 1,  col.max = 4., # scale = F,
        features = unique(top5.markers$gene)) +
  scale_color_distiller(palette = "RdBu", direction = -1) + # YlOrRd
  coord_fixed() +
  theme_cat() +
  theme(
    axis.title = element_blank(),
    #legend.margin = margin(b = -8),
    #legend.position = "top",
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


ggsave("./micro-subset-markersdotplot.pdf",
       ggplot2::last_plot(),width = 7,height = 4)

