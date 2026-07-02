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
load("./data/neuron_subcluster_m1ascore.Rda")

# 主要思路
# 找出每个时间点，Slit2_IN神经元在各个时间点共同的差异功能
#### figureA 通路UPSET venn#########
# 上调基因
load('./figure/Figure_3/B_1dpi/files/up_compareCluster.bp.Rda')
B_1dpi_dfun_up <- xx_BP@compareClusterResult
load('./figure/Figure_3/C_1wpi/files/up_compareCluster.bp.Rda')
C_1wpi_dfun_up <- xx_BP@compareClusterResult
load('./figure/Figure_3/D_3wpi/files/up_compareCluster.bp.Rda')
D_3wpi_dfun_up <- xx_BP@compareClusterResult
load('./figure/Figure_3/E_6wpi/files/up_compareCluster.bp.Rda')
E_6wpi_dfun_up <- xx_BP@compareClusterResult

load('./figure/Figure_3/B_1dpi/files/dn_compareCluster.bp.Rda')
B_1dpi_dfun_dn <- xx_BP@compareClusterResult
load('./figure/Figure_3/C_1wpi/files/dn_compareCluster.bp.Rda')
C_1wpi_dfun_dn <- xx_BP@compareClusterResult
load('./figure/Figure_3/D_3wpi/files/dn_compareCluster.bp.Rda')
D_3wpi_dfun_dn <- xx_BP@compareClusterResult
load('./figure/Figure_3/E_6wpi/files/dn_compareCluster.bp.Rda')
E_6wpi_dfun_dn <- xx_BP@compareClusterResult

# 准备upset图的文件
Slit2_list_up <- list(
  B_1dpi = subset(B_1dpi_dfun_up, Cluster == 'Slit2_IN') %>% 
    .[,'Description'],
  C_1wpi = subset(C_1wpi_dfun_up, Cluster == 'Slit2_IN') %>% 
    .[,'Description'],
  D_3wpi = subset(D_3wpi_dfun_up, Cluster == 'Slit2_IN') %>% 
    .[,'Description'],
  E_6wpi = subset(E_6wpi_dfun_up, Cluster == 'Slit2_IN') %>% 
    .[,'Description']
)

Slit2_list_dn <- list(
  B_1dpi = subset(B_1dpi_dfun_dn, Cluster == 'Slit2_IN') %>% 
    .[,'Description'],
  C_1wpi = subset(C_1wpi_dfun_dn, Cluster == 'Slit2_IN') %>% 
    .[,'Description'],
  D_3wpi = subset(D_3wpi_dfun_dn, Cluster == 'Slit2_IN') %>% 
    .[,'Description'],
  E_6wpi = subset(E_6wpi_dfun_dn, Cluster == 'Slit2_IN') %>% 
    .[,'Description']
)
# devtools::install_github("GuangchuangYu/UpSetR")

library(UpSetR)
groupcol <- c("#E64B35CC", "#4DBBD5CC", "#00A087CC", "#3C5488FF",
              "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "#DC0000FF" ,"#7e6148")

p1 <- upset(fromList(Slit2_list_up),
            query.legend = "bottom",
            nsets = 70,
            #number.angles = 0,
            mb.ratio = c(0.7, 0.3),
            point.size = 3, #point.size更改矩阵中圆圈的大小。 
            line.size = 1, #line.size更改连接矩阵中圆圈的线的大小。
            mainbar.y.label = "intersection size", 
            sets.x.label = "Numbers of Change pathway", 
            text.scale = c(1.5, 1.5, 1.5, 1.5),
            order.by = "freq",
            main.bar.color = "#2a83a2", 
            #sets.bar.color = "#3b7960",
            sets.bar.color = c("#00A087CC", "#3C5488FF","#F39B7FFF", "#4DBBD5CC"),
            queries = list(list(query = intersects, 
                                params = list("B_1dpi","C_1wpi" ,
                                              "D_3wpi","E_6wpi"), active = T,
                                color="#e6b422", query.name = ""))
);p1


library(ggvenn)
ggvenn(Slit2_list_up)
#美化
change_vennplot <- ggvenn(Slit2_list_up,
                          fill_alpha = 0.7,
                          show_percentage = T,
                          stroke_color = "white",
                          fill_color = c("#4DBBD5CC","#00A087CC","#3C5488FF","#F39B7FFF" ),
                          set_name_color = c("#4DBBD5CC","#00A087CC","#3C5488FF","#F39B7FFF"),
                          set_name_size = 3.5,stroke_size = 0.5,
                          text_size = 3);change_vennplot

#将绘图对象转换为可以使用ggplot2的对象
library(ggplotify)
change_upsetplot <- as.ggplot(p1) ;change_upsetplot 

#使用ggimage包嵌图，将结合change_vennplot和upset两种图
library(ggimage)
plot_all <- change_upsetplot + geom_subview(subview = change_vennplot + 
                                theme_void(),x=0.7, y=0.75, w=0.45, h=0.45);plot_all


ggsave(
  file.path( "./figure/Figure_3/UPgene_A_shared_pathway_Upset_venn.pdf"),
  plot = plot_all,
  height = 6,
  width = 10)


# 下调差异基因绘图
p1 <- upset(fromList(Slit2_list_dn),
            query.legend = "bottom",
            nsets = 70,
            #number.angles = 0,
            mb.ratio = c(0.7, 0.3),
            point.size = 3, #point.size更改矩阵中圆圈的大小。 
            line.size = 1, #line.size更改连接矩阵中圆圈的线的大小。
            mainbar.y.label = "intersection size", 
            sets.x.label = "Numbers of Change pathway", 
            text.scale = c(1.5, 1.5, 1.5, 1.5),
            order.by = "freq",
            main.bar.color = "#2a83a2", 
            #sets.bar.color = "#3b7960",
            sets.bar.color = c("#F39B7FFF","#4DBBD5CC", "#3C5488FF","#00A087CC" ),
            queries = list(list(query = intersects, 
                                params = list("B_1dpi","C_1wpi" ,
                                              "D_3wpi","E_6wpi"), active = T,
                                color="#e6b422", query.name = ""))
);p1


library(ggvenn)
ggvenn(Slit2_list_dn)
#美化
change_vennplot <- ggvenn(Slit2_list_dn,
                          fill_alpha = 0.7,
                          show_percentage = T,
                          stroke_color = "white",
                          fill_color = c("#4DBBD5CC","#00A087CC","#3C5488FF","#F39B7FFF" ),
                          set_name_color = c("#4DBBD5CC","#00A087CC","#3C5488FF","#F39B7FFF"),
                          set_name_size = 3.5,stroke_size = 0.5,
                          text_size = 3);change_vennplot

#将绘图对象转换为可以使用ggplot2的对象
library(ggplotify)
change_upsetplot <- as.ggplot(p1) ;change_upsetplot 

#使用ggimage包嵌图，将结合change_vennplot和upset两种图
library(ggimage)
plot_all <- change_upsetplot + geom_subview(subview = change_vennplot + 
                                              theme_void(),x=0.75, y=0.75, w=0.45, h=0.45);plot_all


ggsave(
  file.path( "./figure/Figure_3/Downgene_A_shared_pathway_Upset_venn.pdf"),
  plot = plot_all,
  height = 6,
  width = 10)

####figureb 具体通路##########
shared_up <- intersect(intersect(intersect(Slit2_list_up$B_1dpi, Slit2_list_up$C_1wpi),Slit2_list_up$D_3wpi),
                    Slit2_list_up$E_6wp)

shared_dn <- intersect(intersect(intersect(Slit2_list_dn$B_1dpi, Slit2_list_dn$C_1wpi),Slit2_list_dn$D_3wpi),
                       Slit2_list_dn$E_6wp)


# 这个数据需要保存一下
save(Slit2_list_up,Slit2_list_dn,shared_up,shared_dn,file = './data/Slit2_diff_function.Rda')

# 进行打分

library(msigdbr)
mouse <- msigdbr(species = "Mus musculus")
table(mouse$gs_cat)
mouse_GO_bp = msigdbr(species = "Mus musculus",
                      category = "C5", #GO在C5
                      subcategory = "GO:BP") %>% 
  dplyr::select(gs_name,gene_symbol)

mouse_GO_bp_Set = mouse_GO_bp %>% split(x = .$gene_symbol, f = .$gs_name)

names(mouse_GO_bp_Set) <- names(mouse_GO_bp_Set) %>% 
  gsub('GOBP_','',.) %>% 
  gsub('_',' ',.) %>% 
  str_to_lower()

subneuro_ucell<- AddModuleScore_UCell(sub_seurat_obj1, 
                                      maxRank = 5000,
                                      ncores = 50,
                                      features = mouse_GO_bp_Set, name = "_UCell")

# 提取出打分的值
allpathwy <- grep("_UCell",colnames(subneuro_ucell@meta.data))

table(names(mouse_GO_bp_Set) %in% shared_up)

table(names(mouse_GO_bp_Set) %in% shared_dn)



# 使用msigdbr很多通路没有了，直接使用clusterprofiler来搞
# 先提出需要的通路对应的goid
GO_data <- clusterProfiler:::get_GO_data("org.Mm.eg.db", "ALL", "SYMBOL") 
names(GO_data)
pathtogen <- GO_data$PATHID2EXTID
pathtoname <- GO_data$PATHID2NAME

table(pathtoname %in% shared_up)
goidup <- pathtoname[pathtoname %in% shared_up]
goidupid <- names(goidup)

table(pathtoname %in% shared_dn)
goiddn <- pathtoname[pathtoname %in% shared_dn]
goiddnid <- names(goiddn)

# 再根据goid提取对应通路和基因
pathtogen
upsig <- pathtogen[goidupid]
names(upsig) <- goidup

dnsig <- pathtogen[goiddnid]
names(dnsig) <- goiddn

#先打个下调的分看看
subneuro_ucell_dn<- AddModuleScore_UCell(sub_seurat_obj1, 
                                      maxRank = 1500,
                                      ncores = 50,
                                      features = dnsig, name = "_UCell")

dncordata <- subneuro_ucell_dn@meta.data[]
dncordata <- dncordata[dncordata$celltype == 'Slit2_IN',]
colnames(dncordata)

#简单画相关性
library(corrplot)
M <- dncordata[,c(13,21:45)]
colnames(M)[2:26] <- names(goiddn)

cor<-cor(M)
corrplot(cor,
         method = 'ellipse',#图形形状,椭圆的大小和方向表示相关性强度和方向，可选'circle' (default), 'square', 'ellipse', 'number', 'pie', 'shade' and 'color'
         type='upper',#矩阵填充方式，可选'full'(defult),'lower','upper'
         order = 'AOE',#首先安装相关性的绝对值大小进行排序，然后调整变啊零的顺序，最后根据特征值的大小再次调整排序
         title='',#不设置标题
         tl.col = 'black',#文字标签颜色,默认'red'
         tl.cex=1,#文字标签大小
         addgrid.col = 'black'#网格颜色，默认'black'
)

# m1A和通路单独做相关性
# 循环计算m1A和通路相关性
correlation <- data.frame()
##2.准备数据
data <- M
##3.获取基因列表
genelist <- colnames(M)[2:26]
##4.指定基因
gene <- "m1a_ucell"
genedata <- as.numeric(data[,gene])
pathdata <- M[,-1]
##5.开始for循环
for(i in 1:length(genelist)){
  ## 1.指示
  print(i)
  ## 2.计算
  dd = cor.test(genedata,as.numeric(pathdata[,i]),
                method="spearman")
  ## 3.填充
  correlation[i,1] = gene
  correlation[i,2] = genelist[i]
  correlation[i,3] = dd$estimate
  correlation[i,4] = dd$p.value
}

colnames(correlation) <- c("m1a","pathway","cor","p.value")
## 6.p值矫正
### 校正后的P值会扩大，其更可靠
correlation$padjust = p.adjust(correlation$p.value,method = "BH")



#上调基因pathway分值
subneuro_ucell_up<- AddModuleScore_UCell(sub_seurat_obj1, 
                                         maxRank = 1500,
                                         ncores = 50,
                                         features = upsig, name = "_UCell")

upcordata <- subneuro_ucell_up@meta.data
upcordata <- upcordata[upcordata$celltype == 'Slit2_IN',]

colnames(upcordata)

#简单画相关性
UP <- upcordata[,c(7,13,21:347)]

# 循环计算m1A和通路相关性
correlation_up <- data.frame()
##2.准备数据
data <- UP
##3.获取基因列表
genelist <- colnames(UP)[2:328]
##4.指定基因
gene <- "m1a_ucell"
genedata <- as.numeric(data[,gene])
pathdata <- UP[,-1]
##5.开始for循环
for(i in 1:length(genelist)){
  ## 1.指示
  print(i)
  ## 2.计算
  dd = cor.test(genedata,as.numeric(pathdata[,i]),
                method="spearman")
  ## 3.填充
  correlation_up[i,1] = gene
  correlation_up[i,2] = genelist[i]
  correlation_up[i,3] = dd$estimate
  correlation_up[i,4] = dd$p.value
}

colnames(correlation_up) <- c("m1a","pathway","cor","p.value")
## 6.p值矫正
### 校正后的P值会扩大，其更可靠
correlation_up$padjust = p.adjust(correlation_up$p.value,method = "BH")


colnames(UP)
# 单独提取作图数据
df <- UP[,c('m1a_ucell','condition','regulation of mRNA metabolic process_UCell')]

colnames(df) <- c('m1a_ucell','condition','pathway')
df <- df[df$m1a_ucell!=0,]
pmain <- ggplot(df, aes(x = m1a_ucell, y = pathway,color = condition
                          )) + #color = Species , color = condition
  geom_point(alpha = 0.6) +
  geom_smooth(method = 'lm',#线性回归
              formula = 'y ~ x',
              se=T,#添加置信区间，默认就是T
              lwd=1,#线条宽度
              color = "#9f0000", #拟合曲线颜色
              fill = "lightgrey")+#置信区间颜色
  stat_cor(method='spearman',
           label.x = 0.25, 
           # label.y = 100, 
           # label.sep = "\n",
           size=3) +
  scale_color_manual(values = groupcol)+
  xlab('m1A score')+
  ylab('Regulation of mRNA metabolic process') +
  theme_bw()+
  theme(axis.title = element_text(size = 12, colour = "black"),
        axis.text = element_text(size = 10, colour = "black"))
# color_palette("jco")
# Marginal densities along x axis
xdens <- axis_canvas(pmain, axis = "x") +
  geom_density(data = df, aes(x = m1a_ucell),fill = '#868686',
               alpha = 0.7, size = 0.2) 
  # fill_palette("jco")
# Marginal densities along y axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE) +
  geom_density(data = df, aes(x = pathway,fill = condition), #fill = '#0073c2',
               alpha = 0.4, size = 0.2) + scale_color_manual(values = groupcol)+
  coord_flip() 

p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
ggdraw(p2)

ggsave(
  file.path( "./figure/Figure_3/B_pos_m1A_pathway1_group.pdf"),
  plot = p2,
  height = 8,
  width = 12)


#### 相关性图2
df <- UP[,c('m1a_ucell', 'condition','neuron recognition_UCell')]
colnames(df) <- c('m1a_ucell','condition','pathway')
df <- df[df$m1a_ucell!=0,]
pmain <- ggplot(df, aes(x = m1a_ucell, y = pathway ,color = condition
)) + #color = Species
  geom_point(alpha = 0.6) +
  geom_smooth(method = 'lm',#线性回归
              formula = 'y ~ x',
              se=T, #添加置信区间，默认就是T
              lwd=1, #线条宽度
              color = "#9f0000", #拟合曲线颜色
              fill = "lightgrey")+ #置信区间颜色
  stat_cor(method='spearman',
           label.x = 0.25, 
           # label.y = 100, 
           # label.sep = "\n",
           size=3) +
  scale_color_manual(values = groupcol)+
  xlab('m1A score')+
  ylab('Neuron recognition') +
  theme_bw()+
  theme(axis.title = element_text(size = 12, colour = "black"),
        axis.text = element_text(size = 10, colour = "black"))
#color_palette("jco")
# Marginal densities along x axis
xdens <- axis_canvas(pmain, axis = "x") +
  geom_density(data = df, aes(x = m1a_ucell),fill = '#868686',
               alpha = 0.7, size = 0.2) 
# fill_palette("jco")
# Marginal densities along y axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE) +
  geom_density(data = df, aes(x = pathway,fill = condition),# fill = '#9479ad',
               alpha = 0.4, size = 0.2) + scale_color_manual(values = groupcol)+ 
  coord_flip() 

p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
ggdraw(p2)

ggsave(
  file.path( "./figure/Figure_3/B_neg_m1A_neuron_recognition_group.pdf"),
  plot = p2,
  height = 8,
  width = 12)




# 提取出每个组特有的功能，分别来做相关性(相关性结果不是很好)
{BiocManager::install('VennDetail')
library(VennDetail)
ven <- venndetail(Slit2_list_up)
detail(ven)

eachgroupf <- getSet(ven, subset = c("B_1dpi", "C_1wpi", "D_3wpi", "E_6wpi"))
B_1dpispec <- eachgroupf %>% 
  filter(Subset== 'B_1dpi') %>% 
  .[,'Detail']

table(pathtoname %in% B_1dpispec)
B_1dpispec1 <- pathtoname[pathtoname %in% B_1dpispec]
B_1dpispec1id <- names(B_1dpispec1)

# 再根据goid提取对应通路和基因
B_1sig <- pathtogen[B_1dpispec1id]
names(B_1sig) <- B_1dpispec1

#先打个下调的分看看
B_1dpispec1_up<- AddModuleScore_UCell(sub_seurat_obj1, 
                                         maxRank = 1500,
                                         ncores = 50,
                                         features = B_1sig, name = "_UCell")

B_1cordata <- B_1dpispec1_up@meta.data[]
B_1cordata <- B_1cordata[B_1cordata$celltype == 'Slit2_IN',]
colnames(B_1cordata)

#简单画相关性
library(corrplot)
M <- B_1cordata[,c(13,21:61)]
colnames(M)[2:42] <- names(B_1dpispec1)

cor<-cor(M)
corrplot(cor,
         method = 'ellipse',#图形形状,椭圆的大小和方向表示相关性强度和方向，可选'circle' (default), 'square', 'ellipse', 'number', 'pie', 'shade' and 'color'
         type='upper',#矩阵填充方式，可选'full'(defult),'lower','upper'
         order = 'AOE',#首先安装相关性的绝对值大小进行排序，然后调整变啊零的顺序，最后根据特征值的大小再次调整排序
         title='',#不设置标题
         tl.col = 'black',#文字标签颜色,默认'red'
         tl.cex=1,#文字标签大小
         addgrid.col = 'black'#网格颜色，默认'black'
)

# m1A和通路单独做相关性
# 循环计算m1A和通路相关性
correlation <- data.frame()
##2.准备数据
data <- M
##3.获取基因列表
genelist <- colnames(M)[2:42]
##4.指定基因
gene <- "m1a_ucell"
genedata <- as.numeric(data[,gene])
pathdata <- M[,-1]
##5.开始for循环
for(i in 1:length(genelist)){
  ## 1.指示
  print(i)
  ## 2.计算
  dd = cor.test(genedata,as.numeric(pathdata[,i]),
                method="spearman")
  ## 3.填充
  correlation[i,1] = gene
  correlation[i,2] = genelist[i]
  correlation[i,3] = dd$estimate
  correlation[i,4] = dd$p.value
}

colnames(correlation) <- c("m1a","pathway","cor","p.value")
## 6.p值矫正
### 校正后的P值会扩大，其更可靠
correlation$padjust = p.adjust(correlation$p.value,method = "BH")}

###### 把看起来有相关性的通路取出来基因，进一步做相关性
# 正相关通路
{ # regulation of mRNA metabolic process_UCell
uppathgene <- upsig[['regulation of mRNA metabolic process']]

# 然后提取上述基因的表达量以及m1a_ucell 分组信息等
df1 = FetchData(sub_seurat_obj1, vars=c("condition", "celltype", "m1a_ucell", uppathgene))
df1 <- df1[df1$celltype == 'Slit2_IN',]

df2 <- df1[,-c(1:3)]
usegene <- colSums(df2) > 0.5
df2 <- df2[,usegene]

df <- cbind(df1[,c(1:3)],df2)

# 循环计算m1A和基因相关性
correlation <- data.frame()
##2.准备数据
data <- df
##3.获取基因列表
genelist <- colnames(df)[4:279]
##4.指定基因
gene <- "m1a_ucell"
genedata <- as.numeric(df[,gene])
pathdata <- df[,-c(1:3)]
##5.开始for循环
for(i in 1:length(genelist)){
  ## 1.指示
  print(i)
  ## 2.计算
  dd = cor.test(genedata,as.numeric(pathdata[,i]),
                method="spearman")
  ## 3.填充
  correlation[i,1] = gene
  correlation[i,2] = genelist[i]
  correlation[i,3] = dd$estimate
  correlation[i,4] = dd$p.value
}

colnames(correlation) <- c("m1a","gene","cor","p.value")
## 6.p值矫正
### 校正后的P值会扩大，其更可靠
correlation$padjust = p.adjust(correlation$p.value,method = "BH")

# 分别提取前五和后五的数据来画图
top_bottom <- correlation %>%
  arrange(desc(cor)) %>%
  slice_head(n = 5) %>%
  bind_rows(correlation %>%
              arrange(cor) %>%
              slice_head(n = 5))

top_bottom$gene <- factor(top_bottom$gene,
                          levels = dput(top_bottom$gene)) 
top_bottom$cor <- round(top_bottom$cor,2)

ggplot(top_bottom,
       aes(x=gene,y=cor)) +
  geom_segment(aes(x=gene,xend=gene,y=0,yend=cor),  #绘制火柴的杆儿 
               linetype="solid", #线条类型 实线
               size=0.5, #线条粗细 
               color="black") + #线条颜色 灰色
  geom_hline( #添加水平线
    yintercept = 0.0,  #水平线位置
    linetype="dashed",  #线条类型 虚线
    size=0.5,  
    colour="black") +
  geom_point(aes(color=gene),  #绘制火柴头儿
             #color = col_list,
             size = 7) +   
  geom_text(aes(label = cor ), #添加文字标签
            color = "black", 
            size = 3) +
  scale_y_continuous( #设置y轴
    limits = c(-0.6,0.6),
    # breaks = c(-0.5,0,0.5,1.0),
    # labels = c(-0.5,0,0.5,1.0)
    ) +
  labs(
    x ="",  #设置x轴标题
    y="Spearman correlation coefficient",  #设置y轴标题
    title="Correlations between m1A score and gene expression",    #设置图片主标题
    subtitle = "(Regulation of mRNA metabolic process)"
    )+   #设置图片副标题
  theme_classic() +
  theme(
    plot.title = element_text(size = 10,hjust=0.5), #设置图片主标题的字体大小和设置居中位置
    plot.subtitle = element_text(size = 8,hjust = 0.5), #设置图片副标题的字体大小和设置居中位置
    axis.text.x = element_text(size=10,colour = 'black',angle = 45, hjust=1), #调整x轴刻度字体的大小和位置 
    axis.text.y = element_text(size=10,colour = 'black'),legend.position = "none",
    axis.title.y = element_text(size=10))
    #plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm"
                        #设置绘图的区域边距

ggsave(
  file.path( "./figure/Figure_3/C_correalation_m1a and metabolic pathway gene.pdf"),
  height = 5,
  width = 6)
}

# 负相关通路
{ # regulation of mRNA metabolic process_UCell
  negpathgene <- upsig[['neuron recognition']]
  
  # 然后提取上述基因的表达量以及m1a_ucell 分组信息等
  df1 = FetchData(sub_seurat_obj1, vars=c("condition", "celltype", "m1a_ucell", negpathgene))
  df1 <- df1[df1$celltype == 'Slit2_IN',]
  
  df2 <- df1[,-c(1:3)]
  usegene <- colSums(df2) > 0.5
  df2 <- df2[,usegene]
  
  df <- cbind(df1[,c(1:3)],df2)
  
  # 循环计算m1A和基因相关性
  correlation <- data.frame()
  ##2.准备数据
  data <- df
  ##3.获取基因列表
  genelist <- colnames(df)[4:52]
  ##4.指定基因
  gene <- "m1a_ucell"
  genedata <- as.numeric(df[,gene])
  pathdata <- df[,-c(1:3)]
  ##5.开始for循环
  for(i in 1:length(genelist)){
    ## 1.指示
    print(i)
    ## 2.计算
    dd = cor.test(genedata,as.numeric(pathdata[,i]),
                  method="spearman")
    ## 3.填充
    correlation[i,1] = gene
    correlation[i,2] = genelist[i]
    correlation[i,3] = dd$estimate
    correlation[i,4] = dd$p.value
  }
  
  colnames(correlation) <- c("m1a","gene","cor","p.value")
  ## 6.p值矫正
  ### 校正后的P值会扩大，其更可靠
  correlation$padjust = p.adjust(correlation$p.value,method = "BH")
  
  # 分别提取前五和后五的数据来画图
  top_bottom <- correlation %>%
    arrange(desc(cor)) %>%
    slice_head(n = 5) %>%
    bind_rows(correlation %>%
                arrange(cor) %>%
                slice_head(n = 5))
  
  top_bottom$gene <- factor(top_bottom$gene,
                            levels = dput(top_bottom$gene)) 
  top_bottom$cor <- round(top_bottom$cor,2)
  
  ggplot(top_bottom,
         aes(x=gene,y=cor)) +
    geom_segment(aes(x=gene,xend=gene,y=0,yend=cor),  #绘制火柴的杆儿 
                 linetype="solid", #线条类型 实线
                 size=0.5, #线条粗细 
                 color="black") + #线条颜色 灰色
    geom_hline( #添加水平线
      yintercept = 0.0,  #水平线位置
      linetype="dashed",  #线条类型 虚线
      size=0.5,  
      colour="black") +
    geom_point(aes(color=gene),  #绘制火柴头儿
               #color = col_list,
               size = 7) +   
    geom_text(aes(label = cor ), #添加文字标签
              color = "black", 
              size = 3) +
    # scale_y_continuous( #设置y轴
    #   limits = c(-0.6,0.6),
    #   # breaks = c(-0.5,0,0.5,1.0),
    #   # labels = c(-0.5,0,0.5,1.0)
    # ) +
    labs(
      x ="",  #设置x轴标题
      y="Spearman correlation coefficient",  #设置y轴标题
      title="Correlations between m1A score and gene expression",    #设置图片主标题
      subtitle = "(Neuron recognition)"
    )+   #设置图片副标题
    theme_classic() +
    theme(
      plot.title = element_text(size = 10,hjust=0.5), #设置图片主标题的字体大小和设置居中位置
      plot.subtitle = element_text(size = 8,hjust = 0.5), #设置图片副标题的字体大小和设置居中位置
      axis.text.x = element_text(size=10,colour = 'black',angle = 45, hjust=1), #调整x轴刻度字体的大小和位置 
      axis.text.y = element_text(size=10,colour = 'black'),legend.position = "none",
      axis.title.y = element_text(size=10))
  #plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm"
  #设置绘图的区域边距
  
  ggsave(
    file.path( "./figure/Figure_3/C_correalation_Neuron recognition.pdf"),
    height = 5,
    width = 6)
}