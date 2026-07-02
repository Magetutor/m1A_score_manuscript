# 设置工作目录
setwd('/home/bio/Projects/NC2022/neuro_sub/')
# 清空变量 加载包
# rm(list = ls())
library(Seurat)
# packageVersion('MC.seurat')
library(DESeq2)
library(tidyverse)
library(ggthemes)
library(ggrepel)

# load('./figure/Figure_4_supercell/data/subneuron_MC.seurat.Rda')
# 对应metacell进行打分
# features <- c("Alkbh1", "Alkbh3", "Fto", "Trmt10c", "Trmt6", "Trmt61a", 
#               "Ythdc1", "Ythdf1", "Ythdf2", "Ythdf3")
#signatures <- list()-
# signatures$m1a <- features

DefaultAssay(MC.seurat) <- 'RNA'

MC.seurat <- AddModuleScore_UCell(MC.seurat, features = signatures, name = "_ucell")

#先打个下调的分看看
subneuro_ucell_dn<- AddModuleScore_UCell(MC.seurat, 
                                         maxRank = 1500,
                                         ncores = 50,
                                         features = dnsig, name = "_UCell")

dncordata <- subneuro_ucell_dn@meta.data[]
dncordata <- dncordata[dncordata$celltype == 'Slit2_IN',]
colnames(dncordata)

#简单画相关性
library(corrplot)
M <- dncordata[,c(12:37)]
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
correlation_dn <- data.frame()
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
  correlation_dn[i,1] = gene
  correlation_dn[i,2] = genelist[i]
  correlation_dn[i,3] = dd$estimate
  correlation_dn[i,4] = dd$p.value
}

colnames(correlation_dn) <- c("m1a","pathway","cor","p.value")
## 6.p值矫正
### 校正后的P值会扩大，其更可靠
correlation_dn$padjust = p.adjust(correlation_dn$p.value,method = "BH")


####上调基因pathway分值######
subneuro_ucell_up<- AddModuleScore_UCell(MC.seurat, 
                                         maxRank = 1500,
                                         ncores = 50,
                                         features = upsig, name = "_UCell")

upcordata <- subneuro_ucell_up@meta.data
upcordata <- upcordata[upcordata$celltype == 'Slit2_IN',]

colnames(upcordata)

#简单画相关性
UP <- upcordata[,c(12:339)]

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



####画上下调通路正相关和负相关的相关性图###########
#####先画上调通路正相关前10的棒棒糖图
uptop10 <- correlation_up %>% 
  arrange(desc(cor)) %>%  # 按照相关性系数从大到小排序
  slice(1:10)   # 选择前10个观测
# 修改通路名称
uptop10$pathway <- gsub('_UCell','',uptop10$pathway) %>% 
  str_to_sentence()

p1 <- ggplot(uptop10,aes(cor,reorder(pathway,cor)))+
  geom_point(size = 5,color = "#d62e2d")+ # #3377a9
  geom_segment(aes(x=0,xend=cor,y=pathway,yend=pathway),  #绘制火柴的杆儿 
              linetype="solid", #线条类型 实线
              size=0.5, #线条粗细 
              color="#d62e2d") + #线条颜色 灰色
  scale_x_continuous(limits = c(0,0.5), #设置x轴范围
                     expand = expansion(mult = 0))+ #expansion函数的参数mult: 百分比间距，可以接受一个向量
  theme(panel.background = element_blank(), #删除背景
        # legend.position = "top",
        axis.title = element_text(size = 12, colour = "black"),
        axis.text = element_text(size = 10, colour = "black"),
        panel.grid = element_line (colour = "lightgrey"), #设置网格颜色
        panel.border = element_rect(fill = NA,colour = "black",size = 0.8))+ #设置边框
  ylab("")+ #删除y轴名
  xlab('Correalation between m1A score and pathways')

ggsave(
  file.path( "./figure/Figure_4_supercell/D_uppos_m1A_pathway.pdf"),
  plot = p1,
  height = 4,
  width = 9)

#####先画上调通路正相关前10的棒棒糖图
dntop10 <- correlation %>% 
  arrange(cor) %>%  # 按照相关性系数从大到小排序
  slice(1:10)   # 选择前10个观测
dntop10 <- dntop10[c(1:9),]
# 修改通路名称
dntop10$pathway <- gsub('_UCell','',dntop10$pathway) %>% 
  str_to_sentence()


ggplot(dntop10,aes(cor,reorder(pathway,cor)))+
  geom_point(size = 5,color = "#3377a9")+ # 
  geom_segment(aes(x=0,xend=cor,y=pathway,yend=pathway),  #绘制火柴的杆儿 
               linetype="solid", #线条类型 实线
               size=0.5, #线条粗细 
               color="#3377a9") + #线条颜色 灰色
  scale_x_continuous(limits = c(-0.2,0), #设置x轴范围
                     expand = expansion(mult = 0))+ #expansion函数的参数mult: 百分比间距，可以接受一个向量
  theme(panel.background = element_blank(), #删除背景
        # legend.position = "top",
        axis.title = element_text(size = 12, colour = "black"),
        axis.text = element_text(size = 10, colour = "black"),
        panel.grid = element_line (colour = "lightgrey"), #设置网格颜色
        panel.border = element_rect(fill = NA,colour = "black",size = 0.8))+ #设置边框
  ylab("")+ #删除y轴名
  xlab('Correalation between m1A score and pathways')

ggsave(
  file.path( "./figure/Figure_4_supercell/D_dnneg_m1A_pathway.pdf"),
  plot = last_plot(),
  height = 4,
  width = 7.5)

#### 单独提取代表性通路数据画图 ##########
#####上调通路正相关####
####### regulation of mRNA metabolic process_UCell

# neuron cell-cell adhesion_UCell
df <- upcordata[,c('m1a_ucell','condition','regulation of mRNA metabolic process_UCell')]
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
           label.x = 0.1, 
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
  file.path( "./figure/Figure_4_supercell/E_uppos_m1A_mRNA metabolic_group.pdf"),
  plot = p2,
  height = 5,
  width = 7)

# 通路里面基因相关性
{
  # regulation of mRNA metabolic process_UCell
  uppathgene <- upsig[['regulation of mRNA metabolic process']]
  
  # 然后提取上述基因的表达量以及m1a_ucell 分组信息等
  df1 = FetchData(MC.seurat, vars=c("condition", "celltype", "m1a_ucell", uppathgene))
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
  genelist <- colnames(df)[4:277]
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
    slice_head(n = 11) %>%
    bind_rows(correlation %>%
                arrange(cor) %>%
                slice_head(n = 5))
  top_bottom <- top_bottom[2:11,]
  top_bottom$gene <- factor(top_bottom$gene,
                            levels = dput(top_bottom$gene)) 
  top_bottom$cor <- round(top_bottom$cor,3)
  
  ggplot(top_bottom,
         aes(x=gene,y=cor)) +
    geom_segment(aes(x=gene,xend=gene,y=0,yend=cor),  #绘制火柴的杆儿 
                 linetype="solid", #线条类型 实线
                 size=0.8, #线条粗细 
                 color="black") + #线条颜色 灰色
    geom_hline( #添加水平线
      yintercept = 0.0,  #水平线位置
      linetype="dashed",  #线条类型 虚线
      size=0.5,  
      colour="black") +
    geom_point(aes(color=gene),  #绘制火柴头儿
               #color = col_list,
               size = 7) +   
    geom_text_repel(aes(label = cor ), #添加文字标签 # geom_text_repel(color="grey20",size=3,point.padding = NA)
              color = "black", nudge_y = 0.02,#direction = y,
              size = 3) +
    scale_y_continuous( #设置y轴
      limits = c(0,0.3),
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
      axis.title.y = element_text(size=10)) +
    scale_color_npg()
  #plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm"
  #设置绘图的区域边距
  
  ggsave(
    file.path( "./figure/Figure_4_supercell/E_uppos_m1A_mRNA metabolic_pathway gene.pdf"),
    height = 4,
    width = 5)
}


#####上调通路正相关####
####### neuron cell-cell adhesion_UCell ######
  df <- upcordata[,c('m1a_ucell','condition','neuron cell-cell adhesion_UCell')]
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
             label.x = 0.1, 
             # label.y = 100, 
             # label.sep = "\n",
             size=3) +
    scale_color_manual(values = groupcol)+
    xlab('m1A score')+
    ylab('Neuron cell-cell adhesion_UCell') +
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
    file.path( "./figure/Figure_4_supercell/E_uppos_m1A_neuron cell-cell adhesion_group.pdf"),
    plot = p2,
    height = 5,
    width = 7)
  
  # 通路里面基因相关性
  {
    # regulation of mRNA metabolic process_UCell
    uppathgene <- upsig[['neuron cell-cell adhesion']]
    
    # 然后提取上述基因的表达量以及m1a_ucell 分组信息等
    df1 = FetchData(MC.seurat, vars=c("condition", "celltype", "m1a_ucell", uppathgene))
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
    genelist <- colnames(df)[4:12]
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
      slice_head(n = 8)
    
    top_bottom$gene <- factor(top_bottom$gene,
                              levels = dput(top_bottom$gene)) 
    top_bottom$cor <- round(top_bottom$cor,3)
    
    ggplot(top_bottom,
           aes(x=gene,y=cor)) +
      geom_segment(aes(x=gene,xend=gene,y=0,yend=cor),  #绘制火柴的杆儿 
                   linetype="solid", #线条类型 实线
                   size=0.8, #线条粗细 
                   color="black") + #线条颜色 灰色
      geom_hline( #添加水平线
        yintercept = 0.0,  #水平线位置
        linetype="dashed",  #线条类型 虚线
        size=0.5,  
        colour="black") +
      geom_point(aes(color=gene),  #绘制火柴头儿
                 #color = col_list,
                 size = 7) +   
      geom_text_repel(aes(label = cor ), #添加文字标签 # geom_text_repel(color="grey20",size=3,point.padding = NA)
                      color = "black", nudge_y = 0.011,#direction = y,
                      size = 3) +
      scale_y_continuous( #设置y轴
        limits = c(0,0.25),
        # breaks = c(-0.5,0,0.5,1.0),
        # labels = c(-0.5,0,0.5,1.0)
      ) +
      labs(
        x ="",  #设置x轴标题
        y="Spearman correlation coefficient",  #设置y轴标题
        title="Correlations between m1A score and gene expression",    #设置图片主标题
        subtitle = "(Neuron cell-cell adhesion)"
      )+   #设置图片副标题
      theme_classic() +
      theme(
        plot.title = element_text(size = 10,hjust=0.5), #设置图片主标题的字体大小和设置居中位置
        plot.subtitle = element_text(size = 8,hjust = 0.5), #设置图片副标题的字体大小和设置居中位置
        axis.text.x = element_text(size=10,colour = 'black',angle = 45, hjust=1), #调整x轴刻度字体的大小和位置 
        axis.text.y = element_text(size=10,colour = 'black'),legend.position = "none",
        axis.title.y = element_text(size=10)) +
      scale_color_npg()
    #plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm"
    #设置绘图的区域边距
    
    ggsave(
      file.path( "./figure/Figure_4_supercell/E_uppos_m1A_neuron cell-cell adhesion_pathway gene.pdf"),
      height = 4,
      width = 5)
  }

  #####下调通路负相关####
  ####### aerobic respiration_UCell ######
  df <- dncordata[,c('m1a_ucell','condition','aerobic respiration_UCell')]
  colnames(df) <- c('m1a_ucell','condition','pathway')
  df <- df[df$m1a_ucell!=0,]
  pmain <- ggplot(df, aes(x = m1a_ucell, y = pathway,color = condition
  )) + #color = Species , color = condition
    geom_point(alpha = 0.6) +
    geom_smooth(method = 'lm',#线性回归
                formula = 'y ~ x',
                se=T,#添加置信区间，默认就是T
                lwd=1,#线条宽度
                color = '#0073c2', #拟合曲线颜色
                fill = "lightgrey")+#置信区间颜色
    stat_cor(method='spearman',
             label.x = 0.1, 
             # label.y = 100, 
             # label.sep = "\n",
             size=3) +
    scale_color_manual(values = groupcol)+
    xlab('m1A score')+
    ylab('aerobic respiration') +
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
    file.path( "./figure/Figure_4_supercell/E_dnneg_m1A_aerobic respiration_UCell_group.pdf"),
    plot = p2,
    height = 5,
    width = 7)
  
  # 通路里面基因相关性
  {
    # regulation of mRNA metabolic process_UCell
    dnpathgene <- dnsig[['aerobic respiration']]
    
    # 然后提取上述基因的表达量以及m1a_ucell 分组信息等
    df1 = FetchData(MC.seurat, vars=c("condition", "celltype", "m1a_ucell", dnpathgene))
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
    genelist <- colnames(df)[4:157]
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
    
    # 分别提取后10个数据来画图
    top_bottom <- correlation %>%
      arrange(cor) %>%
      slice_head(n = 10)
    
    top_bottom$gene <- factor(top_bottom$gene,
                              levels = dput(top_bottom$gene)) 
    top_bottom$cor <- round(top_bottom$cor,3)
    
    ggplot(top_bottom,
           aes(x=gene,y=cor)) +
      geom_segment(aes(x=gene,xend=gene,y=0,yend=cor),  #绘制火柴的杆儿 
                   linetype="solid", #线条类型 实线
                   size=0.8, #线条粗细 
                   color="black") + #线条颜色 灰色
      geom_hline( #添加水平线
        yintercept = 0.0,  #水平线位置
        linetype="dashed",  #线条类型 虚线
        size=0.5,  
        colour="black") +
      geom_point(aes(color=gene),  #绘制火柴头儿
                 #color = col_list,
                 size = 7) +   
      geom_text_repel(aes(label = cor ), #添加文字标签 # geom_text_repel(color="grey20",size=3,point.padding = NA)
                      color = "black", nudge_y = -0.01,#direction = y,
                      size = 3) +
      scale_y_continuous( #设置y轴
        limits = c(-0.2,0),
        # breaks = c(-0.5,0,0.5,1.0),
        # labels = c(-0.5,0,0.5,1.0)
      ) +
      labs(
        x ="",  #设置x轴标题
        y="Spearman correlation coefficient",  #设置y轴标题
        title="Correlations between m1A score and gene expression",    #设置图片主标题
        subtitle = "(Aerobic respiration)"
      )+   #设置图片副标题
      theme_classic() +
      theme(
        plot.title = element_text(size = 10,hjust=0.5), #设置图片主标题的字体大小和设置居中位置
        plot.subtitle = element_text(size = 8,hjust = 0.5), #设置图片副标题的字体大小和设置居中位置
        axis.text.x = element_text(size=10,colour = 'black',angle = 45, hjust=1), #调整x轴刻度字体的大小和位置 
        axis.text.y = element_text(size=10,colour = 'black'),legend.position = "none",
        axis.title.y = element_text(size=10)) +
      scale_color_d3('category20c')
    #plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm"
    #设置绘图的区域边距
    
    ggsave(
      file.path( "./figure/Figure_4_supercell/E_dnneg_m1A_aerobic respiration_pathway gene.pdf"),
      height = 4,
      width = 6)
  }

  #####下调通路负相关####
  ####### electron transport chain_UCell ######
  df <- dncordata[,c('m1a_ucell','condition','electron transport chain_UCell')]
  colnames(df) <- c('m1a_ucell','condition','pathway')
  df <- df[df$m1a_ucell!=0,]
  pmain <- ggplot(df, aes(x = m1a_ucell, y = pathway,color = condition
  )) + #color = Species , color = condition
    geom_point(alpha = 0.6) +
    geom_smooth(method = 'lm',#线性回归
                formula = 'y ~ x',
                se=T,#添加置信区间，默认就是T
                lwd=1,#线条宽度
                color = '#0073c2', #拟合曲线颜色
                fill = "lightgrey")+#置信区间颜色
    stat_cor(method='spearman',
             label.x = 0.1, 
             # label.y = 100, 
             # label.sep = "\n",
             size=3) +
    scale_color_manual(values = groupcol)+
    xlab('m1A score')+
    ylab('Electron transport chain') +
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
    file.path( "./figure/Figure_4_supercell/E_dnneg_m1A_electron transport chain_group.pdf"),
    plot = p2,
    height = 5,
    width = 7)
  
  # 通路里面基因相关性
  {
    # regulation of mRNA metabolic process_UCell
    dnpathgene <- dnsig[['electron transport chain']]
    
    # 然后提取上述基因的表达量以及m1a_ucell 分组信息等
    df1 = FetchData(MC.seurat, vars=c("condition", "celltype", "m1a_ucell", dnpathgene))
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
    genelist <- colnames(df)[4:89]
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
    
    # 分别提取后10个数据来画图
    top_bottom <- correlation %>%
      arrange(cor) %>%
      slice_head(n = 10)
    
    top_bottom$gene <- factor(top_bottom$gene,
                              levels = dput(top_bottom$gene)) 
    top_bottom$cor <- round(top_bottom$cor,3)
    
    ggplot(top_bottom,
           aes(x=gene,y=cor)) +
      geom_segment(aes(x=gene,xend=gene,y=0,yend=cor),  #绘制火柴的杆儿 
                   linetype="solid", #线条类型 实线
                   size=0.8, #线条粗细 
                   color="black") + #线条颜色 灰色
      geom_hline( #添加水平线
        yintercept = 0.0,  #水平线位置
        linetype="dashed",  #线条类型 虚线
        size=0.5,  
        colour="black") +
      geom_point(aes(color=gene),  #绘制火柴头儿
                 #color = col_list,
                 size = 7) +   
      geom_text_repel(aes(label = cor ), #添加文字标签 # geom_text_repel(color="grey20",size=3,point.padding = NA)
                      color = "black", nudge_y = -0.01,#direction = y,
                      size = 3) +
      scale_y_continuous( #设置y轴
        limits = c(-0.2,0),
        # breaks = c(-0.5,0,0.5,1.0),
        # labels = c(-0.5,0,0.5,1.0)
      ) +
      labs(
        x ="",  #设置x轴标题
        y="Spearman correlation coefficient",  #设置y轴标题
        title="Correlations between m1A score and gene expression",    #设置图片主标题
        subtitle = "(electron transport chain)"
      )+   #设置图片副标题
      theme_classic() +
      theme(
        plot.title = element_text(size = 10,hjust=0.5), #设置图片主标题的字体大小和设置居中位置
        plot.subtitle = element_text(size = 8,hjust = 0.5), #设置图片副标题的字体大小和设置居中位置
        axis.text.x = element_text(size=10,colour = 'black',angle = 45, hjust=1), #调整x轴刻度字体的大小和位置 
        axis.text.y = element_text(size=10,colour = 'black'),legend.position = "none",
        axis.title.y = element_text(size=10)) +
      scale_color_d3('category20c')
    #plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm"
    #设置绘图的区域边距
    
    ggsave(
      file.path( "./figure/Figure_4_supercell/E_dnneg_m1A_electron transport chain_pathway gene.pdf"),
      height = 4,
      width = 6)
  }
  
  
# 需要保存的数据
save(subneuro_ucell_dn,correlation_dn,
     dncordata,subneuro_ucell_up,
     correlation_up,upcordata,shared_up,shared_dn,upsig,dnsig,
     file = './figure/Figure_4_supercell/data/all_data_used_in_metacell_function.Rda')
