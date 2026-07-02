#### 加载包
library(Seurat)
library(tidyverse)
library(patchwork)
library(AUCell)

# 加载数据
seu <- qs::qread("output/03-2.myeloid.seurat.aucell.qs")

# 提取m1A score和aucell的值
DefaultAssay(seu) <- "AUCell"
dfau <- seu@assays$AUCell@data
dfaut <- t(dfau)

dfm1A <- seu@meta.data[,'m1a_ucell',drop = F]

# 排序一下
dfaut <- dfaut[rownames(dfm1A),]

# cbind合并
dfuse <- cbind(dfm1A,dfaut)

# 循环计算m1A和转录因子相关性
correlation <- data.frame()
##2.准备数据
data <- dfuse
##3.获取基因列表
genelist <- colnames(data)
##4.指定基因
gene <- "m1a_ucell"
genedata <- as.numeric(data[,gene])
##5.开始for循环
for(i in 1:length(genelist)){
  ## 1.指示
  print(i)
  ## 2.计算
  dd = cor.test(genedata,as.numeric(data[,i]),method="spearman")
  ## 3.填充
  correlation[i,1] = gene
  correlation[i,2] = genelist[i]
  correlation[i,3] = dd$estimate
  correlation[i,4] = dd$p.value
}

colnames(correlation) <- c("gene1","gene2","cor","p.value")
## 6.p值矫正
### 校正后的P值会扩大，其更可靠
correlation$padjust = p.adjust(correlation$p.value,method = "BH")

## 7.筛选p值小于0.05，按照相关性系数绝对值选前500个的基因， 数量可以自己定
library(dplyr)
library(tidyr)
cor_data_sig <- correlation %>% 
  filter(padjust < 0.05 & abs(cor) >= 0.25) 
## 横向柱状图展示
#添加上下调分组标签：
dt <- cor_data_sig[-1,]
dt$group <- case_when(dt$cor > 0 ~ 'Positive correlation',
                      dt$cor < 0 ~ 'Negative correlation')
dt <- arrange(dt,desc(cor),)
levels <- dt$gene2
dt$gene2 <- factor(dt$gene2,levels = rev(levels))
p <- ggplot(dt,
            aes(x =cor, y = gene2, fill = group)) + #数据映射 # aes(x = reorder(source, score), y = score))
  geom_col() + #绘制添加条形图
  theme_bw()
mytheme <- theme(
  legend.position = 'none',
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  axis.line.x = element_line(color = 'black',size = 1.1),
  axis.text = element_text(size = 12)
)
p1 <- p + mytheme
  # scale_x_break(c(0.1,0.4),scales =1,space = 0.1)+
  # scale_x_break(c(-0.1,-0.4),scales =1,space = 0.1)
p1

#先根据上下调标签拆分数据框：
up <- dt[which(dt$cor > 0),]
down <- dt[which(dt$cor < 0),]
#添加上调pathway标签：
p2 <- p1 +
  geom_text(data = up,
            aes(x = -0.01, y = gene2, label = gene2),
            size = 5,
            hjust = 1) #标签右对齐
p2
#添加下调pathway标签：
p3 <- p2 +
  geom_text(data = down,
            aes(x = 0.01, y = gene2, label = gene2),
            size = 5,
            hjust = 0) #标签左对齐
p3
#继续调整细节：
p4 <- p3 +
  labs(x = 'Correlation coefficient', y = ' ', title = 'Correlation analysis') + #修改x/y轴标签、标题添加
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14)) #主标题居中、字号调整
p4
p5 <- p4 +
  scale_fill_manual(values = c("#4978bc", "#ea604c"))
p5
p6 <- p5 +
  geom_text(x = -0.705, y = 6.5, label = "Positive correlation", size = 6, color = "#ea604c") +
  geom_text(x = 0.72, y = 1.5, label = "Negative correlation", size = 6, color = "#4978bc")
p6

ggsave('./figures/Figure_2/Figure2G_m1a_tfs_cor.pdf',width = 8,height =7 )

# 做单独的相关性图
c("Bclaf1(+)", "Yy1(+)", "Atf4(+)", "Pura(+)", "Etv6(+)", "Nr3c1(+)", 
  "Sp3(+)", "Zbtb7a(+)", "Chd1(+)", "Nfe2l1(+)", "Ep300(+)", "Phf20(+)", 
  "Bach1(+)", "Kmt2b(+)", "Klf4(+)", "Klf3(+)", "Zfp131(+)", "Elf1(+)", 
  "Egr3(+)", "Klf13(+)", "Nfic(+)", "Churc1(+)", "Fli1(+)", "E2f1(+)", 
  "Usf2(+)", "Hes6(+)")

# 需要去除(+)再来做相关性
dfuse01 <- dfuse
colnames(dfuse01) <- gsub("\\(\\+\\)", "", colnames(dfuse01))

ggscatterhist(
  dfuse01, x = "m1a_ucell", y = "Bclaf1",
  color = "#00AFBB", 
  size = 3, alpha = 0.6,
  palette = c("#00AFBB", "#E7B800", "#FC4E07"),
  margin.params = list(color = "#00AFBB", size = 0.2)
)

library(ggpointdensity)
library(ggplot2)
library(viridis)
library(ggrastr)


p1 <- ggplot(data = dfuse01, mapping = aes(x = m1a_ucell,
                                 y = Bclaf1)) + 
  # stat_density_2d( linemitre = 20)+
   geom_pointdensity() + #密度散点图（geom_pointdensity）
   #geom_point_rast()+
  scale_color_viridis()  +
  geom_smooth(method = lm, se = TRUE) +  ##省略拟合曲线
  stat_cor(method = "spearman")+
  theme_bw()+
  theme(axis.title = element_text(size = 12,
                                    face = "bold", 
                                    vjust = 0.5, 
                                    hjust = 0.5))+
  theme(axis.text = element_text(size = 10,
                                    color = "black", 
                                    vjust = 0.5, 
                                    hjust = 0.5))+
  labs(title = 'Correlation coefficient') +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
    panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        #text=element_text(size=12,  family="serif")
        ) +
  theme(legend.position='none')  ##去除legend

ggsave('./figures/Figure_2/Figure2H_m1a_Bclaf1_cor.pdf',width = 5,height =5 )

# 定义为一个函数

midu <- function(i){
  ggplot(data = dfuse01, mapping = aes(x = m1a_ucell,
                                       y = {{i}})) + 
    # stat_density_2d( linemitre = 20)+
    geom_pointdensity() + #密度散点图（geom_pointdensity）
    #geom_point_rast()+
    scale_color_viridis()  +
    geom_smooth(method = lm, se = TRUE) +  ##省略拟合曲线
    stat_cor(method = "spearman")+
    theme_bw()+
    theme(axis.title = element_text(size = 12,
                                    face = "bold", 
                                    vjust = 0.5, 
                                    hjust = 0.5))+
    theme(axis.text = element_text(size = 10,
                                   color = "black", 
                                   vjust = 0.5, 
                                   hjust = 0.5))+
    labs(title = 'Correlation coefficient') +
    theme(plot.title = element_text(hjust = 0.5, size = 14),
          panel.grid.major=element_line(colour=NA),
          panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_rect(fill = "transparent",colour = NA),
          panel.grid.minor = element_blank(),
          #text=element_text(size=12,  family="serif")
    ) +
    theme(legend.position='none')
}


p2 <- midu(Yy1)
p3 <- midu(Atf4)
p4 <- midu(Pura)
p5 <- midu(Etv6)
p6 <- midu(Nr3c1)

library(cowplot)
plot_grid(p1, p2, p3, p4, p5, p6,
           align = "h",ncol = 3)
ggsave('./figures/Figure_2/Figure2I_m1a_6tf_cor.pdf',width = 10,height =7 )


### 保存数据
qs::qsave(seu, "output/03-2.myeloid.seurat.aucell.qs")
