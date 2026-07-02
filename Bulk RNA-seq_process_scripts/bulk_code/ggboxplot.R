## exprSet
metadata <- data.table::fread("data/STTT/sci_bulkrna_samples.csv",data.table = F)
group <- rep(c("SCI_2d","sham"),each = 3)
exprSet <- spinal_Exp_vst
exprSet <- t(exprSet)
exprSet <- as.data.frame(exprSet)
identical(metadata$sample, rownames(exprSet))
exprSet <- cbind("group" = metadata$group,exprSet)
exprSet <- exprSet %>%
  filter(group %in% c("Sham","SCI_1d","SCI_3d","SCI_7d"))
exprSet$group <- ifelse(exprSet$group %in% "SCI_1d","1 dpi",
                        ifelse(exprSet$group %in% "SCI_3d","3 dpi",
                               ifelse(exprSet$group %in% "SCI_7d","7 dpi","Sham")))

exprSet <- exprSet %>%
  separate(group,c("group","time"),sep = "_")

exprSet$group <- factor(exprSet$group,levels = c("Control","Sham","SCI_0min",
                                                 "SCI_30min","SCI_1h","SCI_6h",
                                                 "SCI_12h","SCI_1d","SCI_3d",
                                                 "SCI_5d","SCI_7d","SCI_14d",
                                                 "SCI_1m","SCI_2m","SCI_3m"))

exprSet$group <- factor(exprSet$group,levels = c("Sham","1 dpi","3 dpi","7 dpi"))
exprSet$time <- factor(exprSet$time,levels = c("0h","2h","4h",
                                               "6h","8h","12h",
                                               "24h","48h"))






diffplot <- function(gene){
ggboxplot(exprSet, x = "group", y = gene, 
         color = "group",
         #palette =c("#E64B357F","#4DBBD57F","#00A0877F","#3C54887F","#F39B7F7F","#8491B47F","#91D1C27F","#DC00007F"),#设置颜色
         add = "jitter",#添加箱线图
         add.params = list(color="black"),#设置箱线图边颜色
         xlab = F, #不显示x轴的标签
         legend = "right"#图例显示在右侧
)+   
  stat_compare_means(method = "t.test",
                     label = "p.signif",size = 7,label.x =1.5,label.y = max(exprSet[,gene])+1, ref.group = "Sham")+
  ggtitle("ourdata")+ 
  scale_y_continuous(limits = c(0, max(exprSet[,gene])+1))+
  theme(plot.title = element_text(size = 20, face = "bold",hjust=0.5),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 16))

}
diffplot("Trmt10c")
m1a <- unlist(markergenes_mice$gene)
genelist <- m1a[m1a %in% colnames(exprSet)]

plist <- list()

for(i in genelist) {
  
  p <- diffplot(i)
  
  plist[[i]] <- p
  
}

plist <- list(p1,p2,p3,p4,p5,p6,p7,p8,p9)
library(patchwork)
m1a_plot <- wrap_plots(p1,p2,p3,p4,p5,p6,p7,p8,p9)
m1a_plot <- plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,p9)

library(export)
## 导成PPT可编辑的格式
graph2ppt(file="output/发表级别图-初版/m1a/10.pptx")


gene <- exprSet[,c("group","time","Jam3")]
library(rstatix)
t <- t_test(group_by(gene, time),  Jam3  ~ group)
tj <- adjust_pvalue(t, method = 'fdr') #p值矫正；
tj
#根据p.adj添加显著性标记符号；
tj <- add_significance(tj, 'p.adj')
tj
lab <- add_xy_position(tj, x = 'group', dodge = 0.65)

ggboxplot(gene, x = "group", y = "Jam3", 
         color = "group",
         #palette =c("#E64B357F","#4DBBD57F","#00A0877F","#3C54887F","#F39B7F7F","#8491B47F","#91D1C27F","#DC00007F"),#设置颜色
         add = "jitter",#添加箱线图
         add.params = list(color="white"),#设置箱线图边颜色
         xlab = F, #不显示x轴的标签
         legend = "right"#图例显示在右侧
) + scale_y_continuous(limits = c(0, max(gene$Jam3)+1))+
  facet_grid(.~time)+
  stat_pvalue_manual(lab, label = 'p.adj.signif', label.size=4, bracket.size=0.5, tip.length = 0.02)
ggsave(file = "ourdata_ssgsea.pdf",width = 18,height = 11)

library(ggplot2)
my_comparisons <- list(
  c("1 dpi", "Sham"),
  c("3 dpi","Sham"),
  c("7 dpi","Sham")
)
#自定义图表主题，对图表主题做精细调整；
library(ggpubr)
library(ggbreak)
p10 <- ggplot(exprSet,aes(x=group,y=Trmt10c))+
  geom_boxplot(aes(fill=group))+
  geom_point()+
  scale_y_continuous(limits = c(0, max(exprSet$Trmt10c)+0.5),breaks = seq(0,max(exprSet$Trmt10c)+0.5, 0.5))+
  scale_y_break(c(0.5,min(exprSet$Trmt10c)),scales =20,space = 0.3)+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",
                      label = "p.format",size = 5)+
  scale_fill_manual(values = c("#4b81b4", "#27a12c","#fda711","#ff4911"))+
  theme_classic()+
  theme(axis.title.y.right = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank())+
  labs(x="",y= "Trmt10c",fill= "time") +
  theme(axis.title.y = element_text(size = 16,colour = "black",face = "italic"),
        axis.text = element_text(size = 15,colour = "black"),
        axis.line.y.right = element_blank(),
        legend.key.size = unit(22, "pt"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15))

ggsave(file = "output/发表级别图-初版/Trmt10c/基因与ssGSEA——boxplot/m1A_Score_boxplot.pdf",width = 6.5,height = 6.8)




