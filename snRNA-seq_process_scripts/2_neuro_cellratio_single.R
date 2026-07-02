# table(scobj$orig.ident)#查看各组细胞数
# prop.table(table(Idents(scobj)))
# table(Idents(scobj), scobj$orig.ident)#各组不同细胞群细胞数
Cellratio <- prop.table(table(sub_seurat_obj1$celltype, sub_seurat_obj1$orig.ident), margin = 2)#计算各组样本不同细胞群比例
Cellratio <- data.frame(Cellratio)
library(reshape2)
cellper <- dcast(Cellratio,Var2~Var1, value.var = "Freq")#长数据转为宽数据
rownames(cellper) <- cellper[,1]
cellper <- cellper[,-1]


###添加分组信息
sample <- c("A-45", "A-46", "A-48", "B-49", "B-50", "B-51", "C-41", "C-42", 
            "C-44", "D-25", "D-26", "D-29", "F-12", "F-13", "F-14")
group <- c("A_Uninj","A_Uninj","A_Uninj","B_1dpi","B_1dpi","B_1dpi",
           "C_1wpi","C_1wpi","C_1wpi","D_3wpi","D_3wpi","D_3wpi",
           "E_6wpi","E_6wpi","E_6wpi")
samples <- data.frame(sample, group)#创建数据框

rownames(samples)=samples$sample
cellper$sample <- samples[rownames(cellper),'sample']#R添加列
cellper$group <- samples[rownames(cellper),'group']#R添加列

###作图展示
pplist = list()
scobj_groups = c("A_Uninj", "B_1dpi", "C_1wpi", 
                 "D_3wpi", "E_6wpi")

library(ggplot2)
library(dplyr)
library(ggpubr)
for(group_ in scobj_groups){
  cellper_  = cellper %>% dplyr::filter(group == group_)#选择一组数据
  cellper_ <- pivot_longer(cellper_, 
               cols = c(1:7),  # 指定需要转换的列范围
               names_to = "celltype",  # 新列的名称，保存原始列的名称
               values_to = "Value")    # 新列的值，保存原始列的值
  
  colnames(cellper_) = c('sample','group',"celltype",'percent')#对选择数据列命名
  cellper_$percent = as.numeric(cellper_$percent)#数值型数据
  cellper_ <- cellper_ %>% group_by(celltype) %>% mutate(upper =  quantile(percent, 0.75), 
                                                      lower = quantile(percent, 0.25),
                                                      mean = mean(percent),
                                                      median = median(percent))#上下分位数
  print(group_)
  cellper_$celltype <- factor(cellper_$celltype,
                              levels = c("Slit2_IN", "Npy_IN", "Cck_EN", "Sox5_EN", "Pde11a_EN", "Tac2_EN", 
                                         "Gal_IN"))
  print(cellper_$median)
  
  pp1 = ggplot(cellper_,aes(x=celltype,y=percent)) + #ggplot作图
    geom_jitter(shape = 21,size=4,aes(fill=celltype),width = 0.25) + 
    stat_summary(fun=mean, geom="point", color="black") +
    theme_cowplot() +
    theme(axis.text = element_text(size = 10),
          axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75),
          axis.title = element_text(size = 12),legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),plot.title = element_text(size = 12,face = 'plain'),legend.position = 'none') + 
    labs(title = group_ ,y='Percentage') +
    geom_errorbar(aes(ymin = lower, ymax = upper),col = "black",width = 0.25) +
    scale_fill_manual(values = col) +
    xlab('')
  
  ###组间t检验分析
  labely = max(cellper_$percent)
  compare_means(percent ~ celltype,  data = cellper_)
  my_comparisons <- list( c("Npy_IN", "Slit2_IN"),
                          c("Cck_EN", "Slit2_IN"),
                          c("Sox5_EN", "Slit2_IN"),
                          c("Pde11a_EN", "Slit2_IN"),
                          c("Tac2_EN", "Slit2_IN"),
                          c("Gal_IN", "Slit2_IN"))
  
  pp1 = pp1 + stat_compare_means(comparisons = my_comparisons, label = "p.signif",
                                                     size = 5,method = "t.test")
  pplist[[group_]] = pp1
}

scobj_groups = c("A_Uninj", "B_1dpi", "C_1wpi", 
                 "D_3wpi", "E_6wpi")
library(cowplot)
plot_grid(pplist[['A_Uninj']],
          pplist[['B_1dpi']],
          pplist[['C_1wpi']],
          pplist[['D_3wpi']],
          pplist[['E_6wpi']],ncol = 3)

ggsave(
  file.path( "./figure/Figure_3/各个时间点中的细胞比例.pdf"),
  height = 10,
  width = 12)
