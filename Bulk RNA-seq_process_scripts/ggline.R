#数据计算好，可以直接用
library(ggpubr)
diffplot <- function(gene){
  ggline(
    spinal_Exp_vst_mean, x = "group", y = gene, color = "#E64B357F")
}

genelist <- m7G_m[m7G_m %in% colnames(spinal_Exp_vst)]

plist <- list()

for(i in genelist) {
  
  p <- diffplot(i)
  
  plist[[i]] <- p
  
}

library("patchwork")
m7G_plot <- wrap_plots(plist,ncol = 3)
ggsave(m7G_plot, file="output/sci/ggline/m7G_plot.pdf", width=23, height=23)

