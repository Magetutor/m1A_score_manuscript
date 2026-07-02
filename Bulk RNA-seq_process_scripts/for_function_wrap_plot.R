group <- rep(c("Control", "Sham", "0min", "30min","1h","6h","12h","1d","3d","5d","7d","14d","1m","2m","3m"),each = 4)
spinal_Exp_vst <- cbind(group,spinal_Exp_vst)
spinal_Exp_vst$group <- factor(spinal_Exp_vst$group, levels = c("Sham", "Control", "0min", "30min","1h","6h","12h","1d","3d","5d","7d","14d","1m","2m","3m"))
diffplot <- function(gene){
  my_comparisons <- list(
    c("Control", "Sham"),
    c("0min", "Sham"), 
    c("30min", "Sham"), 
    c("1h", "Sham"),
    c("6h", "Sham"), 
    c("12h", "Sham"), 
    c("1d", "Sham"), 
    c("3d", "Sham"),
    c("5d", "Sham"),
    c("7d", "Sham"), 
    c("14d", "Sham"), 
    c("1m", "Sham"), 
    c("2m", "Sham"),
    c("3m", "Sham")
  )
  library(ggpubr)
  ggboxplot(
    spinal_Exp_vst, x = "group", y = gene,
    color = "group", palette = c("#E64B357F","#4DBBD57F","#00A0877F","#3C54887F","#F39B7F7F","#8491B47F","#91D1C27F","#DC00007F","#7E61487F","#D31175","#796BAF","#ED7419","#10863B","#E64B323F","#7E616F"),
    add = "jitter"
  )+
    stat_compare_means(comparisons = my_comparisons, method = "t.test")
}
library(cowplot)
genelist <- c("FTO","TRMT6","WTAP","YTHDF2","YTHDF3","YTHDC1")

m5c <- c("NSUN6","NSUN7","DNMT1","DNMT2","DNMT3A",
         "DNMT3B","NSUN4","NSUN2","NSUN3","NOP2",
         "NSUN5","TRDMT1","ALYREF" ,"YBX1","TET2",
         "TET3","ALKBH1")
m7G <- c("METTL1","WDR4","NSUN2","DCP2","DCPS","NUDT10",
         "NUDT11","NUDT16","NUDT3","NUDT4","NUDT4B",
         "AGO2","CYFIP1","EIF4E","EIF4E1B","EIF4E2",
         "EIF4E3","GEMIN5","LARP1","NCBP1","NCBP2",
         "NCBP3","EIF3D","EIF4A1","EIF4G3","IFIT5",
         "LSM1","NCBP2L","SNUPN")
m1a <- c("TRMT10C","TRMT61B","TRMT6","TRMT61A","RRP8",
         "YTHDF1","YTHDF2","YTHDF3","YTHDC1","FTO",
         "ALKBH1","ALKBH3")
m6a <- c("METTL3","METTL14","WTAP","RBM15","RBM15B","VIRMA",
         "KIAA1429","CBLL1","ZC3H13","ALKBH5","FTO","YTHDC1",
         "YTHDC2","YTHDF1","YTHDF2","YTHDF3","IGF2BP1","IGF2BP2",
         "FMR1","LRPPRC","ELAVL1")

genelist <- m1a_m[m1a_m %in% colnames(spinal_Exp_vst)]


plist <- list()

for(i in genelist) {
  
  p <- diffplot(i)
  
  plist[[i]] <- p
  
}

m1a_plot <- wrap_plots(plist,ncol = 3)
ggsave(m1a_plot, file="output/sci/m1a_plot.pdf", width=23, height=20)
