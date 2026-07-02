# 绘制折线图
p <- ggplot(spinal_Exp_vst, aes(x = group, y = "Ncbp3")) +
  geom_line() + 
  geom_smooth()
p
ggline(
  spinal_Exp_vst_mean, x = "group", y = "Ncbp3", color = "#E64B357F") + 
  geom_smooth(method =  "glm")

#小提琴图
library(dplyr)

gene <- spinal_Exp_vst %>% 
  filter(group %in% c("Sham","30min")) %>%
  select(c(group,Ncbp3))

Ncbp3_plot <- ggviolin(gene, x = "group", y = "Ncbp3", 
         color = "group",
         fill="group",
         palette =c("#E64B357F","#4DBBD57F","#00A0877F","#3C54887F","#F39B7F7F","#8491B47F","#91D1C27F","#DC00007F"),#设置颜色
         add = c("boxplot","jitter"),#添加箱线图
         add.params = list(color="white"),#设置箱线图边颜色
         xlab = F, #不显示x轴的标签
         legend = "right"#图例显示在右侧
)+   
  stat_compare_means(method = "t.test",
                     label = "p.signif",size = 7,label.x =1.5,label.y = 11.7, ref.group = "Sham")

ggsave(Ncbp3_plot, file="output/sci/ggviolin/m7G/Ncbp3_plot.pdf", width=6, height=8)




