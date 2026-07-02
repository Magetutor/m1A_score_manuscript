### 1.载入数据
library(readxl)
library(tidyr)
library(dplyr)
blood_sci_fpkm <- read_excel("data/GSE226238_Morrison_et_al_processed_data.xlsx")
blood_sci_fpkm <- as.data.frame(blood_sci_fpkm)  
rownames(blood_sci_fpkm) <- blood_sci_fpkm[,1]
blood_sci_fpkm <- t(blood_sci_fpkm)
blood_sci_fpkm <- as.data.frame(blood_sci_fpkm)
blood_sci_fpkm <- cbind(gene = rownames(blood_sci_fpkm),blood_sci_fpkm)
blood_sci_fpkm <- blood_sci_fpkm %>% 
  tidyr::separate(gene,into = c("gene_id","gene"),sep=">")
blood_sci_fpkm <- blood_sci_fpkm[,-1]




### 基因名称去重(保留最大值法)
### 列转行名一定要去重，因为行名不支持重复
blood_sci_fpkm <- blood_sci_fpkm %>%
  mutate(newcolumn = rowMeans(.[,-1])) %>% 
  arrange(desc(newcolumn)) %>% 
  distinct(gene,.keep_all = T) %>% 
  dplyr::select(-newcolumn)

rownames(blood_sci_fpkm) <- blood_sci_fpkm[,1]
blood_sci_fpkm <- blood_sci_fpkm[,-1]
### 保存数据
save(blood_sci_fpkm,file = "data/blood_sci_fpkm.Rdata")

### 6.行列转置
blood_sci_fpkm <- t(blood_sci_fpkm)
blood_sci_fpkm <- as.data.frame(blood_sci_fpkm)

### 7.添加分组
group <- c(rep(c("acute","3mpi","6mpi","12mpi"),4),rep(c("acute","3mpi","6mpi"),1),rep(c("acute","3mpi","6mpi","12mpi"),1),rep(c("acute","3mpi","6mpi"),1),rep(c("acute","3mpi","6mpi","12mpi"),1),rep(c("acute","3mpi","6mpi"),1),rep(c("acute","3mpi","6mpi","12mpi"),1),rep(c("control"),9))
blood_sci_fpkm <- cbind(group= group,blood_sci_fpkm)
table(blood_sci_fpkm$group)
table(c("METTL3" %in% colnames(blood_sci_fpkm)))
gene <- c("FTO","TRMT6","WTAP","YTHDF2","YTHDF3","YTHDC1")
blood_sci_fpkm$group <- factor(blood_sci_fpkm$group , levels = c("control","acute","3mpi","6mpi","12mpi"))

## steal plot
my_comparisons <- list(
  c("acute", "control"), 
  c("3mpi", "control"), 
  c("6mpi", "control"),
  c("12mpi", "control")
)


## 改写成函数
diffplot <- function(gene){
  my_comparisons <- list(
    c("acute", "control"), 
    c("3mpi", "control"), 
    c("6mpi", "control"),
    c("12mpi", "control")
  )
  library(ggpubr)
  ggboxplot(
    blood_sci_fpkm, x = "group", y = gene,
    color = "group", palette = c("#E64B357F","#4DBBD57F","#00A0877F","#3C54887F","#F39B7F7F"),
    add = "jitter"
  )+
    stat_compare_means(comparisons = my_comparisons, method = "t.test")
}

### AGR3,ESR1,SLC4A10, ALPP,VSIR,PLA2G2F
diffplot("YTHDF2")
diffplot("FTO")

## 多个基因作图查看
## 先把基因提取出来
genelist <- c("FTO","TRMT6","WTAP","YTHDF2","YTHDF3","YTHDC1")
## 再提取表达量，使用名称选取行
data <- blood_sci_fpkm[,c("group",genelist)]
## 用pivot_longer调整数据，数据变长，增加的是行
library(tidyr)
data <- data %>% 
  pivot_longer(cols=-1,
               names_to= "gene",
               values_to = "expression")
## 多基因作图
## 作图
data %>%
  mutate(group = factor(group, levels = c("acute","3mpi","6mpi","12mpi"))) %>%
  ggplot(data = .,aes(x=group,y=expression,fill=group))+
  geom_boxplot()+
  geom_jitter()+
  theme_bw()+
  facet_grid(gene~.)+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test")
## 图片导出
library(export)
## 导成PPT可编辑的格式
graph2ppt(file="output/diffgenboxplot.pptx")
## 其他自己想要的格式
graph2pdf(file="output/diffgenboxplot.pdf")
graph2tif(file="output/diffgenboxplot.tif")
## 导成AI可以编辑的状态
graph2eps(file="output/diffgenboxplot.eps")








