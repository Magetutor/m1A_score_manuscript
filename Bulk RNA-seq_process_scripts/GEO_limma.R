#### 1.公共数据库的读取 ####

### 1.1 GEOquery读取数据 首选
## 加载R包
library(GEOquery)
## 下载数据，GSExxxx
gset = getGEO('GSE32575', destdir=".",getGPL = F)
## 获取ExpressionSet对象，包括的表达矩阵和分组信息
gset=gset[[1]]
## 分组信息-重点关注GSM与可以传给group的信息
pdata=pData(gset)
## 基因表达量-
exprSet=exprs(gset)
## 只要后36个,本次选择的这36个是配对的。
## 所以，别人的芯片我们也不是全要，我们只要适合自己的数据
group_list=c(rep('before',18),rep('after',18))
### 分组变成向量，并且限定leves的顺序
### levels里面，把对照组放在前面(把要比的放在前面)
group <- factor(group,levels = c("before","after"))

## 可以先简单看一下整体样本的表达情况，col为之前给定好的分组信息
boxplot(exprSet[,-c(1:12)],outline=F,notch=T,col=group_list, las=2)

### 1.2 在网站上直接下载数据 






#### 2.芯片数据特有处理方式——芯片数据的表达量不齐，需要人工进行校正 ####
library(limma)
exprSet=normalizeBetweenArrays(exprSet)
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)
## 这步把矩阵转换为数据框很重要
exprSet = as.data.frame(exprSet)


#### 3.数据log化处理判断 ####

## 数值太大的话则需要log化，此处自动判断，输出为exprSet
ex <- exprSet
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

if (LogC) { ex[which(ex <= 0)] <- NaN
exprSet <- log2(ex)
print("log2 transform finished")}else{print("log2 transform not needed")}






#### 4.芯片数据特有处理方式——基因探针名称转化为gene_symbol ####
platformMap <- data.table::fread("data/platformMap.txt")
## 通过GPLXXX 来获取需要安装的包的名称
index = gset@annotation  #获取平台名称
platformDB = paste0(platformMap$bioc_package[grep(index,platformMap$gpl)],".db") #通过平台来获取包的名称
## 第一句是为了增加镜像
if(length(getOption("BioC_mirror"))==0) options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
## 没有这个包就给你装，有就不装
if(!require("illuminaHumanv2.db")) BiocManager::install("illuminaHumanv2.db",update = F,ask = F)
## 加载R包 
library(illuminaHumanv2.db)
## 获取表达矩阵的行名，就是探针名称
probeset <- rownames(exprSet)
## 使用lookup函数，找到探针在illuminaHumanv2.db中的对应基因名称
## 如果分析别的芯片数据，把illuminaHumanv2.db更换即可
SYMBOL <-  annotate::lookUp(probeset,"illuminaHumanv2.db", "SYMBOL")
## 转换为向量
symbol = as.vector(unlist(SYMBOL))

probe2symbol = data.frame(probeset,symbol,stringsAsFactors = F)








#### 5.基因信息合并与基因去重 ####

### 5.1 基因探针与gene_symbol转化
library(dplyr)
library(tibble)
exprSet <- exprSet %>% 
  rownames_to_column(var="probeset") %>% 
  #合并探针的信息
  inner_join(probe2symbol,by="probeset") %>% 
  #去掉多余信息
  select(-probeset) %>% 
  #重新排列
  select(symbol,everything()) %>% 
  #求出平均数(这边的点号代表上一步产出的数据)
  mutate(rowMean =rowMeans(.[grep("GSM", names(.))])) %>% 
  #去除symbol中的NA
  filter(symbol != "NA") %>% 
  #把表达量的平均值按从大到小排序
  arrange(desc(rowMean)) %>% 
  # symbol留下第一个
  distinct(symbol,.keep_all = T) %>% 
  #反向选择去除rowMean这一列
  select(-rowMean) %>% 
  # 列名变成行名
  column_to_rownames(var = "symbol")


#### limma包核心——差异基因分析 ####
### 表达矩阵：经过上述过程处理好的数据框，列名为样本名，行名为基因名
### 分组信息：group_list：factor格式的向量
### 差异分析没有配对
## 1.构建比较矩阵
design=model.matrix(~ group_list)

fit=lmFit(exprSet,design)
fit=eBayes(fit)
allDiff=topTable(fit,adjust='fdr',coef="group_listafter",number=Inf,p.value=0.05) 
### 差异分析，配对--样本有一一对应的关系
##这里配对和非配对的差异分析，只相差一步，也就是是否有配对信息
pairinfo = factor(rep(1:18,2))
design=model.matrix(~ pairinfo+group_list)
fit=lmFit(exprSet,design)
fit=eBayes(fit) 
allDiff_pair=topTable(fit,adjust='BH',coef="group_listafter",number=Inf,p.value=0.05)



#### 配对基因作图展示——循环作图与plot_grid应用 ####
### ggplot用的是清洁数据：行是样本；列为基因；分组信息进行汇总
data_plot = as.data.frame(t(exprSet))
data_plot = data.frame(pairinfo=rep(1:18,2),
                       group=group_list,
                       data_plot,stringsAsFactors = F) #创建数据框的时候字符串会变成因子
## 作图展示，此时还没有有配对信息
library(ggplot2)
ggplot(data_plot, aes(group,CAMKK2,fill=group)) +
  geom_jitter(aes(colour=group), size=2, alpha=0.5)+
  xlab("") +
  ylab(paste("Expression of ","CAMKK2"))+
  theme_classic()+
  theme(legend.position = "none")
## 加上配对信息以后
ggplot(data_plot, aes(group,CAMKK2,fill=group)) +
  geom_boxplot() +
  geom_point(size=2, alpha=0.5) +
  geom_line(aes(group=pairinfo), colour="black", linetype="11") +
  xlab("") +
  ylab(paste("Expression of ","CAMKK2"))+
  theme_classic()+
  theme(legend.position = "none")
## 循环作图后组合
library(dplyr)
library(tibble)
allDiff_arrange <- allDiff_pair %>% 
  rownames_to_column(var="genesymbol") %>% 
  arrange(desc(abs(logFC)))
genes <- allDiff_arrange$genesymbol[1:8]

plotlist <- lapply(genes, function(x){
  data =data.frame(data_plot[,c("pairinfo","group")],gene=data_plot[,x])
  ggplot(data, aes(group,gene,fill=group)) +
    geom_boxplot() +
    geom_point(size=2, alpha=0.5) +
    geom_line(aes(group=pairinfo), colour="black", linetype="11") +
    xlab("") +
    ylab(paste("Expression of ",x))+
    theme_classic()+
    theme(legend.position = "none")
})

library(cowplot)
plot_grid(plotlist=plotlist, ncol=4,labels = LETTERS[1:8])














#### 基因分组信息的质量评估(看不同处理组的样本是否分隔开)：相关性分析、PCA、热图 ####
### 相关性分析 ###

### PCA分析 ###

### 热图 ###


#### 基因表达量数据分析：折线图、箱线图、小提琴图——单个/多个基因长宽转换 ####
### 折线图 ###

### 箱线图/小提琴图——加P值、选定分析方法 ###


#### 差异基因分析：热图、火山图、富集分析(KEGG\GO)、GSEA ####
### 热图 ###

### 火山图——不同类型火山图 ###

### KEGG/GO/pathview ###

### GSEA ###








