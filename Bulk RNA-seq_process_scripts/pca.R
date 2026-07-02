# 获取样本与分组信息
blood_sci_fpkm <- cbind(rownames(blood_sci_fpkm),blood_sci_fpkm)
fpkm <- blood_sci_fpkm
groupdata <- blood_sci_fpkm[,1:2]
colnames(groupdata) <- c("sample","group")
groupdata$group <- factor(groupdata$group,levels=group)

# 表达量过滤
data <- blood_sci_fpkm
data <- data[rowSums(data)>1,]

# 剔除异常值
data <- na.omit(data)

### 将数据进行log（数值过大）
qx <- as.numeric(quantile(data, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

## 开始判断，数值小的则无需log化
if (LogC) { 
  data[which(data <= 0)] <- NaN
  ## 取log2
  dataprSet <- log2(data)
  print("log2 transform finished")
}else{ 
  print("log2 transform not needed")
}


# zscore标准化
data <- scale(data)
data <- t(scale(t(data),center=TRUE,scale=F))

# 主成分分析
pca <- prcomp(t(data),center=FALSE,scale.=FALSE)
# ??prcomp  # 查看帮助文档
# str(pca)

# 获取距离矩阵
pca_mat <- data.frame(pc1=pca$x[,1],pc2=pca$x[,2],sample=groupdata$sample,group=groupdata$group)
library(ggplot2)
p <- ggplot(pca_mat,aes(x=pc1,y=pc2,label=sample,colour=group)) +
  geom_point(size=3) +
  geom_tdatat(colour='black',size=3) +
  geom_hline(yintercept=0,linetype='dotdash',size=0.8,color='grey') +
  geom_vline(xintercept=0,linetype='dotdash',size=0.8,color='grey') +
  theme_bw() + 
  stat_ellipse(aes(fill = group), geom = 'polygon', level = 0.95, alpha = 0.1, show.legend = FALSE)
p
