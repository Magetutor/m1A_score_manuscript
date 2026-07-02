spinal_Exp_vst <- t(spinal_Exp_vst)
spinal_Exp_vst <- as.data.frame(spinal_Exp_vst)
group <- rep(c("Control","Sham","0min","30min","1h","6h","12h",   
              "1d","3d","5d","7d","14d","1m","2m","3m"),each = 4)
spinal_Exp_vst <- cbind(group = group,spinal_Exp_vst)

library(biomaRt)
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
genes <- data.frame(m7G)
genesout=getLDS(attributes=c("hgnc_symbol"),filters="hgnc_symbol",
                values=genes,mart=human,
                attributesL=c("mgi_symbol"),
                martL=mouse,
                uniqueRows=T)
aaa <- as.data.frame(genes[!(genes$m7G %in% genesout$HGNC.symbol),])
m7G_m <- genesout[,2]
m7G_m <- edit(m7G_m)




