library(GEOquery)
library(ggplot2)

# get the data from GEO ( GSE matrix & annotation GPL)
# define destination to save data, SO didn't need load data each time
series = "GSE34670"
platform = "GPL96"
gset = getGEO(series, GSEMatrix =TRUE, AnnotGPL=TRUE, destdir = "Data/")
# if the data have several platforms, choose the certain platform
if (length(gset) > 1) idx <- grep(platform, attr(gset, "names")) else idx <- 1
gset = gset[[idx]]

# define groups for samples
gr = c(rep("cALL_BM", 5) , "cALL_PB" ,rep("cALL_BM", 19), "CL_cALL2" , "CL_697" , "CL_NALM6" , rep("CD10", 6), rep("CD10_pool" , 3) )
# elicit expression matrix from data
ex = exprs(gset)

# log2 transformation
ex = log2(ex + 1 )
exprs(gset) = ex

# principal component analysis (PCA)
# PCA for genes & plot it
pc = prcomp(ex)
pdf("Results/pc.pdf")
plot(pc)
plot(pc$x[,1:2])
dev.off()

# zero mean normalization & compute PCA & plot it
ex.scale = t(scale(t(ex), scale = F))
pc = prcomp(ex.scale)
pdf("Results/pc_scaled.pdf")
plot(pc)
plot(pc$x[,1:2])
dev.off()

# PCA for samples & plot it
pcr = data.frame(pc$r[,1:3], Group=gr)
pdf("Results/PCA_SAmples.pdf")
ggplot(pcr, aes(PC1, PC2, color=Group)) + geom_point(size=3) + theme_bw()
dev.off()
