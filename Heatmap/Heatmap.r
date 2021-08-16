library(GEOquery)
library(pheatmap) 

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

#correlation Heatmap
pdf("Results/corheatmap.pdf", width = 15 , height = 15)
pheatmap(cor(ex), labels_row = gr , labels_col = gr)
dev.off()
