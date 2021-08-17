library(GEOquery)
library(limma)

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

## Differential expression analysis
# turn our groups to factor & assign it as gset$description
gr = factor(gr)
gset$description = gr 

# create a design matrix to show samples are belong to which level
design = model.matrix(~description + 0, gset)
colnames(design) = levels(gr)

# fit a leaner model to data
fit = lmFit(gset, design)

# define the levels which we want make contrast & calculate p-value
cont.matrix = makeContrasts(cALL_BM-CD10, levels=design)
fit2 = contrasts.fit(fit, cont.matrix)
fit2 = eBayes(fit2, 0.01)

# make a top Table / adjust by false discovery rate (fdr)(benjamini-hochberg), default one
tT = topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
# extract a subset of required columns and write the table
tT = subset(tT, select=c("Gene.symbol","Gene.ID","adj.P.Val","logFC"))
write.table(tT, "Results/cALL_CD10.txt", row.names=F, sep="\t", quote = F)

# extract a subset of genes with adjusted_P_Value < 0.05 & log fold change more than 2_fold
ALL.up = subset(tT, logFC > 1 & adj.P.Val < 0.05)
ALL.up.genes = unique(as.character(strsplit2(ALL.up$Gene.symbol, "///")))
write.table(ALL.up.genes, file = "Results/cALL_CD10-UP.txt", quote = F, row.names = F , col.names = F)
 
#extract a subset of genes with adjusted_P_Value < 0.05 & log fold change less than 1/2_fold
ALL.down = subset(tT, logFC <  -1 & adj.P.Val < 0.05)
ALL.down.genes = unique(as.character(strsplit2(ALL.down$Gene.symbol, "///")))
write.table(ALL.down.genes, file = "Results/cALL_CD10-down.txt", quote = F, row.names = F , col.names = F)

