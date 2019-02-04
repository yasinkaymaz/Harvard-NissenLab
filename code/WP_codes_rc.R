.libPaths("/n/home13/yasinkaymaz/biotools/Rlibs/")

source("https://raw.githubusercontent.com/yasinkaymaz/cellRtools/master/main/functions.R")

#Tasic 2018
#ALM data
exp <- read.csv("~/data/mouse_ALM_2018-06-14_exon-matrix.csv",row.names = 1,header = T)
meta <- read.csv("~/data/mouse_ALM_2018-06-14_samples-columns.csv", row.names = 1, header = T)
geneinfo <- read.csv("~/data/mouse_ALM_2018-06-14_genes-rows.csv", row.names = 4, header = T)
identical(rownames(exp), rownames(geneinfo))
rownames(exp) <- geneinfo$gene_symbol
alm.exp <- exp
alm.meta <- meta
rm(exp,meta,geneinfo)

#VISp
exp <- read.csv("~/data/mouse_VISp_2018-06-14_exon-matrix.csv",row.names = 1,header = T)
meta <- read.csv("~/data/mouse_VISp_2018-06-14_samples-columns.csv", row.names = 1, header = T)
geneinfo <- read.csv("~/data/mouse_VISp_2018-06-14_genes-rows.csv", row.names = 4, header = T)
identical(rownames(exp), rownames(geneinfo))
rownames(exp) <- geneinfo$gene_symbol
visp.exp <- exp
visp.meta <- meta
rm(exp,meta)

identical(colnames(visp.meta), colnames(alm.meta))
identical(rownames(visp.exp), rownames(alm.exp))
exp <- cbind(alm.exp, visp.exp)
meta <- rbind(alm.meta, visp.meta)

Tasic2018 <- SeuratWrapper(ExpData = exp, ProjectLabel = "Tasic2018", Normalize = T, NewMeta = meta, scale.only.var = T, PCs = 50, dump.files = F, min.cells = 4, min.genes = 500)
