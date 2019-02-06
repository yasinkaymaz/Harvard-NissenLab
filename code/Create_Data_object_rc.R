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

#Filter genes:
load("~/data/Genes.to.filter.Rdata")

exp <- exp[which(!rownames(exp) %in% genesTofilter),]

Tasic2018 <- SeuratWrapper(ExpData = exp, ProjectLabel = "Tasic2018", Normalize = T, NewMeta = meta, scale.only.var = F, PCs = 50, dump.files = F, min.cells = 4, min.genes = 500)
save(Tasic2018, "~/data/Tasic2018.seurat.Robj")


#Tasic2016
#cd ~/Documents/Harvard_Informatics/Data_Explore/mouse/FullSets/GEO/Tasic2016/
#wget https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE71585&format=file&file=GSE71585%5FClustering%5FResults%2Ecsv%2Egz
#wget https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE71585&format=file&file=GSE71585%5FRefSeq%5FTPM%2Ecsv%2Egz
exp <- read.csv("~/data/GSE71585_RefSeq_counts.csv",row.names = 1,header = T)
meta <- read.csv("~/data/GSE71585_Clustering_Results.csv.gz", header=TRUE,row.names=1)
Tasic2016 <- SeuratWrapper(ExpData = exp, ProjectLabel = "Tasic2016", Normalize = T, scale.only.var = T, PCs = 20, dump.files = T, NewMeta = meta)


g.1 <- head(rownames(Tasic2016@hvg.info), 1000)
g.2 <- head(rownames(Tasic2018@hvg.info), 1000)

genes.use <- unique(c(g.1, g.2))
genes.use <- intersect(genes.use, rownames(Tasic2016@scale.data))
genes.use <- intersect(genes.use, rownames(Tasic2018@scale.data))

multiCCA_list <- list(Tasic2016 , Tasic2018 )
Tasic.combined <- RunCCA(object = Tasic2018, object2 = Tasic2016,genes.use = genes.use, num.cc = 4)

save(Tasic.combined, file="~/data/Tasic.combined.seurat.Robj")
