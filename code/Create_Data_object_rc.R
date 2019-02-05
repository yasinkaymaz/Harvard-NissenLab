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

Tasic2018 <- SeuratWrapper(ExpData = exp, ProjectLabel = "Tasic2018", Normalize = T, NewMeta = meta, scale.only.var = T, PCs = 50, dump.files = F, min.cells = 4, min.genes = 500)


#Tasic2016
#cd ~/Documents/Harvard_Informatics/Data_Explore/mouse/FullSets/GEO/Tasic2016/
#wget https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE71585&format=file&file=GSE71585%5FClustering%5FResults%2Ecsv%2Egz
#wget https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE71585&format=file&file=GSE71585%5FRefSeq%5FTPM%2Ecsv%2Egz
exp <- read.csv("GSE71585_RefSeq_TPM.csv",row.names = 1,header = T)
meta <- read.csv("GSE71585_Clustering_Results.csv.gz", header=TRUE,row.names=1)
Tasic2016 <- SeuratWrapper(ExpData = exp, ProjectLabel = "Tasic2016", Normalize = F, scale.only.var = T, PCs = 20, dump.files = T, NewMeta = meta)