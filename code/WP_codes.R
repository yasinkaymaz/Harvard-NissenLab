
source("https://raw.githubusercontent.com/yasinkaymaz/cellRtools/master/main/functions.R")

#Tasic2016
#cd ~/Documents/Harvard_Informatics/Data_Explore/mouse/FullSets/GEO/Tasic2016/
#wget https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE71585&format=file&file=GSE71585%5FClustering%5FResults%2Ecsv%2Egz
#wget https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE71585&format=file&file=GSE71585%5FRefSeq%5FTPM%2Ecsv%2Egz
exp <- read.csv("GSE71585_RefSeq_TPM.csv",row.names = 1,header = T)
meta <- read.csv("GSE71585_Clustering_Results.csv.gz", header=TRUE,row.names=1)
Tasic2016 <- SeuratWrapper(ExpData = exp, ProjectLabel = "Tasic2016", Normalize = F, scale.only.var = T, PCs = 20, dump.files = T, NewMeta = meta)


#ALM data only
#cd ~/data
#wget http://celltypes.brain-map.org/api/v2/well_known_file_download/694413179
#unzip mouse_ALM_gene_expression_matrices_2018-06-14.zip
here::here()
exp <- read.csv("~/data/mouse_ALM_2018-06-14_exon-matrix.csv",row.names = 1,header = T)
meta <- read.csv("~/data/mouse_ALM_2018-06-14_samples-columns.csv", row.names = 1, header = T)
geneinfo <- read.csv("~/data/mouse_ALM_2018-06-14_genes-rows.csv", row.names = 4, header = T)

identical(rownames(exp), rownames(geneinfo))
rownames(exp) <- geneinfo$gene_symbol

Tasic2018 <- SeuratWrapper(ExpData = exp, ProjectLabel = "Tasic2018", Normalize = T, NewMeta = meta, scale.only.var = T, PCs = 50, dump.files = F)
head(Tasic2018@meta.data)

TSNEPlot(Tasic2018, group.by="brain_subregion")
TSNEPlot(Tasic2018, group.by="subclass")
TSNEPlot(Tasic2018, group.by="class")
TSNEPlot(Tasic2018, group.by="cluster")

#"L5 IT ALM Npw" --> was left out
clusters_of_interest <- c("L5 IT ALM Cbln4 Fezf2",
                          "L5 IT ALM Cpa6 Gpr88",
                          "L5 IT ALM Gkn1 Pcdh19",
                          "L5 IT ALM Lypd1 Gpr88",
                          "L5 IT ALM Pld5",
                          "L5 IT ALM Tmem163 Arhgap25",
                          "L5 IT ALM Tmem163 Dmrtb1",
                          "L5 IT ALM Tnc",
                          "L5 PT ALM Hpgd",
                          "L5 PT ALM Npsr1",
                          "L5 PT ALM Slco2a1")
cells <- rownames(Tasic2018@meta.data[which(Tasic2018@meta.data$cluster %in% clusters_of_interest),])
L5ALM.ITPT <- SubsetData(object = Tasic2018, cells.use = cells,do.center = T,do.scale = T,do.clean = T)

latentVars <- c("nUMI","nGene","percent_ecoli_reads","percent_synth_reads","percent_mt_exon_reads",
                "percent_rrna_reads")

L5ALM.ITPT <- QuickSeurat(L5ALM.ITPT, vars2reg = latentVars)

pdf("clusters.pdf",width = 10,height = 5)
TSNEPlot(L5ALM.ITPT, group.by="brain_subregion")
TSNEPlot(L5ALM.ITPT, group.by="subclass")
TSNEPlot(L5ALM.ITPT, group.by="class")
TSNEPlot(L5ALM.ITPT, group.by="cluster", do.label = T)
TSNEPlot(L5ALM.ITPT, group.by="res.1", do.label = T)
dev.off()

head(L5ALM.ITPT@meta.data)
TSNEPlot(L5ALM.ITPT, group.by="sex")
TSNEPlot(L5ALM.ITPT, group.by="age_days")
TSNEPlot(L5ALM.ITPT, group.by="brain_hemisphere")

L5ALM.ITPT = SetAllIdent(L5ALM.ITPT, id = 'subclass')
allmarkers.subclass <- FindAllMarkers(object = L5ALM.ITPT, test.use = "MAST", latent.vars = latentVars, logfc.threshold = 0.1,min.pct = 0.05,only.pos = T)

save(Tasic2018, L5ALM.ITPT, allmarkers.subclass, file=here::here("data/Tasic2018-reAnalysis.RData"))

