.libPaths("/n/home13/yasinkaymaz/biotools/Rlibs/")

source("https://raw.githubusercontent.com/yasinkaymaz/cellRtools/master/main/functions.R")

#DE between IT and rest of the Gad1 neurons within Tasic2018
load("~/data/Tasic2018.seurat.Robj")
Tasic2018 <- SeuratObj
rm(SeuratObj)

latentVars <- c("nUMI","nGene","percent_ecoli_reads","percent_synth_reads","percent_mt_exon_reads",
                "percent_rrna_reads")
class_of_interest <- c("Glutamatergic")
cells <- rownames(Tasic2018@meta.data[which(Tasic2018@meta.data$class %in% class_of_interest),])
Glutamatergics <- SubsetData(object = Tasic2018, cells.use = cells,do.center = T,do.scale = T,do.clean = T)
#Glutamatergics <- QuickSeurat(Glutamatergics, vars2reg = latentVars, scale.only.var = F)
Glutamatergics <- QuickSeurat(Glutamatergics, scale.only.var = F)
#DE of L5 IT or L5 PT against all other subclasses:
Glutamatergics = SetAllIdent(Glutamatergics, id = 'subclass')
#allmarkers.glu.subclass <- FindAllMarkers(object = Glutamatergics, test.use = "MAST", latent.vars = latentVars, logfc.threshold = 0.1,min.pct = 0.05,only.pos = T)
allmarkers.glu.subclass <- FindAllMarkers(object = Glutamatergics,only.pos=T)

save(Glutamatergics, allmarkers.glu.subclass, file="~/data/Tasic2018-reAnalysis_Glutamatergics.RData")
