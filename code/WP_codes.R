
source("https://raw.githubusercontent.com/yasinkaymaz/cellRtools/master/main/functions.R")

#This was run on the cluster. bigmem
#Tasic2018 <- SeuratWrapper(ExpData = exp, ProjectLabel = "Tasic2018", Normalize = T, NewMeta = meta, scale.only.var = T, PCs = 50, dump.files = F, min.cells = 4, min.genes = 500)
load("~/data/Tasic2018.seurat.Robj")
Tasic2018 <- SeuratObj
rm(SeuratObj)

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
L5ALM.ITPT <- QuickSeurat(L5ALM.ITPT, vars2reg = latentVars, scale.only.var = F)

save(L5ALM.ITPT, file="~/data/L5ALM.ITPT.seurat.Robj")

L5ALM.ITPT = SetAllIdent(L5ALM.ITPT, id = 'subclass')
allmarkers.subclass <- FindAllMarkers(object = L5ALM.ITPT, test.use = "MAST", latent.vars = latentVars, logfc.threshold = 0.1,min.pct = 0.05,only.pos = T)
allmarkers.subclass <- FindAllMarkers(object = L5ALM.ITPT,only.pos = T)

L5ALM.ITPT = SetAllIdent(L5ALM.ITPT, id = 'cluster')
allmarkers.cluster <- FindAllMarkers(object = L5ALM.ITPT, test.use = "MAST", latent.vars = latentVars, logfc.threshold = 0.1,min.pct = 0.05,only.pos = T)
allmarkers.cluster <- FindAllMarkers(object = L5ALM.ITPT,only.pos = T)

save(L5ALM.ITPT, allmarkers.subclass, allmarkers.cluster, file="~/data/Tasic2018-reAnalysis_L5ITPT.RData")
rm(L5ALM.ITPT, allmarkers.subclass, allmarkers.cluster)


save(Tasic2018,
     L5ALM.ITPT,
     allmarkers.subclass,
     allmarkers.cluster,
     file="~/data/Tasic2018-reAnalysis.RData")


DElimma <- function(SeuratObj ){
  library(limma)
  expr <- SeuratObj@data
  # Filter out genes that are 0 for every cell in this cluster
  bad <- which(rowSums(expr) == 0)
  expr <- expr[-bad,]
  meta <- SeuratObj@meta.data[,c("orig.ident", "subclass")]
  meta$subclass <- make.names(meta$subclass)
  head(meta)
  mm <- model.matrix(~0 + subclass, data = meta)

  fit <- lmFit(expr, mm)
  head(coef(fit)) # means in each sample for each gene
  conditions <- colnames(coef(fit))
  contr <- makeContrasts(paste(conditions[1]," - ",conditions[2],sep=""), levels = colnames(coef(fit)))
  tmp <- contrasts.fit(fit, contrasts = contr)
  tmp <- eBayes(tmp)
  topTable(tmp, sort.by = "P", n = 200)

  return(topTable(tmp, p.value = 0.01,lfc = 1,n=20000))
}

detable <- DElimma(L5ALM.ITPT)
DoHeatmap(object = L5ALM.ITPT, use.scaled = T, genes.use = rownames(detable), slim.col.label = TRUE, remove.key = TRUE, group.label.loc="top", group.label.rot=T)

load("~/data/Tasic2018-reAnalysis_L5ITPT.RData")

pdf("output/test.pdf",width = 20,height = 90)
VlnPlot(object = Tasic2018, features.plot = c("Nrgn","Slc17a7","Rgs4","Arpp19","Ptn","Arpp21"), use.raw = TRUE, y.log = TRUE,nCol = 1,x.lab.rot = T)
dev.off()

pdf("output/test2.pdf",width = 20,height = 90)
VlnPlot(object = Tasic2018, features.plot = c("S100a10", "Fam3c", "Deptor", "Cnih3", "Rgs4", "Vstm2l"), use.raw = TRUE, y.log = TRUE,nCol = 1,x.lab.rot = T)
dev.off()

marker.it <- markers.L5IT.vs.Gaba_Glu %>%
  add_column(gene=row.names(markers.L5IT.vs.Gaba_Glu)) %>%
  select(-p_val)%>%
  filter(p_val_adj < 0.01) %>%
  filter(pct.1 > 0.5) %>% filter( ((pct.1 - pct.2)/max(pct.1,pct.2) > 0.5))


pdf("output/test3.pdf",width = 20,height = 90)
VlnPlot(object = Tasic2018, features.plot = c("Crym","S100a10","Hs3st2","Adcyap1"), use.raw = TRUE, y.log = TRUE,nCol = 1,x.lab.rot = T)
dev.off()


save(Tasic2018, file="~/data/Tasic2018.seurat.Robj")
save(Tasic2016, file="~/data/Tasic2016.seurat.Robj")
