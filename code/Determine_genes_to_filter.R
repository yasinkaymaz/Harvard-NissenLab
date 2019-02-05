#Determine gender specific genes within L5 IT/PT cells
source("https://raw.githubusercontent.com/yasinkaymaz/cellRtools/master/main/functions.R")

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

#Find the gene specific ones.
L5ALM.ITPT = SetAllIdent(L5ALM.ITPT, id = 'sex')
genderSpecificDE <- FindAllMarkers(object = L5ALM.ITPT, only.pos = T)
save(genderSpecificDE, file = "~/data/GenderSpecificDE-table.Rdata")

#Visualize
volcanoPlotly(genderSpecificDE)

#Gene filtration criteria (from Tasic et al., 2018 - Methods)
geneinfo <- read.csv("~/data/mouse_VISp_2018-06-14_genes-rows.csv", row.names = 4, header = T)
#1. Filter Gm- genes
gm.genes <- grep(pattern = "^Gm", x = rownames(x = geneinfo$gene_symbol), value = TRUE)
#2. Mitochondrial genes:
mt.genes <- droplevels(geneinfo[which(geneinfo$chromosome == "MT"),]$gene_symbol)
#3. Sex-specific genes:
sex.genes <- genderSpecificDE$gene
#4. Ribosomal genes:
r.genes <- droplevels(geneinfo[grep("ribosomal",x=geneinfo$gene_name,value = F),]$gene_symbol)

# filter genes before creating the seurat object. Rerun the analysis from the begining.
genesTofilter <- c(gm.genes, mt.genes, sex.genes, r.genes)
save(genesTofilter, file="~/data/Genes.to.filter.Rdata")

