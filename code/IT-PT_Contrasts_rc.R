.libPaths("/n/home13/yasinkaymaz/biotools/Rlibs/")

source("https://raw.githubusercontent.com/yasinkaymaz/cellRtools/master/main/functions.R")

#DE between IT and rest of the Gad1 neurons within Tasic2018
load("~/data/Tasic2018.seurat.Robj")
Tasic2018 <- SeuratObj
rm(SeuratObj)



#L5 IT Specific Differential expression
TestGroups <- Tasic2018@meta.data  %>% mutate(TestGroups = if_else( subclass == "L5 IT", "L5 IT", as.character(class) )) %>% select(TestGroups)
Tasic2018@meta.data$TestGroups <- TestGroups$TestGroups
Tasic2018 = SetAllIdent(Tasic2018, id = 'TestGroups')

excludedGroups <- c("Low Quality", "No Class")
#cells <- rownames(Tasic2018@meta.data[which(Tasic2018@meta.data$TestGroups %in% excludedGroups),])
Tasic2018.sub <- SubsetData(object = Tasic2018, ident.remove = excludedGroups, do.center = T,do.scale = T,do.clean = T)

markers.L5IT <- FindAllMarkers(object = Tasic2018.sub, only.pos=T)
markers.L5IT.vs.Gaba <- FindMarkers(object = Tasic2018.sub, ident.1 = "L5 IT", ident.2 = "GABAergic", only.pos = T)
save(markers.L5IT.vs.Gaba, file="~/data/markers.L5IT.vs.Gaba.Rdata")
markers.L5IT.vs.Glu <- FindMarkers(object = Tasic2018.sub, ident.1 = "L5 IT", ident.2 = "Glutamatergic", only.pos = T)
save(markers.L5IT.vs.Glu, file="~/data/markers.L5IT.vs.Glu.Rdata")
markers.L5IT.vs.Gaba_Glu <- FindMarkers(object = Tasic2018.sub, ident.1 = "L5 IT", ident.2 = c("GABAergic","Glutamatergic"), only.pos = T)
save(markers.L5IT.vs.Gaba_Glu, file="~/data/markers.L5IT.vs.Gaba_Glu.Rdata")

load("~/data/markers.L5IT.vs.Gaba_Glu.Rdata")

marker.it <- markers.L5IT.vs.Gaba_Glu %>%
  add_column(gene=row.names(markers.L5IT.vs.Gaba_Glu)) %>%
  select(-p_val)%>%
  filter(p_val_adj < 0.01) %>%
  filter(pct.1 > 0.5) %>% filter( ((pct.1 - pct.2)/max(pct.1,pct.2) > 0.5))

VlnPlot(object = Tasic2018, features.plot = marker.it$gene, use.raw = TRUE, y.log = TRUE,nCol = 1,x.lab.rot = T)

pdf("output/test.pdf",width = 20,height = 270)
VlnPlot(object = Tasic2018, features.plot = marker.it$gene, use.raw = TRUE, y.log = TRUE,nCol = 1,x.lab.rot = T)
dev.off()

Tasic2016 = SetAllIdent(Tasic2016, id = 'broad_type')
pdf("output/test.pdf",width = 20,height = 270)
VlnPlot(object = Tasic2016, features.plot = marker.it$gene, use.raw = TRUE, y.log = TRUE,nCol = 1,x.lab.rot = T)
dev.off()



