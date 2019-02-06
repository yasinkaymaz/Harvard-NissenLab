.libPaths("/n/home13/yasinkaymaz/biotools/Rlibs/")

source("https://raw.githubusercontent.com/yasinkaymaz/cellRtools/master/main/functions.R")

#DE between IT and rest of the Gad1 neurons within Tasic2018
load("~/data/Tasic2018.seurat.Robj")

#L5 PT Specific Differential expression
cells <- rownames(Tasic2018@meta.data[which(Tasic2018@meta.data$sex == "F"),])
Tasic2018.F <- SubsetData(object = Tasic2018, cells.use = cells,do.center = T,do.scale = T,do.clean = T)
save(Tasic2018.F, file="~/data/Tasic2018.F.seurat.Robj")

TestGroups <- Tasic2018.F@meta.data  %>% mutate(TestGroups = if_else( subclass == "L5 PT", "L5 PT", as.character(class) )) %>% select(TestGroups)
Tasic2018.F@meta.data$TestGroups <- TestGroups$TestGroups
Tasic2018.F = SetAllIdent(Tasic2018.F, id = 'TestGroups')

markers.L5PT.vs.Gaba <- FindMarkers(object = Tasic2018.F, ident.1 = "L5 PT", ident.2 = "GABAergic", only.pos = T)
save(markers.L5PT.vs.Gaba, file="~/data/markers.L5PT.vs.Gaba.Rdata")
markers.L5PT.vs.Glu <- FindMarkers(object = Tasic2018.F, ident.1 = "L5 PT", ident.2 = "Glutamatergic", only.pos = T)
save(markers.L5PT.vs.Glu, file="~/data/markers.L5PT.vs.Glu.Rdata")
markers.L5PT.vs.Gaba_Glu <- FindMarkers(object = Tasic2018.F, ident.1 = "L5 PT", ident.2 = c("GABAergic","Glutamatergic"), only.pos = T)
save(markers.L5PT.vs.Gaba_Glu, file="~/data/markers.L5PT.vs.Gaba_Glu.Rdata")


#L5 IT Specific Differential expression
TestGroups <- Tasic2018@meta.data  %>% mutate(TestGroups = if_else( subclass == "L5 IT", "L5 IT", as.character(class) )) %>% select(TestGroups)
Tasic2018@meta.data$TestGroups <- TestGroups$TestGroups
Tasic2018 = SetAllIdent(Tasic2018, id = 'TestGroups')

markers.L5IT.vs.Gaba <- FindMarkers(object = Tasic2018, ident.1 = "L5 IT", ident.2 = "GABAergic", only.pos = T)
save(markers.L5IT.vs.Gaba, file="~/data/markers.L5IT.vs.Gaba.Rdata")
markers.L5IT.vs.Glu <- FindMarkers(object = Tasic2018, ident.1 = "L5 IT", ident.2 = "Glutamatergic", only.pos = T)
save(markers.L5IT.vs.Glu, file="~/data/markers.L5IT.vs.Glu.Rdata")
markers.L5IT.vs.Gaba_Glu <- FindMarkers(object = Tasic2018, ident.1 = "L5 IT", ident.2 = c("GABAergic","Glutamatergic"), only.pos = T)
save(markers.L5IT.vs.Gaba_Glu, file="~/data/markers.L5IT.vs.Gaba_Glu.Rdata")



