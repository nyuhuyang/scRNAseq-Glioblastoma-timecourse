library(SingleR)
library(Seurat)
library(reshape2)
library(pheatmap)
library(kableExtra)
source("../R/Seurat_functions.R")
source("../R/SingleR_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 3.1 Create Singler Object  ==========================================
(load(file = "./data/Glioma_Harmony_20181201.Rda"))
(load(file='./data/GeneSets/ref_GBM_RSEM.RData'))
(load(file='./data/GeneSets/ref_GBMLGG_RSEM.RData'))
length(ref$types)
length(unique(ref$types))
length(unique(ref$main_types))
#pca
TSNEPlot.1(object = Glioma, no.legend = TRUE,
        do.return = TRUE,vector.friendly = F, pt.size = 1,
        do.label = TRUE,label.size = 8, group.by = "ident") + 
        ggtitle("Cluster ID") + 
        theme(plot.title = element_text(hjust = 0.5))
Glioma <- SetAllIdent(Glioma, id = "res.0.6")
singler = CreateSinglerObject(as.matrix(Glioma@data), annot = NULL, 
                              project.name="Glioma",
                              min.genes = 500,technology = "10X", species = "Human", citation = "",
                              ref.list = list(ref),normalize.gene.length = F, variable.genes = "de",
                              fine.tune = F, do.signatures = F, clusters = NULL,
                              numCores = SingleR.numCores/4)

# if singler didn't find all cell labels
length(singler$singler[[1]]$SingleR.single$labels) == ncol(Glioma@data)
singler$meta.data$orig.ident = Glioma@meta.data$orig.ident # the original identities, if not supplied in 'annot'
singler$meta.data$xy = Glioma@dr$tsne@cell.embeddings # the tSNE coordinates
singler$meta.data$clusters = Glioma@ident # the Seurat clusters (if 'clusters' not provided)
save(singler,file="./output/singler_Glioma_20190205.RData")