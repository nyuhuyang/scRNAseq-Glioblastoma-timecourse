# Install devtools from CRAN
install.packages("devtools")
# Use devtools to install hdf5r and loomR from GitHub
install.packages("hdf5r")
devtools::install_github(repo = "mojaveazure/loomR", ref = "develop")

# loomR and Seurat================
library(Seurat) #v2.3.4 main
library(magrittr)
(load(file="data/Glioblastoma_12_20190214.Rda"))
object %<>% SetAllIdent("cell.lines")
table(object@meta.data$cell.lines)
TSNEPlot.1(object,do.label = T)
(cell.lines <- unique(object@meta.data$cell.lines) %>% as.character())
# Convert from Seurat to loom Convert takes and object in 'from', a name of
# a class in 'to', and, for conversions to loom, a filename
for(sample in cell.lines){
        sub_object <- SubsetData(object,ident.use = sample)
        print(sample)
        pfile <- Convert(from = sub_object, to = "anndata", filename = paste0("data/h5ad/",sample,".h5ad"), 
                         display.progress = F)
}
remove(pfile);GC()
