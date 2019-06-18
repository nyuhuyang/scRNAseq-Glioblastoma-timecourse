########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
library(tidyr)
library(kableExtra)
library(magrittr)
library(gplots)
source("../R/Seurat_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#3.1  Compare DE across all major cell types==================
#We would need the data for all clusters, as well the subclusters.
#detect changes in gene expression between young and aged, 
#in the different cell types and subtypes. 
#It will also be interesting to check if there is some subtype enriched in young compared to aged or viceversa. 

# 3.1.1 load data
# Rename ident
(load(file="data/Glioblastoma_12_20190214.Rda"))
df_samples <- readxl::read_excel("doc/190207_scRNAseq_info.xlsx")
colnames(df_samples) <- colnames(df_samples) %>% tolower
sample_n = which(df_samples$tests %in% c("control",paste0("test",1:4)))
df_samples[sample_n,] %>% kable() %>% kable_styling()
table(df_samples$tests);nrow(df_samples)
samples <- df_samples$sample[sample_n]
sample.id <- df_samples$sample.id[sample_n]
conditions <- df_samples$conditions[sample_n]
projects <- df_samples$project[sample_n]
tissues <- df_samples$tissue[sample_n]
tests <- df_samples$tests[sample_n] 

# Doheatmap for Normal / MCL ================
table(object@meta.data$conditions)
object@meta.data$orig.ident = gsub("D21","D28",object@meta.data$orig.ident)
object@meta.data$conditions = gsub("D21","D28",object@meta.data$conditions)
object@meta.data$cell.lines = gsub("GSC","",object@meta.data$cell.lines)
object@meta.data$orig.ident = gsub("-","_",object@meta.data$orig.ident)
object %<>% SetAllIdent(id="conditions")
Split_object <- SplitSeurat(object, split.by = "conditions")
names(Split_object)

dent.1 = c("923","827","923")
dent.2 = c("827","NSC","NSC")

for(n in names(Split_object)){
        Split_object[[n]] %<>% SetAllIdent(id="orig.ident")
        gde.markers <- FindPairMarkers(Split_object[[n]], 
                                       ident.1 = paste(dent.1, n, sep="_"), 
                                       ident.2 = paste(dent.2, n, sep="_"),
                                       only.pos = FALSE,
                                       logfc.threshold = 0.1,min.cells.group =2,
                                       return.thresh = 0.05,
                                       save.path = paste0(path,n,"/"))
        GC();GC()
        all.markers <- FindAllMarkers.UMI(Split_object[[n]],
                                           only.pos = FALSE,
                                           logfc.threshold = 0.1,min.cells.group =2,
                                           return.thresh = 0.05)
        write.csv(all.markers,file = paste0(path,n,"/all_cell_lines_in_",n,".csv"))
        GC();GC()
}
