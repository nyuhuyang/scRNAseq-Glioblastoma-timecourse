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
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#3.1  Compare DE across all major cell types==================
#We would need the data for all clusters, as well the subclusters.
#detect changes in gene expression between young and aged, 
#in the different cell types and subtypes. 
#It will also be interesting to check if there is some subtype enriched in young compared to aged or viceversa. 

# 3.1.1 load data
# Rename ident
(load(file="data/Glioblastoma_V3_12_20190612.Rda"))
df_samples <- readxl::read_excel("doc/190207_scRNAseq_info.xlsx")
colnames(df_samples) <- colnames(df_samples) %>% tolower
sample_n = which(df_samples$tests %in% c("control",paste0("test",1:4)))
df_samples[sample_n,] %>% kable() %>% kable_styling()
table(df_samples$tests);nrow(df_samples)
exp <- AverageExpression(object)

exp_log2 <- as.data.frame(exp[[1]])*log2(exp(1))
write.csv(exp_log2, file = paste0(path,"exp_log2.csv"))
# FindAllmarkers ================
table(object@meta.data$conditions)
object@meta.data$orig.ident = gsub("-","_",object@meta.data$orig.ident)
Idents(object) = "cell.lines"
(cell.lines <- unique(Idents(object)) %>% as.character())
Split_object <- lapply(cell.lines, function(x) subset(object,idents = x))
names(Split_object) = cell.lines

markers_list <- list()
for(i in 1:length(Split_object)){
        Idents(Split_object[[i]]) = "conditions"
        markers_list[[i]] <- FindAllMarkers.UMI(Split_object[[i]],logfc.threshold = 0.1,
                                                  only.pos = FALSE, 
                                                  min.pct = 0.1,return.thresh = 0.05)
}
names(markers_list) = cell.lines
save(markers_list, file = paste0(path,"markers_list.Rda"))

Seurat::DotPlot()