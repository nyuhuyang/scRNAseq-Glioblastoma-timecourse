library(SingleR)
library(Seurat)
library(reshape2)
library(pheatmap)
library(kableExtra)
library(dplyr)
library(tidyr)
source("../R/Seurat_functions.R")
source("../R/SingleR_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

#====== 3.2 SingleR specifications ==========================================
# Step 1: Spearman coefficient
(load(file="data/MCL_Harmony_24_20190128.Rda"))
(load(file="output/singler_MCL_24T_20190128.Rda"))

# if singler didn't find all cell labels
if(length(singler$singler[[1]]$SingleR.single$labels) != ncol(object@data)){
        all.cell = object@cell.names;length(all.cell)
        know.cell = rownames(singler$singler[[1]]$SingleR.single$labels);length(know.cell)
        object = SubsetData(object, cells.use = know.cell)
}
object
##############################
# add singleR label to Seurat
###############################
singlerDF = data.frame("singler1sub" = singler$singler[[1]]$SingleR.single$labels,
                       "singler1main" = singler$singler[[1]]$SingleR.single.main$labels,
                       "orig.ident" = object@meta.data$orig.ident,
                       row.names = rownames(singler$singler[[1]]$SingleR.single$labels))

table(rownames(singlerDF) %in% object@cell.names)

apply(singlerDF,2,function(x) length(unique(x)))
object <- AddMetaData(object = object,
                   metadata = singlerDF)
object <- SetAllIdent(object = object, id = "singler1sub")

##############################
# check the spearman correlation
###############################
#Or by all cell types (showing the top 50 cell types):
jpeg(paste0(path,"/DrawHeatmap_sub1.jpeg"), units="in", width=10, height=7,
     res=600)
print(SingleR.DrawHeatmap(singler$singler[[1]]$SingleR.single, top.n = 50,normalize = F))
dev.off()
jpeg(paste0(path,"/DrawHeatmap_sub1_N.jpeg"), units="in", width=10, height=7,
     res=600)
print(SingleR.DrawHeatmap(singler$singler[[1]]$SingleR.single,top.n = 50,normalize = T))
dev.off()

#Finally, we can also view the labeling as a table compared to the original identities:

kable(table(singler$singler[[1]]$SingleR.single$labels, singler$meta.data$orig.ident)) %>%
        kable_styling()
singler$meta.data$orig.ident %>% table() %>% kable() %>% kable_styling()
object@meta.data$singler1sub %>% table() %>% kable() %>% kable_styling()

##############################
# process color scheme
##############################
singler_colors <- readxl::read_excel("./doc/singler.colors.xlsx")
singler_colors1 = as.vector(singler_colors$singler.color1[!is.na(singler_colors$singler.color1)])
singler_colors1[duplicated(singler_colors1)]
length(singler_colors1)
apply(object@meta.data[,c("singler1sub","singler1main")],
      2,function(x) length(unique(x)))
object@meta.data[,c("singler1sub")] %>% table() %>% kable() %>% kable_styling()
object <- AddMetaColor(object = object, label= "singler1sub", colors = singler_colors1[1:37])
object <- SetAllIdent(object = object, id = "singler1sub")
TSNEPlot.1(object, colors.use = ExtractMetaColor(object),no.legend = F)
##############################
# draw tsne plot
##############################
p3 <- TSNEPlot.1(object = object, do.label = T, group.by = "ident",
                 do.return = TRUE, no.legend = T,
                 colors.use = ExtractMetaColor(object),
                 pt.size = 1,label.size = 3,force = 2)+
  ggtitle("Supervised cell type labeling by Blueprint + Encode + MCL")+
  theme(text = element_text(size=10),							
        plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) 

jpeg(paste0(path,"PlotTsne_sub1.jpeg"), units="in", width=10, height=7,
     res=600)
print(p3)
dev.off()

save(object,file="data/MCL_Harmony_24_20190128.Rda")

##############################
# adjust cell label
##############################
# reduce false positive results (B cells are labeled as MCL in normal samples)
# and false negative results (MCL cells are labeled as B cells in MCL samples)
# singler1main false positive results  ========
table(singlerDF$singler1main, singlerDF$orig.ident) %>% t %>% kable %>% kable_styling()

normal_cells <- singlerDF$orig.ident %in% c("BH","DJ","MD","NZ") %>% rownames(singlerDF)[.]
singlerDF[normal_cells,"singler1main"] = gsub("MCL","B_cells",
                                              singlerDF[normal_cells,"singler1main"])

table(singlerDF$singler1main, singlerDF$orig.ident) %>% t %>% kable %>% kable_styling()

# singler1sub false positive results  =========
table(singlerDF$singler1sub, singlerDF$orig.ident) %>% kable %>% kable_styling()

singlerDF$singler1sub = gsub("MCL:.*$","MCL",singlerDF$singler1sub)
singlerDF[normal_cells,"singler1sub"] = gsub("MCL","B_cells:Memory",
                                              singlerDF[normal_cells,"singler1sub"])

table(singlerDF$singler1sub, singlerDF$orig.ident) %>% kable %>% kable_styling()

# false negative results ======================
object <- SetAllIdent(object, id = "orig.ident")
MCL <- SubsetData(object, ident.remove = c("BH","DJ","MD","NZ"))

FeaturePlot(MCL,features.plot = "CCND1")
CCND1.list <- SplitCells(MCL,split.by = "CCND1")
CCND1 <- SubsetData(object, cells.use = CCND1.list[[1]])
FeaturePlot(CCND1,features.plot = "CCND1")
remove(CCND1);GC()

# singler1main false negative results ======================
table(singlerDF[CCND1.list[[2]],"singler1main"]) %>% kable %>% kable_styling()
singlerDF[CCND1.list[[2]],"singler1main"] = gsub("B_cells","MCL",
                                             singlerDF[CCND1.list[[2]],"singler1main"])
# singler1sub false negative results ======================
table(singlerDF[CCND1.list[[2]],"singler1sub"]) %>% kable %>% kable_styling()
singlerDF[CCND1.list[[2]],"singler1sub"] = gsub("B_cells.*","MCL",
                                                 singlerDF[CCND1.list[[2]],"singler1sub"])
##############################
# process color scheme
##############################
singlerDF$singler1sub %>% table() %>% kable() %>% kable_styling()

singler_colors <- readxl::read_excel("./doc/singler.colors.xlsx")
singler_colors1 = as.vector(singler_colors$singler.color1[!is.na(singler_colors$singler.color1)])
singler_colors1[duplicated(singler_colors1)]
length(singler_colors1)
apply(singlerDF[,c("singler1sub","singler1main")],2,function(x) length(unique(x)))
object@meta.data[,c("singler1sub")] %>% table() %>% kable() %>% kable_styling()

object <- AddMetaData(object = object,metadata = singlerDF)

object <- AddMetaColor(object = object, label= "singler1sub", colors = singler_colors1[1:26])
object <- SetAllIdent(object = object, id = "singler1sub")
TSNEPlot.1(object, colors.use = ExtractMetaColor(object),no.legend = F)

save(object,file="data/MCL_Harmony_24_20190128.Rda")
##############################
# subset Seurat
###############################
# select 1/4 of cell from control
# in Identify_Cell_Types_Manually.R 2.2

table(object@meta.data$orig.ident)
table(object@ident)
object@meta.data$orig.ident = gsub("BH|DJ|MD|NZ","Normal",object@meta.data$orig.ident)
object %<>% SetAllIdent(id = "orig.ident")

df_samples <- readxl::read_excel("doc/190126_scRNAseq_info.xlsx")
colnames(df_samples) <- tolower(colnames(df_samples))

tests <- paste0("test",2:7)
for(test in tests){
        sample_n = which(df_samples$tests %in% test)
        df <- as.data.frame(df_samples[sample_n,])
        samples <- unique(df$sample)
        rownames(df) = samples
        
        print(samples <- c("Normal",df$sample[order(df$tsne)]))
        
        g <- lapply(samples,function(sample) {
                SubsetData(object, ident.use = sample) %>%
                        SetAllIdent(id = "singler1sub") %>%
                        TSNEPlot.1(no.legend = T,do.label =F,label.size=3,size=20,
                                   colors.use = ExtractMetaColor(.),
                                   return.plots =T, label.repel = T,force=2)+
                        ggtitle(sample)+theme(text = element_text(size=20),
                              plot.title = element_text(hjust = 0.5))
        })
        jpeg(paste0(path,test,"_TSNEPlot.jpeg"), units="in", width=10, height=7,
             res=600)
        print(do.call(plot_grid, c(g, nrow = 2)))
        dev.off()
}

###################################
# creat seurat object
###################################
(load(file = 'data/ref_MCL_blue_encode.RData'))
MCL_blue_encode <- CreateSeuratObject(raw.data = ref$data) %>%
  NormalizeData %>% ScaleData

MCL_blue_encode@data = ref$data
ident = ref$main_types
ident[!grepl("B_cells|MCL",ident)] = "others"
names(ident) = MCL_blue_encode@cell.names
MCL_blue_encode@ident = factor(ident, levels = c("others","B_cells","MCL"))
MCL_B_markers <- FindAllMarkers.UMI(MCL_blue_encode, logfc.threshold = 0.1, only.pos = T,
                                   test.use = "MAST")
write.csv(MCL_B_markers, paste0(path,"MCL_B_markers.csv"))
g <- DoHeatmap.1(MCL_blue_encode, MCL_B_markers,Top_n = 50,
                 ident.use = paste("B cells vs. MCL vs. other cell types"),
                 group.label.rot = T,cex.row = 4,remove.key =T)
jpeg(paste0(path,"heatmap_B_MD_others.jpeg"), units="in", width=10, height=7,
     res=600)
print(g)
dev.off()
