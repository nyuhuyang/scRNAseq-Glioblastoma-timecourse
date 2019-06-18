########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(dplyr)
library(tidyr)
library(kableExtra)
library(magrittr)
library(gplots)
library(ggplot2)
library(fgsea)
library(tibble)
library(ggpubr)
library(ggsci)
source("../R/Seurat_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
# 3.1.1 load data
# Load some results from Seurat
load(file = paste0(path,"markers_list.Rda"))
res = markers_list[['NSC']]
table(res$cluster)
head(res)
res = res[order(res["p_val_adj"]),]
head(res, 20)

hallmark <- gmtPathways("../seurat_resources/msigdb/h.all.v6.2.symbols.gmt")
biocarta <- gmtPathways("../seurat_resources/msigdb/c2.cp.biocarta.v6.2.symbols.gmt")
kegg <- gmtPathways("../seurat_resources/msigdb/c2.cp.kegg.v6.2.symbols.gmt")
tft <- gmtPathways("../seurat_resources/msigdb/c3.tft.v6.2.symbols.gmt")
c6 <- gmtPathways("../seurat_resources/msigdb/c6.all.v6.2.symbols.gmt")
GO <- gmtPathways("../seurat_resources/msigdb/c5.all.v6.2.symbols.gmt")
allpathways <- c(hallmark,biocarta,kegg)

hallmark %>% head() %>% lapply(head)
biocarta %>% head() %>% lapply(head)
# Now, run the fgsea algorithm with 1000 permutations:

FgseaPlot <- function(pathways=hallmark, stats=res, nperm=1000,show=NULL,cluster = 1,
                      sample="",pathway.name = "Hallmark"){
        
        res = stats[order(stats["p_val_adj"]),]
        res1 = res[res$cluster == cluster,c("gene","avg_logFC")] %>% deframe

        fgseaRes <- fgsea(pathways=pathways, stats=res1, nperm=nperm)
        fgseaResTidy <- fgseaRes %>% as_tibble() %>% arrange(desc(NES))
        fgseaResTidy %<>% 
                dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
                arrange(padj)
        print(dim(fgseaResTidy))
        if(!is.null(show)){
                (N = nrow(fgseaResTidy)-1)
                fgseaResTidy =fgseaResTidy[c(1:(show/2),(N-show/2):N),]
        }
        path <- paste0("output/",gsub("-","",Sys.Date()),"/")
        if(!dir.exists(path)) dir.create(path, recursive = T)
        jpeg(paste0(path,sample,"_",cluster,"-",pathway.name,".jpeg"), units="in", width=10, height=7,res=600)
        g <- ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
                geom_col(aes(fill=pval<0.05)) +
                guides(fill = guide_legend(reverse = TRUE))+
                coord_flip() +
                labs(x="Pathway", y="Normalized Enrichment Score",
                     title=paste(pathway.name,"pathways in",sample,cluster)) + 
                theme_minimal()
        print(g)
        dev.off()
}

FgseaBarplot <- function(pathways=hallmark, stats=res, nperm=1000,cluster = 1,
                      sample="",pathway.name = "Hallmark", hjust=0.5){
        
        res = stats[order(stats["p_val_adj"]),]
        res1 = res[res$cluster == cluster,c("gene","avg_logFC")] %>% deframe
        
        fgseaRes <- fgsea(pathways=pathways, stats=res1, nperm=nperm)
        print(dim(fgseaRes))
        topPathwaysUp <- fgseaRes[NES > 0][head(order(padj), n=25), pathway]
        topPathwaysDown <- fgseaRes[NES < 0][head(order(padj), n=25), pathway]
        topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
        out_write=fgseaRes[match(topPathways,pathway)]
        colnames(out_write)[1] <- 'Pathways'
        out_write$To_Plot_value <- -log10(out_write$pval)
        out_write$sign <- ifelse(out_write$NES >0,1,-1)
        out_write$To_Plot_value <- out_write$To_Plot_value*out_write$sign
        out_write$sign <- ifelse(out_write$NES >0,"Upregulated", "Downregulated")
        
        p<-ggbarplot(out_write,
                     x = "Pathways",
                     y = "NES",
                     #fill = "sign",           # change fill color by mpg_level
                     color = "white",            # Set bar border colors to white
                     palette = "jco",            # jco journal color palett. see ?ggpar
                     sort.val = "asc",          # Sort the value in descending order
                     sort.by.groups = FALSE,     # Don't sort inside each group
                     x.text.angle = 90,          # Rotate vertically x axis texts
                     ylab = 'Normalized Enrichment Score',
                     legend.title = "padj < 0.25",
                     rotate = TRUE,
                     title = paste(pathway.name,"pathways in",sample,cluster),
                     ggtheme = theme_minimal(base_size = 15))+
                geom_col(aes(fill=padj<0.25))+
                guides(fill = guide_legend(reverse = TRUE))+
                theme(text = element_text(size=12),
                      plot.title = element_text(size = 12,hjust = hjust))
        path <- paste0("output/",gsub("-","",Sys.Date()),"/")
        if(!dir.exists(path)) dir.create(path, recursive = T)
        jpeg(paste0(path,sample,"_",cluster,"-",pathway.name,".jpeg"), units="in", width=10, height=7,res=600)
        print(p)
        dev.off()
}

 
(clusters <- as.character(unique(res$cluster)))
(cell.lines <- as.character(names(markers_list)))
for(cell.line in cell.lines){
        print(cell.line)
        results = markers_list[[cell.line]]
        for(i in 1:length(clusters)) FgseaBarplot(pathways=hallmark, stats=res, nperm=1000,show=NULL,
                                               cluster = clusters[i],sample = cell.line,
                                               pathway.name = "Hallmark")
        for(i in 1:length(clusters)) FgseaBarplot(pathways=allpathways, stats=res, nperm=1000,show=50,
                                               cluster = clusters[i],sample = cell.line,hjust=0,
                                               pathway.name = "Hallmark, biocarta,and KEGG")
}

FgseaDotPlot <- function(pathways=hallmark, stats=results, nperm=1000,
                    sample="",pathway.name = "Hallmark", padj = 0.25, pval=0.05,verbose=T){
        
        (clusters <- as.character(unique(stats$cluster)))
        fgseaRes <- list()
        for(i in 1:length(clusters)){
                res = stats[stats$cluster == clusters[i],]
                res1 = res[order(res["p_val_adj"]),c("gene","avg_logFC")]  %>% deframe
                fgseaRes[[i]] <- fgsea(pathways=pathways, stats=res1, nperm=nperm)
                fgseaRes[[i]] = as.data.frame(fgseaRes[[i]])
                fgseaRes[[i]] = fgseaRes[[i]][,c("pathway","pval","padj","NES")]
                fgseaRes[[i]] = fgseaRes[[i]][fgseaRes[[i]]$pval < pval,]
                fgseaRes[[i]] = fgseaRes[[i]][fgseaRes[[i]]$padj < padj,]

                fgseaRes[[i]]$cluster = clusters[i]
        }
        df_fgseaRes <- data.table::rbindlist(fgseaRes)
        if(verbose) print(round(dim(df_fgseaRes)/length(clusters)))
        plot<- ggballoonplot(df_fgseaRes, x = "cluster", y = "pathway",
                      size = "padj", fill = "NES",
                      size.range = c(1, 5),
                      font.xtickslab=14, 
                      font.ytickslab=round(50/dim(df_fgseaRes)[1]*14),
                      title = paste(pathway.name,"pathways in",sample),
                      xlab = "", ylab = "",
                      font.x =14,font.y= 14,font.main=16) +
                scale_fill_gradientn(colors = pal_gsea()(10))+
                theme(plot.title = element_text(hjust = 0.5))+
                scale_size(breaks=c(0,0.05,0.10,0.15,0.2,0.25),
                           labels=rev(c(0,0.05,0.10,0.15,0.2,0.25)))
        path <- paste0("output/",gsub("-","",Sys.Date()),"/")
        if(!dir.exists(path)) dir.create(path, recursive = T)
        jpeg(paste0(path,"Dotplot_",sample,"_",pathway.name,
                    "_",padj,"_",pval,".jpeg"), units="in", width=10, height=7,res=600)
        print(plot)
        dev.off()
}

cell.line = 'NSC'
results = markers_list[[cell.line]]
FgseaDotPlot(pathways=hallmark, stats=results, nperm=1000,padj = 0.25,pval = 0.05,
             sample = cell.line, pathway.name = "Hallmark")
FgseaDotPlot(pathways=allpathways, stats=results, nperm=1000,
             padj = 0.25,pval = 0.005,
             sample = cell.line, pathway.name = "Hallmark, biocarta,and KEGG")

cell.line = '827'
results = markers_list[[cell.line]]
FgseaDotPlot(pathways=hallmark, stats=results, nperm=1000,padj = 0.25,pval = 0.05,
             sample = cell.line, pathway.name = "Hallmark")
FgseaDotPlot(pathways=allpathways, stats=results, nperm=1000,
             padj = 0.25,pval = 0.15,
             sample = cell.line, pathway.name = "Hallmark, biocarta,and KEGG")

cell.line = '923'
results = markers_list[[cell.line]]
FgseaDotPlot(pathways=hallmark, stats=results, nperm=1000,padj = 0.25,pval = 0.05,
             sample = cell.line, pathway.name = "Hallmark")
FgseaDotPlot(pathways=allpathways, stats=results, nperm=1000,
             padj = 0.15,pval = 0.05,
             sample = cell.line, pathway.name = "Hallmark, biocarta,and KEGG")

