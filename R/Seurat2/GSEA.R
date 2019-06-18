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

# 5.1 load data ==============
(load(file="data/Glioblastoma_12_20190214.Rda"))

##############################
# Split PrepareGSEA
###############################
Split_object <- SplitSeurat(object, split.by = "cell.lines")
names(Split_object)
remove(object);GC()
for(d in c("827","923")){
        print(continuous.label <- paste0(d,"-D",c(14,16,28,0)))
        PrepareGSEA(Split_object[[d]], k = 100,continuous.label = continuous.label)
}

# Run GSEA and generate reports
#'@example ReportGSEA(file = "c5.all.827-D14.827-D16.827-D28.827-D0.Gsea.1552707417126",pos=T,ncol=3)
ReportGSEA <- function(file, pos=T,ncol = 3){
        (gsea_path <- paste("~/gsea_home/output",tolower(format(Sys.Date(), "%b%d")),
                            file,sep ='/'))
        #(pos.xls.path <- list.files(gsea_path,pattern="gsea_report_for_.*pos.*xls"))
        p<-1; if(pos) p <-2
        (pos.xls.path <- list.files(gsea_path,pattern="gsea_report_for_.*xls")[p])
        GSEA_output <- readr::read_delim(paste0(gsea_path,"/",pos.xls.path),"\t", 
                                         escape_double = FALSE, trim_ws = TRUE)
        print(GSEA_output %>% .[-c(2,3,12)] %>% head(50) %>% kable() %>% kable_styling())
        
        (GSEA.plots <- sapply(GSEA_output$NAME[1:9], function(name) {
                paste0("enplot_",name, "_([0-9]+)*\\.png$")}) %>%
                        sapply(function(x) list.files(path = gsea_path, pattern =x)) %>%
                        .[sapply(.,length)>0] %>% #Remove empty elements from list with character(0)
                        paste(gsea_path, ., sep = "/")) 
        CombPngs(GSEA.plots, ncol = ncol)
        path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
        if(!dir.exists(path)) dir.create(path, recursive = T)
        file.rename(paste0(path,list.files(path,pattern=paste0("GSEA.plots_CombPngs"))), 
                    paste0(path,file,"-",p,".jpeg"))
        
}
ReportGSEA(file = "c5.all.827-D14.827-D16.827-D28.827-D0.Gsea.1552707417126",pos=T,ncol=3)
