## options(repos = c(CRAN = "http://cran.rstudio.com"))
##   This uses updated Seurat package 3 - starts with merged counts/demux from step 2

library(Seurat)
library(Matrix)
library(tidyverse)

library(future)

future::plan(strategy = 'multicore', workers = 16)
options(future.globals.maxSize = 20 * 1024 ^ 3)

###########################################
## Testing sc transform           
## 2and3_Diem_Output
## adata <- read_rds("./kb_diem_Output/kb_diem_Seurat.list.rds")

sc <- read_rds("./kb_diem_Output/kb_diem_Seurat.obj.rds")

sc@meta.data$Library <- (gsub("_[ACGT]{6,}","",colnames(sc)))

expNames <- unique(sc@meta.data$Library)
expNames

outFolder= "./2_ST_integrate_Output_v2/"
system(paste0("mkdir -p ",outFolder))
setwd(outFolder)


## 
##Maybe open demuxlet here... or ignore it???
demux <- read_rds("../1_demuxletOutput/demuxlet_all.rds")

demux <- mutate(demux,BARCODE2 = paste0(BATCH,"_",gsub("-.*","",BARCODE)))

selcells <- intersect(demux$BARCODE2,colnames(sc))

demux <- filter(demux,BARCODE2 %in% selcells)
dim(demux)

sc <- sc[,selcells]

aux <- sc@meta.data %>% rownames_to_column(var="BARCODE2") %>% left_join(demux) %>% column_to_rownames("BARCODE2")
stopifnot(identical(rownames(aux),colnames(sc)))
sc@meta.data <- aux

##
sc <- SCTransform(sc, verbose = TRUE)


## These are now standard steps in the Seurat workflow for visualization and clustering

sc <- RunPCA(object = sc, verbose = TRUE, npcs= 100)

######
##setwd("./2and3_Kallisto_Diem_Output/")

pdf("elbowPlot.pdf")
ElbowPlot(sc, ndims = 100)
dev.off()

sc <- RunUMAP(sc, dims = 1:50, verbose = TRUE)

sc <- FindNeighbors(sc, dims = 1:50, verbose = TRUE)

sc <- FindClusters(sc, verbose = TRUE)

##

write_rds(sc,"ST_Integrated.obj.rds")

##
##
pdf("UMAP.pdf")
DimPlot(sc, label = TRUE) + NoLegend()
dev.off()
##

## pdf("UMAP_batch.pdf")
## DimPlot(sc,group.by="BATCH")
## dev.off()

## ##
## pdf("UMAP_MvsF.pdf")
## DimPlot(sc,group.by="TYPE.BEST.GUESS")
## dev.off()

pdf("UMAP_Location_NewCls.pdf")
aa <- FetchData(sc,c("UMAP_1","UMAP_2","seurat_clusters","Location"))
ggplot(aa,aes(UMAP_1,UMAP_2,color=seurat_clusters)) +
    geom_point(size=0.1) +
##    facet_wrap("Location",2,2) +
    theme_bw()
dev.off()


### END- HERE ###
########################################################


