## options(repos = c(CRAN = "http://cran.rstudio.com"))
##   This uses updated Seurat package 3 - starts with merged counts/demux from step 2

library(Seurat)

library(harmony)

library(Matrix)
library(tidyverse)

library(annotables)

library(future)

future::plan(strategy = 'multicore', workers = 16)
options(future.globals.maxSize = 30 * 1024 ^ 3)


###########################################
## Testing sc transform           
## 2and3_Diem_Output
## adata <- read_rds("./kb_diem_Output/kb_diem_Seurat.list.rds")

sc4 <- read_rds("/nfs/rprdata/scilab/cc1-Placenta/analyses/nuc.20200425/2_ST_integrate_Output_v2/ST_Integrated.obj.rds")

sc4 <- subset(sc4, subset = nFeature_RNA > 100)


sc3 <- read_rds("/nfs/rprdata/scilab/novogene/Analyses/Roger_20200218/2_ST_integrate_Output_v2/ST_Integrated.obj.rds")

sc3 <- subset(sc3, subset = nFeature_RNA > 100)


sc2 <- read_rds("/nfs/rprdata/scilab/endometrium/Analysis/20200223/2_ST_integrate_Output_v2/ST_Integrated.obj.rds")

sc2 <- subset(sc2[,sc2$Location=="PVBP"], subset = nFeature_RNA > 100)


sc1 <- read_rds("/nfs/rprdata/scilab/novogene/otherdata/roser/Analysis/20200226/3_scTransferLabel_RoserPrep/ST_Integrated_RVT.obj.rds")

sc1 <- subset(sc1[,sc1$Location!="Blood"], subset = nFeature_RNA > 100)


##geneList <- intersect(intersect(intersect(rownames(sc1),rownames(sc2)),rownames(sc3)),rownames(sc4))

#sc4 <- sc4[geneList,]

#sc3 <- sc3[geneList,!is.na(sc3$FinalName)]

#sc2 <- sc2[geneList,sc2$Location=="PVBP"]

#sc1 <- sc1[geneList,!is.na(sc1$annotation)]

sc1@meta.data$Trimester="1st"
sc2@meta.data$Trimester="2nd"
sc3@meta.data$Trimester="3rd"
sc4@meta.data$Trimester="3rd"

sc4@meta.data$Location="Nucleus"


outFolder="./covid_merge_harmony/"
system(paste0("mkdir -p ", outFolder))
setwd(outFolder)

sc <- merge(sc1,list(sc2,sc3,sc4))

dim(sc)

table(sc$Library)

table(sc$Location) 


## Harmony

DefaultAssay(sc) <- "RNA"

sc <- NormalizeData(sc, verbose=TRUE) 

sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 2000)

sc <- ScaleData(sc, verbose = TRUE) 

sc <- RunPCA(sc,pc.genes = sc@var.genes, npcs = 50, verbose = TRUE)

sc <- RunHarmony(sc,"Location",reduction="pca")

sc <- RunUMAP(sc,reduction = "harmony", dims = 1:30)


###### Cluster

sc <- FindNeighbors(sc, reduction = "harmony", dims = 1:30, verbose = TRUE)

sc <- FindClusters(sc, verbose = TRUE,resolution=0.6)


######   TRANSFER LABELS  #####

scLaborTransAnchors <- FindTransferAnchors(reference = sc, query = sc, dims = 1:30) ## try harmony reduction? Not sure how to do it... 

pred.scLabor <- TransferData(anchorset = scLaborTransAnchors, refdata = sc$FinalName, dims=1:30)

summary(pred.scLabor$predicted.id==sc$FinalName)

table(sc$Trimester, pred.scLabor$prediction.score.max>0.01)

table(sc$Location, pred.scLabor$predicted.id)

sum(pred.scLabor$prediction.score.max>0.01)

sc@meta.data$Pred.Id.scLaborRef=pred.scLabor$predicted.id

sc@meta.data$Pred.Score.scLaborRef=pred.scLabor$prediction.score.max


##  Using Roser ref. 
pred.scRVT <- TransferData(anchorset = scLaborTransAnchors, refdata = sc$final_cluster, dims=1:30)

summary(pred.scRVT$predicted.id==sc$annotation)

table(sc$Trimester, pred.scRVT$prediction.score.max>0.1)

table(sc$Location, pred.scRVT$prediction.score.max>0.01)

sum(pred.scRVT$prediction.score.max>0.01)

sc@meta.data$Pred.Id.scRVT=pred.scRVT$predicted.id

sc@meta.data$Pred.Score.scRVT=pred.scRVT$prediction.score.max


###

write_rds(sc,"sc.Harmonized.NormByLocation.rds")

## Label transfer seems to not be owrking at all.


####
## Remove in a separate script what follows.
## use script to plot..

### END- HERE ###
########################################################


