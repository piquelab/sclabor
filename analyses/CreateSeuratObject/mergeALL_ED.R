## options(repos = c(CRAN = "http://cran.rstudio.com"))
##install.packages("tidyverse")
library("rhdf5")
library(cellrangerRkit)
library(Seurat)
library(Matrix)
library(tidyverse)

## Load the PBMC dataset
basefolder <- "/nfs/rprdata/scilab/novogene/aggr/AGGR"

anFolder<-"/nfs/rprdata/scilab/novogene/Analyses/ED2"

auxUF<-load_cellranger_matrix_h5(basefolder,genome="hg19",barcode_filtered=FALSE)



#a <- do.call(cbind,aux3)

##If not identical we may need to re-order columns. 
#stopifnot(identical(colnames(a),demuxlet$NEW_BARCODE))
libID<-read_csv("LibraryID.csv")
aux80K <- auxUF[,libID$Barcode]

eaux80K<-exprs(aux80K[,])

#a<-do.call(cbind,eaux80K)


afile <- paste0(anFolder,"/80K-all.rds")
write_rds(eaux80K,afile)


sparse.size <- object.size(aux3)
sparse.size

