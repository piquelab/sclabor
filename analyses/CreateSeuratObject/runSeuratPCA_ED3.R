

## options(repos = c(CRAN = "http://cran.rstudio.com"))
##install.packages("tidyverse")
library("rhdf5")
#library(cellrangerRkit)
library(Seurat)
library(Matrix)
library(tidyverse)

## Load other objects, see runScaipSeurat.R
## This script follows after runScaipSeurat.R

selDate <- "2018-06-05"
basefolder <- "/nfs/rprdata/scilab/novogene/Analyses/ED2/"


outFolder <- paste0(basefolder,"outputs1/")
##system(paste0("mkdir -p ",outFolder))
setwd(outFolder)



## Load seurat object, takes ~30 secs to load 
scilab <- read_rds(paste0(outFolder,"novogene_scaled_",selDate,".rds"))

## Perform linear dimensional reduction, takes 1 min
scilab <- RunPCA(scilab, pc.genes = scilab@var.genes, 
                do.print = TRUE, pcs.print = 1:5, genes.print = 5)

## # Examine and visualize PCA results a few different ways
#PrintPCA(scilab, pcs.print = 1:5, genes.print = 5, use.full = FALSE)


## Move to the plot function?
pdf("bowplot.pdf")
PCElbowPlot(object = scilab)
dev.off()


pdf(paste0("vixPCA",selDate,".sng.pdf"),width=12,height=12)
VizPCA(scilab, pcs.use = 1:16,nCol=4,num.genes=10)
dev.off()

pdf(paste0("PCHeatmap",selDate,".sng.pdf"),width=12,height=20)
PCHeatmap(scilab, pc.use = 1:20, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
dev.off()

pdf("PCAplot.pdf")
PCAPlot(object=scilab,dim.1=1,dim.2=2)
dev.off()


##PCAPlot(scilab, dim.1 = 1, dim.2 = 2)
##pbmc <- ProjectPCA(object = pbmc, do.print = FALSE)
##PCHeatmap(scilab, pc.use = 16, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)
##PCHeatmap(scilab, pc.use = 1:12, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)

##JackStrawPlot(object = pbmc, PCs = 1:12)


## Run Non-linear dimensional reduction (tSNE)
## takes a bit of time ~10 min? 
scilab <- RunTSNE(scilab, genes.use=scilab@var.genes, do.fast = TRUE)

##I could save and do some plot here before we take it to the clustering thing. 


## Takes some time.... (separate to another script.)
## Cluster the cells
                                        #
#scilab <- FindClusters(scilab, reduction.type = "pca", genes.use=scilab@var.genes, resolution = 0.7, print.output = 0, save.SNN = TRUE)

scilab<-FindClusters(scilab,reduction.type="pca",dims.use=1:20,,resolution=0.7,print.output=0,save.SNN=TRUE)


scilab<- StashIdent(scilab, save.name = "Clust_genes_RMrem_res0.7")

PrintFindClustersParams(object = scilab)

## Saving at this point
prefix="NovoGene_R3__seurat_genes_Res0.7"
fname = paste0(prefix,selDate,".rds")
write_rds(scilab,fname)
fname




