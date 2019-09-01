## options(repos = c(CRAN = "http://cran.rstudio.com"))
##install.packages("tidyverse")
library("rhdf5")
#library(cellrangerRkit)
library(Seurat)
library(Matrix)
library(tidyverse)

## Load the PBMC dataset

basefolder<-"/nfs/rprdata/scilab/novogene/Analyses/ED2/"

outFolder <- paste0(basefolder,"/outputs1/")
system(paste0("mkdir -p ",outFolder))

#setwd(outFolder)
setwd(basefolder)
##scaip.cv_full <- readr::read_rds("/nfs/rprdata/ALOFT/AL/sc.genotypes/scaip.cv_full.rds")
##sc.cv <- readr::read_tsv("/nfs/rprdata/ALOFT/gencove/SCAIP-all64-merge.txt")

#dfile <- paste0(basefolder,"demux_sel.rds")
#demuxlet <- readr::read_rds(dfile)

afile <- paste0(basefolder,"/80K-all.rds")
all <- readr::read_rds(afile)

## open gene names: (TO DO, make a script to retrieve annotations)
gene_anno <- read_tsv("/nfs/rprdata/our10xData/counts/SCAIP-new/SCAIP1-ctrl/outs/filtered_gene_bc_matrices/hg19/genes.tsv",col_names=c("ensg.id","symbol"))
gene_symbol <- gene_anno$symbol
names(gene_symbol) <- gene_anno$ensg.id


## Symbol is not unique...

##If not identical we may need to re-order columns. 
stopifnot(identical(colnames(all),demuxlet$NEW_BARCODE))

sparse.size <- object.size(all)
sparse.size



#Keep all genes expressed in >= 3 cells, keep all cells with >= 200 genes 
#Perform log-normalization, first scaling each cell to a total of 1e4 molecules
#(as in Macosko et al. Cell 2015) 
#pbmc <- Setup(pbmc, min.cells = 3, min.genes = 200, 
#  do.logNormalize = T, total.expr = 1e4, project = "10X_PBMC_VitA_EtOH")
#setup setting do.scale and do.center to F - this means that we will NOT scale
#genes by default (to speed things up)
scilab <- CreateSeuratObject(all, 
                            min.cells = 3, 
                            min.genes = 200, 
                            project =  "novogene", 
                            do.scale = F, 
                            do.center = F, 
                            names.field = 2,
                            names.delim = "\\_")


#nGene and nUMI are automatically calculated for every object by Seurat. For
#non-UMI data, nUMI represents the sum of the non-normalized values within a
#cell We calculate the percentage of mitochondrial genes here and store it in
#percent.mito using the AddMetaData. The % of UMI mapping to MT-genes is a
#common scRNA-seq QC metric. 
# We may be able to get more MT genes as prefix does not capture everything.

## 
#gene_anno <- read_tsv("/nfs/rprdata/our10xData/counts/SCILAB-new/SCILAB1-ctrl/outs/filtered_gene_bc_matrices/hg19/genes.tsv",col_names=c("ensg.id","symbol"))
gene_anno <- gene_anno %>% filter(ensg.id %in% rownames(scilab@data))
identical(gene_anno$ensg.id,rownames(scilab@data))


mito.genes <- grepl("^MT-", gene_anno$symbol)
percent.mito <- colSums(expm1(scilab@data[mito.genes, ])) / colSums(expm1(scilab@data))

#AddMetaData adds columns to object@data.info, and is a great place to stash QC stats
Code<-read_csv("Code_ALL.csv")

aux <- Code %>% filter(Barcode %in% row.names(scilab@meta.data)) %>% select(-Barcode) %>% as.data.frame

rownames(aux) <-  row.names(scilab@meta.data)





#scilab <- AddMetaData(scilab, percent.mito, "percent.mito")

## prepare metadata from demuxlet
#aux <- (demuxlet %>% select(BARCODE,RD.UNIQ,BATCH,TREAT,EXP,SNG.1ST,type,PRB.DBL,PRB.SNG1) %>% as_data_frame())
#rownames(aux) <- demuxlet$NEW_BARCODE


##aux <- (demuxlet %>% select(BARCODE,RD.UNIQ,BATCH,EXP,SNG.1ST,type,PRB.DBL,PRB.SNG1) %>% as_data_frame()) 
##rownames(aux) <- demuxlet$NEW_BARCODE

scilab <- AddMetaData(scilab, aux)

str(scilab@meta.data)

setwd(outFolder)

pdf("histograms.pdf")
hist(log10(scilab@meta.data$nGene),breaks=100)
#hist(scilab@meta.data$percent.mito,breaks=100)
hist(log10(scilab@meta.data$nUMI),breaks=100)
dev.off()
##VlnPlot(scilab, c("nGene", "nUMI", "percent.mito"), nCol = 3)

#qplot(nUMI,nGene,data=scilab@meta.data,col=TREAT)
#qplot(nUMI,nGene,data=scilab@meta.data,col=BATCH,alpha=0.01)
##I could find the upper limit here for the nGene/nUMI, but importing the PRB.doublet and singlet 

#qplot(nUMI,nGene,data=scilab@meta.data,col=type,log="xy",alpha=0.01)

#We filter out cells that have unique gene counts under 500, and > 5% mitochondrial percentage

sum(scilab@meta.data$percent.mito>0.01)
sum(scilab@meta.data$nUMI<500)
sum(scilab@meta.data$percent.mito>0.05)
sum(scilab@meta.data$nGene<200)

scilab <- SubsetData(scilab, subset.name = "percent.mito", accept.high = 0.01)
##scilab <- SubsetData(scilab, subset.name = "nGene", accept.low = 300)
##scilab <- SubsetData(scilab, subset.name = "PRB.SNG1", accept.low = 0.99)
##scilab <- SubsetData(scilab, subset.name = "numsing", accept.low = 100)
## Only accept singlets form indivs with more than 100 singlets. may need to recalculate...


## maybe filter things with too many genes... ??
##pbmc33k <- SubsetData(pbmc33k, subset.name = "nGene", accept.high = 2500)

## maybe select things with high probability of being SNG. 
## scilab <- SubsetData(scilab, subset.name = "type", ident.use="SNG")

#subset to 16 individuals with >100 cells




### normalize data, very fast
scilabFilt <- NormalizeData(scilab, normalization.method = "LogNormalize", 
                       scale.factor = 10000)

## this step is fast too. 
## Detection of variable genes across the single cells
scilabFilt <- FindVariableGenes(scilabFilt, mean.function = ExpMean, dispersion.function = LogVMR, 
#                           x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5,
                          x.low.cutoff = 0.0125, x.high.cutoff = 5, y.cutoff = 0.5,
                           do.plot=FALSE)

length(scilabFilt@var.genes)

# hist(log10(scilab@hvg.info$gene.mean),breaks=100)
# hist((scilab@hvg.info$gene.dispersion.scaled),breaks=100)
# VariableGenePlot(scilab,do.text=FALSE)
# abline(v=0.01)
# abline(v=5)
# abline(h=0.4)

#  filter out mito and ribo (MT,RP) genes from scilab@var.genes
rpg<-grepl("^RP",gene_anno$symbol)
gene_anno2<-gene_anno[!rpg,]
mtg<-grepl("^MT",gene_anno2$symbol)
gene_anno3<-gene_anno2[!mtg,]
nvg<-intersect(scilabFilt@var.genes, gene_anno3$ensg.id)
length(scilabFilt@var.genes)
length(nvg)

scilabFilt@var.genes<-nvg



var_genes <- gene_anno %>% filter(ensg.id %in% scilabFilt@var.genes) %>% select(symbol) %>% unlist()

## This step I think it is slow. It takes hours
## Scaling the data and removing unwanted sources of variation
##scilab <- ScaleData(scilab, vars.to.regress = c("nUMI", "percent.mito"))
#I'm keeping the basic thing, but not sure I need all this. 
# maybe we should take out also the exp and indiv effect?
## Removing exp and indiv, EXP includes batch + treatment, SNG.1ST includes individual effect.
## Takes longer, so increased block.size to see if everything can be done in one block. 
## Not sure what may be the best model.use

#scilabFilt <- ScaleData(scilabFilt, vars.to.regress = c("nUMI","SNG.1ST","TREAT"),
                   #                   model.use = "negbinom",
#                   block.size=30000)


hv.genes<-head(rownames(scilabFilt@hvg.info),1000)
scilabFilt<-ScaleData(scilabFilt,genes.use=hv.genes,display.progress=FALSE,vars.to.regress=c("nUMI"),do.par=TRUE,num.cores=14)


## Save here
fname <- paste0("novogene_scaled_",Sys.Date(),".rds")
fname
write_rds(scilabFilt,fname)
##scilab = read_rds(paste0("scaip_20180111_v0.rds"))

## May be try if doing variable gene selection after makes a difference. 
## It does not. 
### END- HERE ###
########################################################


