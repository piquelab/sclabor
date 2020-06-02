library(Seurat)

##library(tidyverse)

##library(Matrix)

##library(data.table)

library(scater)

library(slingshot)

load("/nfs/rprdata/scilab/novogene/Analyses/ED6/PrimaryAnalysis/CreateSeuratObject/outputs1/scilabObject-UMAP.Rdata")

gene_anno <- readr::read_tsv("/nfs/rprdata/our10xData/counts/SCAIP-new/SCAIP1-ctrl/outs/filtered_gene_bc_matrices/hg19/genes.tsv",col_names=c("ensg.id","symbol"))
gene_symbol <- gene_anno$symbol                                                  
names(gene_symbol) <- gene_anno$ensg.id

scilab <- UpdateSeuratObject(scilab)

newNames = c("B_cells"="B-cell"
            ,"ColumnCyto"="npiCTB"
            ,"Cytotroph"="CTB"
            ,"Decidual"="Decidual"
            ,"Dendr_Macro_A"="Macrophage-1"
            ,"Dendr_Macro_B"="Macrophage-2"
            ,"Endometrial"="Endometrial"
            ,"Endothelial"="LED"
            ,"EVT"="EVT"
            ,"Fibroblasts"="Fibroblast"
            ,"HSC"="HSC"
            ,"Monocytes"="Monocyte"
            ,"(Myeloid)Progenitor"="Stromal-3"
            ,"NK_cells"="NK-cell"
            ,"Stromal_A"="Stromal-1"
            ,"Stromal_B"="Stromal-2"
            ,"Synciotrophoblasts"="STB"
            ,"Tcells_activated"="T-cell-activated"
            ,"Tcells_resting"="T-cell-resting")


scilab@meta.data$FinalNames <- newNames[scilab@meta.data$NewClsName]


Idents(scilab) = "FinalNames"
##subset(pbmc, subset = replicate == "rep2")


troph = subset(scilab, idents = c("STB", "EVT", "CTB", "npiCTB"))


## Maybe not needed or embed from PC to UMAP.
sce <- Seurat::as.SingleCellExperiment(troph)

sce <- slingshot(sce, clusterLabels = 'FinalNames', reducedDim = 'PCA')

ssds <- SlingshotDataSet(sce)

##c1 <- getCurves(ssds)



sds <- slingshot(Embeddings(troph,"umap") , clusterLabels = troph$FinalNames,
                                  start.clus = "CTB", stretch = 0)


##save(troph,sds,sds2,file="Slingshot.Rdata")


group.colors<-c(B_cells="red",ColumnCyto="limegreen",Cytotroph="purple4",Decidual="sienna4",Dendr_Macro_A="hotpink",Dendr_Macro_B="tomato",Endometrial="maroon",Endothelial="yellow4",EVT="plum2",Fibroblasts="navy",HSC="magenta",Monocytes="lightslateblue",Progenitor="peru",NK_cells="royalblue",Stromal_A="turquoise4",Stromal_B="lightpink3",Synciotrophoblasts="yellowgreen",Tcells_activated="seagreen",Tcells_resting="powderblue")

names(group.colors) <- newNames 


pdf("./SlingshotTrajectory.pdf")
plot(reducedDim(sds), col = group.colors[troph$FinalNames], pch = 16, cex = 0.5)
lines(sds, lwd = 2, type = 'lineages', col = 'black')
##lines(sds, lwd = 2, col = 'black',)
dev.off()


##res.0.7 , start 3 blue, or 15 yellow


sds2 <- slingshot(Embeddings(troph,"umap") , clusterLabels = troph$res.0.7,
                                  start.clus = 3, stretch = 0)





mycol = RColorBrewer::brewer.pal(8, "Set1")
names(mycol)= unique(troph$res.0.7)
mycol

table(troph$res.0.7,troph$FinalNames)

pdf("./SlingshotTrajectory.v2.pdf")
plot(reducedDim(sds2), col = mycol[troph$res.0.7], pch = 16, cex = 0.5)
lines(sds2, lwd = 2, type = 'lineages', col = 'black')
##lines(sds, lwd = 2, col = 'black',)
dev.off()

pdf("./SlingshotTrajectory.v3.pdf")
plot(reducedDim(sds2), col = group.colors[troph$FinalNames], pch = 16, cex = 0.5, xlim=c(-10,6),ylim=c(3,15))
lines(sds2, lwd = 2, type = 'lineages', col = '#000000AA')
##lines(sds2, lwd = 2, col = 'black',)
dev.off()


