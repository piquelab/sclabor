library(Seurat)

library(tidyverse)

library(Matrix)

library(data.table)

library(DoubletFinder)

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



sl <- NormalizeData(scilab)

sl <- ScaleData(sl)

sl <- FindVariableFeatures(sl, selection.method = "vst", nfeatures = 2000)

sl <- RunPCA(sl)

sl <- RunUMAP(sl, dims = 1:20)

sl

## ## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
## sweep.res <- paramSweep_v3(sl, PCs = 1:20, sct = FALSE)

## sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)

## bcmvn <- find.pK(sweep.stats)

## ## so that is the pK????? 0.005

## bcmvn$pK[which.max(bcmvn$BCmetric)]


## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------


homotypic.prop <- modelHomotypic(sl@meta.data$NewClsName)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.01*ncol(sl))  ## Assuming 1.0% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
##sl <- doubletFinder_v3(sl, PCs = 1:20, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
sl <- doubletFinder_v3(sl, PCs = 1:20, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)

##save(sl,file="doubletFinder.sl3.Rdata")

load("doubletFinder.sl3.Rdata")

##seu_kidney <- doubletFinder_v3(seu_kidney, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)

table(sl@meta.data$NewClsName,sl@meta.data$DF.classifications_0.25_0.005_700)

table(sl@meta.data$DF.classifications_0.25_0.005_700)


umapdf <- cbind(as.data.frame(scilab@reductions$umap@cell.embeddings),
                sl@meta.data) 

umapdf$FinalName <- newNames[umapdf$NewClsName]

umapdf$Location <- factor(umapdf$Location)

levels(umapdf$Location)=c("BP","PV","CAM")

table(umapdf$Location)


##umapdf %>% filter(prediction.score.max>0.001) %>% select(predicted.id) %>% table
gg<- ggplot(umapdf %>% arrange(desc(DF.classifications_0.25_0.005_700))
            ## %>% filter(DF.classifications_0.25_0.005_1558=="Doublet")
           ,aes(x = UMAP_1, y = UMAP_2, col = DF.classifications_0.25_0.005_700)) + 
            geom_point( size = 0.25)  +
    scale_color_brewer(type="qual",palette=8) + 
            ##scale_colour_gradientn(colors = c('#ffffcc','#ffeda0','#fed976','#feb24c','#fd8d3c','#fc4e2a','#e31a1c','#b10026','#800026')) + 
#            facet_grid(Group ~ Location) +
##                facet_grid( Location ~ .) +
     guides(colour = guide_legend(override.aes = list(size=10))) +
##    labs(col = "Auto labels") + 
##            ggtitle(paste(s,gene_symbol[s])) +
            theme_bw()

ggsave("plotDoubletFinder.UMAP.pdf",gg,width=8,height=6)


tt <- umapdf %>% select(FinalName,DF=DF.classifications_0.25_0.005_700) %>% table()

freq <- tt[,1]/rowSums(tt) *100

freq <- data.frame(FinalName=rownames(tt),Perc=tt[,1]/rowSums(tt) *100)


group.colors<-c(B_cells="red",ColumnCyto="limegreen",Cytotroph="purple4",Decidual="sienna4",Dendr_Macro_A="hotpink",Dendr_Macro_B="tomato",Endometrial="maroon",Endothelial="yellow4",EVT="plum2",Fibroblasts="navy",HSC="magenta",Monocytes="lightslateblue",Progenitor="peru",NK_cells="royalblue",Stromal_A="turquoise4",Stromal_B="lightpink3",Synciotrophoblasts="yellowgreen",Tcells_activated="seagreen",Tcells_resting="powderblue")
names(group.colors) <- newNames 


cat("# Doublets %:",sum(tt[,1])/sum(tt)*100,"\n")

gg <- ggplot(freq,aes(x=FinalName,y=Perc, fill = FinalName)) + geom_bar(stat = "identity") +
        scale_fill_manual(values=group.colors) +
    geom_hline(yintercept=sum(tt[,1])/sum(tt)*100,lty=3) +
    coord_flip() +
        labs(y = "% Doublets in cell-type", x = "Cell-type") +
    theme_bw()

ggsave("doubletBarplot.pdf",gg,width=5,height=7)
