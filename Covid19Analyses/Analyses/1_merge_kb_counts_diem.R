## options(repos = c(CRAN = "http://cran.rstudio.com"))
##install.packages("tidyverse")
## library("rhdf5")

library(Seurat)

library(Matrix)

library(tidyverse)

library(diem)

## Load the datasets
##basefolder <- "/nfs/rprdata/scilab/cc1-Placenta/counts/ALL-recount-6.17.19-preRNA/"
##basefolder <- "/nfs/rprdata/scilab/endometrium/kallisto2/bus/"
##expNames <- dir(basefolder,"^s*")
basefolder <- "/nfs/rprdata/scilab/cc1-Placenta/kallisto2/bus/"
expNames <- dir(basefolder,"^CC[78]-*")


expNames

##subfolder <- "/outs/filtered_gene_bc_matrices/hg19/"
##subfolder <- "/outs/raw_gene_bc_matrices/hg19/"
##subfolder <- "/outs/"
subfolder <- ""
folders <- paste0(basefolder,expNames,subfolder)
ind <- file.info(folders)$isdir
ind[is.na(ind)]<- FALSE
folders <- folders[ind]
expNames <- expNames[ind]
names(folders) <- expNames
folders

##aa2 <- ReadH5AD(file = "../kallisto2/bus/HPL20289_R_C_1/counts_filtered/adata.h5ad",verbose=TRUE)

adata <- map(expNames,function(ii){
    ##
    expPrefix = ii;
    folders[expNames[1]]
    fName= paste0(folders[ii],"/counts_unfiltered/adata.h5ad")
    cat("#Loading ",fName,"... \n")
    adata = ReadH5AD(file = fName,verbose=TRUE)
    cat(dim(adata),"\n")
    ## May need to rename rownames or split...
    sce <- create_SCE(adata@assays$RNA@data)
    cat(dim(sce),"\n")
    ## Remove debris...
    sce <- diem(sce)
    ##
    sc = convert_to_seurat(sce)
    cat("#Final: ",dim(sc),"\n")
    adata <- adata[,colnames(sc)]
    ## I could keep the background debris model, or also return cells not debris
    ##  adata <- SCTransform(adata, verbose = TRUE)
    adata
})


names(adata) <- expNames


system("mkdir -p kb_diem_Output/")
write_rds(adata,"kb_diem_Output/kb_diem_Seurat.list.rds")

sparse.size <- object.size(adata)
sparse.size


sc <- merge(adata[[1]], y = adata[-1], project = "KallistoDiem",add.cell.ids=expNames)

dim(sc)

system("mkdir -p kb_diem_Output/")
write_rds(sc,"kb_diem_Output/kb_diem_Seurat.obj.rds")


##############################
