library(tidyverse)
library(knitr)
library(DESeq2)
library(annotables)
library(qqman)


basefolder <- "/nfs/rprdata/scilab/novogene/Analyses/ED3/"

outFolder <- paste0(basefolder,"/outputs_DEseq_by_Location-byCluster2/")
system(paste0("mkdir -p ",outFolder))
setwd(outFolder)



#convert annos to gene names
gene_anno <- read_tsv("/nfs/rprdata/our10xData/counts/SCAIP-new/SCAIP1-ctrl/outs/filtered_gene_bc_matrices/hg19/genes.tsv",col_names=c("ensg.id","symbol"))
gene_symbol <- gene_anno$symbol
names(gene_symbol) <- gene_anno$ensg.id


afile <- paste0(basefolder,"../ED2/80K-all.rds")
all <- readr::read_rds(afile)



clusterfile <- paste0("../../ED2/MergedClusters/newMergedClustersData.rds")
clusters <- readr::read_rds(clusterfile)

clusters$NewClsName <- as.character(clusters$NewClsName)


cSize <- clusters %>% dplyr::count(Location,NewClsName) %>%
    group_by(NewClsName) %>%
    summarize(minc=min(n)) %>%
    filter(minc>100) %>% arrange(-minc)
cSize


CellTypes=cSize$NewClsName
CellTypes

library(qqman)
cat("Celltype","\t","CAM_v_DB","\t","PV_v_DB","\t","PV_v_CAM","\n",sep="",file="QuickFDR10.txt",append=FALSE)
for (i in CellTypes[c(1,2,3,4,6,7,8,9,10)]){  #exlude (4) Progenitors, bec. 0 cells in some locations...
##resList <- lapply(CellTypes,function(i){
    ##    
    ## subset cell-types that I want to use.
    clusters_i  <- filter(clusters,NewClsName %in% c(i))
    dim(clusters_i)    
    all_i <- all[,clusters_i$NEW_BARCODE]
    ##
    ## Double check columns are matching ok
    ##stopifnot(identical(demuxlet$NEW_BARCODE,colnames(all)))
    ##
    ## Building a design matrix to get aggregate
    ## RPR: I'm changing this to aggregate across locations. 
    ## EDS:  Now, aggregating across groups
    
    bti <- clusters_i %>% transmute(bti=paste(Location,Indiv,sep="_")) %>% unlist %>% factor
    ##
    X <- model.matrix( ~ 0 + bti)
    qr.X <- qr(X)
    qr.X$rank
    dim(X)
    ##To aggregate the data $\bm X' \bm y_j$
    YtX <- all_i %*% X
    ##    
    YtX <- as.matrix(YtX)
    dim(YtX)
    ##
    ## Running DEseq2 on the aggregate
    ## ========================================================
    ##
    bti2 <- gsub("bti","",colnames(YtX))
    colnames(YtX) <- bti2
    ##
    cmat<-YtX
    ##
    anno <- tibble(ensgene=rownames(cmat),rs=rowSums(cmat),nz=rowSums(cmat>0)) %>%
        inner_join(grch37) %>%
        filter(chr %in% 1:22,rs>10,nz>3) ## keep only autosomal
    ##
    table(anno$chr)
    ##
    genesel <- (unique(anno$ensgene))
    cmat <- cmat[genesel,]
    dim(cmat)
    ##
    ##Create sample table
    cn<-colnames(cmat)
    x<-strsplit(cn,"_")
    ##
    cvt <- data.frame(matrix(unlist(x), nrow=length(x), byrow=T),stringsAsFactors=FALSE)
    colnames(cvt)<-c("Location","Indiv")
    ##
    #  ****         cvt$Group = relevel(factor(cvt$Group),"TL")
    ##
    ##
    ## Running DEseq2 on the aggregate
    ## ========================================================
    ##
    dds <- DESeqDataSetFromMatrix(cmat,cvt, ~Indiv+Location)    
    dds <- DESeq(dds,parallel=TRUE)
    ##
    ## save DEseq object. 
    ##fname <- paste0("scilab_group3_deseq_",Sys.Date(),".rds")
    ##fname
    ##write_rds(dds,fname)
    ##
    ##
    ## Parse the results 
    ## ========================================================
    ##
    ## DB v W
    resDW <- results(dds)
    ##
    resDW<-results(dds,contrast=c("Location","DB","W"))
    resDW$Gene<-gene_symbol[rownames(resDW)]
    write.table(resDW,paste0("CL_",gsub("/","Or",i),"_W_v_DB_ALL",Sys.Date(),".txt"))
    fname=paste0("CL_",gsub("/","Or",i),"W_v_DB_QQ.pdf")
    pdf(fname,width=12,height=12)
    qqman::qq(resDW$pvalue,main=fname,cex.main=3)
    dev.off()
    ##
    ##
    ## DB v P
    resDP <- results(dds)
    resDP <-results(dds,contrast=c("Location","DB","P"))
    resDP$Gene<-gene_symbol[rownames(resDP)]
    write.table(resDP,paste0("CL_",gsub("/","Or",i),"_P_v_DB_ALL",Sys.Date(),".txt"))
    fname=paste0("CL_",gsub("/","Or",i),"P_v_DB_QQ.pdf")
    pdf(fname,width=12,height=12)
    qqman::qq(resDP$pvalue,main=fname,cex.main=3)
    dev.off()
    ##
    ## W v P
    resWP <- results(dds)
    resWP <-results(dds,contrast=c("Location","W","P"))
    resWP$Gene<-gene_symbol[rownames(resWP)]
    write.table(resWP,paste0("CL_",gsub("/","Or",i),"_P_v_W_ALL",Sys.Date(),".txt"))
fname=paste0("CL_",gsub("/","Or",i),"W_v_P_QQ.pdf")
    pdf(fname,width=12,height=12)
    qqman::qq(resWP$pvalue,main=fname,cex.main=3)
    dev.off()

    
    cat(i,"\t",sum(resDW$padj<0.1,na.rm=TRUE),"\t",sum(resDP$padj<0.1,na.rm=TRUE),"\t",sum(resWP$padj<0.1,na.rm=TRUE),"\n",sep="")
    cat(i,"\t",sum(resDW$padj<0.1,na.rm=TRUE),"\t",sum(resDP$padj<0.1,na.rm=TRUE),"\t",sum(resWP$padj<0.1,na.rm=TRUE),"\n",sep="",file="QuickFDR10.txt",append=TRUE)
    ##
##    list(resTL=resTL,resPL=resPL)
    ##    
}

# Parse results - compare DE genes per group per 3 classes



#  ****************************************************************
  # Parse res files to get # of genes with padj <.1, FC>1
# Parse results - compare DE genes per group per 3 classes

cat("Celltype","\t","PV v DB","\t","PV v CAM","\t","CAM v DB","\n",sep="",file="DEgenes_padj_0.1_FC_1.txt",append=FALSE)

date="2018-11-07"

for (i in CellTypes[c(1:4,6:10)]){

#i=CellTypes[1]
  resPD<-read.table(paste0("CL_",gsub("/","Or",i),"_P_v_DB_ALL",date,".txt"))
  resPW<-read.table(paste0("CL_",gsub("/","Or",i),"_P_v_W_ALL",date,".txt"))
  resWD<-read.table(paste0("CL_",gsub("/","Or",i),"_W_v_DB_ALL",date,".txt"))

pt=0.1

PDs<-nrow(resPD %>% filter(padj<pt,abs(log2FoldChange)>1))
PWs<-nrow(resPW %>% filter(padj<pt,abs(log2FoldChange)>1))
WDs<-nrow(resWD %>% filter(padj<pt,abs(log2FoldChange)>1))

cat(i,"\t",PDs,"\t",PWs,"\t",WDs,"\n",sep="",file="DEgenes_padj_0.1_FC_1.txt",append=TRUE)

}



 #  assemble complete file with DESeq results
elist<-read.table("ENSG.list",header=FALSE)
colnames(elist)<-"ENSG"

CellTypes<-gsub("/","Or",CellTypes)


comp<-c("W_v_DB","P_v_DB","P_v_W")
for (i in comp) {
  cmp<-comp[i]
  for (j in CellTypes[c(1,2,3,4,6,7,8,9,10)]){

    #i=comp[1]
    #j=CellTypes[1]


    
    x<-read.table(paste0("CL_",j,"_",i,"_ALL2018-11-07.txt"),header=TRUE)
    n<-colnames(x)
    n<-paste(i,"_",j,"_",n,sep="")
    colnames(x)<-n

    x$ENSG=rownames(x)
    
    elist<-left_join(elist,x,by="ENSG")

  }
