library(tidyverse)
library(knitr)
library(DESeq2)
library(annotables)
library(qqman)



basefolder <- "/nfs/rprdata/scilab/novogene/Analyses/ED3"


outFolder <- paste0(basefolder,"/outputs_DESeq_ConditionsByCluster_3Comps/")
system(paste0("mkdir -p ",outFolder))
setwd(outFolder)

#convert annos to gene names
gene_anno <- read_tsv("/nfs/rprdata/our10xData/counts/SCAIP-new/SCAIP1-ctrl/outs/filtered_gene_bc_matrices/hg19/genes.tsv",col_names=c("ensg.id","symbol"))
gene_symbol <- gene_anno$symbol
names(gene_symbol) <- gene_anno$ensg.id

afile <- paste0(basefolder,"/../ED2/80K-all.rds")
all <- readr::read_rds(afile)

clusterfile <- paste0(basefolder,"/../ED2/MergedClusters/newMergedClustersData.rds")
clusters <- readr::read_rds(clusterfile)

clusters$NewClsName <- as.character(clusters$NewClsName)

#only use clusters with at least 100 per Condition
cSize <- clusters %>% dplyr::count(Group,NewClsName) %>%
    group_by(NewClsName) %>%
    summarize(minc=min(n)) %>%
    filter(minc>100) %>% arrange(-minc)
cSize


CellTypes=cSize$NewClsName
CellTypes

cat("Celltype","\t","TIL_v_PTL","\t","TNL_v_TIL","\t","PTL_v_TNL","\n",sep="",file="QuickFDR10.txt",append=FALSE)
for (i in CellTypes) {
    clusters_i  <- filter(clusters,NewClsName %in% c(i))
    dim(clusters_i)    
    all_i <- all[,clusters_i$NEW_BARCODE]
    ## Building a design matrix to get aggregate
    ## RPR: I'm changing this to aggregate across locations. 
    bti <- clusters_i %>% transmute(bti=paste(Group,Indiv,sep="_")) %>% unlist %>% factor
    ##
    X <- model.matrix( ~ 0 + bti)
    qr.X <- qr(X)
    qr.X$rank
    dim(X)
    YtX <- all_i %*% X
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
    colnames(cvt)<-c("Group","Indiv")
    ##
    cvt$Group = relevel(factor(cvt$Group),"TL")
    ##
    ##
    ## Running DEseq2 on the aggregate
    ## ========================================================
    ##
    dds <- DESeqDataSetFromMatrix(cmat,cvt, ~ Group)    
    dds <- DESeq(dds,parallel=TRUE)
    ##
    ## save DEseq object. 
    fname <- paste0("ConditionsBetweenCluster_",CellTypes[i],"__dds_",Sys.Date(),".rds")
    fname
    write_rds(dds,fname)
    ##
    ##
    ## Parse the results 
    ## ========================================================
    ##
    ## TL v PL
    resPL <- results(dds)
    ##
    resPL<-results(dds,contrast=c("Group","TL","PL"))
    resPL$Gene<-gene_symbol[rownames(resPL)]
    write.table(resPL,paste0("CL_",gsub("/","Or",i),"_TL_v_PL_ALL",Sys.Date(),".txt"))
    ##
    ##
    ## TL v TN
    resTL <- results(dds)
    resTL <-results(dds,contrast=c("Group","TN","TL"))
    resTL$Gene<-gene_symbol[rownames(resTL)]
    write.table(resTL,paste0("CL_",gsub("/","Or",i),"_TN_v_TL_ALL",Sys.Date(),".txt"))
    ##
    ## ****  (add now)  TN v PL  ****
    resPN <- results(dds)
    resPN <-results(dds,contrast=c("Group","TN","PL"))
    resPN$Gene<-gene_symbol[rownames(resPN)]
    write.table(resPN,paste0("CL_",gsub("/","Or",i),"_TN_v_PL_ALL",Sys.Date(),".txt"))

    #filter !na, padj<.1, |FC|>1
    nPL<-nrow(data.frame(resPL@listData) %>% filter(padj<.1,abs(log2FoldChange)>1))
    nTL<-nrow(data.frame(resTL@listData) %>% filter(padj<.1,abs(log2FoldChange)>1))
    nPN<-nrow(data.frame(resPN@listData) %>% filter(padj<.1,abs(log2FoldChange)>1))
    

    
    cat(i,"\t",sum(resPL$padj<0.1,na.rm=TRUE),"\t",sum(resTL$padj<0.1,na.rm=TRUE),"\n",sep="")
    cat(i,"\t",nPL,"\t",nTL,"\t",nPN,"\n",sep="",file="QuickFDR10.txt",append=TRUE)
    ##
##    list(resTL=resTL,resPL=resPL)
    ##    
}

# Parse results - compare DE genes per group per 3 classes

cat("Celltype","\t","TIL_v_PTL [ TP ]","\t","TNL_v_TIL [ NL ]","\t","PTL_v_TNL [PN ]","TP_NL","\t","TP_PN","\t","NL_PN","\t","TP_NL_PN","\n",sep="",file="Intersections.txt",append=FALSE)

date="2018-10-15"
library(venneuler)

for (i in CellTypes){


#i=CellTypes[1]
  resTP<-read.table(paste0("CL_",gsub("/","Or",i),"_TL_v_PL_ALL",date,".txt"))
  resNL<-read.table(paste0("CL_",gsub("/","Or",i),"_TN_v_TL_ALL",date,".txt"))
  resPN<-read.table(paste0("CL_",gsub("/","Or",i),"_TN_v_PL_ALL",date,".txt"))

pt=0.1

resTPs<-resTP %>% filter(padj<pt) %>% select(Gene)
resNLs<-resNL %>% filter(padj<pt) %>% select(Gene)
resPNs<-resPN %>% filter(padj<pt) %>% select(Gene)
resTPs<-resTPs$Gene
resNLs<-resNLs$Gene
resPNs<-resPNs$Gene

TP<-length(resTPs)
NL<-length(resNLs)
PN<-length(resPNs)
TP_NL<-length(intersect(resTPs,resNLs))
TP_PN<-length(intersect(resPNs,resTPs))
NL_PN<-length(intersect(resNLs,resPNs))
TP_NL_PN<-length(intersect(intersect(resNLs,resTPs),resPNs))


vl<-c("TP"=TP,"NL"=NL,"PN"=PN,"TP&NL"=TP_NL,"TP&PN"=TP_PN,"NL&PN"=NL_PN,"TP&NL&PN"=TP_NL_PN)

cat(i,"\t",TP,"\t",NL,"\t",PN,"\t",TP_NL,"\t",TP_PN,"\t",NL_PN,"\t",TP_NL_PN,"\n",sep="",file="Intersections.txt",append=TRUE)


  
}

pdf("Venn6.pdf")
v<-venneuler(vl)
plot(v)
dev.off()
