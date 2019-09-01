## options(repos = c(CRAN = "http://cran.rstudio.com"))
##install.packages("tidyverse")
##library("rhdf5")
##library(cellrangerRkit)
library(Seurat)
library(Matrix)
library(tidyverse)

## This script follows after runScaipSeuratPCA.R We may want to take out the
## clustering part to another script.


selDate <- "2018-06-06"
basefolder <- "/nfs/rprdata/scilab/novogene/Analyses/ED2/"


outFolder <- paste0(basefolder,"outputs1/")
##system(paste0("mkdir -p ",outFolder))
setwd(outFolder)

#prefix="SCAIP4_seurat_genes_Res0.2_rmRiboMito"
fname = "NovoGene_R3__seurat_genes_Res0.72018-06-05.rds"
#Load seurat object, takes ~30 secs to load 
scilab <- read_rds(fname)

####################################################### 
gene_anno <- read_tsv("/nfs/rprdata/our10xData/counts/SCAIP-new/SCAIP1-ctrl/outs/filtered_gene_bc_matrices/hg19/genes.tsv",col_names=c("ensg.id","symbol"))
gene_symbol <- gene_anno$symbol                                                  

names(gene_symbol) <- gene_anno$ensg.id


##   *********    add MergedNewClusters to scilab object
setwd("../MergedClusters")
newCls<-read.table("ClusterConversion.txt",header=TRUE,row.names=1)
aux<-data.frame(row.names=rownames(scilab@meta.data),OldCls=scilab@meta.data$res.0.7)
aux$NewCls<-newCls[as.character(aux$OldCls),"NewClusterID"]
aux$NewClsName<-newCls[as.character(aux$OldCls),"NewClustername"]
aux<-aux[,2:3]
scilab<-AddMetaData(scilab,aux)

##      compute markers between clusters
##
scilab<-SetAllIdent(scilab,id="NewClsName")



# Run UMAP
scilab <- RunUMAP(scilab, reduction.use = "pca", dims.use = 1:20)

save("scilabObject-UMAP.Rdata")



#  compare with markers from Tsang paper
# compare with markers - paper #2
selsymb<-c("DKK1","IGFBP1","PRL","CDH5","CD34","ICAM1","PLVAP","CNN1","MYH11","ECM1","FIBIN","FMOD","THY1","COL1A1","VIM","CD14","CSF1R","CD53","AIF1","CD52","CD83"
,"CD4","CD86","CD163","CD209","HBB","ALAS2","HBA1","HBG1","MMP11","PAPPA2","HLA-G","PARP1","CGA","CYP19A1","GH2","CD2","CD3G","CD247","GZMA")

#selsymb<-c("HLA-DQB1","FSTL3","PAPPA","CD34","HLA-F","BCL2","NOS1","ACE")
##selsymb<-c("CD90","CD146","CD166","CD44","CD73","CD105","CD45","HLA-DR","CD19")
#selsymb<-t12$symbol[1:9]
selgene <- (gene_anno %>% filter(symbol %in% selsymb) %>% select(ensg.id) %>% unlist)

##q<-gene_anno %>% filter(ensg.id %in% selgene)
##selgene2<-q[match(selsymb,q$symbol),]
##selgene<-selgene2$ensg.id

##selsymb="XIST"
##selgene <- (gene_anno %>% filter(symbol %in% selsymb) %>% select(ensg.id) %>% unlist) 

pdf("Novogenemarkers3_Roger.pdf",width=48,height=32)
#selgene <- (gene_anno %>% filter(symbol=="TMSB4") %>% select(ensg.id) %>% unlist)

FeaturePlot(scilab, features.plot = selgene, cols.use = c("grey", "blue"),
            reduction.use = "tsne")

dev.off()

#    get computed markers
m1<-read.table("../outputs1/NovoGene_R3__seurat_genes_Res0.7markersB_2018-06-06.txt",header=TRUE)

m2<-m1 %>% top_n(5,-p_val_adj)

m3<-m2 %>% group_by(cluster) %>% top_n(10,avg_logFC)

m3$symbol<-gene_symbol[m3$gene]

M10_13<-FindMarkers(scilab,ident.1="10",ident.2="13",min.pct=.25)

# get m1
load("outputs1/Markers")



m1$gene<-gene_symbol[m1$gene]
#m2<-m1 %>% top_n(10,-p_val_adj)

m2<-m1 %>% group_by(cluster) %>% top_n(2,avg_logFC)

# have to get first 15, because extra digits added to rownames to make distinct...
m1$ENSG<-substring(rownames(m1),1,15)
m10<-m1 %>% group_by(cluster) %>% filter(avg_logFC>1.0) %>% top_n(-10,p_val)
##############################################################################################

#####   Features plot - ADD gene names    #####
library(rlist)
selsymb<-c("TSIX","XIST","HDHD1","KDM6A","PRAME","DDX3X","KDM5C","EIF1AX","SMC1A")
#selsymb<-c("SPP1","F13A1","PLTP","FCGBP","CD14","HLA-DRA","HLA-DPA1","HLA-DPB1","HLA-DQA1") 
                                        #selsymb<-t10$symbol[1:9]
selsymb<-m2$gene
selgene <- (gene_anno %>% filter(symbol %in% selsymb) %>% select(ensg.id) %>% unlist)
#selsymb<-m2$gene


q<-gene_anno %>% filter(ensg.id %in% selgene)
selgene2<-q[match(selsymb,q$symbol),]
selgene<-selgene2$ensg.id

pdf("TopMarkersExpression.pdf",height=12,width=12)
#par(mfrow=c(2,3),oma=0,0,2,0)
#<-FeaturePlot(scilab,pt.size=1.5,features.plot=selgene[1],cols.use=c("grey","blue"),reduction.use="tsne",do.return=TRUE)
#p-rep(t,length(selgene))

p<-list()

for (i in 1:length(selgene)){
  t   <- FeaeturePlot2(scilab,pt.size=1.5,features.plot=selgene[i],cols.use=c("grey","blue"),reduction.use="tsne",do.return=TRUE)
t   <-t[[1]]    +ggtitle(selsymb  [i])+theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank())
p[[i]]<-t
}
plot_grid(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],p[[10]],p[[11]],p[[12]],p[[13]],p[[14]],p[[15]],p[[16]],p[[17]],p[[18]],p[[19]],p[[20]],p[[21]],p[[22]],p[[23]],p[[24]],p[[25]],p[[26]],p[[27]],p[[28]],p[[29]],p[[30]],p[[31]],p[[32]],p[[33]],p[[34]],p[[35]],p[[36]],p[[37]],p[[38]],ncol=8,nrow=5)
dev.off()

##   ***********  heatmap of marker genes, top 5 per cluster   *****
pdf("Heatmap_MarkerGenes_10PerCluster.pdf",width=12,height=12)
DoHeatmap(scilab,genes.use = m10e$ENSG,slim.col.label = TRUE,remove.key = TRUE,group.cex=5,group.label.rot=TRUE)
dev.off()

###############################################################################


# 3x3 Heatmap by DE genes, by Group

mydf <- cbind(as.data.frame(scilab@dr$umap@cell.embeddings),scilab@meta.data)

mydf$Indiv <- gsub("[WDP][B]?", "", mydf$LibraryID)
mydf$NEW_BARCODE <- rownames(mydf)

write_rds(mydf,"newMergedClustersData.rds")

#  Downsample, assign correct
 # get 3000 from each Group,Location
m_PL_DB<-mydf %>% filter(Group=="PL",Location=="DB") %>%  sample_n(3000,replace=TRUE)
m_PL_P<-mydf %>% filter(Group=="PL",Location=="P") %>%  sample_n(3000,replace=TRUE)
m_PL_W<-mydf %>% filter(Group=="PL",Location=="W") %>%  sample_n(3000,replace=TRUE)
m_TL_DB<-mydf %>% filter(Group=="TL",Location=="DB") %>%  sample_n(3000,replace=TRUE)
m_TL_P<-mydf %>% filter(Group=="TL",Location=="P") %>%  sample_n(3000,replace=TRUE)
m_TL_W<-mydf %>% filter(Group=="TL",Location=="W") %>%  sample_n(3000,replace=TRUE)
m_TN_DB<-mydf %>% filter(Group=="TN",Location=="DB") %>%  sample_n(3000,replace=TRUE)
m_TN_P<-mydf %>% filter(Group=="TN",Location=="P") %>%  sample_n(3000,replace=TRUE)
m_TN_W<-mydf %>% filter(Group=="TN",Location=="W") %>%  sample_n(3000,replace=TRUE)

mydf_3000<-rbind(m_PL_DB,m_PL_P,m_PL_W,m_TL_DB,m_TL_P,m_TL_W,m_TN_DB,m_TN_P,m_TN_W)


# change some cluster names
cls<-levels(mydf_3000$NewClsName)
levels(mydf_3000$NewClsName)<-c(cls[1:4],"Dendr_Macro_A","Dendr_Macro_B",cls[7:12],"Progenitor",cls[14],"Stromal_A","Stromal_B",cls[17],"Tcells_activated","Tcells_resting")

#specify more distinct colors
cls<-levels(mydf_3000$NewClsName)
group.colors<-c(B_cells="red",ColumnCyto="limegreen",Cytotroph="purple4",Decidual="sienna4",Dendr_Macro_A="hotpink",Dendr_Macro_B="tomato",Endometrial="maroon",Endothelial="yellow4",EVT="plum2",Fibroblasts="navy",HSC="magenta",Monocytes="lightslateblue",Progenitor="peru",NK_cells="royalblue",Stromal_A="turquoise4",Stromal_B="lightpink3",Synciotrophoblasts="yellowgreen",Tcells_activated="seagreen",Tcells_resting="powderblue")

mydf_3000$color.codes=group.colors[mydf_3000$NewClsName]
mydf_3000s<-mydf_3000[order(mydf_3000$NewClsName),]

gg <- ggplot(mydf_3000s,aes(x=UMAP1,y=UMAP2,col=NewClsName)) + geom_point(alpha=1,size=0.5) + guides(colour = guide_legend(override.aes = list(size=4)))+scale_color_manual(values=unique(as.character(mydf_3000s$color.codes)))
ggsave("Tissue_Full4.pdf",gg,device="pdf",width=12,height=12)

f_3000f<-mydf_3000s
mydf_3000f$Group<-factor(mydf_3000f$Group)
mydf_3000f$Location<-factor(mydf_3000f$Location)

levels(mydf_3000f$Group)<-c("PTL","TIL","TNL")
levels(mydf_3000f$Location)<-c("DB","PV","CAM")

mydf_3000f$Group2<-factor(mydf_3000f$Group,levels=c("TNL","TIL","PTL"))
mydf_3000f$Location2<-factor(mydf_3000f$Location,levels=c("DB","PV","CAM"))


gg <- ggplot(mydf_3000f,aes(x=UMAP1,y=UMAP2,col=NewClsName)) + geom_point(alpha=1,size=0.5) + facet_grid(Group2 ~ Location2)+ guides(colour = guide_legend(override.aes = list(size=4)))+scale_color_manual(values=unique(as.character(mydf_3000f$color.codes)))+theme_bw()
ggsave("Downsampled3000_facet4.pdf",gg,device="pdf",width=12,height=12)

gg <- ggplot(mydf,aes(x=tSNE_1,y=tSNE_2,col=NewClsName)) +
    geom_point(alpha=1,size=0.5,show.legend=F) +
    facet_grid(Location ~ .)

ggsave("Location_facet.pdf",gg,device="pdf",width=6,height=12)










gg <- ggplot(mydf,aes(x=tSNE_1,y=tSNE_2,col=NewClsName)) + geom_point(alpha=1,size=0.5) + guides(colour = guide_legend(override.aes = list(size=4)))
ggsave("Tissue_Full_legend2.pdf",gg,device="pdf",width=12,height=12)


gg <- ggplot(mydf,aes(x=tSNE_1,y=tSNE_2,col=NewClsName)) + geom_point(alpha=1,size=0.5) + facet_grid(Group ~ Location) 
ggsave("Group_Location_facet.pdf",gg,device="pdf",width=12,height=12)

gg <- ggplot(mydf,aes(x=tSNE_1,y=tSNE_2,col=NewClsName)) +
    geom_point(alpha=1,size=0.5,show.legend=F) +
    facet_grid(Location ~ .) 

ggsave("Location_facet.pdf",gg,device="pdf",width=6,height=12)

gg <- ggplot(mydf,aes(NewClsName,fill=Location)) +
    geom_bar(aes(y=..prop..,group=Location),position = position_dodge2(preserve = "single")) +
    coord_flip()
##+
##    scale_y_continuous(labels=percent_format())

ggsave("Location_Barplot.pdf",gg,device="pdf",width=5,height=7)


gg <- ggplot(mydf,aes(x=tSNE_1,y=tSNE_2,col=NewClsName)) +
    geom_point(alpha=1,size=0.5,show.legend=F) +
     facet_grid(Group ~ .) 

ggsave("Group_facet.pdf",gg,device="pdf",width=6,height=12)


gg <- ggplot(mydf %>% filter(!Group %in% "PL"),aes(NewClsName,fill=Group)) +
    geom_bar(aes(y=..prop..,group=Group),position = position_dodge2(preserve = "single")) +
    coord_flip()
##+
##    scale_y_continuous(labels=percent_format())

ggsave("Labor_Barplot.pdf",gg,device="pdf",width=5,height=7)



gg <- ggplot(mydf %>% filter(!Group %in% "TN"),aes(NewClsName,fill=Group)) +
    geom_bar(aes(y=..prop..,group=Group),position = position_dodge2(preserve = "single")) +
    coord_flip()
##+
##    scale_y_continuous(labels=percent_format())

ggsave("Preterm_Barplot.pdf",gg,device="pdf",width=5,height=7)


gg <- ggplot(mydf,aes(NewClsName,fill=Group)) +      facet_grid(Group ~ .) 
    geom_bar(aes(y=..prop..,group=Group),position = position_dodge2(preserve = "single")) +
    coord_flip()
##+
##    scale_y_continuous(labels=percent_format())

ggsave("Group_Barplot.pdf",gg,device="pdf",width=5,height=7)


gg <- ggplot(mydf,aes(NewClsName,fill=Group)) +     
    geom_bar(aes(y=..prop..,group=Group),position = position_dodge2(preserve = "single")) +
    coord_flip() +
    facet_grid(. ~ Location ) 

ggsave("NEWGroup_Barplot.pdf",gg,device="pdf",width=10,height=7)

##mydf %>% filter(!Group %in% "PL")
gg <- ggplot(mydf %>% filter(!Group %in% "PL"),aes(NewClsName,fill=Group)) +     
    geom_bar(aes(y=..prop..,group=Group),position = position_dodge2(preserve = "single")) +
    coord_flip() +
    facet_grid(. ~ Location ) 
ggsave("NEWLabor_Barplot.pdf",gg,device="pdf",width=10,height=7)





### I need a better way to do this testing. 

mytG <- mydf %>% group_by(Indiv,Group,NewClsName) %>%
    summarize(count=n(),.drop=FALSE) %>%
    mutate(proportion=count/sum(count))
mytG$Group = relevel(factor(mytG$Group),"TL")

myfun <- function(dat){
    mysum=summary(lm(proportion ~ Group  ,data=dat))
    mysum$coefficients[-(1:2),4]
}


res <- mytG %>% group_by(NewClsName) %>%
    do(mod=myfun(.)) 
    
minc <- mydf %>% count(Group,NewClsName) %>%
    group_by(NewClsName) %>%
    summarize(minc=min(n)) %>%
    filter(minc>100)


# Show expression of one gene
gsel_symb= "IL8"
gsel <- names(which(gene_symbol==gsel_symb))
gene.expr <- (scilab@data[gsel,])
identical(names(gene.expr),row.names(mydf))
mydf$gene.expr <- gene.expr

gg <- ggplot(mydf,aes(x=tSNE_1,y=tSNE_2,col=gene.expr)) +   geom_point(alpha=1,size=0.5) +  scale_colour_gradient(low="#e0ecf4",high="#6e016b") + facet_grid(Group ~ .)
#  scale_color_gradientn(colors=colorRampPalette(brewer.pal(8, "BuPu"))(20)[-(1:2)]) +
ggsave(paste0("TSNE_ggplot_",gsel_symb,".pdf"),gg,device="pdf",width=6,height=12)


  #Heat map of DE genes per cluster (p.adj<.05)
#TL_v_PL<-read.table("TL_v_PL.DEgenes")
#TN_v_TL<-read.table("TN_v_TL.DEgenes")

#library(tibble)
#gs<-as.data.frame(gene_symbol)
#gss<-gs %>% rownames_to_column("ENS")

#q<-as.character(TL_v_PL$V1)
#TL_v_PL.markers<-filter(gss,gene_symbol %in% q)$ENS
#TL_v_PL.markers.scale<-TL_v_PL.markers[TL_v_PL.markers %in% rownames(scilab@scale.data)]


#q<-as.character(TN_v_TL$V1)
#TN_v_TL.markers<-filter(gss,gene_symbol %in% q)$ENS
#TN_v_TL.markers.scale<-TN_v_TL.markers[TN_v_TL.markers %in% rownames(scilab@scale.data)]

#scilab<-SetAllIdent(scilab,id="res.0.7")
#Scilab0<-SubsetData(scilab,ident.use=0)
#Scilab0<-SetAllIdent(Scilab0,id="Group")
#Scilab0<-SubsetData(Scilab0,ident.use=c("TL","PL"))
#Scilab0@ident<-factor(Scilab0@ident,levels=(c("PL","TL")))

## New heatmap-Roger
library(pheatmap)

allmarkers <- unique(c(TL_v_PL.markers,TL_v_PL.markers))

anno_col <- scilab@meta.data %>% rownames_to_column() %>% 
  filter(Location=="W",Group!="TN") %>%
  select(rowname,res.0.7,Group) %>%
  arrange(res.0.7,Group) 
anno_col <- filter(anno_col,res.0.7 %in% names(which(table(anno_col$res.0.7)>200))) %>%
  as.data.frame() %>% column_to_rownames()


omat <- as.matrix(scilab@data[allmarkers,rownames(anno_col)])
rmax <- apply(omat,1,max)
omat <- omat / rmax

prefix <- "novogene";
pdf(paste0(prefix,"_",selDate,"_heatmap.pdf"), width=24,height=24)
pheatmap(omat,cluster_cols=FALSE,annotation_col=anno_col,scale="none",
         show_colnames=FALSE)
dev.off()


## New heatmap #2 - organize by DE Genes
xda<-data.frame(matrix(ncol = 9, nrow = 0))
colnames(xda)<-c("ENS","baseMean","log2FC","lfcSE","stat","pvalue","padj","Gene","CL")
for (i in 0:18){
  if (i==9) {next}
  x<-read.table(paste0("CL_",i,"_TL_v_PL2018-06-08.txt"),header=FALSE,skip=1)
  xd<-data.frame(x)
  xd$CL<-i
  xda<-rbind(xda,xd)
  #colnames(xd)<-c("ENS","baseMean","log2FC","lfcSE","stat","pvalue","padj","Gene","CL")
}
colnames(xda)<-c("ENS","baseMean","log2FC","lfcSE","stat","pvalue","padj","Gene","CL")                

xda<-xda[xda$padj<.1,]
#(q calculated below...)
xda<-filter(xda,CL %in% q)
mg<-unique(xda$ENS)

library(pheatmap)

#allmarkers <- unique(c(TL_v_PL.markers,TL_v_PL.markers))
q<-names(which(table(anno_col$res.0.7)>200))



anno_col <- scilab@meta.data %>% rownames_to_column() %>% 
  filter(Location=="P",Group!="TN") %>%
  select(rowname,res.0.7,Group) %>%
  arrange(res.0.7,Group) 
anno_col <- filter(anno_col,res.0.7 %in% names(which(table(anno_col$res.0.7)>200))) %>%
  as.data.frame() %>% column_to_rownames()


omat <- as.matrix(scilab@data[mg,rownames(anno_col)])
omat<-log(omat+1)
rmax <- apply(omat,1,max)
omat <- omat / rmax

prefix <- "novogene";
pdf(paste0(prefix,"_",selDate,"_heatmap_P_TL_v_PL.pdf"), width=24,height=24)
pheatmap(omat,cluster_cols=FALSE,cluster_rows=FALSE,annotation_col=anno_col,scale="none",
         show_colnames=FALSE)
dev.off()
                
## New heatmap #3   TN v TL - organize by DE Genes
xda<-data.frame(matrix(ncol = 9, nrow = 0))
colnames(xda)<-c("ENS","baseMean","log2FC","lfcSE","stat","pvalue","padj","Gene","CL")
for (i in 0:18){
  if (i==4) {next}
  if (i==12) {next}
  x<-read.table(paste0("CL_",i,"_TN_v_TL2018-06-08.txt"),header=FALSE,skip=1)
  xd<-data.frame(x)
  xd$CL<-i
  xda<-rbind(xda,xd)
  #colnames(xd)<-c("ENS","baseMean","log2FC","lfcSE","stat","pvalue","padj","Gene","CL")
}
colnames(xda)<-c("ENS","baseMean","log2FC","lfcSE","stat","pvalue","padj","Gene","CL")                

anno_col <- scilab@meta.data %>% rownames_to_column() %>% 
  filter(Location=="P",Group!="PL") %>%
  select(rowname,res.0.7,Group) %>%
  arrange(res.0.7,Group) 
anno_col <- filter(anno_col,res.0.7 %in% names(which(table(anno_col$res.0.7)>200))) %>%
  as.data.frame() %>% column_to_rownames()

q<-names(which(table(anno_col$res.0.7)>200))

xda<-xda[xda$padj<.01,]
#(q calculated below...)
xda<-filter(xda,CL %in% q)
mg<-unique(xda$ENS)

library(pheatmap)

allmarkers <- mg


omat <- as.matrix(scilab@data[mg,rownames(anno_col)])
omat<-log(omat+1)
rmax <- apply(omat,1,max)
omat <- omat / rmax

prefix <- "novogene";
pdf(paste0(prefix,"_",selDate,"_heatmap_P_TN_v_TL.pdf"), width=24,height=24)
pheatmap(omat,cluster_cols=FALSE,cluster_rows=FALSE,annotation_col=anno_col,scale="none",
         show_colnames=FALSE)
dev.off()

##     xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
FeaeturePlot2<-function (object, features.plot, min.cutoff = NA, max.cutoff = NA, 
    dim.1 = 1, dim.2 = 2, cells.use = NULL, pt.size = 1, cols.use = c("yellow", 
        "red"), pch.use = 16, overlay = FALSE, do.hover = FALSE, 
    data.hover = "ident", do.identify = FALSE, reduction.use = "tsne", 
    use.imputed = FALSE, nCol = NULL, no.axes = FALSE, no.legend = TRUE, 
    coord.fixed = FALSE, dark.theme = FALSE, do.return = FALSE, 
    vector.friendly = FALSE) 
{
    cells.use <- SetIfNull(x = cells.use, default = colnames(x = object@data))
    if (is.null(x = nCol)) {
        nCol <- 2
        if (length(x = features.plot) == 1) {
            nCol <- 1
        }
        if (length(x = features.plot) > 6) {
            nCol <- 3
        }
        if (length(x = features.plot) > 9) {
            nCol <- 4
        }
    }
    num.row <- floor(x = length(x = features.plot)/nCol - 1e-05) + 
        1
    if (overlay | do.hover) {
        num.row <- 1
        nCol <- 1
    }
    par(mfrow = c(num.row, nCol))
    dim.code <- GetDimReduction(object = object, reduction.type = reduction.use, 
        slot = "key")
    dim.codes <- paste0(dim.code, c(dim.1, dim.2))
    data.plot <- as.data.frame(GetCellEmbeddings(object = object, 
        reduction.type = reduction.use, dims.use = c(dim.1, dim.2), 
        cells.use = cells.use))
    x1 <- paste0(dim.code, dim.1)
    x2 <- paste0(dim.code, dim.2)
    data.plot$x <- data.plot[, x1]
    data.plot$y <- data.plot[, x2]
    data.plot$pt.size <- pt.size
    names(x = data.plot) <- c("x", "y")
    data.use <- t(x = FetchData(object = object, vars.all = features.plot, 
        cells.use = cells.use, use.imputed = use.imputed))
    min.cutoff <- mapply(FUN = function(cutoff, feature) {
        ifelse(test = is.na(x = cutoff), yes = min(data.use[feature, 
            ]), no = cutoff)
    }, cutoff = min.cutoff, feature = features.plot)
    max.cutoff <- mapply(FUN = function(cutoff, feature) {
        ifelse(test = is.na(x = cutoff), yes = max(data.use[feature, 
            ]), no = cutoff)
    }, cutoff = max.cutoff, feature = features.plot)
    check_lengths = unique(x = vapply(X = list(features.plot, 
        min.cutoff, max.cutoff), FUN = length, FUN.VALUE = numeric(length = 1)))
    if (length(x = check_lengths) != 1) {
        stop("There must be the same number of minimum and maximum cuttoffs as there are features")
    }
    if (overlay) {
        pList <- list(BlendPlot(data.use = data.use, features.plot = features.plot, 
            data.plot = data.plot, pt.size = pt.size, pch.use = pch.use, 
            cols.use = cols.use, dim.codes = dim.codes, min.cutoff = min.cutoff, 
            max.cutoff = max.cutoff, coord.fixed = coord.fixed, 
            no.axes = no.axes, no.legend = no.legend, dark.theme = dark.theme))
    }
    else {
        pList <- mapply(FUN = SingleFeaturePlot, feature = features.plot, 
            min.cutoff = min.cutoff, max.cutoff = max.cutoff, 
            coord.fixed = coord.fixed, MoreArgs = list(data.use = data.use, 
                data.plot = data.plot, pt.size = pt.size, pch.use = pch.use, 
                cols.use = cols.use, dim.codes = dim.codes, no.axes = no.axes, 
                no.legend = no.legend, dark.theme = dark.theme, 
                vector.friendly = vector.friendly), SIMPLIFY = FALSE)
    }
    if (do.hover) {
        if (length(x = pList) != 1) {
            stop("'do.hover' only works on a single feature or an overlayed FeaturePlot")
        }
        if (is.null(x = data.hover)) {
            features.info <- NULL
        }
        else {
            features.info <- FetchData(object = object, vars.all = data.hover)
        }
        return(HoverLocator(plot = pList[[1]], data.plot = data.plot, 
            features.info = features.info, dark.theme = dark.theme, 
            title = features.plot))
    }
    else if (do.identify) {
        if (length(x = pList) != 1) {
            stop("'do.identify' only works on a single feature or an overlayed FeaturePlot")
        }
        return(FeatureLocator(plot = pList[[1]], data.plot = data.plot, 
            dark.theme = dark.theme))
    }
    else {
        #print(x = cowplot::plot_grid(plotlist = pList, ncol = nCol))
    }
    ResetPar()
    if (do.return) {
        return(pList)
    }
}

e<-environment(FeaturePlot)
environment(FeaeturePlot2)<-e
