## options(repos = c(CRAN = "http://cran.rstudio.com"))
##   This uses updated Seurat package 3 - starts with merged counts/demux from step 2

library(Seurat)
library(Matrix)
library(tidyverse)

library(annotables)

###########################################
## Testing sc transform           
## 2and3_Diem_Output
## adata <- read_rds("./kb_diem_Output/kb_diem_Seurat.list.rds")


##sc <- read_rds("/nfs/rprdata/scilab/novogene/Analyses/Roger_20200218/2_ST_integrate_Output_v2/ST_Integrated.obj.rds")

##sc <- read_rds("./sc.integrated.NormByLocation.rds")

sc <- read_rds("./sc.integrated.NormByLocation.Harmony.rds")

outFolder="./covid_plots_Harmony_merge/"
system(paste0("mkdir -p ", outFolder))
setwd(outFolder)


##table(sc$FinalName)

table(sc$sclabor.tlabel)

table(sc$Location)

table(sc$Location,sc$sclabor.tlabel)

##

##source("../theme_black.R")
source("../../scLaborColors.R")

##sc2 <- sc[,!is.na(sc$FinalName)]

sc2 <- sc[,sc$Location!="Blood"]

##sc2 <- sc[,!is.na(sc$Location)]


sc2@meta.data$Pred.Id.scLaborRef <- factor(sc2$sclabor.tlabel,levels=c("CTB","EVT","npiCTB","STB","Endometrial","Decidual","Fibroblast","Stromal-1","Stromal-2","Stromal-3","HSC","LED","Monocyte","Macrophage-1","Macrophage-2","B-cell","T-cell-activated","T-cell-resting","NK-cell"))

##sc2$SeuratClusters <- Idents(sc2)

sc2$Location2 = paste(sc2$Trimester,"trim.",sc2$Location)

table(sc2$Location2)


sc2@meta.data$LocationShort <- recode_factor(sc2@meta.data$Location2,
                              "1st trim. Decidua"="1DP",
                              "1st trim. Placenta"="1DP",
                              "2nd trim. PVBP"="2DP",
                              "3rd trim. BP"="3DP",
                              "3rd trim. PV"="3DP",
                              "3rd trim. Nucleus"="3Nuc",
                              "3rd trim. CAM"="3CAM")


sc2@meta.data$Location2 <- recode_factor(sc2@meta.data$Location2,
                              "1st trim. Decidua"="1st trim. Decidua + Placenta",
                              "1st trim. Placenta"="1st trim. Decidua + Placenta",
                              "2nd trim. PVBP"="2nd trim. Decidua + Placenta",
                              "3rd trim. BP"="3rd trim. Decidua + Placenta",
                              "3rd trim. PV"="3rd trim. Decidua + Placenta",
                              "3rd trim. Nucleus"="3rd trim. Dec. + Placenta Nuclei",
                              "3rd trim. CAM"="3rd trim. Chorioamniotic membranes")




table(sc2$Location2,sc2$Pred.Id.scLaborRef)

table(sc2$LocationShort,sc2$Pred.Id.scLaborRef)


### 

Idents(sc2) <- "Pred.Id.scLaborRef"

DefaultAssay(sc2) <- "RNA"

group.colors["STB"]="olivedrab1"

pdf("UMAP_Location.sclabor.tl.v2.pdf",width=12,height=5)
aa <- FetchData(sc2,c("UMAP_1","UMAP_2","FinalName","Location2","Group","Pred.Id.scLaborRef")) 
p1 <- ggplot(aa,aes(UMAP_1,UMAP_2,color=Pred.Id.scLaborRef)) +
    geom_point(size=0.1) +
    scale_color_manual(values=group.colors) +
    guides(colour = guide_legend(override.aes = list(size=5),title="Cell type")) +
    facet_grid(.~Location2) +
    theme_bw()
p1
##    theme_black()
dev.off()



grch38nodup <- grch38[!duplicated(grch38$ensgene),]
anno <- tibble(kbid=rownames(sc2)) %>% mutate(ensgene=gsub("\\..*","",kbid)) %>% left_join(grch38nodup)


## SARS-CoV-2  ACE2 TMPRSS2
##CD209/DC-SIGN AXL TYRO3 Zika
## CMV "NRP2","PDGFRA" cell.com/cell/fulltext/S0092-8674(18)30796-7?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867418307967%3Fshowall%3Dtrue

##"(Microbial infection)" "coronavirus" site:genecards.org
## "(Microbial infection)"  "receptor for" virus  site:genecards.org 

##receptors <- tibble(symbol=c("ACE2","TMPRSS2","AXL","CD209","TYRO3","NRP2","PDGFRA"),virus=c("SARS-CoV-2","SARS-CoV-2","Zika","Zika","Zika","CMV","CMV"))
##aa$Location <- recode_factor(aa$Location,C="Cervix",M="Myometrium",PVBP="Placenta & Basal Plate",MP="Myometrium + Placenta Accr.")

## CD209 could be CoV


receptors <- tibble(symbol=c("ACE2","TMPRSS2","FURIN","CTSL","BSG","SIGLEC1","NRP2","PDGFRA","AXL","CD209",
                             "ADAM17","TYRO3","CLEC4M","ANPEP","DPP4","CALM1","C1QBP","CD46","SLAMF1"),
                    virus=c("SARS-CoV-2","SARS-CoV-2","SARS-CoV-Alt","SARS-CoV-Alt","SARS-CoV-Alt","SARS-CoV-Alt","CMV","CMV","ZIKV","ZIKV/Dengue",
                            "ACE-Modifier","Ebola/Zika","229E-CoV","229E-CoV","MERS-CoV","Rubella","Rubella","Measles","Measles"),
                    vpanel=c(rep("A",10),rep("B",9)))

selGenes <- left_join(receptors,anno) %>% select(kbid,symbol,virus,vpanel)
selGenes

##features = selGenes$kbid



a2 <- FetchData(sc2,c("UMAP_1","UMAP_2","FinalName","Location2","LocationShort","Group","Pred.Id.scLaborRef")) 

myscale = 1/colSums(sc2@assays$RNA@counts)*1000000

rec <- map_dfr(1:nrow(selGenes),function(ii){
 a2$symbol=selGenes$symbol[ii]
 a2$Expression=sc2@assays$RNA@counts[selGenes$kbid[ii],]*myscale
## a2$Expression=sc2@assays$RNA@scale.data[selGenes$kbid[ii],]
 a2
})

dim(rec)

sum_rec <- rec %>% group_by(symbol,LocationShort,Pred.Id.scLaborRef) %>% summarize(Prop=mean(Expression>0),Expr=mean(Expression[Expression>0]))

dim(sum_rec)

head(sum_rec)

sum_rec <- sum_rec %>% left_join(selGenes) 

##sum_rec <- mutate(sum_rec,virus = fct_relevel(virus,c("SARS-CoV-2","ARS-CoV-Alt","CMV","CMV","ZIKV/Dengue","229E-CoV","MERS-CoV","Rubella","Measles","Ebola/Zika"))
                                         
pdf("scLabor.DotPlotViralReceptors.ALL.pdf",width=18,height=6)
ggplot(sum_rec,aes(x=LocationShort,y=Pred.Id.scLaborRef,color=Expr,size=Prop)) +
    geom_point() +
##    scale_size(trans='log10',breaks=c(0,0.0001,0.001,0.01,0.1,0.2,0.4,0.6,0.8,1.0)) +
    scale_size_area(breaks=c(0.001,0.01,0.1,0.2,0.4,0.8,1.0),max_size=6,na.value=0) +
##    scale_size(breaks=c(0.001,0.01,0.1,0.2,0.4,0.8,1.0),range=c(0.1,6),limits=c(0,1)) +
##    scale_color_continuous(trans='log10') +
##    scale_color_distiller(trans='log10', palette = "Spectral") +
    scale_color_distiller(trans='log10', palette = "RdYlBu",direction = -1,na.value=0,limits=c(10,NA)) +
##    scale_color_gradient(trans='log10', low = "Blue", high="red", na.value=0) +
    facet_grid( .~ virus + symbol) +
##   RotatedAxis() +
    theme_classic() +
    labs(x="Tissue",y="Cell type",size="Prop.",color="Expr.(TPM)") + 
    ## theme(legend.position="none") +
    theme(axis.text.x = element_text(angle = 45,hjust=1)) 
dev.off()


pdf("Figure2.pdf",width=14,height=6)
##sum_rec %>% filter(grepl("SARS",virus)) %>%
sum_rec %>% filter(vpanel=="A") %>%
    ##mutate(virus = factor(virus,c("SARS-CoV-2","ARS-CoV-Alt","CMV","CMV","ZIKV/Dengue","229E-CoV","MERS-CoV","Rubella","Measles","Ebola/Zika")) %>%
    mutate(virus = factor(virus,c("SARS-CoV-2","SARS-CoV-Alt","CMV","ZIKV","ZIKV/Dengue"))) %>%
    ggplot(aes(x=LocationShort,y=Pred.Id.scLaborRef,color=Expr,size=Prop)) +
    geom_point() +
##    scale_size(trans='log10',breaks=c(0,0.0001,0.001,0.01,0.1,0.2,0.4,0.6,0.8,1.0)) +
    scale_size_area(breaks=c(0.001,0.01,0.1,0.2,0.4,0.8,1.0),max_size=6,na.value=0) +
##    scale_size(breaks=c(0.001,0.01,0.1,0.2,0.4,0.8,1.0),range=c(0.1,6),limits=c(0,1)) +
##    scale_color_continuous(trans='log10') +
##    scale_color_distiller(trans='log10', palette = "Spectral") +
    scale_color_distiller(trans='log10', palette = "RdYlBu",direction = -1,na.value=0,limits=c(10,NA)) +
##    scale_color_gradient(trans='log10', low = "Blue", high="red", na.value=0) +
    facet_grid( .~ virus + symbol) +
##   RotatedAxis() +
    theme_classic() +
    labs(x="Tissue",y="Cell type",size="Prop.",color="Expr.(TPM)") +
    scale_y_discrete(limits = rev(levels(sum_rec$Pred.Id.scLaborRef))) +
    ## theme(legend.position="none") +
    theme(axis.text.x = element_text(angle = 45,hjust=1)) 
dev.off()

pdf("FigureS2.pdf",width=12,height=6)
##sum_rec %>% filter(grepl("SARS",virus)) %>%
sum_rec %>% filter(vpanel=="B") %>%
    ##mutate(virus = factor(virus,c("SARS-CoV-2","ARS-CoV-Alt","CMV","CMV","ZIKV/Dengue","229E-CoV","MERS-CoV","Rubella","Measles","Ebola/Zika")) %>%
    mutate(virus = factor(virus,c("ACE-Modifier","229E-CoV","MERS-CoV","Rubella","Measles","Ebola/Zika"))) %>%
    ggplot(aes(x=LocationShort,y=Pred.Id.scLaborRef,color=Expr,size=Prop)) +
    geom_point() +
##    scale_size(trans='log10',breaks=c(0,0.0001,0.001,0.01,0.1,0.2,0.4,0.6,0.8,1.0)) +
    scale_size_area(breaks=c(0.001,0.01,0.1,0.2,0.4,0.8,1.0),max_size=6,na.value=0) +
##    scale_size(breaks=c(0.001,0.01,0.1,0.2,0.4,0.8,1.0),range=c(0.1,6),limits=c(0,1)) +
##    scale_color_continuous(trans='log10') +
##    scale_color_distiller(trans='log10', palette = "Spectral") +
    scale_color_distiller(trans='log10', palette = "RdYlBu",direction = -1,na.value=0,limits=c(10,NA)) +
##    scale_color_gradient(trans='log10', low = "Blue", high="red", na.value=0) +
    facet_grid( .~ virus + symbol) +
##   RotatedAxis() +
    theme_classic() +
    labs(x="Tissue",y="Cell type",size="Prop.",color="Expr.(TPM)") +
    scale_y_discrete(limits = rev(levels(sum_rec$Pred.Id.scLaborRef))) +
    ## theme(legend.position="none") +
    theme(axis.text.x = element_text(angle = 45,hjust=1)) 
dev.off()




ind = selGenes$kbid[selGenes$virus=="SARS-CoV-2"]
aa$ExprL <- colSums(sc2@assays$RNA@counts[ind,]>0)==2



#pdf("Seurat_Ridge_plot_Covid19.pdf",width=20,height=8)
#RidgePlot(sc2, features = selGenes$kbid, ncol = 9)
#dev.off()


##library(ggrastr)

pdf("UMAP_Location.sclabor.GeneExpr.Covid19.v2.pdf",width=12,height=4)
ind = selGenes$kbid[selGenes$virus=="SARS-CoV-2"]
aa$ExprL <- colSums(sc2@assays$RNA@counts[ind,]>0)==2
p2 <- ggplot(arrange(aa,ExprL),aes(UMAP_1,UMAP_2,color=ExprL,size=ExprL*5+0.1)) +
    geom_point(size=0.1) +
##    geom_point_rast() + ## ggrastr
##    scale_color_manual(values=c("#33333322","white")) +
    scale_color_manual(values=c("#DDDDDD22","#AA0000")) +
    facet_grid(.~Location2) + 
    theme_bw() + theme(legend.position="none")
p2
dev.off()


library(patchwork)


pdf("Figure1.bw.horz.pdf",width=14,height=7)
p1 / p2 + plot_layout(guides = 'collect') +
    plot_annotation(tag_levels = 'A')
dev.off()



##Try rotate

table(aa$Pred.Id.scLaborRef,aa$ExprL)

tableS1 <- sc2@meta.data %>%
    select(Location2,LocationShort,Library,nCount_RNA,nFeature_RNA,PREG.BEST.GUESS) %>%
    group_by(Location2,LocationShort,Library)  %>%
    summarize(Num.Pooled=length(unique(PREG.BEST.GUESS)),Num.Cells=n(),Avg.UMI=mean(nCount_RNA),Avg.Genes=mean(nFeature_RNA))

tableS1

write_tsv(tableS1,"tableS1.tsv.xls")


sc2@meta.data %>% mutate(preg.id=paste0(PREG.BEST.GUESS,Fetus)) %>% select(LocationShort,preg.id) %>%
    group_by(LocationShort) %>%
    summarize(Num.Preg=length(unique(preg.id)))

sc2@meta.data %>% select(LocationShort,Library) %>%
    group_by(LocationShort) %>%
    summarize(Num.Preg=length(unique(Library)))


##    select(Location2,LocationShort,Library,nCount_RNA,nFeature_RNA,PREG.BEST.GUESS) %>%
##    group_by(Location2,LocationShort,Library)  %>%
##    summarize(Num.Pooled=length(unique(PREG.BEST.GUESS)),Num.Cells=n(),Avg.UMI=mean(nCount_RNA),Avg.Genes=mean(nFeature_RNA))



### END- HERE ###
########################################################


