## options(repos = c(CRAN = "http://cran.rstudio.com"))
##   This uses updated Seurat package 3 - starts with merged counts/demux from step 2

library(Seurat)
library(Matrix)
library(tidyverse)

library(future)

future::plan(strategy = 'multicore', workers = 16)
options(future.globals.maxSize = 20 * 1024 ^ 3)

###########################################

sc <- read_rds("./2_ST_integrate_Output_v2/ST_Integrated.obj.rds")

scref <- read_rds("/nfs/rprdata/scilab/novogene/Analyses/Roger_20200218/2_ST_integrate_Output_v2/ST_Integrated.obj.rds")

outFolder="./3_scTransferLabel_scLabor/"
system(paste0("mkdir -p ",outFolder))
setwd(outFolder)

selgene <- intersect(rownames(scref),rownames(sc))

length(selgene)

DefaultAssay(scref) <- "RNA"

scref <- scref[selgene,!is.na(scref$FinalName)]

scref <- SCTransform(scref, verbose = TRUE)

DefaultAssay(sc) <- "RNA"

query <- sc[selgene,]

query <- SCTransform(query, verbose = TRUE)


## Try reducing the number of PCs if it fails.  

anchors <- FindTransferAnchors(reference = scref, query = query, dims = 1:30)

predictions <- TransferData(anchorset = anchors, refdata = scref$FinalName, dims = 1:30)

table(query$seurat_clusters,predictions$predicted.id)

##summary(predictions$predicted.id==sc$FinalName)

mean(predictions$prediction.score.max>0.2)

stopifnot(identical(rownames(predictions),rownames(sc@meta.data)))

sc@meta.data$sclabor.tlabel <- predictions$predicted.id

table(sc@meta.data$sclabor.tlabel)

##sc@meta.data$sclabor.tlabel[predictions$prediction.score.max < 0.5] = "unk"

write_rds(sc,"ST_Integrated.scLabor.obj.rds")

##
theme_black = function(base_size = 12, base_family = "") {
    theme_grey(base_size = base_size, base_family = base_family) %+replace%
        theme(
            ## Specify axis options
            axis.line = element_blank(),
            axis.text.x = element_text(size = base_size*0.8, color = "white", lineheight = 0.9),
            axis.text.y = element_text(size = base_size*0.8, color = "white", lineheight = 0.9),
            axis.ticks = element_line(color = "white", size  =  0.2),
            axis.title.x = element_text(size = base_size, color = "white", margin = margin(0, 10, 0, 0)),
            axis.title.y = element_text(size = base_size, color = "white", angle = 90, margin = margin(0, 10, 0, 0)),
            axis.ticks.length = unit(0.3, "lines"),
            ## Specify legend options
            legend.background = element_rect(color = NA, fill = "black"),
            legend.key = element_rect(color = "white",  fill = "black"),
            legend.key.size = unit(1.2, "lines"),
            legend.key.height = NULL,
            legend.key.width = NULL,
            legend.text = element_text(size = base_size*0.8, color = "white"),
            legend.title = element_text(size = base_size*0.8, face = "bold", hjust = 0, color = "white"),
            legend.position = "right",
            legend.text.align = NULL,
            legend.title.align = NULL,
            legend.direction = "vertical",
            legend.box = NULL,
            ## Specify panel options
            panel.background = element_rect(fill = "black", color  =  NA),
            panel.border = element_rect(fill = NA, color = "white"),
            panel.grid.major = element_line(color = "grey35"),
            panel.grid.minor = element_line(color = "grey20"),
            panel.margin = unit(0.5, "lines"),
            ## Specify facetting options
            strip.background = element_rect(fill = "grey30", color = "grey10"),
            strip.text.x = element_text(size = base_size*0.8, color = "white"),
            strip.text.y = element_text(size = base_size*0.8, color = "white",angle = -90),
            ## Specify plot options
            plot.background = element_rect(color = "black", fill = "black"),
            plot.title = element_text(size = base_size*1.2, color = "white"),
            plot.margin = unit(rep(1, 4), "lines")
        )
}


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

group.colors<-c(B_cells="red",ColumnCyto="limegreen",Cytotroph="purple4",Decidual="sienna4",Dendr_Macro_A="hotpink",Dendr_Macro_B="tomato",Endometrial="maroon",Endothelial="yellow4",EVT="plum2",Fibroblasts="navy",HSC="magenta",Monocytes="lightslateblue",Progenitor="peru",NK_cells="royalblue",Stromal_A="turquoise4",Stromal_B="lightpink3",Synciotrophoblasts="yellowgreen",Tcells_activated="seagreen",Tcells_resting="powderblue")

names(group.colors) <- newNames 

pdf("UMAP_sclabor.tl.v2.pdf",width=10,height=8)
aa <- FetchData(sc,c("UMAP_1","UMAP_2","FinalName","Location","sclabor.tlabel"))
ggplot(aa,aes(UMAP_1,UMAP_2,color=sclabor.tlabel)) +
    geom_point(size=0.5) +
    scale_color_manual(values=group.colors) +
##    facet_grid(Location ~ Group) 
##     theme_bw()
    theme_black()
dev.off()

pdf("UMAP_Location.sclabor.tl.v2.pdf",width=10,height=8)
##aa <- FetchData(sc,c("UMAP_1","UMAP_2","FinalName","Location","Group","Pred.Id"))
ggplot(aa,aes(UMAP_1,UMAP_2,color=sclabor.tlabel)) +
    geom_point(size=0.1) +
    scale_color_manual(values=group.colors) +
##    facet_wrap("Location",nrow=2) + 
    theme_black()
dev.off()


## To-do add heatmap and alluvial plot with the cluster correspondance. 

X <- model.matrix(~ 0 + seurat_clusters,data=sc@meta.data)
Y = t(X) %*% as.matrix(predictions[,-c(1,ncol(predictions))])
Ybar = Y * 1/colSums(X)

rownames(Ybar) <- gsub("seurat_clusters","SC_",rownames(Ybar))
colnames(Ybar) <- gsub("prediction.score.","",colnames(Ybar))

Ybar[1:5,1:3]


library(pheatmap)

pdf("pheatmap.pdf",width=8,height=6)
pheatmap(t(Ybar),cluster_rows=FALSE,cluster_cols=FALSE,scale="none")
dev.off()

### END- HERE ###
########################################################


