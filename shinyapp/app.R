##library(Seurat)
library(Matrix)
library(tidyverse)

library(data.table)
library(ggridges)
##library(DT)

##save(gene_anno,gene_symbol,gene.expr.mat,mydf,file="./data/data.Rdata")


## loads scilab
load("data/data.Rdata")
## See if opens faster as h5


##conditions-all.txt 
deg_table <- fread("cat data/conditions-all.txt | tr ' ' '\t'")

##deg_table <- deg_table %>% mutate_if(is.numeric, scales::scientific, 3)
deg_table <- deg_table %>% mutate_if(is.numeric, signif, 3) %>% filter(padj<0.5)

newNames = c("B_cells"="B-cell"
                         ,"ColumnCyto"="npiCTB"
                         ,"Cytotroph"="CTB"
                         ,"Decidual"="Decidual"
                         ,"Dendr_Macro_A"="Macrophage-1"
                         ,"Dendr_Macro_B"="Macrophage-2"
                         ,"Endometrial"="Endometrial"
                         ,"Endothelial"="Endothelial"
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

deg_table$Cluster=newNames[deg_table$Cluster]

##locations-all.txt 
Ldeg_table <- fread("cat data/locations-all.txt | tr ' ' '\t'")

##deg_table <- deg_table %>% mutate_if(is.numeric, scales::scientific, 3)
Ldeg_table <- Ldeg_table %>% mutate_if(is.numeric, signif, 3) %>% filter(padj<0.5)

Ldeg_table$Cluster=newNames[Ldeg_table$Cluster]

Ldeg_table$Compare=gsub("DB","BP",Ldeg_table$Compare)

Ldeg_table$Compare=gsub("_P","_PV",Ldeg_table$Compare)

Ldeg_table$Compare=gsub("W_","CAM_",Ldeg_table$Compare)

Ldeg_table$Compare=gsub("_W","_CAM",Ldeg_table$Compare)

#threshold baseMean
Ldeg_table <- Ldeg_table %>% filter(baseMean>50)
deg_table <- deg_table %>% filter(baseMean>50)



# Marker table
Marker_table <- fread("cat data/MarkersAll.txt | tr ' ' '\t'")
Marker_table<-Marker_table %>% mutate_if(is.numeric,signif,3) %>% filter(p_val_adj<0.1)
Marker_table<-Marker_table[c(1,7,8)]
newNames = c("B_cells"="B-cell"
                         ,"ColumnCyto"="npiCTB"
                         ,"Cytotroph"="CTB"
                         ,"Decidual"="Decidual"
                         ,"Dendritic/Macrophage1"="Macrophage-1"
                         ,"Dendritic/Macrophage2"="Macrophage-2"
                         ,"Endometrial"="Endometrial"
                         ,"Endothelial"="Endothelial"
                         ,"EVT"="EVT"
                         ,"Fibroblasts"="Fibroblast"
                         ,"HSC"="HSC"
                         ,"Monocytes"="Monocyte"
                         ,"(Myeloid)Progenitor"="Stromal-3"
                         ,"NK_cells"="NK-cell"
                         ,"Stromal1"="Stromal-1"
                         ,"Stromal2"="Stromal-2"
                         ,"Synciotrophoblasts"="STB"
                         ,"Tcells_activated"="T-cell-activated"
                         ,"Tcells_resting"="T-cell-resting")

Marker_table$cluster=newNames[Marker_table$cluster]


MClusters<-unique(Marker_table$cluster)
##

newNames = c("B_cells"="B-cell"
                         ,"ColumnCyto"="npiCTB"
                         ,"Cytotroph"="CTB"
                         ,"Decidual"="Decidual"
                         ,"Dendr_Macro_A"="Macrophage-1"
                         ,"Dendr_Macro_B"="Macrophage-2"
                         ,"Endometrial"="Endometrial"
                         ,"Endothelial"="Endothelial"
                         ,"EVT"="EVT"
                         ,"Fibroblasts"="Fibroblast"
                         ,"HSC"="HSC"
                         ,"Monocytes"="Monocyte"
                         ,"Progenitor"="Stromal-3"
                         ,"NK_cells"="NK-cell"
                         ,"Stromal_A"="Stromal-1"
                         ,"Stromal_B"="Stromal-2"
                         ,"Synciotrophoblasts"="STB"
                         ,"Tcells_activated"="T-cell-activated"
                         ,"Tcells_resting"="T-cell-resting")

mydf$NewClsName<-newNames[mydf$NewClsName]

mydf$NewClsName<-factor(mydf$NewClsName)

levels(mydf$Location)=c("BP","PV","CAM")

cls <- levels(mydf$NewClsName)


## this is because we can then do +scale_fill_celltype() in any ggplot2 plot to color appropriatelly. 
scale_fill_celltype <- function(...){
    ggplot2:::manual_scale('fill', values = group.colors,...)
}


location.colors  = c(BP="#a6d854", PV="#fc8d62", CAM="#8da0cb")
condition.colors = c(TNL="#377eb8",TIL="#4daf4a", PTL="#e41a1c")

scale_fill_location <-  function(...){
    ggplot2:::manual_scale('fill', values = location.colors,...)
}
scale_fill_condition <-  function(...){
    ggplot2:::manual_scale('fill', values = condition.colors,...)
}


gene_anno <- gene_anno %>% filter(ensg.id %in% rownames(gene.expr.mat))

gene_symbol <- gene_anno$symbol                                     
names(gene_symbol) <- gene_anno$ensg.id

gene_ensg <- gene_anno$ensg.id
names(gene_ensg) <- gene_anno$symbol


theme_black = function(base_size = 12, base_family = "") {
  
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    
    theme(
      # Specify axis options
      axis.line = element_blank(),  
      axis.text.x = element_text(size = base_size*0.8, color = "white", lineheight = 0.9),  
      axis.text.y = element_text(size = base_size*0.8, color = "white", lineheight = 0.9),  
      axis.ticks = element_line(color = "white", size  =  0.2),  
      axis.title.x = element_text(size = base_size, color = "white", margin = margin(0, 10, 0, 0)),  
      axis.title.y = element_text(size = base_size, color = "white", angle = 90, margin = margin(0, 10, 0, 0)),  
      axis.ticks.length = unit(0.3, "lines"),   
      # Specify legend options
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
      # Specify panel options
      panel.background = element_rect(fill = "black", color  =  NA),  
      panel.border = element_rect(fill = NA, color = "white"),  
      panel.grid.major = element_line(color = "grey35"),  
      panel.grid.minor = element_line(color = "grey20"),  
      panel.margin = unit(0.5, "lines"),   
      # Specify facetting options
      strip.background = element_rect(fill = "grey30", color = "grey10"),  
      strip.text.x = element_text(size = base_size*0.8, color = "white"),  
      strip.text.y = element_text(size = base_size*0.8, color = "white",angle = -90),  
      # Specify plot options
      plot.background = element_rect(color = "black", fill = "black"),  
      plot.title = element_text(size = base_size*1.2, color = "white"),  
      plot.margin = unit(rep(1, 4), "lines")
      
    )
  
}


#  dynamic ui example:https://shiny.rstudio.com/gallery/dynamic-ui.html

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(##theme = "bootstrap.css",
    
    titlePanel("Placenta Gene Expression Atlas at Parturition"),
    
    headerPanel(
        img(src = 'Clusters2.png',width="1400px")##, height = '100px', width = '100px')
    ),
    mainPanel(
           column(8,plotOutput("mainPlot",width = "100%")),
           column(4,plotOutput("sidePlot",width = "100%")),
           column(12,
                tabsetPanel(id="tabs",
                    tabPanel('DE by Condition',DT::dataTableOutput("mytable")),
                    tabPanel('DE by Location',DT::dataTableOutput("Lmytable")),
                    tabPanel('Select Gene',selectInput("gene","Gene:",choices=names(gene_ensg)),hr(),helpText("Select Gene to Display")),
                    tabPanel('Marker Genes',
                             titlePanel("Select marker genes per cluster"),
                               fluidRow(
                                       column(3, wellPanel(
                                                           selectInput("Select_Cluster", "Select Cluster",MClusters)
                                                 )
                                                         )),

                                       column(3, wellPanel(
                                                            selectInput("Select_Gene","Select Gene","")
                                                         ))
                )
           )
       )
    )
  )
# Define server logic required to draw a histogram
server <- function(input, output,session) {

  output$mytable = DT::renderDataTable(##{deg_table},##{marker_table},
                                       ##                                         container = sketch,   # layout (from helper.R)
                           DT:::datatable(deg_table,
                                       selection = "single", # only select 1 row at a time
                                       escape = FALSE,       # Allow HTML in table
##                                       rownames = FALSE,
                                       options = list(
                                         lengthMenu = list(c(10, 25, 50, -1), c(10, 25, 50, "All")),
                                         search = list(regex = TRUE, caseInsensitive = TRUE),
                                         order = list(list(10, 'asc'), list(5, 'desc'))
                                       ),
                                       filter = list(position = 'top', clear = TRUE)
                                       ) ##%>% formatStyle(columns=colnames(deg_table),color='white',background = 'black',target = 'row')
                           
        ) 

output$Lmytable = DT::renderDataTable(##{deg_table},##{marker_table},
                                       ##                                         container = sketch,   # layout (from helper.R)
                           DT:::datatable(Ldeg_table,
                                       selection = "single", # only select 1 row at a time
                                       escape = FALSE,       # Allow HTML in table
##                                       rownames = FALSE,
                                       options = list(
                                         lengthMenu = list(c(10, 25, 50, -1), c(10, 25, 50, "All")),
                                         search = list(regex = TRUE, caseInsensitive = TRUE),
                                         order = list(list(10, 'asc'), list(5, 'desc'))
                                       ),
                                       filter = list(position = 'top', clear = TRUE)

                                       ) ##%>% formatStyle(columns=colnames(deg_table),color='white',background = 'black',target = 'row')
                          
     ) 

    
    output$mainPlot <- renderPlot({
        ##        s = as.character(marker_table$ENSG[input$mytable_rows_selected])
        if (input$tabs=="DE by Condition") {
            s = as.character(deg_table$ENSG[input$mytable_rows_selected])}
        else if (input$tabs=="DE by Location") {s = as.character(Ldeg_table$ENSG[input$Lmytable_rows_selected])}
        else if (input$tabs=="Marker Genes") {s = names(which(gene_symbol==input$Select_Gene))}
        else { s=names(which(gene_symbol==input$gene))}
        print(s)
        if(length(s)==0){
            s="ENSG00000109320";
        }
        gene.expr <- gene.expr.mat[s,rownames(mydf)]
        ##gene.expr <- gene.expr.mat[input$gsel,rownames(mydf)]
        stopifnot(identical(names(gene.expr),row.names(mydf)))
        mydf$gene.expr <- gene.expr        
        mydf %>% arrange(gene.expr) %>% 
            ggplot(aes(x = UMAP1, y = UMAP2, col = gene.expr)) + 
            geom_point(alpha = 1, size = 0.25)  +  
            scale_colour_gradient(low="#e0ecf4",high="#6e016b") + 
            ##scale_colour_gradientn(colors = c('#ffffcc','#ffeda0','#fed976','#feb24c','#fd8d3c','#fc4e2a','#e31a1c','#b10026','#800026')) + 
            facet_grid(Group ~ Location) + 
            ggtitle(paste(s,gene_symbol[s])) +
            theme_bw()
    })
    
    output$sidePlot <- renderPlot({
        if (input$tabs=="DE by Condition") {
            s = as.character(deg_table$ENSG[input$mytable_rows_selected])}
        else if (input$tabs=="DE by Location") {s = as.character(Ldeg_table$ENSG[input$Lmytable_rows_selected])}
        else if (input$tabs=="Marker Genes") {s = names(which(gene_symbol==input$Select_Gene))}
        else { s=names(which(gene_symbol==input$gene))}
        
        #s = as.character(deg_table$ENSG[input$mytable_rows_selected])
        if (input$tabs=="DE by Condition") {
            mycell = as.character(deg_table$Cluster[input$mytable_rows_selected])}
        else if (input$tabs=="DE by Location") {mycell = as.character(Ldeg_table$Cluster[input$Lmytable_rows_selected])}
        #else if (input$tabs=="Marker Genes") {s = names(which(gene_symbol==input$Select_Gene))}
        else { mycell="Clusters"}

#        mycell = as.character(deg_table$Cluster[input$mytable_rows_selected])
        print(s)
        if(length(s)==0){
            s="ENSG00000109320";
            mycell="Macrophage-1";
        }
  if (mycell!="Clusters") {        
        gene.expr <- gene.expr.mat[s,rownames(mydf)]
        ##gene.expr <- gene.expr.mat[input$gsel,rownames(mydf)]
#        stopifnot(identical(names(gene.expr),row.names(mydf)))
        mydf$gene.expr <- gene.expr
        ## Option to clip small values
        ##mydf$gene.expr[gene.expr==0] <- (min(gene.expr[gene.expr>0]))
        mydf3 <- mydf %>% dplyr::filter(NewClsName %in% mycell) %>% add_count(Location,NewClsName,gene.expr>0) %>% dplyr::filter(n>200)            
        gg <- ggplot(mydf3,aes(y=log2(gene.expr),
                               x=reorder(Location, desc(Location)),
                               col=reorder(Group,desc(Group)) )) +
            geom_point(alpha=0.1,position = position_jitterdodge()) +
            geom_boxplot(outlier.colour=NA,fill=NA,na.rm=FALSE,notch=TRUE) +
            labs(x="Location",y=expression(log[2](Expr.)),col="Condition") +
            scale_fill_condition() +
            scale_color_manual(values=condition.colors) + 
            guides(colour = guide_legend(reverse=T)) +
            coord_flip() + 
            theme_bw() +
            ggtitle(paste(s,gene_symbol[s]))
        gg
       }
  else {
           gene.expr <- gene.expr.mat[s,rownames(mydf)]
        ##gene.expr <- gene.expr.mat[input$gsel,rownames(mydf)]
        #stopifnot(identical(names(gene.expr),row.names(mydf)))
           mydf$gene.expr <- gene.expr
          gg <- ggplot(mydf %>% dplyr::filter(gene.expr>0.0),aes(x=log2(gene.expr),y=reorder(NewClsName, gene.expr, FUN = mean),alpha=0.8,fill=Location)) +
         geom_density_ridges(scale = 2, size = 0.25, rel_min_height = 0.01) +
         labs(x=expression(log[2](Expr.)),y="Cell-type",fill="Location",title=paste(s," ",gene_symbol[s])) +
         scale_fill_location() +
         theme_ridges()
      gg           
      }
    })

    observe({

        x<-input$Select_Cluster


#        output$Select_Gene_List <- updateSelectInput({
              updateSelectInput(session, "Select_Gene",
                   #label = paste("Select input label", length(x)),
                   choices = Marker_table %>% filter(cluster==x) %>% select(gene) %>% slice(1:n())
                   
                   #selected = head(choices, 1)
                   )
                                                    
 #             })
      })
    
}

# Run the application 
shinyApp(ui = ui, server = server)
