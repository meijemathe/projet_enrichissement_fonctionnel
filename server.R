# Projet d'annotation fonctionnel avec R Shiny 

# Auteurs :

# MATHE Meije : meije.mathe@univ-rouen.fr
# PETY Solene : solene.pety@etu.univ-rouen.fr
# LETERRIER Bryce : bryce.leterrier@univ-rouen.fr
# OLLIVIER Louis : louis.ollivier@etu.univ-rouen.fr

# M2.1 BIMS - Univ. Rouen Normandie 

# 2021 - 2022

#####################################################################
# Initialisation des packages 

if (!require('shiny', quietly = T)) install.packages('shiny');
if (!require('shinydashboard', quietly = T)) install.packages('shinydashboard');
if (!require('shinybusy', quietly = T)) install.packages('shinybusy');
if (!require('shinycssloaders', quietly = T)) install.packages('shinycssloaders');
if (!require('shinyjs', quietly = T)) install.packages('shinyjs');
if (!require('shinyalert', quietly = T)) install.packages('shinyalert');
if (!require('BiocManager', quietly = T)) install.packages('BiocManager');
if (!require('biomaRt', quietly = T)) BiocManager::install('biomaRt');
if (!require('biomaRt', quietly = T)) install.packages('biomartr');
if (!require('biomaRt', quietly = T)) BiocManager::install('pathview');
if (!require('biomaRt', quietly = T)) BiocManager::install('clusterProfiler');
if (!require('stringr', quietly = T)) install.packages('stringr');
if (!require('stringi', quietly = T)) install.packages('stringi');
if (!require('DT', quietly = T)) install.packages('DT');
if (!require('plotly', quietly = T)) install.packages('plotly');
if (!require('htmlwidgets', quietly = T)) install.packages('htmlwidgets');
if (!require('biomaRt', quietly = T)) BiocManager::install('reactomePA');


library(pathview)
library(reactomePA)
library(clusterProfiler)
library(shiny)
library(shinydashboard)
library(shinybusy)
library(shinycssloaders)
library(shinyjs)
library(shinyalert)
library(biomaRt)
library(biomartr)
library(stringr)
library(stringi)
library(DT)
library(plotly)
library(htmlwidgets)

#####################################################################
#to run for each R session for connexion to ensembl with biomart
#allow R to try other ensembl mirrors site if connexion doesn't work
httr::set_config(httr::config(ssl_verifypeer = FALSE))

############################################################## Biomart, convertion ID en fonction organisme sélectionné###########
name_short <- paste(sapply(strsplit(listGenomes(db="ensembl"), "_"), "[", 1)
                    %>%stri_extract_first_regex(".{1}"),
                    sapply(strsplit(listGenomes(db="ensembl"), "_"), "[", 2),
                    sep="")%>%
        sort()%>%
        unique()

mart<- useMart("ensembl")
dataset<-biomaRt::listDatasets(mart)
name_short2=name_short[grep(paste(sapply(strsplit(dataset$dataset, "_"), "[", 1),collapse="|"), name_short)]
dataset=dataset[grep(paste(sapply(strsplit(dataset$dataset, "_"), "[", 1),collapse="|"), name_short2),]

#ui
# Adaptation liste en fonction deux listes récupérées par Biomart qui sont de tailles différentes. Conversion dans server.R

liste<-listGenomes(db="ensembl") %>%
        str_extract("[^_]+_[^_]+") %>%
        gsub(pattern = "_",replacement = " ") %>%
        unique() %>%
        str_to_title()

frame_names=as.data.frame(cbind(name_short,liste))
frame_names=frame_names[order(frame_names$liste),]

###Récupère dataframe avec id dataset à utiliser, description espèce et version génome
mart<- useMart("ensembl")
dataset<-biomaRt::listDatasets(mart)

frame_names=frame_names[grep(paste(sapply(strsplit(dataset$dataset, "_"), "[", 1),collapse="|"), frame_names$name_short),]

### Préparation liste pour menu déroulant, enlever tiret bas, unique, majuscule début mot
choices=setNames(frame_names$liste,frame_names$liste)
############################################################################################################################""

# Debut
function(input, output) {
        
        output$choices<- renderUI({
                selectInput("select_organism", label = "Organism of interest",
                            choices=choices,
                            selected = 1
                )
        })
        
        shinyjs::disable("start")
        observeEvent(input$file, {
          req(input$file)
          #test extension
          ext <- tools::file_ext(input$file$datapath)
          validate(need(ext == "csv", "Please upload a csv file"))
          # Vector with good columns names
          good <-c("GeneName", "ID", "baseMean", "log2FC", "pval", "padj")
          # Get data from the uploaded file
          df <- read.csv(input$file$datapath, sep = ";")
          #Check data format, colnames and GeneID
          if(ncol(df)!=6)
          {
            shinyalert("Column Error","Uploaded Data has the wrong number of columns\n Need 6 columns : GeneName, ID, baseMean, log2FC, pval, padj",
                       type="error")
            returnValue()
          }

          else if(identical(tolower(good),tolower(colnames(df)))==FALSE)
          {
            shinyalert("Column names Error",paste0("Wrong format for column ",which(colnames(df) != good),", it must be : ",good[which(colnames(df) != good)],"\n",collapse=""),type = "error")
            returnValue()
          }
          else {
                enable("start")
          }
        })
        # When the analysis is started
        observeEvent(input$start, {
                
                req(input$file)
                # Get data from the uploaded file
                data <- reactive({
                        req(input$file)
                        df <- read.csv(input$file$datapath, sep = ";")
                        df["log2padj"] <- -log2(df["padj"])
                        if(startsWith(df$ID[1], "ENS")){
                                selected<-grep(paste(stri_extract_first_regex(input$select_organism,".{1}"),unlist(strsplit(input$select_organism," "))[2],sep=""),dataset$dataset,ignore.case = TRUE,value=TRUE)
                                mart=useDataset("mmusculus_gene_ensembl",mart=mart)   
                                genes <- getBM(filters = "ensembl_gene_id",
                                               attributes = c("ensembl_gene_id","entrezgene_id"),
                                               values = df$ID, 
                                               mart = mart)
                                #changer ordre colonnes data frame
                                frame<-merge(df, genes, by.x ="ID", by.y = "ensembl_gene_id")
                                frame<-frame[,c("GeneName", "ID", "entrezgene_id","baseMean", "log2FC", "pval", "padj")]
                                colnames(frame)<-c("GeneName", "ID", "EntrezID","baseMean", "log2FC", "pval", "padj")
                                frame
                        }
                        else{
                             df
                        }
                })
                
                # Data view
                output$table <- DT::renderDataTable({
                        data2()
                })
                # Create the action when the curser hovers on the plot
                addHoverBehavior <- "function(el, x){
                        el.on('plotly_hover', function(data){
                                var infotext = data.points.map(function(d){
                                        console.log(d)
                                        return (d.data.name[d.pointNumber]+': x= '+d.x+', y= '+d.y.toPrecision(3));
                                });
                                console.log(infotext)
                                Shiny.onInputChange('hover_data', infotext)
                        })
                }"
                renderPlotly2 <- function (expr, env = parent.frame(), quoted = FALSE){
                        if (!quoted) {
                                expr <- substitute(expr)
                        }
                        shinyRenderWidget(expr, plotlyOutput, env, quoted = TRUE)
                }
                # Generate the volcano plot
                volcanoPlot <- reactive({
                        df <- data()
                        df["DEG"] <- ifelse(df["log2FC"] >= input$FC & df["padj"] <= input$pvalue, "UP", 
                                        ifelse(df["log2FC"] <= -input$FC & df["padj"] <= input$pvalue,"DOWN",
                                        "NO")
                                )
                        mycolors <- c("blue", "orange", "grey")
                        names(mycolors) <- c("DOWN", "UP", "NO")
                        plot_ly(data = df, x = df$log2FC, y = df$log2padj, text = df$GeneName, mode = "markers", color = df$DEG, colors = mycolors) %>%
                          layout(xaxis = list(title = 'log2(FC)'),
                                 yaxis = list(title = 'log2(padj)'),
                                 title = list(text = paste("Volcano Plot ", input$select_organism), x = 0))                        
                })
                output$volcano <- renderPlotly2({
                        p <- volcanoPlot()
                        as_widget(p) %>% 
                          onRender(addHoverBehavior) %>%
                          config(modeBarButtons = list(list("zoomIn2d"), list("zoomOut2d"), list("select2d"), list("resetScale2d"), list("toImage")))
                        
                })
                # Print genes name when the curser hovers on the points
                output$hover <- renderText({
                        input$hover_data
                })
                MAPlot <- reactive({
                        df <- data()
                        df["DEG"] <- ifelse(df["log2FC"] >= input$FC & df["padj"] <= input$pvalue, "UP", 
                                            ifelse(df["log2FC"] <= -input$FC & df["padj"] <= input$pvalue,"DOWN",
                                                   "NO")
                        )
                        mycolors <- c("blue", "orange", "grey")
                        names(mycolors) <- c("DOWN", "UP", "NO")
                        plot_ly(data = df, x = log2(df$baseMean), y = df$log2FC, text = df$GeneName, mode = "markers", color = df$DEG, colors = mycolors) %>%
                                layout(xaxis = list(title = "log2(baseMean)"),
                                       yaxis = list(title = "log2(FC)"),
                                       title = list(text = paste("MA Plot ", input$select_organism), x = 0))                        
                })
                output$MA <- renderPlotly2({
                        p <- MAPlot()
                        as_widget(p) %>% 
                          onRender(addHoverBehavior) %>% 
                          config(modeBarButtons = list(list("zoomIn2d"), list("zoomOut2d"), list("select2d"), list("resetScale2d"), list("toImage")))
                })
        })
}
