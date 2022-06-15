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

options(timeout = 2000)
if (!require('tibble', quietly = T)) install.packages('tibble');
if (!require('shiny', quietly = T)) install.packages('shiny');
if (!require('shinydashboard', quietly = T)) install.packages('shinydashboard');
if (!require('shinybusy', quietly = T)) install.packages('shinybusy');
if (!require('shinycssloaders', quietly = T)) install.packages('shinycssloaders');
if (!require('shinyjs', quietly = T)) install.packages('shinyjs');
if (!require('shinyalert', quietly = T)) install.packages('shinyalert');
if (!require('BiocManager', quietly = T)) install.packages('BiocManager');
if (!require('Biostrings', quietly = T)) BiocManager::install('Biostrings');
if (!require('biomaRt', quietly = T)) BiocManager::install('biomaRt');
if (!require('biomartr', quietly = T)) install.packages('biomartr');
if (!require('stringr', quietly = T)) install.packages('stringr');
if (!require('stringi', quietly = T)) install.packages('stringi');
if (!require('DT', quietly = T)) install.packages('DT');
if (!require('plotly', quietly = T)) install.packages('plotly');
if (!require('htmlwidgets', quietly = T)) install.packages('htmlwidgets');
if (!require('clusterProfiler', quietly = T)) BiocManager::install('clusterProfiler');
if (!require('pathview', quietly = T)) BiocManager::install('pathview');
if (!require('ReactomePA', quietly = T)) BiocManager::install('ReactomePA');
if (!require('AnnotationHub', quietly = T)) BiocManager::install('AnnotationHub');
if (!require('enrichplot', quietly = T)) BiocManager::install("enrichplot");
if (!require('topGO', quietly = T)) BiocManager::install("topGO");
if (!require('ggplot2', quietly = T)) install.packages("ggplot2");
if (!require('DOSE', quietly = T)) BiocManager::install("DOSE");

library(tibble)
library(shiny)
library(shinydashboard)
library(shinybusy)
library(shinycssloaders)
library(shinyjs)
library(shinyalert)
library(biomaRt)
library(Biostrings)
library(biomartr)
library(stringr)
library(stringi)
library(DT)
library(plotly)
library(htmlwidgets)
library(clusterProfiler)
library(pathview)
library(ReactomePA)
library(AnnotationHub)
library(enrichplot)
library(topGO)
library(ggplot2)
library(DOSE)

R.utils::setOption("clusterProfiler.download.method","auto")

#####################################################################
#to run for each R session for connexion to ensembl with biomart
#allow R to try other ensembl mirrors site if connexion doesn't work
httr::set_config(httr::config(ssl_verifypeer = FALSE))

####################################################### COMMANDES A UTILISER #######
# validate(need(dim(frame)[1]>0, "No ID found :\n Are you sure you have selected the right organism ?"))
############################################################## Preparation de la liste des organismes disponibles ###########

ah<-AnnotationHub()
liste_all <- query(ah, c("OrgDb", "maintainer@bioconductor.org"))
data_annot=as.data.frame(cbind(liste_all$species,liste_all$title))
data_db=data_annot[grep(".*db.*", data_annot$V2),]
data_db_ordered=data_db[order(data_db$V2),]
data_db_ordered$V2<-gsub(".sqlite","",data_db_ordered$V2)
choices=setNames(data_db_ordered$V1,data_db_ordered$V1)
rm(data_db,ah,data_annot,liste_all)

#####################################################################
# Usefull functions
#####################################################################
download <- function(title,plotname){
        downloadHandler(
                filename = function() { title },
                content = function(file) {
                        ggsave(file, plot = plotname, device = "png")
                        }
                )
}

############################################################################################################################

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
                # load organism library
                organism <- reactive({
                        req(input$select_organism)
                        data_db_ordered[which(data_db_ordered$V1==input$select_organism),2]
                })
                req(organism())
                BiocManager::install(organism(), update = F)
                # Get data from the uploaded file
                data_prev <- reactive({
                        req(input$file)
                        df <- read.csv(input$file$datapath, sep = ";")
                        df["log2padj"] <- -log2(df["padj"])
                        
                        # df <- df[df$padj < input$pvalue,]
                        
                        if(startsWith(df$ID[1], "ENS")){
                                df
                        }
                        else{
                                shinyalert("Wrong ID format !","Need Ensembl ID", type="error")
                                return(NULL)
                        }
                })
                
                data <- reactive({
                        req(data_prev())
                        df <- data_prev()
                        # View(df)
                        # print(sum(df$padj < input$pvalue & (df$log2FC < -input$FC | df$log2FC > input$FC)))
                        df <- df[df$padj < input$pvalue & (df$log2FC < -input$FC | df$log2FC > input$FC),]
                        df
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
                        req(data_prev())
                        df <- data_prev()
                        df["DEG"] <- ifelse(df["log2FC"] >= input$FC & df["padj"] <= input$pvalue, "UP", 
                                        ifelse(df["log2FC"] <= -input$FC & df["padj"] <= input$pvalue,"DOWN",
                                        "NO")
                                )
                        mycolors <- c("blue", "orange", "grey")
                        names(mycolors) <- c("DOWN", "UP", "NO")
                        plot_ly(data = df, x = df$log2FC, y = df$log2padj, text = df$GeneName, mode = "markers", type = "scatter", color = df$DEG, colors = mycolors) %>%
                          layout(xaxis = list(title = 'log2(FC)'),
                                 yaxis = list(title = 'log2(padj)'),
                                 title = list(text = paste("Volcano Plot ", input$select_organism), x = 0))
                        
                })
                output$volcano <- renderPlotly2({
                        req(volcanoPlot())
                        p <- volcanoPlot()
                        as_widget(p) %>% 
                          onRender(addHoverBehavior) %>%
                          config(modeBarButtons = list(list("zoomIn2d"), list("zoomOut2d"), list("select2d"), list("resetScale2d"), list("toImage")))
                        
                })
                # Print genes name when the curser hovers on the points
                output$hover <- renderText({
                        input$hover_data
                })
          
                # Coord selected box in graph
                coord <- reactive({
                  req(volcanoPlot())

                  event_register(volcanoPlot(), 'plotly_brushed')
                  event_data(event = 'plotly_brushed')
                })

                # Data view
                output$table <- DT::renderDataTable({
                  req(data())
                  df <- data()

                  if(is.null(coord()))
                  {
                    return(df)
                  }
                  else
                  {
                    x <- coord()$x
                    y <- coord()$y

                    return(df[ which(df$log2FC >= x[1] & df$log2FC <= x[2] &
                        df$log2padj >= y[1] & df$log2padj <= y[2]), ])
                  }
                },
                extensions = 'Buttons',
                options = list(
                        fixedColumns = TRUE,
                        autoWidth = FALSE,
                        ordering = TRUE,
                        scrollX = TRUE,
                        dom = 'Bfrtip',
                        buttons = c('csv', 'excel')),
                class = "display"
                )
                MAPlot <- reactive({
                        req(data_prev())
                        df <- data_prev()
                        df["DEG"] <- ifelse(df["log2FC"] >= input$FC & df["padj"] <= input$pvalue, "UP", 
                                            ifelse(df["log2FC"] <= -input$FC & df["padj"] <= input$pvalue,"DOWN",
                                                   "NO")
                        )
                        mycolors <- c("blue", "orange", "grey")
                        names(mycolors) <- c("DOWN", "UP", "NO")
                        plot_ly(data = df, x = log2(df$baseMean), y = df$log2FC, text = df$GeneName, mode = "markers", type = "scatter", color = df$DEG, colors = mycolors) %>%
                                layout(xaxis = list(title = "log2(baseMean)"),
                                       yaxis = list(title = "log2(FC)"),
                                       title = list(text = paste("MA Plot ", input$select_organism), x = 0))
                        
                })
                output$MA <- renderPlotly2({
                        req(MAPlot)
                        p <- MAPlot()
                        as_widget(p) %>% 
                          onRender(addHoverBehavior) %>%
                          config(modeBarButtons = list(list("zoomIn2d"), list("zoomOut2d"), list("select2d"), list("resetScale2d"), list("toImage")))

                })
#####################################################################################
##                        Onglet GO Term Enrichment                                ##
#####################################################################################
                organism <- reactive({
                        req(input$select_organism)
                        data_db_ordered[which(data_db_ordered$V1==input$select_organism),2]
                })
                data_go <- reactive({
                        req(data())
                        if(input$go_filter == 'DEG+'){
                               return(data()[data()$log2FC > 0,])
                        }
                        else if(input$go_filter == 'DEG-'){
                                return(data()[data()$log2FC < 0,])
                        }
                        else {
                                return(data())
                        }
                        
                })
                gene_list <- reactive({
                        req(data_prev())
                        return(get_gene_list(data_prev()))
                })
                ego <- reactive({
                        req(data_go(), organism(), input$go_ontology, input$go_pvalue)
                        return(get_ego(data_go(), organism(), input$go_ontology, input$go_pvalue))
                })
                output$table_go_ora <- DT::renderDataTable({
                  req(ego())
                  df <- as.data.frame(ego()) 
                  df <-df[ , -which(names(df) %in% c("geneID"))]
                  Links <- paste0('<a href="https://www.ebi.ac.uk/QuickGO/GTerm?id=',df$ID,'" target="_blank">GO link</a>')
                  df <- df %>% add_column(Links, .after ="ID")
                  return(df)
                  
                },
                extensions = 'Buttons',
                rownames = FALSE,
                escape = FALSE,
                options = list(
                  fixedColumns = TRUE,
                  autoWidth = FALSE,
                  ordering = TRUE,
                  scrollX = TRUE,
                  dom = 'Bfrtip',
                  buttons = c('csv', 'excel')),
                class = "display"
                )
                gsego <- reactive({
                        req(gene_list(), organism(), input$go_ontology, input$go_pvalue)
                        x <- get_gsego(gene_list(), organism(), input$go_ontology, input$go_pvalue)
                        validate(need(expr = (! isEmpty(as.data.frame(x))), message = "No gene found for this p-value cutoff !"))
                        return(x)
                })
                output$table_go_gsea <- DT::renderDataTable({
                  req(gsego())
                  df <- as.data.frame(gsego()) 
                  df <-df[ , -which(names(df) %in% c("core_enrichment"))]
                  Links <- paste0('<a href="https://www.ebi.ac.uk/QuickGO/GTerm?id=',df$ID,'" target="_blank">GO link</a>')
                  df <- df %>% add_column(Links, .after ="ID")
                  return(df)
                },
                extensions = 'Buttons',
                rownames = FALSE,
                escape = FALSE,
                options = list(
                  fixedColumns = TRUE,
                  autoWidth = FALSE,
                  ordering = TRUE,
                  scrollX = TRUE,
                  dom = 'Bfrtip',
                  buttons = c('csv', 'excel')),
                class = "display"
                )
                #barplot ORA
                GO_ORA_barplot_input <- reactive({
                        req(ego())
                        barplot(ego(), showCategory = 10)
                })
                output$GO_ORA_barplot <- renderPlot({
                        print(GO_ORA_barplot_input())
                })
                #goplot ORA
                GO_ORA_goplot_input <- reactive({
                        req(ego())
                        goplot(ego(), showCategory = 10)
                })
                output$GO_ORA_goplot <- renderPlot({
                        req(GO_ORA_barplot_input())
                        print(GO_ORA_goplot_input())
                })
                #dotplot GSEA
                GO_GSEA_dotplot_input <- reactive({
                        req(gsego())
                        dotplot(gsego(), showCategory=10, split=".sign", font.size = 7) + 
                                facet_grid(.~.sign)
                })
                output$GO_GSEA_dotplot <- renderPlot({
                        req(GO_GSEA_dotplot_input())
                        print(GO_GSEA_dotplot_input())
                })
                #plot GSEA
                output$GO<- renderUI({
                        req(gsego())
                        table = as.data.frame(gsego())
                        choices=setNames(1:nrow(table),table$Description)
                        selectInput("select_GO", label = "GO term of interest for GSEA Plot",
                                    choices=choices,
                                    selected = 1
                        )
                })
                GO_GSEA_plot_input <- reactive({
                        req(gsego())
                        gseaplot(gsego(), by = "all", title = gsego()$Description[as.numeric(input$select_GO)], geneSetID = as.numeric(input$select_GO))
                })
                output$GO_GSEA_plot <- renderPlot({
                        req(GO_GSEA_plot_input())
                        print(GO_GSEA_plot_input())
                })
                #boutons downloads
                output$download_go_barplot <- download(
                        paste('barplotORA.png', sep=''),
                        GO_ORA_barplot_input()
                        )
                output$download_go_dotplot <- download(
                        paste('dotplotORA.png', sep=''),
                        GO_ORA_goplot_input()
                )
                output$download_go_gseaplot <- download(
                        paste('ploGSEA.png', sep=''),
                        GO_GSEA_dotplot_input()
                )
                output$download_go_goplot <- download(
                        paste('goplotGSEA.png', sep=''),
                        GO_GSEA_plot_input()
                )
                
#####################################################################################
##                        Onglet Pathway Enrichment                                ##
#####################################################################################
                # get data depending on the DEG selection
                data_path <- reactive({
                  req(data())
                  if(input$path_filter == 'DEG+'){
                    return(data()[data()$log2FC > 0,])
                  }
                  else if(input$path_filter == 'DEG-'){
                    return(data()[data()$log2FC < 0,])
                  }
                  else {
                    return(data())
                  }
                  
                })
                kegg_gene_list <- reactive({
                        req(data_path())
                        return(get_kegg_gene_list(data_path(), organism()))
                })
                ekk <- reactive({
                        req(kegg_gene_list())
                        return(get_ekk(kegg_gene_list()[[1]], input$path_pvalue, organism()))
                })
                output$table_ekk <- DT::renderDataTable({
                  req(ekk())
                  df <- as.data.frame(ekk())
                  df <-df[ , -which(names(df) %in% c("geneID"))]
                  Links <- paste0('<a href="https://www.genome.jp/entry/',df$ID,'" target="_blank">KEGG link</a>')
                  df <- df %>% add_column(Links, .after ="ID")
                  return(df)
                },
                extensions = 'Buttons',
                rownames = FALSE,
                escape = FALSE,
                options = list(
                  fixedColumns = TRUE,
                  autoWidth = FALSE,
                  ordering = TRUE,
                  scrollX = TRUE,
                  dom = 'Bfrtip',
                  buttons = c('csv', 'excel')),
                class = "display"
                )
                gsekk <- reactive({
                        req(kegg_gene_list())
                        x <- get_gsekk(kegg_gene_list()[[2]], input$path_pvalue, organism())
                        validate(need(expr = (! isEmpty(as.data.frame(x@result))), message = "No gene found for this p-value cutoff !"))
                        return(x)
                })
                output$table_gsekk <- DT::renderDataTable({
                  req(gsekk())
                  df <- as.data.frame(gsekk()@result)
                  df <-df[ , -which(names(df) %in% c("core_enrichment"))]
                  Links <- paste0('<a href="https://www.genome.jp/entry/',df$ID,'" target="_blank">KEGG link</a>')
                  df <- df %>% add_column(Links, .after ="ID")
                  return(df)
                },
                extensions = 'Buttons',
                rownames = FALSE,
                escape = FALSE,
                options = list(
                  fixedColumns = TRUE,
                  autoWidth = FALSE,
                  ordering = TRUE,
                  scrollX = TRUE,
                  dom = 'Bfrtip',
                  buttons = c('csv', 'excel')),
                class = "display"
                )
                #barplot ORA
                path_barplot_input <- reactive({
                        req(ekk())
                        barplot(ekk(), showCategory = 10)
                })
                output$path_barplot <- renderPlot({
                        req(path_barplot_input())
                        print(path_barplot_input())
                })
                output$pathway<- renderUI({
                  req(gsekk())
                  choices=setNames(1:nrow(as.data.frame(gsekk())),gsekk()$Description)
                  selectInput("select_path", label = "Pathway of interest for GSEA Plot",
                              choices=choices,
                              selected = 1
                  )
                })
                #plot GSEA
                path_gseaplot_input <- reactive({
                        req(gsekk())
                        gseaplot(gsekk(), by = "all", title = gsekk()$Description[as.numeric(input$select_path)], geneSetID = as.numeric(input$select_path))
                })
                output$path_gseaplot <- renderPlot({
                        req(path_gseaplot_input())
                        print(path_gseaplot_input())
                })
                #dotplot GSEA
                path_dotplot_input <- reactive({
                        req(gsekk())
                        dotplot(gsekk(), showCategory = 10, title = "Enriched Pathways" , split=".sign") + 
                                facet_grid(.~.sign)
                })
                output$path_dotplot <- renderPlot({
                        req(path_dotplot_input())
                        print(path_dotplot_input())
                })
                #kegg_gene_list gsea
                output$pathplot_list<- renderUI({
                        req(gsekk())
                        choices=setNames(gsekk()$ID,gsekk()$Description)
                        print(choices)
                        selectInput("select_path2", label = "Pathway of interest for Pathway plot",
                                    choices=choices,
                                    selected = 1
                        )
                })
                path_pathplot_input <- reactive({
                        req(input$select_path2)
                        pathview(gene.data=kegg_gene_list()[[2]], pathway.id=gsekk()[input$select_path2]$ID, species = db_to_organism(organism()), kegg.dir = 'www')
                        #print(paste(gsekk()[1]$ID, ".pathview.png", sep = ""))
                })
                output$path_pathplot <- renderImage({
                        req(path_pathplot_input())
                        print(path_pathplot_input())
                        list(src = paste(gsekk()[input$select_path2]$ID, ".pathview.png", sep = ""),
                            alt = "No pathview image found.",
                            width = '20%')
                        #knitr::include_graphics(paste(gsekk()[2]$ID, ".pathview.png", sep = ""))
                }, deleteFile = FALSE)
                #boutons downloads
                output$download_path_barplot <- download(
                        paste('path_barplotORA.png', sep=''),
                        path_barplot_input()
                )
                output$download_path_dotplot <- download(
                        paste('path_dotplotGSEA.png', sep=''),
                        path_gseaplot_input()
                )
                output$download_path_gseaplot <- download(
                        paste('path_plotGSEA.png', sep=''),
                        path_dotplot_input()
                )
                #output$download_path_goplot <- downloadHandler(
                #        filename = paste('path_goplotGSEA.png', sep=''),
                #         contentType = "png",
                #         content = function(file){
                #                 file.copy(path_pathplot_input(), file)
                #                 }
                # )
                
#####################################################################################
##                        Onglet Domain Enrichment                                ##
#####################################################################################
                # get data based on the DEG selection
                data_domain <- reactive({
                        req(data())
                        if(input$domain_filter == 'DEG+'){
                                return(data()[data()$log2FC > 0,])
                        }
                        else if(input$domain_filter == 'DEG-'){
                                return(data()[data()$log2FC < 0,])
                        }
                        else {
                                return(data())
                        }
                        
                })
                # Domain enrichment results datatable
                domains_ORA_results <- reactive({
                        req(data_domain(), organism())
                        return(get_table_ORA_domains(data_domain(), organism(), input$domain_pvalue))
                })
                output$domain_ORA_datatable <- DT::renderDataTable({
                        req(domains_ORA_results())
                        df <- as.data.frame(domains_ORA_results())
                        Links <- paste0('<a href="http://www.ebi.ac.uk/interpro/entry/InterPro/',df$interproID,'" target="_blank">Interpro link</a>')
                        df <- df %>% add_column(Links, .after ="interproID")
                        return(df)
                        
                },
                extensions = 'Buttons',
                rownames = FALSE,
                escape = FALSE,
                options = list(
                        fixedColumns = TRUE,
                        autoWidth = FALSE,
                        ordering = TRUE,
                        scrollX = TRUE,
                        dom = 'Bfrtip',
                        buttons = c('csv', 'excel')),
                class = "display"
                )
                # barplot
                domain_barplot_input <- reactive({
                        req(domains_ORA_results())
                        return(domains_ORA_barplot(domains_ORA_results()))
                })
                output$domain_barplot <- renderPlotly({
                        req(domain_barplot_input())
                        return(ggplotly(domain_barplot_input()))
                })
                # dotplot
                domain_dotplot_input <- reactive({
                        req(domains_ORA_results())
                        return(domains_ORA_dotplot(domains_ORA_results()))
                })
                output$domain_dotplot <- renderPlotly({
                        req(domain_dotplot_input())
                        return(ggplotly(domain_dotplot_input()))
                })
        })
        observe({
                req(input$select_path2)
                image.path = file.path(getwd(), paste0(input$select_path2, ".pathview.png"))
                print(image.path)
                if(file.exists(image.path)) file.show(image.path)
                else showModal(modalDialog(h4("The file does not exists."), easyClose = T, footer = modalButton("Ok")))
        }) %>% bindEvent(input$png)
}

