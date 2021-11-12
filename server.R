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
if (!require('DT', quietly = T)) install.packages('DT');
if (!require('plotly', quietly = T)) install.packages('plotly');
if (!require('htmlwidgets', quietly = T)) install.packages('htmlwidgets');
if (!require('shinyalert', quietly = T)) install.packages('shinyalert');
if (!require('shinyjs', quietly = T)) install.packages('shinyjs');

library(shiny)
library(shinydashboard)
library(DT)
library(plotly)
library(htmlwidgets)
library(shinyalert)
library(shinyjs)
#####################################################################

# Debut
function(input, output) {
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
                else if(identical(good,colnames(df))==FALSE)
                {
                        shinyalert("Column names Error",paste0("Wrong format for column ",which(colnames(df) != good),", it must be : ",good[which(colnames(df) != good)],"\n",collapse=""),type = "error")
                        returnValue()
                } 
        })
        # When the analysis is started
        observeEvent(input$start, {
                req(input$file)
                ext <- tools::file_ext(input$file$datapath)
               
                validate(need(ext == "csv", "Please upload a csv file"))
                # Get data from the uploaded file
                data <- reactive({
                        req(input$file)
                        df <- read.csv(input$file$datapath, sep = ";")
                        df["log2padj"] <- -log2(df["padj"])
                        df
                })
                # Data view
                output$table <- DT::renderDataTable({
                        data()
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
                                 title = "Volcano Plot")
                        
                })
                output$volcano <- renderPlotly2({
                        p <- volcanoPlot()
                        as_widget(p) %>% onRender(addHoverBehavior)
                        
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
                                       title = "MA Plot")
                        
                })
                output$MA <- renderPlotly2({
                        p <- MAPlot()
                        as_widget(p) %>% onRender(addHoverBehavior)
                })
        })
}
