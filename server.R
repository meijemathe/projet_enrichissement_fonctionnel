#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggplot2)

server <- function(input, output){
    
    observeEvent(input$go, {

        # Data input managment
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
        
        # Volcano plot managment
        volcanoPlot <- reactive({
            df <- data()
            
            df["DEG"] <- ifelse(df["log2FC"] >= input$foldChangeIn & df["padj"] <= input$pvalueIn, "DEG+", 
                                ifelse(df["log2FC"] <= -input$foldChangeIn & df["padj"] <= input$pvalueIn,"DEG-",
                                       "Non regulated")
            )
            
            View(df)
            
            ggplot(data = df, aes(x=log2FC, y=log2padj, col=DEG)) + geom_point()
        })
        
        # Volcano plot view
        output$volcano <- renderPlot({
            volcanoPlot()
        }) 
        
        # Volcano plot download
        output$downloadVolcano <- downloadHandler(
            filename = "volcano.png",
            content = function(file) {
                ggsave(file,volcanoPlot())
            }
        )
    })
}


