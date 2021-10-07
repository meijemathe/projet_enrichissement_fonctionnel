library(shiny)
library(plotly)
library(shinydashboard)
library(htmlwidgets)

function(input, output) {
        # When the analysis is started
        observeEvent(input$start, {
                # Get data from the uploaded file
                data <- reactive({
                        read.csv(file = input$file$datapath, sep = ";", header = TRUE)
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
                output$volcano <- renderPlotly2({
                        x <- data()$log2FC
                        y <- data()$pval
                        d <- data.frame(x,y)
                        print(d$y)
                        d$genes <- data()$GeneName
                        d$diffexpressed <- "NO"
                        d$diffexpressed[d$x > input$FC & d$y < input$pvalue] <- "UP"
                        d$diffexpressed[d$x < -input$FC & d$y < input$pvalue] <- "DOWN"
                        mycolors <- c("blue", "red", "black")
                        names(mycolors) <- c("DOWN", "UP", "NO")
                        p <- plot_ly(data = d, x = x, y = -log2(y), text = d$genes, mode = "markers", color = d$diffexpressed, colors = mycolors)
                        as_widget(p) %>% onRender(addHoverBehavior)
                })
                # Print genes name when the curser hovers on the points
                output$hover <- renderText({
                        input$hover_data
                })
                output$MA <- renderPlotly2({
                        x <- data()$baseMean
                        y <- data()$pval
                        d <- data.frame(x,y)
                        print(d$y)
                        d$genes <- data()$GeneName
                        d$diffexpressed <- "NO"
                        d$diffexpressed[d$x > input$FC & d$y < input$pvalue] <- "UP"
                        d$diffexpressed[d$x < -input$FC & d$y < input$pvalue] <- "DOWN"
                        mycolors <- c("blue", "red", "black")
                        names(mycolors) <- c("DOWN", "UP", "NO")
                        p <- plot_ly(data = d, x = -log2(y), y = x, text = d$genes, mode = "markers", color = d$diffexpressed, colors = mycolors)
                        as_widget(p) %>% onRender(addHoverBehavior)
                })
        })
}
