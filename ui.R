## M2.1 BIMS / Promo 2021-2023
## LETERRIR Bryce
## Universite de Rouen

library(shiny)
library(shinydashboard)
# library(DT)

# Define UI for application that draws a histogram
ui <- dashboardPage(
    dashboardHeader(
        title = "Enrichissement fonctionnel",
        titleWidth = 500
    ),
    
    ## ---
    ## SIDEBAR
    ## ---
    
    dashboardSidebar(            
        fileInput("file", label = h2("File input")),
        selectInput(inputId = "organism",label = "Organism",choices = c("Arabidopsis thaliana", "Homo sapiens")),
        
        actionButton(inputId = "go",label = "Run Application"),
        sidebarMenu(
            menuItem("Whole Data Inspection",tabName = "WDInspection"),
            menuItem("GO Term Enrichment",tabName = "GOEnrich"),
            menuItem("Pathway Enrichment",tabName = "PathEnrich"),
            menuItem("Protein Domain Enrichment",tabName = "PDEnrich"),
            menuItem("About",tabName = "About")
        )
    ),
    
    ## ---
    ## BODY
    ## ---
    
    dashboardBody(
        tabItems(
            tabItem(tabName = "WDInspection",
                    
                    # Box valcano plot
                    box(status = "warning", solidHeader = TRUE, collapsible = FALSE, title = "Volcano Plot",
                        plotOutput(outputId = "volcano"),
                    ),
                    
                    # Box input options
                    box(status = "warning", solidHeader = TRUE, collapsible = FALSE, title = "Options",
                        sliderInput(inputId = "pvalueIn",
                                    label = "p-value cutoff from input", 
                                    value = 0.05, 
                                    min = 0, max = 1, 
                                    step = 0.01
                        ),
                        
                        sliderInput(inputId = "foldChangeIn",label = "log2 FoldChange cutoff from input", 
                                    value = 1, 
                                    min = 0, max = 5, 
                                    step = 0.5
                        ), 
                        
                        downloadButton(outputId = "downloadVolcano",
                                     label = "Download Volcano Plot", 
                                     icon = icon("fas fa-download")
                        ), 
                        
                        actionButton(inputId = "downloadMA",
                                     label = "Download MA Plot", 
                                     icon = icon("fas fa-download")
                        )
                    ),
                    
                    # Box MA plot
                    box(status = "warning", solidHeader = TRUE, collapsible = FALSE, title = "MA Plot",
                        
                    ),
                    
                    # Box output options
                    box(status = "warning", solidHeader = TRUE, collapsible = FALSE, title = "Set up for subsequent analyses",
                        sliderInput(inputId = "pvalueOut",
                                    label = "p-value cutoff for subsequent analyses", 
                                    value = 0.05, 
                                    min = 0, max = 1, 
                                    step = 0.01,)
                    ),
                    
                    DT::dataTableOutput(outputId = "table")
            ),
            
            tabItem(tabName = "GOEnrich",
                    box(status = "warning", solidHeader = TRUE, collapsible = FALSE, title = "Volcano plot",width = 12,
                        plotOutput(outputId = "show4")
                    )
            ),
            
            tabItem(tabName = "PathEnrich",
                    # Box colnames
                    box(status = "warning", solidHeader = TRUE, collapsible = FALSE, title = "Column names",
                        verbatimTextOutput(outputId = "colNames"),
                        verbatimTextOutput(outputId = "selectedOrganism")
                    ),
            ),
            
            tabItem(tabName = "PDEnrich",
                   
            ),
            tabItem(tabName = "About",
                    
            )
        )
    )
)




