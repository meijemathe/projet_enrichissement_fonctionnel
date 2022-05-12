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
if (!require('shiny', quietly = T)) install.packages('shiny');
if (!require('shinydashboard', quietly = T)) install.packages('shinydashboard');
if (!require('shinybusy', quietly = T)) install.packages('shinybusy');
if (!require('shinycssloaders', quietly = T)) install.packages('shinycssloaders');
if (!require('shinyjs', quietly = T)) install.packages('shinyjs');
if (!require('shinyalert', quietly = T)) install.packages('shinyalert');
if (!require('BiocManager', quietly = T)) install.packages('BiocManager');
if (!require('biomaRt', quietly = T)) BiocManager::install('biomaRt');
if (!require('Biostrings', quietly = T)) BiocManager::install('Biostrings');
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

#####################################################################


box2 <- function(...){
        box(
                status = "danger",
                solidHeader = TRUE,
                width = 12,
                ...
        )
}


ui = dashboardPage(
        skin = "red",
        dashboardHeader(title = "EnF'R", titleWidth = 350),
        # Sidebar content
        dashboardSidebar(
                #width sidebar
                width = 350,
                sidebarMenu(
                        # File selection
                        fileInput("file", multiple = FALSE, label = "File input", accept = c(
                                "text/csv",
                                "text/comma-separated-values,text/plain",
                                ".csv")),
                        # verbatimTextOutput("file"),
                        # Select box : organism selection
                        uiOutput("choices"),
                        # verbatimTextOutput("select_organism"),
                        shinyjs::useShinyjs(),
                        # Start Button
                        actionButton("start", label = "Start analysis"),
                        hr(),
                        # Menu 
                        menuItem("Whole Data Inspection", tabName = "whole_data_inspection", icon = icon("database")),
                        menuItem("GO Term Enrichment", tabName = "GO_enrichment", icon = icon("sitemap")),
                        menuItem("Pathway Enrichment", tabName = "path_enrichment", icon = icon("chart-pie")),
                        menuItem("Protein Domain Enrichment", tabName = "prot_enrichment", icon = icon("th")),
                        menuItem("About", tabName = "about", icon = icon("question-circle"))
                )
        ),
        # Body content
        dashboardBody(
                tabItems(
                        # First tab content : Whole data inspection
                        tabItem(tabName = "whole_data_inspection",
                                conditionalPanel('input.start',{
                                        fluidRow(
                                                # Box volcano plot
                                                column(width = 6,
                                                       box2(
                                                               title = "Volcano plot",
                                                               plotlyOutput("volcano")%>% withSpinner(color="#0dc5c1")
                                                       )
                                                ),
                                                column(width = 6,
                                                       box2(
                                                               title = "MA plot",
                                                               plotlyOutput("MA")%>% withSpinner(color="#0dc5c1")
                                                       )
                                                )
                                        )
                                }),
                                fixedRow(
                                        # Box selection of log2FC and p-value
                                        box2(
                                                title = "OPTIONS",
                                                sliderInput(
                                                        inputId = "pvalue", 
                                                        label = "adjusted p-value cutoff from input :",
                                                        min = 0, 
                                                        max = 1,
                                                        value = 1),
                                                sliderInput(
                                                        inputId = "FC", 
                                                        label = "log2 Fold-Change cutoff from input :", 
                                                        min = 0, 
                                                        max = 5, 
                                                        value = 0)
                                        ),
                                ),
                                # add_busy_spinner(spin = "fading-circle", position = 'full-page'),
                                add_busy_gif(src = 'loader.gif', position = 'full-page'),
                                DT::dataTableOutput(outputId = "table")
                        ),
                        
                        # Second tab content : GO Terms Enrichment
                        tabItem(tabName = "GO_enrichment",
                                fluidRow(
                                        column(
                                                width = 6,
                                                box2(
                                                        title = "Analysis method", 
                                                        radioButtons(inputId = "go_analysis_method", 
                                                                     label = NULL, 
                                                                     choiceValues = c("ORA","GSEA"), 
                                                                     choiceNames = c("Over representation analysis (ORA)","Gene set enrichment analysis (GSEA)")
                                                        ),
                                                        hr(),
                                                        radioButtons("go_filter",
                                                                     label = NULL,
                                                                     choiceValue = c("DEG+", "DEG-", "both"),
                                                                     choiceNames = c("Over expressed DEG only", "Under expressed DEG only", "Both")
                                                        )
                                                )
                                        ),
                                        column(
                                                width = 6,
                                                box2(
                                                        title = "Gene ontology settings", 
                                                        selectInput("go_ontology", 
                                                                    label = NULL, 
                                                                    choices = c("Biological process", "Molecular function", "Cellular component")
                                                        ),
                                                        br(),
                                                        radioButtons("go_level",
                                                                     label = NULL,
                                                                     choiceValue = c("all_level", "one_level"),
                                                                     choiceNames = c("All-level GO terms", "One-level GO terms")
                                                        ),
                                                        conditionalPanel("input.go_level == 'one_level'", 
                                                                         sliderInput("go_level_selected",
                                                                                     label = NULL,
                                                                                     min = 1, max = 10, value = 3
                                                                         )
                                                        )
                                                )
                                        )
                                ),
                                conditionalPanel(
                                        "input.go_analysis_method == 'ORA'",
                                        fluidRow(
                                                column(
                                                        width = 6,
                                                        box2(
                                                                title = "Barplot",
                                                                plotOutput("GO_ORA_barplot")%>% withSpinner(color="#0dc5c1")
                                                        )
                                                ),
                                                column(
                                                        width = 6,
                                                        box2(
                                                                title = "GO plot",
                                                                plotOutput("GO_ORA_goplot")%>% withSpinner(color="#0dc5c1")
                                                        )
                                                )
                                                
                                        )
                                ),
                                conditionalPanel(
                                        "input.go_analysis_method == 'GSEA'",
                                        fluidRow(
                                                column(
                                                        width = 6,
                                                        box2(
                                                                title = "Dotplot",
                                                                plotOutput("GO_GSEA_dotplot")%>% withSpinner(color="#0dc5c1")
                                                        )
                                                ),
                                                column(
                                                        width = 6,
                                                        box2(
                                                                title = "GSEA plot",
                                                                plotOutput("GO_GSEA_plot")%>% withSpinner(color="#0dc5c1")
                                                        )
                                                )
                                        )
                                ),
                                fluidRow(
                                        box2(
                                                title = "Plot settings",
                                                sliderInput("go_pvalue",
                                                            min = 0, 
                                                            max = 1, 
                                                            value = 0.05, 
                                                            label = "Select a adjusted p-value cutoff"
                                                )
                                        )
                                ),
                                fluidRow(
                                        box2(
                                                title = "Download",
                                                conditionalPanel(
                                                        "input.go_analysis_method == 'ORA'",
                                                        downloadButton("download_go_barplot","Barplot"),
                                                        downloadButton("download_go_dotplot", "Dotplot")
                                                ),
                                                conditionalPanel(
                                                        "input.go_analysis_method == 'GSEA'",
                                                        downloadButton("download_go_gseaplot", "GSEA plot"),
                                                        downloadButton("download_go_goplot", "GO plot")
                                                )
                                        )
                                ),
                                dataTableOutput("go_datatable")
                        ),
                        
                        # Third tab content : Pathways Enrichment
                        tabItem(tabName = "path_enrichment",
                                fluidRow(
                                        column(
                                                width = 6,
                                                box2(
                                                        title = "Analysis method", 
                                                        radioButtons(inputId = "path_analysis_method", 
                                                                     label = NULL, 
                                                                     choiceValues = c("ORA","GSEA"), 
                                                                     choiceNames = c("Over representation analysis (ORA)","Gene set enrichment analysis (GSEA)")
                                                        ),
                                                        hr(),
                                                        radioButtons("path_filter",
                                                                     label = NULL,
                                                                     choiceValue = c("DEG+", "DEG-", "both"),
                                                                     choiceNames = c("Over expressed DEG only", "Under expressed DEG only", "Both")
                                                        )
                                                )
                                        ),
                                        column(
                                                width = 6,
                                                box2(
                                                        title = "Database", 
                                                        radioButtons("path_database",
                                                                     label = NULL,
                                                                     choiceValue = c("kegg", "reactome"),
                                                                     choiceNames = c("KEGG", "Reactome")
                                                        )
                                                )
                                        )
                                ),
                                conditionalPanel(
                                        "input.path_analysis_method == 'ORA'",
                                        box2(
                                                title = "Barplot",
                                                plotOutput("path_barplot")%>% withSpinner(color="#0dc5c1")
                                        )
                                        
                                ),
                                conditionalPanel(
                                        "input.path_analysis_method == 'GSEA'",
                                        fluidRow(
                                                column(
                                                        width = 6,
                                                        box2(
                                                                title = "GSEA plot",
                                                                plotOutput("path_gseaplot")%>% withSpinner(color="#0dc5c1")
                                                        )
                                                ),
                                                column(
                                                        width = 6,
                                                        box2(
                                                                title = "Dotplot",
                                                                plotOutput("path_dotplot")%>% withSpinner(color="#0dc5c1")
                                                        )
                                                )
                                        ),
                                        box2(
                                                title = "Pathway plot",
                                                plotOutput("path_pathplot")%>% withSpinner(color="#0dc5c1")
                                        )
                                ),
                                box2(
                                        title = "Plot settings",
                                        sliderInput("path_pvalue",
                                                    min = 0, 
                                                    max = 1, 
                                                    value = 0.05, 
                                                    label = "Select a adjusted p-value cutoff"
                                        )
                                ),
                                box2(
                                        title = "Download",
                                        conditionalPanel(
                                                "input.path_analysis_method == 'ORA'",
                                                downloadButton("download_path_barplot","Barplot"),
                                        ),
                                        conditionalPanel(
                                                "input.path_analysis_method == 'GSEA'",
                                                downloadButton("download_path_dotplot", "Dotplot"),
                                                downloadButton("download_path_gseaplot", "GSEA plot"),
                                                downloadButton("download_path_goplot", "Path plot")
                                        )
                                ),
                                dataTableOutput("path_datatable")
                        ),
                        tabItem("prot_enrichment",
                                fluidRow(
                                        column(
                                                width = 6,
                                                box2(
                                                        title = "Analysis method", 
                                                        radioButtons(inputId = "domain_analysis_method", 
                                                                     label = NULL, 
                                                                     choiceValues = c("ORA","GSEA"), 
                                                                     choiceNames = c("Over representation analysis (ORA)","Gene set enrichment analysis (GSEA)")
                                                        ),
                                                        hr(),
                                                        radioButtons("domain_filter",
                                                                     label = NULL,
                                                                     choiceValue = c("DEG+", "DEG-", "both"),
                                                                     choiceNames = c("Over expressed DEG only", "Under expressed DEG only", "Both")
                                                        )
                                                )
                                        ),
                                        column(
                                                width = 6,
                                                box2(
                                                        title = "Plot settings",
                                                        sliderInput("domain_pvalue",
                                                                    min = 0, 
                                                                    max = 1, 
                                                                    value = 0.05, 
                                                                    label = "Select a adjusted p-value cutoff"
                                                        )
                                                )
                                        )
                                ),
                                box2(
                                        title = "Domain barplot plot",
                                        plotlyOutput("domain_barplot")%>% withSpinner(color="#0dc5c1")
                                ),
                                box2(
                                        title = "Domain barplot plot",
                                        plotlyOutput("domain_barplot2")%>% withSpinner(color="#0dc5c1")
                                ),
                                box2(
                                        title = "Domain dotplot plot",
                                        plotlyOutput("domain_dotplot")%>% withSpinner(color="#0dc5c1")
                                ),
                                dataTableOutput("domain_datatable")
                        ),
                        tabItem(tabName = "about",
                                includeMarkdown('README.md')
                        )
                )
        )
)

        
        
        
        
        
        
