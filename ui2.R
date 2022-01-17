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

#####################################################################


box2 <- function(...){
        box(
                status = "primary",
                solidHeader = TRUE,
                width = 12,
                ...
        )
}




dashboardPage(
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
                                                box(
                                                        title = "Volcano plot",
                                                        solidHeader = TRUE,
                                                        status = "primary",
                                                        plotlyOutput("volcano")%>% withSpinner(color="#0dc5c1")
                                                ),
                                                box(
                                                        title = "MA plot",
                                                        solidHeader = TRUE,
                                                        status = "primary",
                                                        plotlyOutput("MA")%>% withSpinner(color="#0dc5c1")
                                                )
                                        )
                                }),
                                fixedRow(
                                        # Box selection of log2FC and p-value
                                        box(
                                                title = "OPTIONS",
                                                solidHeader = TRUE,
                                                status = "primary",
                                                sliderInput(
                                                        inputId = "pvalue", 
                                                        label = "fitted p-value cutoff from input :",
                                                        min = 0, 
                                                        max = 1,
                                                        value = 0.05),
                                                sliderInput(
                                                        inputId = "FC", 
                                                        label = "log2 Fold-Change cutoff from input :", 
                                                        min = 0, 
                                                        max = 5, 
                                                        value = 1)
                                        ),
                                        box(
                                                title = "SET UP FOR SUBSEQUENT ANALYSES",
                                                solidHeader = TRUE,
                                                status = "primary",
                                                sliderInput(
                                                        inputId = "FC", 
                                                        label = "fitted p-value cutoff for subsequent analyses :", 
                                                        min = 0, 
                                                        max = 1, 
                                                        value = 0.05)
                                        ),
                                ),
                                add_busy_spinner(spin = "fading-circle", position = 'full-page'),
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
                                                                plotlyOutput("go_barplot")%>% withSpinner(color="#0dc5c1")
                                                        )
                                                ),
                                                column(
                                                        width = 6,
                                                        box2(
                                                                title = "GO plot",
                                                                plotlyOutput("go_gseaplot")%>% withSpinner(color="#0dc5c1")
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
                                                                title = "GSEA plot",
                                                                plotlyOutput("GO_gseaplot")%>% withSpinner(color="#0dc5c1")
                                                        )
                                                ),
                                                column(
                                                        width = 6,
                                                        box2(
                                                                title = "Dotplot",
                                                                plotlyOutput("go_dotplot")%>% withSpinner(color="#0dc5c1")
                                                        )
                                                )
                                        )
                                ),
                                fluidRow(
                                        box(
                                                title = "Plot settings",
                                                solidHeader = TRUE,
                                                status = "primary",
                                                width = 12,
                                                sliderInput("go_pvalue",
                                                            min = 0, 
                                                            max = 1, 
                                                            value = 0.05, 
                                                            label = "Select a adjusted p-value cutoff"
                                                )
                                        )
                                ),
                                fluidRow(
                                        box(
                                                title = "Download",
                                                solidHeader = TRUE,
                                                status = "primary",
                                                width = 12,
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
                                )
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
                                                plotlyOutput("path_barplot")%>% withSpinner(color="#0dc5c1")
                                        )
                                        
                                ),
                                conditionalPanel(
                                        "input.path_analysis_method == 'GSEA'",
                                        fluidRow(
                                                column(
                                                        width = 6,
                                                        box2(
                                                                title = "GSEA plot",
                                                                plotlyOutput("path_gseaplot")%>% withSpinner(color="#0dc5c1")
                                                        )
                                                ),
                                                column(
                                                        width = 6,
                                                        box2(
                                                                title = "Dotplot",
                                                                plotlyOutput("path_dotplot")%>% withSpinner(color="#0dc5c1")
                                                        )
                                                )
                                        ),
                                        box2(
                                                title = "Pathway plot",
                                                plotlyOutput("path_pathplot")%>% withSpinner(color="#0dc5c1")
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
                                )
                        )
                )
        )
)
        
        
        
        
        
        
        