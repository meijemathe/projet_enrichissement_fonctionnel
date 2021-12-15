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

dashboardPage(
        skin = "black",
        dashboardHeader(title = "EnF'R"),
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
                        useShinyalert(),
                        # Start Button
                        actionButton("start", label = "Start analysis"),
                        hr(),
                        # Menu 
                        menuItem("Whole Data Inspection", tabName = "whole_data_inspection", icon = icon("database")),
                        menuItem("GO Term Enrichment", tabName = "GO_enrichment", icon = icon("sitemap")),
                        menuItem("Pathway Enrichment", tabName = "path_enrichment", icon = icon("chart-pie")),
                        menuItem("Protein Domain Enrichment", tabName = "prot_enrichment", icon = icon("th")),
                        menuItem("About", tabName = "about", icon = icon("th"))
                )
        ),
        # Body content
        dashboardBody(
                tabItems(
                        # First tab content : Whole data inspection
                        tabItem(tabName = "whole_data_inspection",
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
                                ),
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
                                add_busy_spinner(spin = "fading-circle"),
                                DT::dataTableOutput(outputId = "table")
                        ),
                        # Second tab content : GO Terms Enrichment
                        tabItem(tabName = "GO_enrichment"
                        )
                )
        )
)