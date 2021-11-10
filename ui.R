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

library(shiny)
library(shinydashboard)
library(DT)
library(plotly)
library(htmlwidgets)
#####################################################################

dashboardPage(
        skin = "black",
        dashboardHeader(title = "EnF'R"),
        # Sidebar content
        dashboardSidebar(
                sidebarMenu(
                        # File selection
                        fileInput("file", label = "File input", accept = c(
                                "text/csv",
                                "text/comma-separated-values,text/plain",
                                ".csv")),
                        # verbatimTextOutput("file"),
                        # Select box : organism selection
                        selectInput("select_organism", label = "Organism of interest",
                                choices = list("Homo Sapiens" = 1, "Mus musculus" = 2, "Arabidopsis Thaliana" = 3),
                                selected = 1
                        ),
                        # verbatimTextOutput("select_organism"),
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
                                                plotlyOutput("volcano")
                                        ),
                                        box(
                                                title = "MA plot",
                                                solidHeader = TRUE,
                                                status = "primary",
                                                plotlyOutput("MA")
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
                                                        label = "p-value cutoff from input :",
                                                        min = 0, 
                                                        max = 1,
                                                        value = 0.05),
                                                sliderInput(
                                                        inputId = "FC", 
                                                        label = "log2 Fold-Change cutoff from input :", 
                                                        min = 0, 
                                                        max = 5, 
                                                        value = 1),
                                                downloadButton(
                                                        outputId = "download_volcano", 
                                                        label = "Download volcano plot", 
                                                ),
                                                downloadButton(
                                                        outputId = "download_MA", 
                                                        label = "Download MA plot", 
                                                )
                                        ),
                                        box(
                                                title = "SET UP FOR SUBSEQUENT ANALYSES",
                                                solidHeader = TRUE,
                                                status = "primary",
                                                sliderInput(
                                                        inputId = "FC", 
                                                        label = "p-value cutoff for subsequent analyses :", 
                                                        min = 0, 
                                                        max = 1, 
                                                        value = 0.05)
                                        ),
                                ),
                                DT::dataTableOutput(outputId = "table")
                        ),
                        # Second tab content : GO Terms Enrichment
                        tabItem(tabName = "GO_enrichment"
                        )
                )
        )
)