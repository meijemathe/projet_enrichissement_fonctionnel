library(shinydashboard)

dashboardPage(
        skin = "black",
        dashboardHeader(title = "Application"),
                # Sidebar content
                dashboardSidebar(
                        sidebarMenu(
                                # File selection
                                fileInput("file", label = "File input"),
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
                                # Box vo
                                box(
                                        title = "Volcano plot",
                                        solidHeader = TRUE,
                                        status = "primary",
                                        plotlyOutput("volcano")
                                ),
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
                                        title = "MA plot",
                                        solidHeader = TRUE,
                                        status = "primary",
                                        plotlyOutput("MA")
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
                                )
                        ),
                        # Second tab content : GO Terms Enrichment
                        tabItem(tabName = "GO_enrichment"
                        )
                )
        )
)