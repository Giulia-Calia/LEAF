#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(DT)
library(readr)
library(aweSOM)
library(bslib)


# load dataset 
dset <- read_tsv("/home/giulia/Workspace/PhytoPhD/period_abroad/functional_analysis/unify_scaled_comp_choice_ALL_dset_func_analysis.tsv")
dset_to_subset <- read_tsv("/home/giulia/Workspace/PhytoPhD/period_abroad/functional_analysis/predicted_eff_to_subset_shiny.tsv")
features_colnames <- c("sequence length","signal peptide","transmembrane domain","aa in tr domain","first 60 aa",
                       "prob N-in","MobiDB-lite", 
                       "CAMP_PHOSPHO_SITE", "PKC_PHOSPHO_SITE", "CK2_PHOSPHO_SITE", "MYRISTYL", "AMIDATION", "ASN_GLYCOSYLATION","L=0", 
                       "LEUCINE_ZIPPER","MICROBODIES_CTER", "RGD", "TYR_PHOSPHO_SITE_1", "TYR_PHOSPHO_SITE_2",
                       "CLUMP0","CLUMP3","CLUMP5","CLUMP6","CLUMP7","CLUMP8", "bin1", "bin2", "bin3", "bin4")

features <- as.matrix(dset[, features_colnames])
# Define UI for application that draws a histogram
ui <- fluidPage(
  # Define app theme
  theme = bslib::bs_theme(bootswatch = "flatly"),  
  # App title
  titlePanel("Self-Organizing Map"),
  
  # Show the dataset 
  # DT::dataTableOutput("dataset"),
  
  # Create tabset panles for different king of SOMaps
  tabsetPanel(
    tabPanel("SOM-Cloudplot",
            # give possible som label-colors
             sidebarLayout(
               sidebarPanel(selectInput(
               inputId = "clouds_color",
               label = "Choose how to color the Cloudplot",
               choices = c("origin",	"ncs",	"PCA_group"),
               selected = "ncs"
               ),
               downloadButton("download")),
               # Show SOM
               mainPanel(htmlOutput("cloud_som"))
               )),
    tabPanel("SOM-Barplot",
             # give possible som label-bars
             sidebarLayout(
               sidebarPanel(selectInput(
                 inputId = "bars",
                 label = "Choose which features to plot",
                 choices = features_colnames,
                 selected = c("sequence length","signal peptide","transmembrane domain","aa in tr domain","first 60 aa",
                              "prob N-in","MobiDB-lite"),
                 multiple = TRUE,
                 selectize = TRUE
               ),
               textInput("text", label = "Insert output file name", value = "Enter text..."),
               downloadButton("download1")),
               # Show SOM
               mainPanel(htmlOutput("bar_som"))
               )),
    tabPanel("SOM-Pieplot",
             sidebarLayout(
               sidebarPanel(selectInput(
                 inputId = "pie_label",
                 label = "Choose which features to plot",
                 choices = c("spec_id", "16S_group", "model", "vote"),
                 selected = c("spec_id")
               ),
               downloadButton("download2")),
               # Show SOM
               mainPanel(htmlOutput("pie_som"))
             )),
    tabPanel("Download subDataset",             
             sidebarLayout(
                sidebarPanel(selectInput(
                  inputId = "cells",
                  label = "Choose which SOM_cell to download",
                  choices = c(1:64),
                  selected = 1,
                  selectize = TRUE
                ),
                downloadButton("download3")),
                mainPanel(DTOutput("dataset"))
       )
  ))
)

  
# Define server logic required to draw a histogram
server <- function(input, output) {
  
  
  ### RNG Seed (for reproducibility)
  set.seed(42)
  ### Initialization (PCA grid)
  init <- somInit(features, 8, 8)
  data.som <- kohonen::som(
    features, 
    grid = kohonen::somgrid(8, 8, "hexagonal"),
    rlen = 250,
    alpha = c(0.05, 0.01),
    radius = c(2.65, -2.65),
    dist.fcts = "euclidean",
    init = init
  )
  
  ### COULDPLOT
  cloud_som_plot <- reactive({aweSOMplot(
    som = data.som, 
    type = "Cloud",
    data = dset,
    variables = c(input$clouds_color))})
  
  output$cloud_som <- renderUI({cloud_som_plot()})
  output$download <- downloadHandler(
    filename = function() {paste("LEAF_SOM_clouds_", input$clouds_color, '.html', sep='')},
    content = function(file) {
      saveWidget(cloud_som_plot(), file)
    })
  ### BARPLOT
  bar_som_plot <- reactive({aweSOMplot(
    som = data.som, 
    type = "Barplot", 
    data = features,
    variables = c(input$bars),
    values = "prototype")})
  
  output$bar_som <- renderUI({bar_som_plot()})
  output$download1 <- downloadHandler(
    filename = function() {paste(input$text, '.html', sep='')},
    content = function(file) {
      saveWidget(bar_som_plot(), file)
    })
  ### PIEPLOT
  pie_som_plot <- reactive({aweSOMplot(
    som = data.som, 
    type = "Pie", 
    data = dset, 
    variables = input$pie_label)})
  
  output$pie_som <- renderUI({pie_som_plot()})
  output$download2 <- downloadHandler(
    filename = function() {paste("LEAF_SOM_pie_", input$pie_label, '.html', sep='')},
    content = function(file) {
      saveWidget(pie_som_plot(), file)
    })
  ### Add column to dset with the cell assignment for each protein
  dset_to_subset["cell"] <- data.som[["unit.classif"]]
  cell_to_download = reactive({input$cells})
  output$dataset = DT::renderDataTable({
    subset(dset_to_subset, cell == cell_to_download())
  })
  output$download3 <- downloadHandler(
    filename = function() {paste("LEAF_SOM_cell", input$cells, "seq.tsv", sep="\t")},
    content = function(file) {
      write_tsv(subset(dset_to_subset, cell == cell_to_download()), file)
    }
  )

}

# Run the application 
shinyApp(ui = ui, server = server)
