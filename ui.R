library(shiny)

shinyUI(fluidPage(
                  
                  fluidRow(
                            h3("Pathway and Gene Search"),
                            h6("Welcome! Please select a SUM cell line and a KEGG pathway to see a pathway diagram
                               enriched by a Cellecta RNAi screen."),
                            uiOutput("sumlineSelect"),
                    h6("To search for a particular pathway, press the backspace
                       key and type either the pathway name or KEGG pathway ID."),
                    selectizeInput("pathway", label = "Pathway:", choices = read.csv("./www/kegg.txt", header = FALSE)[1]),
                    radioButtons("enrich", label = "Enrichment Dataset:", choices = c("shRNA Screen", "Gene Expression")),
                    h6("For information on a gene in the selected cell line, use the search bar below."),
                    textInput("gene", label = "Gene:"),
                    fluidRow(column(4, dataTableOutput("geneInfo")))
                  ),
                  fluidRow(
                    column(8, dataTableOutput("pathwayGenes"))
                  ),
                  fluidRow(
                    imageOutput("pathway")
                  )
                    )
  
)