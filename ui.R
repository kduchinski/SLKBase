library(shiny)

shinyUI(fluidPage(
  # navbarPage("SLKBase",
             # tabPanel("Welcome",
             #          titlePanel("The SUM Breast Cancer Cell Line Knowledge Base"),
             #          sidebarLayout(
             #            sidebarPanel(
             #              img(src = 'cropped-sextant.png'),
             #              h3("Sextant Oncology"),
             #              h4("The SUM breast cancer cell line
             #                 knowledge base"),
             #              h6("Welcome to the first edition of the SUM breast cancer cell line Knowledge Base. 
             #                 Several years ago, I and my laboratory developed a series of breast cancer cell 
             #                 lines from biopsy specimens obtained from patients being treated at the University 
             #                 of Michigan Comprehensive Cancer Center.")
             #            ),
             #            mainPanel(
             #              br(),
             #              br(),
             #              tags$blockquote(h3("This is your one-stop resource for knowledge, data,
             #                 and information relevant to the SUM series of human
             #                 breast cancer cell lines, developed in the Ethier
             #                 laboratory while at the University of Michigan.")),
             #              br(),
             #              h4("Here is what you'll find in these pages:"),
             #              h5("A brief description of how each cell line was derived,
             #                 including information on the patient from which the cells
             #                 were isolated"),
             #              h5("A bibliography that contains every paper ever published that
             #                 contains data obtained with that particular cell line"),
             #              h5("Links to a number of KEGG pathways where all of the genomic
             #                 data obtained with each cell line is displayed in a pathway-specific
             #                 manner. This is also your gateway to viewing all of the data on
             #                 copy number, mutation status, expression level, and hit score for
             #                 each annotated gene in the pathway. There will also be a narrated
             #                 summary of the results for each pathway.")
             #            )
             #          )),
             # tabPanel("Pathway Visualization",
                      fluidRow(
                            h3("Pathway and Gene Search"),
                            h6("Welcome! Please select a SUM cell line and a KEGG pathway to see a pathway diagram
                               enriched by a Cellecta RNAi screen."),
                            uiOutput("sumlineSelect"),
                            h6("To search for a particular pathway, press the backspace
                               key and type either the pathway name or KEGG pathway ID."),
                            selectizeInput("pathway", label = "Pathway:", choices = read.csv("./www/kegg.txt", header = FALSE)[1]),
                            h6("For information on a gene in the selected cell line, use the search bar below."),
                            textInput("gene", label = "Gene:"),
                            fluidRow(column(4, dataTableOutput("geneInfo")))
                        ),
                        fluidRow(
                          imageOutput("pathway")
                        )
                       )
             # navbarMenu("Cell lines",
             #            tabPanel("SUM149"),
             #            tabPanel("SUM185"),
             #            tabPanel("SUM190"),
             #            tabPanel("SUM225"),
             #            tabPanel("SUM229"),
             #            tabPanel("SUM44"),
             #            tabPanel("SUM52")
             #            )
             # )
  
)