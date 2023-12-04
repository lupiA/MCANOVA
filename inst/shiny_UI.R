library(shiny)
library(data.table)
library(ggplot2)


# Define the UI
fluidPage(
  titlePanel("SNP Portability App"),
  sidebarLayout(
    sidebarPanel(
      selectInput("dataset", "Select Genotype Array:",
                  # choices = c("Calls SNPs", "Imputed SNPs")
                  choices = c("Calls SNPs")
      ),
      # h6("*imputed selection will take longer to load and run."),
      br(),
      uiOutput("ancestry_dropdown"),
      br(),
      radioButtons("input_range", "Marker Input:",
                   choices = c("Single Marker", "Range of Markers (within chromosome)", "Comma-separated List of SNP RS IDs","Single Gene")
      ),
      br(),
      conditionalPanel(
        condition = "input.input_range == 'Single Marker'",
        radioButtons("input_type", "Input Type:",
                     choices = c("SNP RS ID", "Base Pair Position & Chromosome")
        )
      ),
      conditionalPanel(
        condition = "input.input_range == 'Range of Markers (within chromosome)'",
        radioButtons("input_type2", "Input Type:",
                     choices = c("SNP RS ID", "Base Pair Position & Chromosome")
        )
      ),
      conditionalPanel(
        condition = "input.input_type == 'SNP RS ID' && input.input_range == 'Single Marker'",
        textInput("rs_id", "Enter SNP ID (e.g., rs4422948):")
      ),
      conditionalPanel(
        condition = "input.input_type == 'Base Pair Position & Chromosome' && input.input_range == 'Single Marker'",
        numericInput("chromosome_input", "Enter Chromosome:", value = 1),
        numericInput("bp_position", "Enter Base Pair Position:", value = 1)
      ),
      conditionalPanel(
        condition = "input.input_type2 == 'SNP RS ID' && input.input_range == 'Range of Markers (within chromosome)'",
        textInput("start_snp_id", "Enter Start SNP ID:"),
        textInput("end_snp_id", "Enter End SNP ID:")
      ),
      conditionalPanel(
        condition = "input.input_type2 == 'Base Pair Position & Chromosome' && input.input_range == 'Range of Markers (within chromosome)'",
        numericInput("chromosome_input2", "Enter Chromosome:", value = 1),
        numericInput("start_bp_position", "Enter Start Base Pair Position:", value = 1),
        numericInput("end_bp_position", "Enter End Base Pair Position:", value = 1)
      ),
      conditionalPanel(
        condition = "input.input_range == 'Comma-separated List of SNP RS IDs'",
        textInput("snp_list", "Enter SNP List (e.g., rs4422948,rs6686302):")
      ),
      conditionalPanel(
        condition = "input.input_range == 'Single Gene'",
        textInput("gene_name", "Enter Gene Name (e.g., ABCG2):")
      )
    ),
    mainPanel(
      fluidRow(width = 10,
               column(10,
                      verbatimTextOutput("error_message")
               )
      ),
      fluidRow(width = 10, 
               h3("Portability Information"),
               tableOutput("table1")),
      fluidRow(
        width = 10, 
        plotOutput("histogram")
      ),
      fluidRow(width = 10, 
               h3("Additional Information"),
               downloadButton("download_table", "Download Table"),
               tableOutput("table2")
      )
    )
  )
)
