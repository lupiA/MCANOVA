#' Run the PGS portability shiny app
#' 
#' @description Calls the R shiny app
#'
#' @export
PGS_portability_app <- function() {
  library(shiny)
  library(data.table)
  library(ggplot2)


  # Define the UI
  ui <- fluidPage(
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
                     choices = c("Single Marker", "Range of Markers (within chromosome)",
                     "Comma-separated List of SNP RS IDs","Single Gene")
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
          condition = "input.input_type == 'Base Pair Position & Chromosome' && 
              input.input_range == 'Single Marker'",
          numericInput("chromosome_input", "Enter Chromosome:", value = 1),
          numericInput("bp_position", "Enter Base Pair Position:", value = 1)
        ),
        conditionalPanel(
          condition = "input.input_type2 == 'SNP RS ID' && input.input_range == 
              'Range of Markers (within chromosome)'",
          textInput("start_snp_id", "Enter Start SNP ID:"),
          textInput("end_snp_id", "Enter End SNP ID:")
        ),
        conditionalPanel(
          condition = "input.input_type2 == 'Base Pair Position & Chromosome' && 
              input.input_range == 'Range of Markers (within chromosome)'",
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
  
  # Define the server
  server <- function(input, output, session) {
    ancestry_choices <- reactive({
      if (input$dataset == "Calls SNPs") {
        c("African","Caribbean","East Asian","South Asian")
      }
    })
    
    observe({
      output$ancestry_dropdown <- renderUI({
        selectInput("ancestry", "Select Ancestry:", choices = ancestry_choices())
      })
    })
    
    dat <- reactive({
      path <- system.file("dat", package = "MCANOVA")
      if (input$dataset == "Calls SNPs") {
        if (input$ancestry == "African"){
          fread(paste0(path, "AF_map.csv"), data.table = FALSE)
        } else if (input$ancestry == "Caribbean") {
          fread(paste0(path, "CR_map.csv"), data.table = FALSE)
        } else if (input$ancestry == "East Asian") {
          fread(paste0(path, "EA_map.csv"), data.table = FALSE)
        } else if (input$ancestry == "South Asian") {
          fread(paste0(path, "SA_map.csv"), data.table = FALSE)
        }
      }
    }) 
    
    ancestry_label <- reactive({
      if(input$ancestry=="African"){
        "AF"
      } else if(input$ancestry=="Caribbean"){
        "CR"
      } else if(input$ancestry=="East Asian"){
        "EA"
      } else if(input$ancestry=="South Asian"){
        "SA"
      }
    })
    
    filtered_snp_data <- reactive({
      snps <- dat()
      colnames(snps)[3] <- "Relative Accuracy (RA)"
      snp_data <- NULL
      
      if (input$input_range == "Single Marker") {
        
        if (input$input_type == "SNP RS ID") {
          snp_data <- snps[which(snps$SNP == input$rs_id & snps$Cohort == input$ancestry),
                           c("SNP","Relative Accuracy (RA)", paste0("Corr. EU=>", ancestry_label()), "Corr. EU=>EU", "Chromosome", "BP position", "Allele", paste0("Corr. EU=>", ancestry_label(), " S.E."), "Corr. EU=>EU S.E.","Gene")]
        } 
        
        else if (input$input_type == "Base Pair Position & Chromosome") {
          snps_chr <- snps[which(snps$Chromosome == input$chromosome_input),]
          snp_data <- snps_chr[which(snps_chr$`BP position` == input$bp_position & snps_chr$Cohort == input$ancestry),
                               c("SNP","Relative Accuracy (RA)", paste0("Corr. EU=>", ancestry_label()), "Corr. EU=>EU", "Chromosome", "BP position", "Allele", paste0("Corr. EU=>", ancestry_label(), " S.E."), "Corr. EU=>EU S.E.","Gene")]
        }
        
      }  else if (input$input_range == "Range of Markers (within chromosome)") {
        
        if (input$input_type2 == "SNP RS ID") {
          snps <- snps[which(snps$Cohort == input$ancestry), ]
          start_snp_idx_1 <- which(snps$SNP == input$start_snp_id)
          end_snp_idx_1 <- which(snps$SNP == input$end_snp_id)
          
          if (length(start_snp_idx_1) > 0 && length(end_snp_idx_1) > 0 &&
              snps$Chromosome[start_snp_idx_1] == snps$Chromosome[end_snp_idx_1]) {
            snp_data <- snps[start_snp_idx_1:end_snp_idx_1,
                             c("SNP","Relative Accuracy (RA)",paste0("Corr. EU=>", ancestry_label()), "Corr. EU=>EU", "Chromosome", "BP position", "Allele", paste0("Corr. EU=>", ancestry_label(), " S.E."), "Corr. EU=>EU S.E.","Gene")]
          }
          
        }  else if (input$input_type2 == "Base Pair Position & Chromosome") {
          snps <- snps[which(snps$Cohort == input$ancestry), ]
          snps <- snps[which(snps$Chromosome == input$chromosome_input2), ]
          
          if(is.null(input$start_bp_position) || is.null(input$end_bp_position) || is.na(input$start_bp_position) || is.na(input$end_bp_position)){
            snp_data <- NULL
          } else if(length(input$start_bp_position:input$end_bp_position) > 0 && input$start_bp_position < input$end_bp_position && 
                                      input$start_bp_position > 0 && nrow(snps)>0){
            snp_data <- snps[which(snps$`BP position` %in% input$start_bp_position:input$end_bp_position),
                             c("SNP","Relative Accuracy (RA)", paste0("Corr. EU=>", ancestry_label()), "Corr. EU=>EU", "Chromosome", "BP position", "Allele", paste0("Corr. EU=>", ancestry_label(), " S.E."), "Corr. EU=>EU S.E.","Gene")]
          }
        }
      }  else if (input$input_range == "Comma-separated List of SNP RS IDs") {
        snps <- snps[which(snps$Cohort == input$ancestry), ]
        snps_in_list <- unlist(strsplit(as.character(input$snp_list), split = ","))
        
        if (length(snps_in_list) > 0) {
          snp_data <- snps[which(snps$SNP %in% snps_in_list),
                           c("SNP","Relative Accuracy (RA)", paste0("Corr. EU=>", ancestry_label()), "Corr. EU=>EU", "Chromosome", "BP position", "Allele", paste0("Corr. EU=>", ancestry_label(), " S.E."), "Corr. EU=>EU S.E.","Gene")]
        }
        
      }  else if (input$input_range == "Single Gene") {
        snps <- snps[which(snps$Cohort == input$ancestry), ]
        snps_in_list <- snps$SNP[which(snps$Gene == input$gene_name)]
        
        if (length(snps_in_list) > 0) {
          snp_data <- snps[which(snps$SNP %in% snps_in_list),
                           c("SNP","Relative Accuracy (RA)", paste0("Corr. EU=>", ancestry_label()), "Corr. EU=>EU", "Chromosome", "BP position", "Allele", paste0("Corr. EU=>", ancestry_label(), " S.E."), "Corr. EU=>EU S.E.","Gene")]
        }
      }
      
      if (is.null(snp_data) || nrow(snp_data) == 0) {
        snp_data <- NULL
      }
      snp_data
    })
    
    filtered_hist_data <- reactive({
      hist <- dat()
      hist_data <- hist[which(hist$Cohort == input$ancestry), c("RA","Cohort")]
      if (nrow(hist_data) == 0) {
        hist_data <- NULL
      }
      if (!is.null(hist_data) && length(which(is.na(hist_data$RA))) > 0) {
        hist_data <- hist_data[-which(is.na(hist_data$RA)),]
      }
      hist_data
    })
    
    output$table1 <- renderTable({
      if (!is.null(input$rs_id) && !is.null(input$ancestry) && !is.null(filtered_snp_data())) {
        if (!is.null(filtered_snp_data())) {
          averages <- colMeans(filtered_snp_data()[, c("Relative Accuracy (RA)", paste0("Corr. EU=>", ancestry_label()), "Corr. EU=>EU")], na.rm = TRUE)
          if (input$input_range == "Range of Markers (within chromosome)") {
            data.frame("Metric" = names(averages), "Average" = averages)
          } else if (input$input_range == "Single Marker") {
            data.frame("Metric" = names(averages), "Value" = averages)
          } else if (input$input_range == "Comma-separated List of SNP RS IDs") {
            data.frame("Metric" = names(averages), "Average" = averages)
          } else if (input$input_range == "Single Gene") {
            data.frame("Metric" = names(averages), "Average" = averages)
          }
        }
      }
    })
    
    output$table2 <- renderTable({
      if (!is.null(input$rs_id) && !is.null(input$ancestry) && !is.null(filtered_snp_data())) {
        if (!is.null(filtered_snp_data())) {
          modified_table <- filtered_snp_data()[, c("SNP", "Chromosome", "BP position", 
                                                    "Allele", "Relative Accuracy (RA)", paste0("Corr. EU=>", ancestry_label()), paste0("Corr. EU=>", ancestry_label(), " S.E."), "Corr. EU=>EU", "Corr. EU=>EU S.E.","Gene")]
          data.frame(modified_table,check.names = FALSE)
        }
      }
    })
    
    output$histogram <- renderPlot({
      if (input$input_range == "Single Marker" && !is.null(input$rs_id) && !is.null(input$ancestry) && !is.null(filtered_hist_data())) {
        ra_values <- as.numeric(filtered_hist_data()$RA)
        snp_ra <- as.numeric(filtered_snp_data()$`Relative Accuracy (RA)`)
        
        snp_percentile <- ecdf(ra_values)(snp_ra) * 100
        
        hist(ra_values, main = "Relative Accuracy Histogram (input in red)\n\ ", xlab = "RA Value", ylab = "Frequency", col = "lightblue1", breaks = 75)
        abline(v = snp_ra, col = "red3", lwd = 3)
        
        mtext(paste("SNP RA Percentile: ", round(snp_percentile, 2), "%\n\ ",sep=""), side = 3, col = "black")
      } else if (input$input_range == "Range of Markers (within chromosome)" || input$input_range == "Comma-separated List of SNP RS IDs" || 
                 input$input_range == "Single Gene" && !is.null(input$rs_id) && !is.null(input$ancestry) && !is.null(filtered_hist_data())) {
        ra_values <- as.numeric(filtered_hist_data()$RA)
        snp_ra <- as.numeric(filtered_snp_data()$`Relative Accuracy (RA)`)
        
        snp_percentile <- ecdf(ra_values)(snp_ra) * 100
        
        hist(ra_values, main = "Relative Accuracy Histogram (inputs in red)\n\ ", xlab = "RA Value", ylab = "Frequency", col = "lightblue1", breaks = 75)
        if (length(table(snp_ra)) > 10) {
          if(input$dataset=="Calls SNPs"){
            points(x = snp_ra, y = rep(12000, length(snp_ra)), col = "red3", lwd = .6, pch = 4)
          } else if(input$dataset=="Imputed SNPs"){
            points(x = snp_ra, y = rep(150000, length(snp_ra)), col = "red3", lwd = .6, pch = 4)
          }
        } else {
          abline(v = snp_ra, col = "red3", lwd = 3)
        }
      }
    })
    
    output$download_table <- downloadHandler(
      filename = function() {
        paste("table_", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        if (!is.null(input$rs_id) && !is.null(input$ancestry) && !is.null(filtered_snp_data())) {
          filtered_table <- filtered_snp_data()[, c("SNP", "Chromosome", "BP position", "Allele", "Relative Accuracy (RA)", paste0("Corr. EU=>", ancestry_label()), paste0("Corr. EU=>", ancestry_label(), " S.E."), "Corr. EU=>EU", "Corr. EU=>EU S.E.","Gene")]
          data.frame(filtered_table,check.names = FALSE)
          write.csv(filtered_table, file, row.names = FALSE)
        }
      }
    )
    
    output$error_message <- renderText({
      if (is.null(filtered_snp_data())) {
        "Enter valid inputs. Range of entries must be within-chromosome,\n\ and a gene name, SNP RS IDs, or single BP position entry must be\n in the UK Biobank array."
      } else {
        return(NULL)
      }
    })
  }
  
  # Run
  shinyApp(ui = ui, server = server)
}
