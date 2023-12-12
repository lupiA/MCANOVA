#' Run the PGS portability shiny app
#' 
#' @description Calls the R shiny app
#'
#' @export
PGS_portability_app <- function() {
  
  library(shiny)
  library(ggplot2)
  
  # Define the UI
  ui <- fluidPage(
    titlePanel("SNP Portability App"),
    sidebarLayout(
      sidebarPanel(
        selectInput("ancestry", "Target Ancestry:",
                    choices = c("African", "Caribbean", "East Asian", "South Asian")
        ),
        br(),
        radioButtons("input_range", "Marker Input:",
                     choices = c("Range of Markers (within chromosome)",
                                 "Comma-separated List of SNP RS IDs", "Single Gene", "Single Marker")
        ),
        br(),
        conditionalPanel(
          condition = "input.input_range == 'Single Marker'",
          radioButtons("input_type", "Input Type:",
                       choices = c("Base Pair Position & Chromosome", "SNP RS ID")
          )
        ),
        conditionalPanel(
          condition = "input.input_range == 'Range of Markers (within chromosome)'",
          radioButtons("input_type2", "Input Type:",
                       choices = c("Base Pair Position & Chromosome", "SNP RS ID")
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
          numericInput("chromosome_input2", "Enter Chromosome:", value = 6),
          numericInput("start_bp_position", "Enter Start Base Pair Position (e.g., the MHC region):", value = 28477797),
          numericInput("end_bp_position", "Enter End Base Pair Position:", value = 33448354)
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
                 h3("Table 1: Portability Information"),
                 tableOutput("table1")),
        fluidRow(
          width = 10, 
          plotOutput("histogram")
        ),
        fluidRow(width = 10, 
                 h3("Table 2: Additional Information"),
                 downloadButton("download_table", "Download Table"),
                 tableOutput("table2")
        )
      )
    )
  )
  
  # Define the server
  server <- function(input, output, session) {
    
    dat <- reactive({
      dat <- MCANOVA::MAP
    })
    
    ancestry_label <- reactive({
      if(input$ancestry == "African"){
        "AF"
      } else if(input$ancestry == "Caribbean"){
        "CR"
      } else if(input$ancestry == "East Asian"){
        "EA"
      } else if(input$ancestry == "South Asian"){
        "SA"
      }
    })
    
    filtered_snp_data <- reactive({
      snps <- dat()
      colnames(snps)[6:10] <- c("Relative Accuracy", paste0("Corr. EU=>", ancestry_label()), "Corr. EU=>EU", paste0("Corr. EU=>", ancestry_label(), " S.E."), "Corr. EU=>EU S.E.")
      snp_data <- NULL
      
      if (input$input_range == "Single Marker") {
        
        if (input$input_type == "SNP RS ID") {
          snp_data <- snps[which(snps$SNP == input$rs_id & snps$Cohort == input$ancestry),
                           c("SNP","Relative Accuracy", paste0("Corr. EU=>", ancestry_label()), "Corr. EU=>EU", "Chromosome", "BP position", "Allele", paste0("Corr. EU=>", ancestry_label(), " S.E."), "Corr. EU=>EU S.E.","Gene")]
        } 
        
        else if (input$input_type == "Base Pair Position & Chromosome") {
          snps_chr <- snps[which(snps$Chromosome == input$chromosome_input),]
          snp_data <- snps_chr[which(snps_chr$`BP position` == input$bp_position & snps_chr$Cohort == input$ancestry),
                            c("SNP","Relative Accuracy", paste0("Corr. EU=>", ancestry_label()), "Corr. EU=>EU", "Chromosome", "BP position", "Allele", paste0("Corr. EU=>", ancestry_label(), " S.E."), "Corr. EU=>EU S.E.","Gene")]
        }
        
      }  else if (input$input_range == "Range of Markers (within chromosome)") {
        
        if (input$input_type2 == "SNP RS ID") {
          snps <- snps[which(snps$Cohort == input$ancestry), ]
          start_snp_idx_1 <- which(snps$SNP == input$start_snp_id)
          end_snp_idx_1 <- which(snps$SNP == input$end_snp_id)
          
          if (length(start_snp_idx_1) > 0 && length(end_snp_idx_1) > 0 &&
              snps$Chromosome[start_snp_idx_1] == snps$Chromosome[end_snp_idx_1]) {
            snp_data <- snps[start_snp_idx_1:end_snp_idx_1,
                             c("SNP","Relative Accuracy",paste0("Corr. EU=>", ancestry_label()), "Corr. EU=>EU", "Chromosome", "BP position", "Allele", paste0("Corr. EU=>", ancestry_label(), " S.E."), "Corr. EU=>EU S.E.","Gene")]
            
          }
          
        }  else if (input$input_type2 == "Base Pair Position & Chromosome") {
          snps <- snps[which(snps$Cohort == input$ancestry), ]
          snps <- snps[which(snps$Chromosome == input$chromosome_input2), ]
          
          if(is.null(input$start_bp_position) || is.null(input$end_bp_position) || is.na(input$start_bp_position) || is.na(input$end_bp_position)){
            snp_data <- NULL
          } else if(length(input$start_bp_position:input$end_bp_position) > 0 && input$start_bp_position < input$end_bp_position && 
                    input$start_bp_position > 0 && nrow(snps)>0){
            snp_data <- snps[which(snps$`BP position` %in% input$start_bp_position:input$end_bp_position),
                             c("SNP","Relative Accuracy", paste0("Corr. EU=>", ancestry_label()), "Corr. EU=>EU", "Chromosome", "BP position", "Allele", paste0("Corr. EU=>", ancestry_label(), " S.E."), "Corr. EU=>EU S.E.","Gene")]
          }
        }
      }  else if (input$input_range == "Comma-separated List of SNP RS IDs") {
        snps <- snps[which(snps$Cohort == input$ancestry), ]
        snps_in_list <- unlist(strsplit(as.character(input$snp_list), split = ","))
        
        if (length(snps_in_list) > 0) {
          snp_data <- snps[which(snps$SNP %in% snps_in_list),
                           c("SNP","Relative Accuracy", paste0("Corr. EU=>", ancestry_label()), "Corr. EU=>EU", "Chromosome", "BP position", "Allele", paste0("Corr. EU=>", ancestry_label(), " S.E."), "Corr. EU=>EU S.E.","Gene")]
        }
        
      }  else if (input$input_range == "Single Gene") {
        snps <- snps[which(snps$Cohort == input$ancestry), ]
        snps_in_list <- snps$SNP[which(snps$Gene == input$gene_name)]
        
        if (length(snps_in_list) > 0) {
          snp_data <- snps[which(snps$SNP %in% snps_in_list),
                           c("SNP","Relative Accuracy", paste0("Corr. EU=>", ancestry_label()), "Corr. EU=>EU", "Chromosome", "BP position", "Allele", paste0("Corr. EU=>", ancestry_label(), " S.E."), "Corr. EU=>EU S.E.","Gene")]
        }
      }
      
      if (is.null(snp_data) || nrow(snp_data) == 0) {
        snp_data <- NULL
      }
      snp_data
    })
    
    filtered_hist_data <- reactive({
      hist <- dat()
      colnames(hist)[6:10] <- c("Relative Accuracy", paste0("Corr. EU=>", ancestry_label()), "Corr. EU=>EU", paste0("Corr. EU=>", ancestry_label(), " S.E."), "Corr. EU=>EU S.E.")
      hist_data <- hist[which(hist$Cohort == input$ancestry), c("Relative Accuracy","Cohort")]
      if (nrow(hist_data) == 0) {
        hist_data <- NULL
      }
      if (!is.null(hist_data) && length(which(is.na(hist_data$`Relative Accuracy`))) > 0) {
        hist_data <- hist_data[-which(is.na(hist_data$`Relative Accuracy`)),]
      }
      hist_data
    })
    
    output$table1 <- renderTable({
      if (!is.null(input$rs_id) && !is.null(input$ancestry) && !is.null(filtered_snp_data())) {
        if (!is.null(filtered_snp_data())) {
          averages <- colMeans(filtered_snp_data()[, c("Relative Accuracy", paste0("Corr. EU=>", ancestry_label()), "Corr. EU=>EU")], na.rm = TRUE)
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
                                                    "Allele", "Relative Accuracy", paste0("Corr. EU=>", ancestry_label()), paste0("Corr. EU=>", ancestry_label(), " S.E."), "Corr. EU=>EU", "Corr. EU=>EU S.E.","Gene")]
          data.frame(modified_table,check.names = FALSE)
        }
      }
    })
    
    output$histogram <- renderPlot({
      theme.ggplot <- theme(legend.position = "none",
                        axis.text = element_text(size = 20), axis.title = element_text(size=25), strip.text = element_text(size=20),
                        panel.background = element_rect(fill = "aliceblue", colour = "grey",
                        linewidth = 2, linetype = "solid"), panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid', colour = "grey"), 
                        axis.title.y = element_text(vjust=1, margin = margin(t=0, r=5, b=0, l=5)),
                        axis.title.x = element_text(vjust = -.05,margin = margin(t=10, r=0, b=5, l=0)),
                        panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid', colour = "grey"),
                        plot.title = element_text(hjust = .5, size=25))
      if (input$input_range == "Single Marker" && !is.null(input$rs_id) && !is.null(input$ancestry) && !is.null(filtered_hist_data())) {
        ra_values <- as.numeric(filtered_hist_data()$`Relative Accuracy`)
        snp_ra <- as.numeric(filtered_snp_data()$`Relative Accuracy`)
        snp_percentile <- ecdf(ra_values)(snp_ra) * 100
        
        ggplot() +
          geom_histogram(aes(x = ra_values), fill = "deepskyblue2", color = "midnightblue", binwidth = 0.02) +
          geom_vline(xintercept = snp_ra, color = "darkred", linewidth = 1) +
          labs(title = "Relative Accuracy Histogram (input in red)\n", x = "Relative Accuracy Value", y = "Frequency") +
          annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 2, 
                   label = paste("SNP Relative Accuracy Percentile: ", round(snp_percentile, 2), "%"), color = "black", size = 6) +
          theme.ggplot
        
      } else if (input$input_range == "Range of Markers (within chromosome)" || input$input_range == "Comma-separated List of SNP RS IDs" || 
                 input$input_range == "Single Gene" && !is.null(input$ancestry) && !is.null(filtered_hist_data())) {
        ra_values <- as.numeric(filtered_hist_data()$`Relative Accuracy`)
        snp_ra <- as.numeric(filtered_snp_data()$`Relative Accuracy`)
        
        p1 <- ggplot()
        if (length(table(snp_ra)) > 10) {
          p1 +
            geom_histogram(aes(x = ra_values), fill = "deepskyblue2", color = "midnightblue", binwidth = 0.02) +
            geom_histogram(aes(x = snp_ra), fill = "darkred", alpha = .4, color = "darkred", binwidth = .12) +
            labs(title = "Relative Accuracy Histogram (inputs in red)\n", x = "Relative Accuracy Value", y = "Frequency") +
            theme.ggplot
        } else {
          p1 +
            geom_histogram(aes(x = ra_values), fill = "deepskyblue2", color = "midnightblue", binwidth = 0.02) +
            geom_vline(xintercept = snp_ra, color = "darkred", linewidth = 1) +
            labs(title = "Relative Accuracy Histogram (inputs in red)\n", x = "Relative Accuracy Value", y = "Frequency") +
            theme.ggplot
        }
      }
    })
    
    output$download_table <- downloadHandler(
      filename = function() {
        paste("table_", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        if (!is.null(input$rs_id) && !is.null(input$ancestry) && !is.null(filtered_snp_data())) {
          filtered_table <- filtered_snp_data()[, c("SNP", "Chromosome", "BP position", "Allele", "Relative Accuracy", paste0("Corr. EU=>", ancestry_label()), paste0("Corr. EU=>", ancestry_label(), " S.E."), "Corr. EU=>EU", "Corr. EU=>EU S.E.","Gene")]
          data.frame(filtered_table, check.names = FALSE)
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
