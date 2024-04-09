#' Run the PGS portability shiny app
#' 
#' @description Calls the R shiny app
#'
#' @export
PGS_portability_app <- function() {
  
  library(shiny)
  library(ggplot2)
  library(MCANOVA)
  
  # Define the UI
  ui <- fluidPage(
  titlePanel("PGS Portability App"),
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
        textInput("rs_id", "Enter SNP ID (e.g., rs4422948):",
                  value = "rs4422948")
      ),
      conditionalPanel(
        condition = "input.input_type == 'Base Pair Position & Chromosome' && 
              input.input_range == 'Single Marker'",
        numericInput("chromosome_input", "Enter Chromosome:", value = 1),
        numericInput("bp_position", "Enter Base Pair Position (in bp):", value = 835499)
      ),
      conditionalPanel(
        condition = "input.input_type2 == 'SNP RS ID' && input.input_range == 
              'Range of Markers (within chromosome)'",
        textInput("start_snp_id", "Enter Start SNP ID:",
                  value = "rs35062161"),
        textInput("end_snp_id", "Enter End SNP ID:",
                  value = "rs2074468")
      ),
      conditionalPanel(
        condition = "input.input_type2 == 'Base Pair Position & Chromosome' && 
              input.input_range == 'Range of Markers (within chromosome)'",
        numericInput("chromosome_input2", "Enter Chromosome:", value = 6),
        numericInput("start_bp_position", "Enter Start Base Pair Position in Mbp (e.g., part of the MHC region):",
                     value = 29.75),
        numericInput("end_bp_position", "Enter End Base Pair Position in Mbp:",
                     value = 30.25)
      ),
      conditionalPanel(
        condition = "input.input_range == 'Comma-separated List of SNP RS IDs'",
        textInput("snp_list", "Enter SNP List (e.g., rs4422948,rs6686302):",
                  value = c("rs4422948,rs6686302"))
      ),
      conditionalPanel(
        condition = "input.input_range == 'Single Gene'",
        textInput("gene_name", "Enter Gene Name (e.g., ABCG2):",
                  value = "ABCG2")
      ),
      br(),br(),br(),br(),br(),
      radioButtons("dataset.input", "SNP Set (Note, the HapMap set is not recommended.):",
                   choices = c("UK Biobank Arrays", "HapMap Variants")
      ),
    ),
    mainPanel(
      fluidRow(width = 10,
               column(10,
                      htmlOutput("error_message")
               )
      ),
      fluidRow(width = 10,
               column(10,
                      verbatimTextOutput("error_message2")
               )
      ),
      fluidRow(width = 9, 
               h3("Table 1:"),
               h5("The PGS portability estimates averaged across the selected input SNPs (relative 
                  accuracy, cross-ancestry R-squared [EU to target ancestry], and within-ancestry 
                  R-squared [EU to EU])."),
               tableOutput("table1")
      ),
      fluidRow(width = 10, 
        plotOutput("histogram"),
        h3("Figure 1:"),
        h5("Relative accuracy distribution. The genome-wide relative accuracy distribution is in blue for the 
           selected target ancestry group. The relative accuracy distribution for the subset of selected variants is
           shown in purple. The number of variants entering into the subset (if applicable) is noted in the 
           upper-right corner. The x-axis has been capped at a RA of 2.5 for plotting purposes.")
      ),
      fluidRow(width = 9, 
               h3("Table 2:"),
               h5("Relative accuracy map for each of the selected variants and target ancestry group. 
                  Each row is a SNP and the columns contain SNP information (rsID, chromosome, 
                  base pair position, reference allele), PGS portability measures (relative 
                  accuracy, cross-ancestry R-squared [EU to target ancestry], and within-ancestry 
                  R-squared [EU to EU]), corresponding SEs and MC error for each portability measure, 
                  and the annotated gene name if available."),
               downloadButton("download_table", "Download Table"),
               tableOutput("table2")
      )
    )
  )
)

# Define the server
server <- function(input, output, session) {
  
  dat <- reactive({
    if(input$dataset.input == "UK Biobank Arrays"){
      dat <- MCANOVA::MAP_UKB
    } else{
      dat <- MCANOVA::MAP_HAPMAP
    }
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
  
  start_bp_input <- reactive({
    if(input$input_type2 == 'Base Pair Position & Chromosome' && input$input_range == 'Range of Markers (within chromosome)'){
      input$start_bp_position * 1e6
    }
  })
  
  end_bp_input <- reactive({
    if(input$input_type2 == 'Base Pair Position & Chromosome' && input$input_range == 'Range of Markers (within chromosome)'){
      input$end_bp_position * 1e6
    }
  })
  
  filtered_snp_data <- reactive({
    snps <- dat()
    snps <- snps[snps$Cohort == input$ancestry, ]
    if(input$ancestry %in% c("African", "East Asian")){
      snps <- snps[-(which(is.na(snps[,7]))), ]
    }
    colnames(snps)[6:15] <- c("Relative Accuracy", 
                              paste0("Rsq. EU\u2192", ancestry_label()), 
                              "Rsq. EU\u2192EU", 
                              "MC Error RA", 
                              "S.E. RA", 
                              paste0("MC Error Rsq. EU\u2192", ancestry_label()), 
                              paste0("S.E. Rsq. EU\u2192", ancestry_label()), 
                              "MC Error Rsq. EU\u2192EU",
                              "S.E. Rsq. EU\u2192EU",
                              "Target Ancestry")
    
    snp_data <- NULL
    
    if (input$input_range == "Single Marker") {
      
      if (input$input_type == "SNP RS ID") {
        snp_data <- snps[which(snps$SNP == input$rs_id),
                         c("Target Ancestry", "SNP","Relative Accuracy", paste0("Rsq. EU\u2192", ancestry_label()), 
                           "Rsq. EU\u2192EU", "Chromosome", "BP position", "Allele", "Gene",
                           "MC Error RA", 
                           "S.E. RA", 
                           paste0("MC Error Rsq. EU\u2192", ancestry_label()), 
                           paste0("S.E. Rsq. EU\u2192", ancestry_label()), 
                           "MC Error Rsq. EU\u2192EU",
                           "S.E. Rsq. EU\u2192EU")]
      } 
      
      else if (input$input_type == "Base Pair Position & Chromosome") {
        snps_chr <- snps[which(snps$Chromosome == input$chromosome_input),]
        if(length(which(snps_chr$`BP position` == input$bp_position))>0){
          snp_data <- snps_chr[which(snps_chr$`BP position` == input$bp_position),
                            c("Target Ancestry", "SNP","Relative Accuracy", paste0("Rsq. EU\u2192", ancestry_label()), 
                              "Rsq. EU\u2192EU", "Chromosome", "BP position", "Allele", "Gene",
                              "MC Error RA", 
                              "S.E. RA", 
                              paste0("MC Error Rsq. EU\u2192", ancestry_label()), 
                              paste0("S.E. Rsq. EU\u2192", ancestry_label()), 
                              "MC Error Rsq. EU\u2192EU",
                              "S.E. Rsq. EU\u2192EU")]
        } else {
          lb <- max(which(snps_chr$`BP position` <= input$bp_position))
          ub <- min(which(snps_chr$`BP position` >= input$bp_position))
  
          if(length(lb) == 1 & length(ub) == 1){
            tmp_data <- snps_chr[lb:ub,
                            c("Target Ancestry", "SNP","Relative Accuracy", paste0("Rsq. EU\u2192", ancestry_label()), 
                              "Rsq. EU\u2192EU", "Chromosome", "BP position", "Allele", "Gene",
                              "MC Error RA", 
                              "S.E. RA", 
                              paste0("MC Error Rsq. EU\u2192", ancestry_label()), 
                              paste0("S.E. Rsq. EU\u2192", ancestry_label()), 
                              "MC Error Rsq. EU\u2192EU",
                              "S.E. Rsq. EU\u2192EU")]
            snp_data <- tmp_data[1,]
            snp_data[, "SNP"] <- paste0(tmp_data[1, "SNP"], "; ", tmp_data[2, "SNP"])
            snp_data[, "Allele"] <- paste0(tmp_data[1, "Allele"], "; ", tmp_data[2, "Allele"])
            snp_data[, "BP position"] <- paste0(tmp_data[1, "BP position"], "; ", tmp_data[2, "BP position"])
            snp_data[, "Gene"] <- paste0(tmp_data[1, "Gene"], "; ", tmp_data[2, "Gene"])
            
            snp_data[, "Relative Accuracy"] <- mean(tmp_data[, "Relative Accuracy"],na.rm=T)
            snp_data[, paste0("Rsq. EU\u2192", ancestry_label())] <- mean(tmp_data[, paste0("Rsq. EU\u2192", ancestry_label())],na.rm=T)
            snp_data[, "Rsq. EU\u2192EU"] <- mean(tmp_data[, "Rsq. EU\u2192EU"],na.rm=T)

            snp_data[, "S.E. RA"] <- mean(tmp_data[, "S.E. RA"],na.rm=T)
            snp_data[, paste0("S.E. Rsq. EU\u2192", ancestry_label())] <- mean(tmp_data[, paste0("S.E. Rsq. EU\u2192", ancestry_label())],na.rm=T)
            snp_data[, "S.E. Rsq. EU\u2192EU"] <- mean(tmp_data[, "S.E. Rsq. EU\u2192EU"],na.rm=T)
            snp_data[, "MC Error RA"] <- mean(tmp_data[, "MC Error RA"],na.rm=T)
            snp_data[, paste0("MC Error Rsq. EU\u2192", ancestry_label())] <- mean(tmp_data[, paste0("MC Error Rsq. EU\u2192", ancestry_label())],na.rm=T)
            snp_data[, "MC Error Rsq. EU\u2192EU"] <- mean(tmp_data[, "MC Error Rsq. EU\u2192EU"],na.rm=T)
          }

        }        
        
      }
      
    }  else if (input$input_range == "Range of Markers (within chromosome)") {
      
      if (input$input_type2 == "SNP RS ID") {
        start_snp_idx_1 <- which(snps$SNP == input$start_snp_id)
        end_snp_idx_1 <- which(snps$SNP == input$end_snp_id)
        
        if (length(start_snp_idx_1) > 0 && length(end_snp_idx_1) > 0 &&
            snps$Chromosome[start_snp_idx_1] == snps$Chromosome[end_snp_idx_1]) {
          snp_data <- snps[start_snp_idx_1:end_snp_idx_1,
                           c("Target Ancestry", "SNP","Relative Accuracy", paste0("Rsq. EU\u2192", ancestry_label()), 
                              "Rsq. EU\u2192EU", "Chromosome", "BP position", "Allele", "Gene",
                              "MC Error RA", 
                              "S.E. RA", 
                              paste0("MC Error Rsq. EU\u2192", ancestry_label()), 
                              paste0("S.E. Rsq. EU\u2192", ancestry_label()), 
                              "MC Error Rsq. EU\u2192EU",
                              "S.E. Rsq. EU\u2192EU")]
          
        }
        
      }  else if (input$input_type2 == "Base Pair Position & Chromosome") {
        snps <- snps[which(snps$Chromosome == input$chromosome_input2), ]
        
        if(is.null(start_bp_input()) || is.null(end_bp_input()) || is.na(start_bp_input()) || is.na(end_bp_input())){
          snp_data <- NULL
        } else if(length(start_bp_input():end_bp_input()) > 0 && start_bp_input() < end_bp_input() && 
                  start_bp_input() > 0 && nrow(snps)>0){
          snp_data <- snps[which(snps$`BP position` %in% start_bp_input():end_bp_input()),
                           c("Target Ancestry", "SNP","Relative Accuracy", paste0("Rsq. EU\u2192", ancestry_label()), 
                              "Rsq. EU\u2192EU", "Chromosome", "BP position", "Allele", "Gene",
                              "MC Error RA", 
                              "S.E. RA", 
                              paste0("MC Error Rsq. EU\u2192", ancestry_label()), 
                              paste0("S.E. Rsq. EU\u2192", ancestry_label()), 
                              "MC Error Rsq. EU\u2192EU",
                              "S.E. Rsq. EU\u2192EU")]
        }
      }
    }  else if (input$input_range == "Comma-separated List of SNP RS IDs") {
      snps_in_list <- unlist(strsplit(as.character(input$snp_list), split = ","))
      
      if (length(snps_in_list) > 0) {
        snp_data <- snps[which(snps$SNP %in% snps_in_list),
                         c("Target Ancestry", "SNP","Relative Accuracy", paste0("Rsq. EU\u2192", ancestry_label()), 
                              "Rsq. EU\u2192EU", "Chromosome", "BP position", "Allele", "Gene",
                              "MC Error RA", 
                              "S.E. RA", 
                              paste0("MC Error Rsq. EU\u2192", ancestry_label()), 
                              paste0("S.E. Rsq. EU\u2192", ancestry_label()), 
                              "MC Error Rsq. EU\u2192EU",
                              "S.E. Rsq. EU\u2192EU")]
      }
      
    }  else if (input$input_range == "Single Gene") {
      snps_in_list <- snps$SNP[which(snps$Gene == input$gene_name)]
      
      if (length(snps_in_list) > 0) {
        snp_data <- snps[which(snps$SNP %in% snps_in_list),
                         c("Target Ancestry", "SNP","Relative Accuracy", paste0("Rsq. EU\u2192", ancestry_label()), 
                              "Rsq. EU\u2192EU", "Chromosome", "BP position", "Allele", "Gene",
                              "MC Error RA", 
                              "S.E. RA", 
                              paste0("MC Error Rsq. EU\u2192", ancestry_label()), 
                              paste0("S.E. Rsq. EU\u2192", ancestry_label()), 
                              "MC Error Rsq. EU\u2192EU",
                              "S.E. Rsq. EU\u2192EU")]
      }
    }
    
    if (is.null(snp_data) || nrow(snp_data) == 0) {
      snp_data <- NULL
    }
    snp_data
  })
  
  filtered_hist_data <- reactive({
    hist <- dat()
    hist <- hist[hist$Cohort == input$ancestry, ]
    if(input$ancestry %in% c("African", "East Asian")){
      hist <- hist[-(which(is.na(hist[,7]))), ]
    }
  colnames(hist)[6:14] <- c("Relative Accuracy", 
                        paste0("Rsq. EU\u2192", ancestry_label()), 
                        "Rsq. EU\u2192EU", 
                        "MC Error RA", 
                        "S.E. RA", 
                        paste0("MC Error Rsq. EU\u2192", ancestry_label()), 
                        paste0("S.E. Rsq. EU\u2192", ancestry_label()), 
                        "MC Error Rsq. EU\u2192EU",
                        "S.E. Rsq. EU\u2192EU")
    hist_data <- hist[ , c("Relative Accuracy","Cohort")]
    if (nrow(hist_data) == 0) {
      hist_data <- NULL
    }
    if (!is.null(hist_data) && length(which(is.na(hist_data$`Relative Accuracy`))) > 0) {
      hist_data <- hist_data[-which(is.na(hist_data$`Relative Accuracy`)), ]
    }
    hist_data
  })
  
  output$table1 <- renderTable({
    if (!is.null(input$rs_id) && !is.null(input$ancestry) && !is.null(filtered_snp_data())) {
      if (!is.null(filtered_snp_data())) {
        averages <- colMeans(filtered_snp_data()[ , c("Relative Accuracy", paste0("Rsq. EU\u2192", ancestry_label()), "Rsq. EU\u2192EU")], na.rm = TRUE)
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
        modified_table <- filtered_snp_data()[, c("Target Ancestry", "SNP", "Chromosome", "BP position",  
                                                  "Relative Accuracy", 
                                                  paste0("Rsq. EU\u2192", ancestry_label()), 
                                                  "Rsq. EU\u2192EU",
                                                  "MC Error RA", 
                                                  "S.E. RA", 
                                                  paste0("MC Error Rsq. EU\u2192", ancestry_label()), 
                                                  paste0("S.E. Rsq. EU\u2192", ancestry_label()), 
                                                  "MC Error Rsq. EU\u2192EU",
                                                  "S.E. Rsq. EU\u2192EU",
                                                  "Allele", "Gene")]

        colnames(modified_table)[3]="Chr"
        data.frame(modified_table,check.names = FALSE)
      }
    }
  })
  
  output$histogram <- renderPlot({
    if(input$dataset.input == "UK Biobank Arrays"){
      theme.ggplot <- theme(legend.position = "none",
                        axis.text = element_text(size = 17), axis.title = element_text(size=19), strip.text = element_text(size=20),
                        panel.background = element_rect(fill = "aliceblue", colour = "grey",
                        linewidth = 2, linetype = "solid"), panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid', colour = "grey"), 
                        axis.title.y = element_text(vjust=1, margin = margin(t=0, r=5, b=0, l=5)),
                        axis.title.x = element_text(vjust = -.05,margin = margin(t=10, r=0, b=5, l=0)),
                        panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid', colour = "grey"),
                        plot.title = element_blank())
      if (input$input_range == "Single Marker" && !is.null(input$rs_id) && !is.null(input$ancestry) && !is.null(filtered_hist_data())) {
        ra_values <- as.numeric(filtered_hist_data()$`Relative Accuracy`)
        snp_ra <- as.numeric(filtered_snp_data()$`Relative Accuracy`)
        snp_percentile <- ecdf(ra_values)(snp_ra) * 100
        
        label_text <- paste0(round(snp_percentile, 2), "%")
        ggplot() +
          geom_histogram(aes(x = ra_values), fill = "skyblue", color = "midnightblue", binwidth = 0.02) +
          geom_vline(xintercept = snp_ra, color = "white", linewidth = 1.7) +
          geom_vline(xintercept = snp_ra, color = "#2E0854", linewidth = 1.3) +
          labs(x = "Relative Accuracy", y = "Frequency") +
          annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.4, 
                   label = paste0("Relative Accuracy\nPercentile: ", label_text, sep = ""), 
                   color = "#2E0854", size = 6, fontface = 2) +
          xlim = c(0,2.5)+
          theme.ggplot
        
      } else if (input$input_range == "Range of Markers (within chromosome)" || input$input_range == "Comma-separated List of SNP RS IDs" || 
                 input$input_range == "Single Gene" && !is.null(input$ancestry) && !is.null(filtered_hist_data())) {
        ra_values <- as.numeric(filtered_hist_data()$`Relative Accuracy`)
        snp_ra <- as.numeric(filtered_snp_data()$`Relative Accuracy`)
        
        p1 <- ggplot()
        if(input$ancestry == "East Asian"){
          y_val = 10000
        } else{
          y_val = 20000
        }
        if (length(table(snp_ra)) > 10) {
          p1 +
            geom_histogram(aes(x = ra_values), fill = "skyblue", color = "midnightblue", binwidth = 0.02) +
            geom_violin(aes(x = snp_ra, y = y_val), fill = "#2E0854", alpha = .35, color = "#2E0854", width = 7500, linewidth = 1.2) +
            annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.4,
                     label = paste0("Number of SNPs\nin Input: ", length(snp_ra), sep=""), color = "#2E0854", size = 6, fontface = 2) +
            labs(x = "Relative Accuracy", y = "Frequency") +
            xlim = c(0,2.5)+
            theme.ggplot
        } else {
          p1 +
            geom_histogram(aes(x = ra_values), fill = "skyblue", color = "midnightblue", binwidth = 0.02) +
            geom_vline(xintercept = snp_ra, color = "white", linewidth = 1.7) +
            geom_vline(xintercept = snp_ra, color = "#2E0854", linewidth = 1.3) +
            labs(x = "Relative Accuracy", y = "Frequency") +
            xlim = c(0,2.5)+
            theme.ggplot
        }
      }
    } else{
      theme.ggplot <- theme(legend.position = "none",
                        axis.text = element_text(size = 17), axis.title = element_text(size=19), strip.text = element_text(size=20),
                        panel.background = element_rect(fill = "aliceblue", colour = "grey",
                        linewidth = 2, linetype = "solid"), panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid', colour = "grey"), 
                        axis.title.y = element_text(vjust=1, margin = margin(t=0, r=5, b=0, l=5)),
                        axis.title.x = element_text(vjust = -.05,margin = margin(t=10, r=0, b=5, l=0)),
                        panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid', colour = "grey"),
                        plot.title = element_blank())
      if (input$input_range == "Single Marker" && !is.null(input$rs_id) && !is.null(input$ancestry) && !is.null(filtered_hist_data())) {
        ra_values <- as.numeric(filtered_hist_data()$`Relative Accuracy`)
        snp_ra <- as.numeric(filtered_snp_data()$`Relative Accuracy`)
        snp_percentile <- ecdf(ra_values)(snp_ra) * 100
        
        label_text <- paste0(round(snp_percentile, 2), "%")
        ggplot() +
          geom_histogram(aes(x = ra_values), fill = "skyblue", color = "midnightblue", binwidth = 0.1) +
          geom_vline(xintercept = snp_ra, color = "white", linewidth = 1.7) +
          geom_vline(xintercept = snp_ra, color = "#2E0854", linewidth = 1.3) +
          labs(x = "Relative Accuracy", y = "Frequency") +
          xlim = c(0,2.5)+
          annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.4, 
                   label = paste0("Relative Accuracy\nPercentile: ", label_text, sep = ""), 
                   color = "#2E0854", size = 6, fontface = 2) +
          theme.ggplot
        
      } else if (input$input_range == "Range of Markers (within chromosome)" || input$input_range == "Comma-separated List of SNP RS IDs" || 
                 input$input_range == "Single Gene" && !is.null(input$ancestry) && !is.null(filtered_hist_data())) {
        ra_values <- as.numeric(filtered_hist_data()$`Relative Accuracy`)
        snp_ra <- as.numeric(filtered_snp_data()$`Relative Accuracy`)
        
        p1 <- ggplot()
        if(input$ancestry == "East Asian"){
          y_val = 10000
        } else{
          y_val = 20000
        }
        if (length(table(snp_ra)) > 10) {
          p1 +
            geom_histogram(aes(x = ra_values), fill = "skyblue", color = "midnightblue", binwidth = 0.1) +
            geom_violin(aes(x = snp_ra, y = y_val), fill = "#2E0854", alpha = .35, color = "#2E0854", width = 7500, linewidth = 1.2) +
            annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.4,
                     label = paste0("Number of SNPs\nin Input: ", length(snp_ra), sep=""), color = "#2E0854", size = 6, fontface = 2) +
            labs(x = "Relative Accuracy", y = "Frequency") +
            xlim = c(0,2.5)+
            theme.ggplot
        } else {
          p1 +
            geom_histogram(aes(x = ra_values), fill = "skyblue", color = "midnightblue", binwidth = 0.1) +
            geom_vline(xintercept = snp_ra, color = "white", linewidth = 1.7) +
            geom_vline(xintercept = snp_ra, color = "#2E0854", linewidth = 1.3) +
            labs(x = "Relative Accuracy", y = "Frequency") +
            xlim = c(0,2.5)+
            theme.ggplot
        }
      }
    }
  })
  
  output$download_table <- downloadHandler(
    filename = function() {
      paste("table_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      if (!is.null(input$rs_id) && !is.null(input$ancestry) && !is.null(filtered_snp_data())) {
        filtered_table <- filtered_snp_data()[, c("Target Ancestry", "SNP", "Chromosome", "BP position",  
                                                  "Relative Accuracy", 
                                                  paste0("Rsq. EU\u2192", ancestry_label()), 
                                                  "Rsq. EU\u2192EU",
                                                  "MC Error RA", 
                                                  "S.E. RA", 
                                                  paste0("MC Error Rsq. EU\u2192", ancestry_label()), 
                                                  paste0("S.E. Rsq. EU\u2192", ancestry_label()), 
                                                  "MC Error Rsq. EU\u2192EU",
                                                  "S.E. Rsq. EU\u2192EU",
                                                  "Allele", "Gene")]
        data.frame(filtered_table, check.names = FALSE)
        colnames(filtered_table)[c(1,4,5)] <- c("Ancestry", "BP_position", "RA")
        colnames(filtered_table)[6:7] <- c("Rsq.across", "Rsq.within")
        colnames(filtered_table)[8:13] <- c("MC_error_RA", "SE_RA","MC_error_Rsq.across","SE_Rsq.across","MC_error_Rsq.within","SE_Rsq.within")
        write.csv(filtered_table, file, row.names = FALSE)
      }
    }
  )
  
  output$error_message <- renderText({
    filtered_data <- filtered_snp_data()
    
    if (is.null(filtered_data)) {
      return(
        HTML(paste("<div style='white-space: pre-wrap;'><font color=\"#EB5406\" size=\"4\"><b>",
                   "Invalid inputs. Note that if an input is monomorphic for the\n",
                   "selected ancestry group it is not a valid input for that group.\n",
                   "</b></font><br></div>")
        )
      )
    } else {
      return(NULL)
    }
  })
  
  output$error_message2 <- renderText({
    
    filtered_data <- filtered_snp_data()
    
    if (is.null(filtered_data)) {
      
      if (input$input_range == "Range of Markers (within chromosome)" && input$input_type2 == "SNP RS ID") {
        return("Error: Note that the range of entries must be within-chromosome and SNP RS\nIDs must be in the UK Biobank genotyping array and are case sensitive.")
      }
      
      if (input$input_range == "Range of Markers (within chromosome)" && input$input_type2 == "Base Pair Position & Chromosome") {
        return("Error: Note that the range of entries must be within-chromosome and\nin Mbp units.")
      }
      
      if (input$input_range == "Comma-separated List of SNP RS IDs") {
        return("Error: Note that the comma-separated list of SNP RS IDs should not\nhave any spaces and only valid SNP RS IDs are allowed.")
      }
      
      if (input$input_range == "Single Gene") {
        return("Error: Note that gene name must be annotated in the UK Biobank array\nand is case sensitive.")
      }
      
      if (input$input_range == "Single Marker" && input$input_type == "SNP RS ID") {
        return("Error: Note that the SNP RS ID must be in the UK Biobank genotyping\narray and is case sensitive.")
      }
      
      if (input$input_range == "Single Marker" && input$input_type == "Base Pair Position & Chromosome") {
        return("Error: Note that a single BP position entry must be within range\n and in base pair units.")
      }
      
    } else {
      return(NULL)
    }
  })

}
  
  # Run
  shinyApp(ui = ui, server = server)
                             
}
