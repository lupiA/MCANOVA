#' Run the PGS portability shiny app
#' 
#' @description Calls the R shiny app
#'
#' @export
PGS_portability_app <- function() {
  library(shiny)
  shinyApp(
    ui = sourceSystem(file.path("inst", "shiny_UI.R"), local = TRUE)$value,
    server = sourceSystem(file.path("inst", "shiny_server.R"), local = TRUE)$value
  )
}
