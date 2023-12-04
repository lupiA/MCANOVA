#' Run the PGS portability shiny app
#' 
#' @description Calls the R shiny app
#'
#' @export
PGS_portability_app <- function() {
  shinyApp(
    ui = source("shinyApp/shiny_UI.R", local = TRUE)$value,
    server = source("shinyApp/shiny_server.R", local = TRUE)$value
  )
}