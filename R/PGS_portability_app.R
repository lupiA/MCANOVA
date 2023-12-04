#' Run the PGS portability shiny app
#' 
#' @description Calls the R shiny app
#'
#' @export
PGS_portability_app <- function() {
  shinyApp(
    ui = source("inst/shiny_UI.R", local = TRUE)$value,
    server = source("inst/shiny_server.R", local = TRUE)$value
  )
}
