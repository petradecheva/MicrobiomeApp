#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(
  fileInput("upload", NULL, buttonLabel = "Upload...", multiple = TRUE),
  tableOutput("files")
)
server <- function(input, output, session) {
  output$files <- renderTable(input$upload)
  observe({
    if (is.null(input$upload)) return()
    str1 <- getwd()
    print(str1)
    str2 <- "\\www"
    dp <- paste(str1,str2, sep = "")
    file.copy(input$upload, dp, overwrite = FALSE)
  })
  
    #if (is.null(input$upload)) return()
    # str1 <- wd
    # str2 <- "/Plots/www"
    # dp <- paste(str1,str2, sep = "")
    # file.copy(input$upload, dp, overwrite = TRUE)
}
# Run the application 
shinyApp(ui = ui, server = server)
