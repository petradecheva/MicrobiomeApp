
#! Please make a "www" folder in the working directory first !

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(
  fileInput("upload", NULL, buttonLabel = "Upload...", multiple = TRUE),
  tableOutput("files")
)
server <- function(input, output, session) {
  
  observe({
    if (is.null(input$upload)) return()
    str1 <- getwd()
    print(str1)
    str2 <- "/www"
    dp <- paste(str1,str2, sep = "")
    # file.create(dp)
    file.copy(input$upload, dp, overwrite = FALSE)
    output$files <- renderTable("1.tsv")
  })
  
  
    #if (is.null(input$upload)) return()
    # str1 <- wd
    # str2 <- "/Plots/www"
    # dp <- paste(str1,str2, sep = "")
    # file.copy(input$upload, dp, overwrite = TRUE)
}
shinyApp(ui = ui, server = server)
