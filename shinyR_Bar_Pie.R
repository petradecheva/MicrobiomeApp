library(shiny)
library(shinythemes)
library(tidyr)
library(tidyverse)
library(readxl)
library(ggtext)
library(RColorBrewer)
library(dplyr)
library(scales)
library(ggplot2)
library(gridExtra)
server <- function(session, input, output) {
    data <- reactive({
      req(input$probe1, input$probe2, input$probe3)
      metadata <- read_excel(input$probe1$datapath, na="NA") %>%
        select(sample_id, disease_stat) %>%
        drop_na(disease_stat)
      otu_counts <- read_tsv(input$probe2$datapath) %>%
        select(Group, starts_with("Otu")) %>%
        rename(sample_id = Group) %>%
        pivot_longer(-sample_id, names_to="otu", values_to = "count")
      taxonomy <- read_tsv(input$probe3$datapath) %>%
        select("OTU", "Taxonomy") %>%
        rename_all(tolower) %>%
        mutate(taxonomy = str_replace_all(taxonomy, "\\(\\d+\\)", ""),
               taxonomy = str_replace(taxonomy, ";$", "")) %>%
        separate(taxonomy,
                 into=c("kingdom", "phylum", "class", "order", "family", "genus"),
                 sep=";")
      
      otu_rel_abund <- inner_join(metadata, otu_counts, by="sample_id") %>%
        inner_join(., taxonomy, by="otu") %>%
        group_by(sample_id) %>%
        mutate(rel_abund = count / sum(count)) %>%
        ungroup() %>%
        select(-count) %>%
        pivot_longer(c("kingdom", "phylum", "class", "order", "family", "genus", "otu"),
                     names_to="level",
                     values_to="taxon") %>%
        mutate(disease_stat = factor(disease_stat,
                                     levels=c(unique(metadata$disease_stat))))
    })

    observe({
      req(data())
      updateSelectInput(session, "case", choices = unique(data()$disease_stat), selected = "DiarrhealControl")
      updateSelectInput(session, "class", choices = unique(data()$level), selected = "phylum")
    })
   data2 <- reactive({
      req(input$pool, input$class)
      taxon_rel_abund <- data() %>%
        filter(level==input$class) %>%
        group_by(disease_stat, sample_id, taxon) %>%
        summarize(rel_abund = sum(rel_abund), .groups="drop") %>%
        group_by(disease_stat, taxon) %>%
        summarize(mean_rel_abund = 100*mean(rel_abund), .groups="drop") %>%
        mutate(taxon = str_replace(taxon,
                                   "(.*)_unclassified", "Unclassified *\\1*"),
               taxon = str_replace(taxon,
                                   "^(\\S*)$", "*\\1*"))
      taxon_pool <- taxon_rel_abund %>%
        group_by(taxon) %>%
        summarize(pool = max(mean_rel_abund) < input$pool,
                  mean = mean(mean_rel_abund),
                  .groups="drop")
      p <- inner_join(taxon_rel_abund, taxon_pool, by="taxon") %>%
        mutate(taxon = if_else(pool, "Other", taxon)) %>%
        group_by(disease_stat, taxon) %>%
        summarize(mean_rel_abund = sum(mean_rel_abund),
                  mean = min(mean),
                  .groups="drop") %>%
        mutate(taxon = factor(taxon),
               taxon = fct_reorder(taxon, mean, .desc=TRUE),
               taxon = fct_shift(taxon, n=1)) 
        if(input$panel == FALSE){
          p %>% filter(disease_stat == input$case)
        }else{
          p
        }
  })
  
  myfilename<-reactive(paste("plot-",input$case,"-", input$class,"-", input$plot,".", input$format, sep = ""))
  barchart <- function(){
    if(input$panel){
      samples <- unique(data2()$disease_stat)
      # print(unique(data2()$disease_stat))
      ptlist <- list()
      # print(var, "var")
      for( i in 1:length(samples)) {
       plist <- filter(data2(), disease_stat == samples[[i]])
       # print("i")
       # print(i)
       # print("plist")
       # print(plist)
       # print("data2()")
       # print(data2())
        pl <- ggplot(plist, aes(x=disease_stat, y=mean_rel_abund, fill=taxon)) +
              geom_col() +
              scale_y_continuous(expand=c(0, 0)) +
              labs(x=plist$disease_stat,
                   y="Mean Relative Abundance (%)",
                   fill = input$class) +
              theme_classic() +
              theme(axis.text.x = element_blank(),
                    legend.text = element_markdown(),
                    legend.key.size = unit(10, "pt"))
            
         ptlist[[i]] <- pl 
         # print("ptlist")
         # print(ptlist)
      }
      do.call("grid.arrange", c(ptlist, ncol = length(samples)))  
    } else{
      ggplot(data2(), aes(x=disease_stat, y=mean_rel_abund, fill=taxon)) +
      geom_col() +
      scale_y_continuous(expand=c(0, 0)) +
      # facet_wrap(.~ data2()$disease_stat, scales = "free_x")+
      labs(x=input$case,
           y="Mean Relative Abundance (%)",
           fill = input$class) +
      theme_classic() +
      theme(axis.text.x = element_blank(),
            legend.text = element_markdown(),
            legend.key.size = unit(10, "pt"))
    }
      
  }
  piechart <- function() {
    if(input$panel){
      samples <- unique(data2()$disease_stat)
      # print(unique(data2()$disease_stat))
      ptlist <- list()
      # print(var, "var")
      for( i in 1:length(samples)) {
        plist <- filter(data2(), disease_stat == samples[[i]])
        # print("i")
        # print(i)
        # print("plist")
        # print(plist)
        # print("data2()")
        # print(data2())
        pl <- ggplot(plist, aes(x="", y=mean_rel_abund, fill=taxon)) +
          geom_bar(stat="identity", width=1) +
          labs(x=NULL,
               y=plist$disease_stat,
               fill = input$class) +
          coord_polar("y", start=0) +
          # facet_grid(.~ data2()$disease_stat)+
          theme_classic() +
          theme(axis.line = element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank())
        
        ptlist[[i]] <- pl 
        # print("ptlist")
        # print(ptlist)
      }
      do.call("grid.arrange", c(ptlist, ncol = length(samples)))  
      
    } else{
      ggplot(data2(), aes(x="", y=mean_rel_abund, fill=taxon)) +
        geom_bar(stat="identity", width=1) +
        labs(x=NULL,
             y=input$case,
             fill = input$class) +
        coord_polar("y", start=0) +
        # facet_grid(.~ data2()$disease_stat)+
        theme_classic() +
        theme(axis.line = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank())
    }
    
  }
  output$plot1 <- renderPlot({
    if(input$plot == "bar"){
      barchart()

    }
    else if(input$plot == "pie") {
      piechart()

    }
  })
  output$down <- downloadHandler(
      filename = myfilename,
      content = function(file){
        if(input$format == "pdf")
          pdf(file)
        else if(input$format == "png")
          png(file)
        if(input$plot == "bar"){
          print(
              barchart()
            )
        }
        else if(input$plot == "pie"){
          print(
             piechart()
          )

        }
        dev.off()
      }

    )
 }

ui <- fluidPage(
  headerPanel('Microbiome'),
  sidebarPanel(
    fileInput("probe1", "Metadata .XLSX File", accept = "xlsx", buttonLabel = "Browse"),
    fileInput("probe2", "Otu_counts (subsample.shared) .TSV File", accept = "tsv", buttonLabel = "Browse"),
    fileInput("probe3", "Taxonomy .TSV File", accept = "tsv", buttonLabel = "Browse"),
    sliderInput(inputId = "pool", label = strong("Treshold"),
                min = 0, max = 100, value = 3),
    selectInput(inputId = "case", label = strong("Sample"),
                choices = ""),
    checkboxInput("panel", strong("Panel view with every sample"), value = F),
    selectInput(inputId = "plot", label = strong("Plot output"),
                choices = list("bar chart" = "bar","pie chart" = "pie")),
    selectInput(inputId = "class", label = strong("Taxonomic level"),
                choices = "", selected = "phylum"),
    radioButtons(inputId = "format", label ="File type", choices = list("png", "pdf")),
    downloadButton(outputId = "down", label = "Download")
  ),
  mainPanel(
    plotOutput('plot1')
  )
)
shinyApp(ui = ui, server = server)

