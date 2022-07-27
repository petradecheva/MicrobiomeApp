library(shiny)
library(shinythemes)
library(tidyr)
library(tidyverse)
library(readxl)
library(ggtext)
library(RColorBrewer)
library(dplyr)

metadata <- read_tsv("https://raw.githubusercontent.com/riffomonas/minimalR-raw_data/master/schubert.metadata.tsv", na="NA") %>%
  select(sample_id, disease_stat) %>%
  drop_na(disease_stat)


otu_counts <- read_tsv("https://raw.githubusercontent.com/riffomonas/minimalR-raw_data/master/schubert.subsample.shared") %>%
  select(Group, starts_with("Otu")) %>%
  rename(sample_id = Group) %>%
  pivot_longer(-sample_id, names_to="otu", values_to = "count")

taxonomy <- read_tsv("https://raw.githubusercontent.com/milievaiv/shinyexperiments/main/schubert.cons.taxonomy.tsv") %>%
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
                               levels=c("NonDiarrhealControl",
                                        "DiarrhealControl",
                                        "Case")))

server <- function(input, output) {
  myfilename<-reactive(paste("plot-",input$case,"-", input$class,"-", input$plot,".", input$format, sep = ""))
  data <- reactive({
    metadata <- read_tsv(input$probe1, na="NA") %>%
      select(sample_id, disease_stat) %>%
      drop_na(disease_stat)
    
    
    otu_counts <- read_tsv(input$probe2) %>%
      select(Group, starts_with("Otu")) %>%
      rename(sample_id = Group) %>%
      pivot_longer(-sample_id, names_to="otu", values_to = "count")
    
    taxonomy <- read_tsv(input$probe3) %>%
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
                                   levels=c("NonDiarrhealControl",
                                            "DiarrhealControl",
                                            "Case")))
    req(input$case,input$pool, input$class)
    taxon_rel_abund <- otu_rel_abund %>%
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
             taxon = fct_shift(taxon, n=1)) %>%
      filter(disease_stat == input$case)
  })

  barchart <- function(){
      ggplot(data(), aes(x=disease_stat, y=mean_rel_abund, fill=taxon)) +
      geom_col() +
      scale_x_discrete(breaks=c("NonDiarrhealControl",
                                "DiarrhealControl",
                                "Case"),
                       labels=c("Healthy",
                                "Diarrhea,<br>*C. difficile*<br>negative",
                                "Diarrhea,<br>*C. difficile*<br>positive")) +
      scale_y_continuous(expand=c(0, 0)) +
      labs(x=NULL,
           y="Mean Relative Abundance (%)") +
      theme_classic() +
      theme(axis.text.x = element_markdown(),
            legend.text = element_markdown(),
            legend.key.size = unit(10, "pt"))
  }
  piechart <- function(){
      ggplot(data(), aes(x="", y=mean_rel_abund, fill=taxon)) +
      geom_bar(stat="identity", width=1) +
      # geom_text(aes(label = round(mean_rel_abund, digits = 0)),
      #           position = position_stack(vjust = 0.5))+
      coord_polar("y", start=0) +

      labs(x=NULL,
           y="Mean Relative Abundance (%)") +
      theme_classic()
  }
  output$moreControls <- renderUI({
    tagList(
        textInput(inputId = "probe1", label = strong("enter raw .tsv metadata link")),
        textInput(inputId = "probe2", label = strong("enter raw .tsv otu counts data link")),
        textInput(inputId = "probe3", label = strong("enter raw .tsv taxonomy data link")),
        sliderInput(inputId = "pool", label = strong("choose value for grouping"),
                    min = 0, max = 100, value = 3),
        selectInput(inputId = "case", label = strong("Case"),
                    choices = unique(otu_rel_abund$disease_stat),
                    selected = "DiarrhealControl"),
        selectInput(inputId = "plot", label = strong("choose plot output"),
                    choices = list("bar chart" = "bar","pie chart" = "pie")),
        selectInput(inputId = "class", label = strong("choose taxonomic level"),
                    choices = unique(otu_rel_abund$level), selected = "phylum"),
        radioButtons(inputId = "format", label ="select file type", choices = list("png", "pdf")),
        downloadButton(outputId = "down", label = "Download")
    )
  })

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
  sidebarPanel(
    uiOutput("moreControls")
  ),
  mainPanel(
    plotOutput('plot1'),
  )
)
  



shinyApp(ui = ui, server = server)

