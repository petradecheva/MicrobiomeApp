library(shiny)
library(shinythemes)
library(tidyr)
library(tidyverse)
library(readxl)
library(ggtext)
library(RColorBrewer)
library(dplyr)
library(bslib)
library(shinydashboard)
library(tools)
library(gridExtra)
library(cowplot)

tid <<- 0 
data_list <<- list()
pl_list <<- list()
ora_list <<- list()
ui <- fluidPage(
  theme = bs_theme(version = 4, bootswatch = "minty"),
  navbarPage("Microbiome", tabsetPanel(id = "tabs",
                                       tabPanel(
                                         title = "Home",
                                         value = "home",
                                         HTML(r"(<div class="row">
                                              <div class="col-sm-5" role="main">
                                              <h3>Generate plots</h3>
                                              <div class="form-group shiny-input-container">
                                              <label class="control-label" id="probe1-label" for="probe1">Metadata .XLSX File</label>
                                              <div class="input-group">
                                              <label class="input-group-btn input-group-prepend">
                                              <span class="btn btn-default btn-file">
                                              Browse
                                              <input id="probe1" name="probe1" type="file" style="position: absolute !important; top: -99999px !important; left: -99999px !important;" accept="xlsx" class="shiny-bound-input">
                                              </span>
                                              </label>
                                              <input type="text" class="form-control" placeholder="No file selected" readonly="readonly">
                                              </div>
                                              <div id="probe1_progress" class="progress active shiny-file-input-progress">
                                              <div class="progress-bar"></div>
                                              </div>
                                              </div>
                                              <div class="form-group shiny-input-container">
                                              <label class="control-label" id="probe2-label" for="probe2">Otu_counts (subsample.shared) .TSV File</label>
                                              <div class="input-group">
                                              <label class="input-group-btn input-group-prepend">
                                              <span class="btn btn-default btn-file">
                                              Browse
                                              <input id="probe2" name="probe2" type="file" style="position: absolute !important; top: -99999px !important; left: -99999px !important;" accept="tsv" class="shiny-bound-input">
                                              </span>
                                              </label>
                                              <input type="text" class="form-control" placeholder="No file selected" readonly="readonly">
                                              </div>
                                              <div id="probe2_progress" class="progress active shiny-file-input-progress">
                                              <div class="progress-bar"></div>
                                              </div>
                                              </div>
                                              <div class="form-group shiny-input-container">
                                              <label class="control-label" id="probe3-label" for="probe3">Taxonomy .TSV File</label>
                                              <div class="input-group">
                                              <label class="input-group-btn input-group-prepend">
                                              <span class="btn btn-default btn-file">
                                              Browse
                                              <input id="probe3" name="probe3" type="file" style="position: absolute !important; top: -99999px !important; left: -99999px !important;" accept="tsv" class="shiny-bound-input">
                                              </span>
                                              </label>
                                              <input type="text" class="form-control" placeholder="No file selected" readonly="readonly">
                                              </div>
                                              <div id="probe3_progress" class="progress active shiny-file-input-progress">
                                              <div class="progress-bar"></div>
                                              </div>
                                              </div>
                                              <div class="form-group shiny-input-container">
                                              <label class="control-label" id="caption-label" for="caption">Name of the plot</label>
                                              <input id="caption" type="text" class="form-control shiny-bound-input" value="Example: plot1">
                                              </div>
                                              <pre class="shiny-text-output noplaceholder shiny-bound-output" id="value" aria-live="polite"></pre>
                                              <button id="add" type="button" class="btn btn-default action-button shiny-bound-input">
                                              <i class="fa fa-plus-circle" role="presentation" aria-label="plus-circle icon"></i>
                                              Add
                                              </button>
                                              </div>
                                              <div class="col-sm-4">
                                              <br>
                                              <form class="well" role="complementary">
                                              <div class="row">
                                              <div class="col-sm-12">
                                              <span class="help-block">Read me:</span>
                                              </div>
                                              </div>
                                              <div class="row">
                                              <div class="col-sm-12">
                                              <span class="help-block">The required tables need to have information about:...</span>
                                              </div>
                                              </div>
                                              </form>
                                              </div>
                                              <br>
                                              </div>)")
                                       )
  ))
  
)

server <- function(input, output, session) {
  shinyInput <- function(name, id) paste(name, id, sep = "_")
  rv <- reactiveValues(counter = 0L)
  
  #go to the new created tab
  observeEvent(input$add, {
    rv$counter <- rv$counter + 1L
    updateTabsetPanel(session, "tabs", shinyInput("new_tab", rv$counter))
  }, ignoreInit = TRUE)
  
  #creating tab with plot after uploading the files
  observeEvent(input$add, {
    tid <<- tid + 1
    inputs <- reactiveValues(input1 = input$probe1, input2 = input$probe2, input3 = input$probe3)
    # print(inputs)
    # print(paste0("This is input1: ", inputs$input1))
    #data_list <<- append(data_list,
    
    #data processing
    metadata <- read_excel(inputs$input1$datapath, na="NA") %>%
      select(sample_id, disease_stat) %>%
      drop_na(disease_stat)
    
    otu_counts <- read_tsv(inputs$input2$datapath) %>%
      select(Group, starts_with("Otu")) %>%
      rename(sample_id = Group) %>%
      pivot_longer(-sample_id, names_to="otu", values_to = "count")
    
    taxonomy <- read_tsv(inputs$input3$datapath) %>%
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
    
    taxon_rel_abund <- otu_rel_abund %>%
      filter(level=="phylum") %>%
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
      summarize(pool = max(mean_rel_abund) < 3, 
                mean = mean(mean_rel_abund),
                .groups="drop")
    
    df <- inner_join(taxon_rel_abund, taxon_pool, by="taxon") %>%
      mutate(taxon = if_else(pool, "Other", taxon)) %>%
      group_by(disease_stat, taxon) %>%
      summarize(mean_rel_abund = sum(mean_rel_abund),
                mean = min(mean),
                .groups="drop") %>%
      mutate(taxon = factor(taxon),
             taxon = fct_reorder(taxon, mean, .desc=TRUE),
             taxon = fct_shift(taxon, n=1))
    #print(c(df))
    #print(df$taxon)
    data_list[[tid]] <<- as.data.frame(df)
    ora_list[[tid]] <<- as.data.frame(otu_rel_abund)
    
    #creating plot
    pl_list[[tid]] <- ggplot(data_list[[tid]], aes(x="", y=mean_rel_abund, fill=taxon)) +
      geom_bar(stat="identity", width=1) +
      # facet_grid(.~ dfg$disease_stat)+
      theme_classic() +
      theme(axis.line = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank())
    
    #print(paste("Before append:", data_list[[tid]]$disease_stat ))
    #output$plot <- render....
    
    #finally creating tab after data processing and plot creation
    appendTab(inputId = "tabs", tabPanel(title = input$caption, value = tid, 
                                         headerPanel('Microbiome'),
                                         sidebarLayout(
                                           sidebarPanel(selectInput(inputId = shinyInput("case", rv$counter), label = strong("Sample"),
                                                                    choices = unique(data_list[[tid]]$disease_stat),
                                                                    selected = "DiarrhealControl"),
                                                        sliderInput(inputId = shinyInput("pool", rv$counter), label = strong("Treshold"),
                                                                    min = 0, max = 100, value = 3),
                                                        checkboxInput(inputId = shinyInput("panel", rv$counter), strong("Panel view with every sample"), value = F),
                                                        selectInput(inputId = shinyInput("plot", rv$counter), label = strong("Plot output"),
                                                                    choices = list("bar chart" = "bar","pie chart" = "pie")),
                                                        selectInput(inputId = shinyInput("class", rv$counter), label = strong("Taxonomic level"),
                                                                    choices = unique(ora_list[[tid]]$level), selected = "phylum"),
                                                        # radioButtons(inputId = shinyInput("format", rv$counter), label ="select file type", choices = list("png", "pdf")),
                                                        # downloadButton(outputId = "down", label = "Download")
                                                        actionButton(shinyInput("remove_btn", rv$counter), "Remove", icon = icon("minus-circle"))),
                                           mainPanel(
                                             plotOutput(paste("pl_list", tid, sep = "_")),
                                             
                                           )
                                         )
                                         
    ), select = TRUE)

    
    #print(pl_list[[tid]])
    #print(paste("After append", data_list[[tid]]$disease_stat))
    ##########################
    #print(inputs)
    
  })
  
  ## REACTIVITY TO ARRANGE TAB NAMES:
  current.tab <- eventReactive(input$tabs, {
    # don't accidentally remove main tab:
    if (!identical(input$tabs, "home")) {
      input$tabs
    } else {
      NULL
    }
  })
  
  
  customRender <- function(x)
  {
    taxon_rel_abund <- ora_list[[tid]] %>%
      filter(level==input[[paste("class", x, sep = "_")]]) %>%
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
      summarize(pool = max(mean_rel_abund) < input[[paste("pool", x, sep = "_")]], 
                mean = mean(mean_rel_abund),
                .groups="drop")
    dfg <- inner_join(taxon_rel_abund, taxon_pool, by="taxon") %>%
      mutate(taxon = if_else(pool, "Other", taxon)) %>%
      group_by(disease_stat, taxon) %>%
      summarize(mean_rel_abund = sum(mean_rel_abund),
                mean = min(mean),
                .groups="drop") %>%
      mutate(taxon = factor(taxon),
             taxon = fct_reorder(taxon, mean, .desc=TRUE),
             taxon = fct_shift(taxon, n=1))
      if(input[[paste("panel", x, sep = "_")]]==TRUE){
        dfg
        print("panel")
        print(paste("disease stat:", dfg$disease_stat))
      }else {
        dfg <- dfg %>% filter(disease_stat == input[[paste("case", x, sep = "_")]])
        print("not panel")
        print(paste("disease_stat:", dfg$disease_stat))
        
      }
# filter(disease_stat == input[[paste("case", x, sep = "_")]])
    #print(paste("this is dfg disease_stat:", dfg$disease_stat))
    v2 <- rainbow(length(dfg$taxon))
    names(v2) <- unique(dfg$taxon)
    print(length(names(v2)))
    print(paste(dfg$taxon, " ", round(dfg$mean_rel_abund, 1),"%"))
    if (length(names(v2)) == 1)
    {
      # v2["Other"] = "#FFFFFF"
    }
    else
    {
      v2["Other"] = "#808080"

    }
    #print(pl_list[[x]])
    #overwriting the first created plot on the tab with the newly created one based on the new data
    barchart <- function(){
      
      if(input[[paste("panel", x, sep = "_")]]){
        samples <- unique(dfg$disease_stat)
        ptlist <- list()
        for( i in 1:length(samples)) {
          plist <- filter(dfg, disease_stat == samples[[i]])
          v2 <- rainbow(length(plist$taxon))
          names(v2) <- unique(plist$taxon)
          print(length(names(v2)))
          if (length(names(v2)) == 1)
          {
            v2["Other"] = "#FFFFFF"
          }
          else
          {
            v2["Other"] = "#808080"
          }
          pl <- ggplot(plist, aes(x=disease_stat, y=mean_rel_abund, fill=taxon)) +
            geom_col() +
            scale_y_continuous(expand=c(0, 0)) +
            labs(x=plist$disease_stat,
                 y="Mean Relative Abundance (%)",
                 fill = input[[paste("class", x, sep = "_")]]) +
            scale_fill_manual(values = v2, labels=paste(gsub('\\*', '', plist$taxon), str_replace_all(paste(round(plist$mean_rel_abund, 1),"%"), " ", "")))+
            theme_classic() +
            theme(axis.text.x = element_blank(),
                  legend.text = element_markdown(),
                  legend.key.size = unit(18, "pt")) 
          
          ptlist[[i]] <- pl 
        }
        # plot_grid(ptlist, ncol = 2)
        # do.call("grid.arrange", c(ptlist, ncol = 2))  
        do.call("plot_grid", c(ptlist, ncol = 2))
        # grid.arrange(plist, ncol=2)
        # plist + plot_layout(ncol = 2)
        
      } else{
        v2 <- rainbow(length(dfg$taxon))
        names(v2) <- unique(dfg$taxon)
        print(length(names(v2)))
        print(paste(dfg$taxon, " ", round(dfg$mean_rel_abund, 1),"%"))
        if (length(names(v2)) == 1)
        {
          v2["Other"] = "#FFFFFF"
        }
        else
        {
          v2["Other"] = "#808080"
        }
        ggplot(dfg, aes(x=disease_stat, y=mean_rel_abund, fill=taxon)) +
          geom_col() +
          scale_y_continuous(expand=c(0, 0)) +
          labs(x=input$case,
               y="Mean Relative Abundance (%)",
               fill = input$class) +
          theme_classic() +
          theme(axis.text.x = element_blank(),
                legend.text = element_markdown(),
                legend.key.size = unit(18, "pt")) +
          scale_fill_manual(values = v2, labels=paste(gsub('\\*', '', dfg$taxon), str_replace_all(paste(round(dfg$mean_rel_abund, 1),"%"), " ", "")))
      }
      
    }
    piechart <- function() {
      if(input[[paste("panel", x, sep = "_")]]){
        samples <- unique(dfg$disease_stat)
        ptlist <- list()
        for( i in 1:length(samples)) {
          plist <- filter(dfg, disease_stat == samples[[i]])
          v2 <- rainbow(length(plist$taxon))
          names(v2) <- unique(plist$taxon)
          print(length(names(v2)))
          if (length(names(v2)) == 1)
          {
            v2["Other"] = "#FFFFFF"
          }
          else
          {
            v2["Other"] = "#808080"
          }
          pl <- ggplot(plist, aes(x="", y=mean_rel_abund, fill=taxon)) +
            geom_bar(stat="identity", width=1) +
            labs(x=NULL,
                 y=plist$disease_stat,
                 fill = input[[paste("class", x, sep = "_")]]) +
            coord_polar("y", start=0) +
            scale_fill_manual(values = v2, labels=paste(gsub('\\*', '', plist$taxon), str_replace_all(paste(round(plist$mean_rel_abund, 1),"%"), " ", "")))+
            theme_classic() +
            theme(axis.line = element_blank(),
                  axis.text = element_blank(),
                  axis.ticks = element_blank(),
                  legend.text = element_markdown(),
                  legend.key.size = unit(18, "pt"))
          
          # renderBillboarder({billboarder() %>% bb_piechart(plist) %>% bb_legend(position = 'right')})
          
          
          ptlist[[i]] <- pl
        }
        do.call("plot_grid", c(ptlist, ncol = 2))
        
      } else{
        v2 <- rainbow(length(dfg$taxon))
        names(v2) <- unique(dfg$taxon)
        print(length(names(v2)))
        print(paste(dfg$taxon, " ", round(dfg$mean_rel_abund, 1),"%"))
        if (length(names(v2)) == 1)
        {
          # v2["Other"] = "#FFFFFF"
        }
        else
        {
          v2["Other"] = "#808080"
        }
        ggplot(dfg, aes(x="", y=mean_rel_abund, fill=taxon)) +
          geom_bar(stat="identity", width=1) +
          labs(x=NULL,
               y=input[[paste("case", x, sep = "_")]],
               fill = input[[paste("class", x, sep="_")]]) +
          coord_polar("y", start=0) +
          scale_fill_manual(values = v2, labels=paste(gsub('\\*', '', dfg$taxon), str_replace_all(paste(round(dfg$mean_rel_abund, 1),"%"), " ", "")))+
          # facet_grid(.~ data2()$disease_stat)+
          theme_classic() +
          theme(axis.line = element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                legend.text = element_markdown(),
                legend.key.size = unit(18, "pt"))
        # renderBillboarder({billboarder() %>% bb_piechart(data2()) %>% bb_legend(position = 'right')})
        
      }
      
    }
  
   
    pl_list[[x]] <<- 
      if(input[[paste("plot", x, sep = "_")]] == "bar"){
        barchart()
        
      }
      else if(input[[paste("plot", x, sep = "_")]] == "pie") {
        piechart()
        
      }
  
    
    output[[paste("pl_list", x, sep = "_")]] <- renderPlot({ pl_list[[x]] })
  }
  
  #observe event for when a choice from a filter is selected
  observe({
    if (rv$counter > 0L) {
      lapply(seq(rv$counter), function(x) {
        observeEvent(input[[paste("case", x, sep = "_")]],{
          customRender(x)
        })
        observeEvent(input[[paste("plot", x, sep = "_")]],{
          customRender(x)
        })
        observeEvent(input[[paste("panel", x, sep = "_")]],{
          customRender(x)
        })
        observeEvent(input[[paste("pool", x, sep = "_")]],{
          customRender(x)
        })
        observeEvent(input[[paste("class", x, sep = "_")]],{
          customRender(x)
        })
      })
    }
  })
  
  #removing a tab
  observe({
    if (rv$counter > 0L) {
      lapply(seq(rv$counter), function(x) {
        observeEvent(input[[paste("remove_btn", x, sep = "_")]], {
          #print(paste0("This is x: ",x))
          removeTab(inputId = "tabs", target = current.tab())
          #print(paste0("Removing: ", input[[paste("remove_btn", x, sep = "_")]]))
        })
      })
    }
  })
  
  output$text <- renderText({paste0("You are viewing tab \"", input$tabs, "\"", " Rv$counter is: ", rv$counter)})
  
}

shinyApp(ui, server)

