#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(DT)
library(dplyr)
library(tidyverse)


# Define UI for application that draws a histogram
ui <- fluidPage(
  # Title
  titlePanel("Translation of gene names"),
  # Selection of experiments and genes
  sidebarLayout(
    sidebarPanel(width = 6,
      selectizeInput("group", 
                     label = "Genes to translate:",
                     choices = c(""), 
                     multiple = T,
                     width = '900px',
                     options = list(delimiter = " ", create = T)),
      actionButton("wb","Translate into WormBase ID", style="width: 210px"),
      actionButton("symbol","Translate into symbol", style="width: 210px"),
      actionButton("whole","Whole translation", style="width: 210px"),
      actionButton("clear","Clear", style="width: 210px")
    ),
  
  mainPanel(width = 4,
      DT::dataTableOutput("table")
    )
  )
)
#_________________________________________________
# Define server logic required to draw a histogram
server <- function(input, output,session) {
  
  rv <- reactiveValues(kind_of_translation = '', list = NULL)
  
  # Main table  
  output$table <- renderDT({
    
    c_elegans_annotations <- vroom::vroom("Caenorhabditis_elegans.WBcel235.cdna_ncrna.all.fa.fasta_headers",col_names = c("transcript_id","type","chromosome","gene_id","gene_biotype","transcript_biotype","symbol","description"),progress = T)
    c_elegans_annotations <- c_elegans_annotations[!(c_elegans_annotations$`transcript_biotype`=="ncRNA"),]
    c_elegans_annotations <- c_elegans_annotations %>% dplyr::rename(gene_id_version=gene_id) %>% dplyr::mutate(gene_id=gsub("^(.*)\\.\\d+","\\1",gene_id_version))
    c_elegans_annotations <- c_elegans_annotations %>% dplyr::rename(transcript_id_version=transcript_id) %>% dplyr::mutate(transcript_id=gsub("^(.*)\\.\\d+","\\1",transcript_id_version))
    c_elegans_annotations <- c_elegans_annotations %>% dplyr::rename(transcript_id_version2=transcript_id) %>% dplyr::mutate(transcript_id=gsub("[A-Za-z]$","",transcript_id_version2))
    
    rv$list <- input$group
    print(rv$list)
    
    translated_table <- c_elegans_annotations %>% 
      filter(symbol %in% rv$list | gene_id %in% rv$list | transcript_id %in% rv$list | transcript_id_version2 %in% rv$list)
    
    translated_table <- unique(translated_table, by = "gene_id")
    # translated_table <- translated_table[,c(7,9,11)]
    
    if (rv$kind_of_translation == 'wb')
    {table <- translated_table[,9]
    colnames(table) = ''}
    else if (rv$kind_of_translation == 'symbol')
    {table <- translated_table[,7]
    colnames(table) = ''}
    else if (rv$kind_of_translation == '')
    {table <- translated_table[,c(7,9,11)]}
    
    return(table %>%
             DT::datatable(rownames = FALSE,
                           selection = list(mode = "single",target = 'column'),
                           extensions = c('Scroller','Buttons'), #"Buttons','Select'
                           options = list(deferRender = TRUE,
                                          scrollY = 350,
                                          scroller = TRUE,
                                          dom = 'Blfrtip',
                                          rowId = 0,
                                          selection = 'none',
                                          buttons = list(list(extend = 'copy', title = NULL)),
                                          scrollX = TRUE,
                                          columnDefs = list(list(className = 'dt-center', targets = '_all'))
                                          )
                           )
           )

   
  }, server = FALSE)

  observeEvent(input$clear, {
    updateSelectizeInput(session, "group", choices = c(""))
    rv$list <- NULL
    rv$kind_of_translation <- ''
  })
  observeEvent(input$wb, {
    rv$kind_of_translation <- 'wb'
  })
  observeEvent(input$symbol, {
    rv$kind_of_translation <- 'symbol'
  })
  observeEvent(input$whole, {
    rv$kind_of_translation <- ''
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
