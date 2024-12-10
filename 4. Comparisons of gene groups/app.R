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
library(VennDiagram)


# Define UI for application that draws a histogram
ui <- fluidPage(
  # Title
  titlePanel("Comparisons of gene groups"),
  # Selection of experiments and genes
  fluidRow(
    column(1, style = "margin-top: 25px;",
           actionButton("clear1","Clear group 1")
    ),
    column(2,
           selectizeInput("group1", 
                          label = "Group 1:",
                          choices = c(""), 
                          multiple = T,
                          width = '500px',
                          options = list(delimiter = " ", create = T)),
    ),
    column(1, style = "margin-top: 25px;",
           actionButton("clear2","Clear group 2")
    ),
    column(2,
           selectizeInput("group2", 
                          label = "Group 2:",
                          choices = c(""), 
                          multiple = T,
                          width = '500px',
                          options = list(delimiter = " ", create = T))
    ),
    column(3,
          DT::dataTableOutput("table")
    ),
    column(3,
           plotOutput("Venn")
    )
  )
)
#_________________________________________________
# Define server logic required to draw a histogram
server <- function(input, output,session) {
  
  rv <- reactiveValues(list1 = NULL, list2 = NULL, number1 = NULL, number2 = NULL, number3 = NULL)
  
  # Main table  
  output$table <- renderDT({
    
    rv$list1 <- input$group1
    rv$number1 <- length(rv$list1)
    rv$list2 <- input$group2
    rv$number2 <- length(rv$list2)
    
    mutual_list <- intersect(rv$list1, rv$list2)
    rv$number3 <- length(mutual_list)
    table <- data.frame(mutual_list)
    colnames(table) <- ''
    
    return(table %>%
             DT::datatable(rownames = FALSE,
                           selection = list(mode = "single",target = 'column'),
                           extensions = c('Scroller','Buttons'), #"Buttons','Select'
                           options = list(deferRender = TRUE,
                                          scrollY = 300,
                                          scroller = TRUE,
                                          dom = 'Blfrtip',
                                          rowId = 0,
                                          selection = 'none',
                                          buttons = list(list(extend = 'copy', title = NULL)),
                                          lengthMenu = c(10, 50, 100, nrow(table)),
                                          scrollX = TRUE,
                                          columnDefs = list(list(className = 'dt-center', targets = '_all'))
                           )
             )
    )
    
    
  }, server = FALSE)
  
  output$Venn <- renderPlot({

    venn <- draw.pairwise.venn(area1      = rv$number1,
                               area2      = rv$number2,
                               cross.area = rv$number3,
                               category   = c("Group 1", "Group 2"),
                               fill = c("cornsilk3", "bisque4"),
                               alpha = 0.5,
                               cex = 1.5,
                               cat.pos = 1,
                               cat.dist = c(0.02, 0.02))
    
    print(venn)

  })
  

  observeEvent(input$clear1, {
    updateSelectizeInput(session, "group1", choices = c(""))
    rv$list1 <- NULL
  })
  observeEvent(input$clear2, {
    updateSelectizeInput(session, "group2", choices = c(""))
    rv$list2 <- NULL
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
