#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(readxl)
library(dplyr)
library(tidyverse)
library(ggplot2)

options(shiny.maxRequestSize = 30*1024^2)

# Define UI for application that draws a histogram
ui <- fluidPage(
    # Application title
    titlePanel("qPCR analysis"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
           fileInput("file", "Choose file:", multiple = FALSE, accept = c(".xls", ".xlsx")),
           # tags$hr(),
           selectInput("reference", "Reference gene:", choices = ""),
           selectInput("control", "Control condition:", choices = ""),
           selectInput("target", "Target gene to show:", choices = "")
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("plot")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output, session){
  
  rv <- reactiveValues()
  qpcr_results <- data_frame()
   
  observeEvent(input$file, {
    qpcr_results <- read_excel(input$file$datapath, sheet = "Results")
    colnames(qpcr_results) <- qpcr_results[51,]
    qpcr_results <- qpcr_results[-c(1:51),]
    qpcr_results <- dplyr::slice(qpcr_results, 1:(n()-5))
    qpcr_results <- qpcr_results[,c(4,5,15)]
    qpcr_results <- qpcr_results[!(qpcr_results$`Target Name`=="-RT"),]
    qpcr_results$CT <- as.numeric(qpcr_results$CT)
    
    qpcr_results2 <- qpcr_results %>% 
      group_by(`Target Name`, `Sample Name`) %>% 
      summarise(ct_mean = mean(CT))
    
    qpcr_results2 <- qpcr_results2 %>% 
      mutate(ct_2 = 2^-ct_mean)

    rv$qpcr_results3 <- split(qpcr_results2,qpcr_results2$`Target Name`)

    samples <- unique(qpcr_results$`Sample Name`)
    rv$targets <- unique(qpcr_results$`Target Name`)
    conditions <- unique(gsub(' [0-9]*','',samples))
    
    updateSelectInput(session, "reference", "Reference gene:", choices = rv$targets)
    updateSelectInput(session, "control", "Control condition:", choices = conditions)
  })
  
  observeEvent(input$control, {
    rv$control2 <- input$control
  })
  
  observeEvent(input$reference, {
    rv$control <- input$reference
    targets2 <- rv$targets[rv$targets != rv$control]
    updateSelectInput(session, "target", "Target gene to show:", choices = targets2)
  })
  
  observeEvent(input$target, {
    rv$target <- input$target
  })
  
  output$plot <- renderPlot({
    
    qpcr_results3 <- rv$qpcr_results3
    control <- rv$control
    control2 <- rv$control2
    target <- rv$target
    
    for (i in 1:length(qpcr_results3)){
      ifelse(qpcr_results3[[i]]$`Target Name` == control,rv$control_df <- qpcr_results3[[i]],print(''))
    }
    
    for (i in 1:length(qpcr_results3)){
      qpcr_results3[[i]]$delta_ct <- qpcr_results3[[i]]$ct_2/rv$control_df$ct_2
      qpcr_results3[[i]]$group <- gsub(' [0-9]*','',qpcr_results3[[i]]$`Sample Name`)
      
      qpcr_results3[[i]] <- qpcr_results3[[i]] %>% 
        group_by(group) %>% 
        mutate(delta_mean = mean(delta_ct))
      
      for (j in 1:nrow(qpcr_results3[[i]])){
        if (qpcr_results3[[i]]$group[j] == control2){
          control_delta_mean <- qpcr_results3[[i]]$delta_mean[j]
        }
      }
      
      qpcr_results3[[i]] <- qpcr_results3[[i]] %>%
        mutate(delta_delta = delta_ct/control_delta_mean)
      
      qpcr_results3[[i]] <- qpcr_results3[[i]] %>%
        group_by(group) %>%
        mutate(delta_delta_mean = mean(delta_delta)) %>%
        mutate(delta_delta_sd = sd(delta_delta))
    }
    
    for (i in 1:length(qpcr_results3)){
      # print(qpcr_results3[[i]])
      qpcr_results3[[i]] <- unique(qpcr_results3[[i]][,c(1,6,9,10)])
    }
    
    for (i in 1:length(qpcr_results3)){
      ifelse(qpcr_results3[[i]]$`Target Name` == target,rv$final_result <- qpcr_results3[[i]], print(''))
    }

    final_result <- rv$final_result
    # print(final_result)
    plot <- ggplot(final_result) +
        geom_bar(aes(x=factor(rv$final_result$group, level=c("pusty", "nspc", "N2", "ADZ31", "ADZ74", "ADZ83", "NSPC", "tm", "rtt5", "NSPC/tm", "NSPC/rtt5", "NSPC(1-10)/tm", "74_1", "83_1", "74_2", "83_2", "74_3", "83_3", "74_4", "83_4", "74_5", "83_5", "bez","atx","fndc","larp")), y=final_result$delta_delta_mean), stat="identity", width = 0.8, color="black",fill="pink", alpha=0.7) + ggtitle(final_result[1,1]) +
        geom_errorbar(aes(x=final_result$group, ymin=final_result$delta_delta_mean-final_result$delta_delta_sd, ymax=final_result$delta_delta_mean+final_result$delta_delta_sd), width=0.1, colour="black", alpha=0.9, size=0.5) + xlab("") + ylab("Expression level") +   geom_text(aes(x=factor(final_result$group, level=c("pusty", "nspc", "N2", "ADZ31", "ADZ74", "ADZ83","NSPC", "tm", "rtt5", "NSPC/tm", "NSPC/rtt5", "NSPC(1-10)/tm", "74_1", "83_1", "74_2", "83_2", "74_3", "83_3", "74_4", "83_4", "74_5", "83_5", "bez","atx","fndc","larp")), y=final_result$delta_delta_mean, label = round(final_result$delta_delta_mean, digits = 2), hjust = 1.8, vjust = -0.5)) + theme_light() + theme(axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 10), text = element_text(size = 15)) 
    show(plot)
  })

  }

# Run the application 
shinyApp(ui = ui, server = server)
