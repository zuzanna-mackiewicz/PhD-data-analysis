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
library(readxl)
library(ggplot2)
library(ggrepel)
library(tidyverse)


final_eggs <- read.csv("~/datasets/final_eggs.csv")
final_tent <- read.csv("~/datasets/final_tent.csv")
final_gld2 <- read.csv("~/datasets/final_gld2.csv")
final_gld4 <- read.csv("~/datasets/final_gld4.csv")
final_gld4_tm <- read.csv("~/datasets/final_gld4_tm.csv")
final_males_herm <- read.csv("~/datasets/final_males_herm.csv")
final_males <- read.csv("~/datasets/final_males.csv")
final_rnai_fndc3 <- read.csv("~/datasets/final_rnai_fndc3.csv")
final_fndc3 <- read.csv("~/datasets/final_fndc3.csv")
final_atx <- read.csv("~/datasets/final_atx.csv")
final_larp <- read.csv("~/datasets/final_larp.csv")
final_nspc <- read.csv("~/datasets/final_nspc.csv")

final_N2_NSPC <- read.csv("~/datasets/final_N2_NSPC.csv")
final_tm_NSPCtm <- read.csv("~/datasets/final_tm_NSPCtm.csv")
final_daf_NSPCdaf <- read.csv("~/datasets/final_daf_NSPCdaf.csv")
final_N2_tm <- read.csv("~/datasets/final_N2_tm.csv")
final_N2_daf <- read.csv("~/datasets/final_N2_daf.csv")
final_NSPC_NSPCtm <- read.csv("~/datasets/final_NSPC_NSPCtm.csv")
final_NSPC_NSPCdaf <- read.csv("~/datasets/final_NSPC_NSPCdaf.csv")
final_N2_NSPCtm <- read.csv("~/datasets/final_N2_NSPCtm.csv")
final_N2_NSPCdaf <- read.csv("~/datasets/final_N2_NSPCdaf.csv")
final_rnai <- read.csv("~/datasets/final_rnai.csv")

final_74_83_1 <- read.csv("~/datasets/final_74_83_1.csv")
final_74_83_2 <- read.csv("~/datasets/final_74_83_2.csv")
final_74_83_3 <- read.csv("~/datasets/inal_74_83_3.csv")
final_74_83_1_2_3 <- read.csv("~/datasets/final_74_83_1_2_3.csv")
final_N2_83_2 <- read.csv("~/datasets/final_83_N2_2.csv")
final_N2_83_3 <- read.csv("~/datasets/final_83_N2_3.csv")
final_N2_83_2_3 <- read.csv("~/datasets/final_N2_83_2_3.csv")
final_N2_74_2 <- read.csv("~/datasets/final_74_N2_2.csv")
final_N2_74_3 <- read.csv("~/datasets/final_74_N2_3.csv")
final_N2_74_2_3 <- read.csv("~/datasets/final_N2_74_2_3.csv")
final_male <- read.csv("~/datasets/final_male_herm.csv")


# Define UI for application that draws a histogram
ui <- fluidPage(
  # Title
  titlePanel("Visualisation of gene groups"),
  # Selection of experiments and genes
  fluidRow(
    column(6,
           selectInput("experiment", "Experiment:",
                       choices=c("",
                                 "Nanopore: WT vs tent-5",
                                 "Nanopore: WT vs tent-5 (males)",
                                 "Nanopore: WT vs tent-5 (eggs)",
                                 "Nanopore: WT vs gld-2", 
                                 "Nanopore: WT vs gld-4", 
                                 "Nanopore: WT vs gld-4/tent-5",
                                 "Nanopore: WT hermaphrodites vs males",
                                 "Nanopore: WT vs larp-5 (RNAi)",
                                 "Nanopore: WT vs atx-2 (RNAi)",
                                 "Nanopore: WT vs C34F6.10 (RNAi)",
                                 "Nanopore: WT vs C34F6.10 (KO)",
                                 "Nanopore: WT vs nspc",
                                 "RNA-seq: WT vs nspc",
                                 "RNA-seq: tent-5 vs nspc/tent-5",
                                 "RNA-seq: daf-16 vs nspc/daf-16",
                                 "RNA-seq: WT vs tent-5",
                                 "RNA-seq: WT vs daf-16",
                                 "RNA-seq: nspc vs nspc/tent-5",
                                 "RNA-seq: nspc vs nspc/daf-16",
                                 "RNA-seq: WT vs nspc/tent-5",
                                 "RNA-seq: WT vs nspc/daf-16",
                                 "RNA-seq: WT vs nspc (RNAi)",
                                 "RNA-seq: WT hermaphrodites vs males",
                                 "RNA-seq: inactive miniSOG vs active miniSOG (1)",
                                 "RNA-seq: inactive miniSOG vs active miniSOG (2)",
                                 "RNA-seq: inactive miniSOG vs active miniSOG (3)", 
                                 "RNA-seq: inactive miniSOG vs active miniSOG (1+2+3)", 
                                 "RNA-seq: WT vs active miniSOG (2)",
                                 "RNA-seq: WT vs active miniSOG (3)",
                                 "RNA-seq: WT vs inactive miniSOG (2)",
                                 "RNA-seq: WT vs inactive miniSOG (3)"), 
                       width = '400px',selected = NULL),
    ),
    column(6,
           fileInput("file", "Choose file with gene groups:", multiple = FALSE, accept = c(".xls", ".xlsx")),
    )
  ),
  # White space   
  fluidRow(
    column(12, style="height:30px")
  ),
  # Main table and plot
  fluidRow(
    column(6,
           plotOutput("MAplot")
    ),
    column(1, style = "margin-top: 25px",
           actionButton("previous_button","< Previous group", style="width: 150px")
    ),
    column(1, style = "margin-top: 25px",
           actionButton("next_button","Next group >", style="width: 150px")
    ),
    column(4,
           selectInput("selection", "Select group from the list:",
                       choices=""))
  ),    
  # White space
  fluidRow(
    column(12,style="height:10px")
  )
)
#_________________________________________________
# Define server logic required to draw a histogram
server <- function(input, output,session) {
  
  rv <- reactiveValues()
  
  rv$j=1
  rv$list_of_names <- list()
  rv$list <- list()
  
  observeEvent(input$file, {
    rv$df_of_groups <- read_excel(input$file$datapath)
    df_of_groups <- rv$df_of_groups
    
    for (i in 1:ncol(df_of_groups))
    {
        selected_cell <- df_of_groups[,i]
        selected_cell <- drop_na(selected_cell)
        colnames(selected_cell) <- "cell"
        group_of_genes <- selected_cell$cell
        
        if(input$experiment == "Nanopore: WT vs tent-5") {
          data_for_MA <- final_tent %>% mutate(log2FoldChange = log2(fold_change)) %>% dplyr::arrange(significance)
        } else if(input$experiment == "Nanopore: WT vs tent-5 (males)") {
          data_for_MA <- final_males %>% mutate(log2FoldChange = log2(fold_change)) %>% dplyr::arrange(significance)
        } else if(input$experiment == "Nanopore: WT vs tent-5 (eggs)") {
          data_for_MA <- final_eggs %>% mutate(log2FoldChange = log2(fold_change)) %>% dplyr::arrange(significance)
        } else if(input$experiment == "Nanopore: WT vs gld-2") {
          data_for_MA <- final_gld2 %>% mutate(log2FoldChange = log2(fold_change)) %>% dplyr::arrange(significance)
        } else if(input$experiment == "Nanopore: WT vs gld-4") {
          data_for_MA <- final_gld4 %>% mutate(log2FoldChange = log2(fold_change)) %>% dplyr::arrange(significance)
        } else if(input$experiment == "Nanopore: WT vs gld-4/tent-5") {
          data_for_MA <- final_gld4_tm %>% mutate(log2FoldChange = log2(fold_change)) %>% dplyr::arrange(significance)
        } else if(input$experiment == "Nanopore: WT hermaphrodites vs males") {
          data_for_MA <- final_males_herm %>% mutate(log2FoldChange = log2(fold_change)) %>% dplyr::arrange(significance)
        } else if(input$experiment == "Nanopore: WT vs larp-5 (RNAi)") {
          data_for_MA <- final_larp %>% mutate(log2FoldChange = log2(fold_change)) %>% dplyr::arrange(significance)
        } else if(input$experiment == "Nanopore: WT vs atx-2 (RNAi)") {
          data_for_MA <- final_atx %>% mutate(log2FoldChange = log2(fold_change)) %>% dplyr::arrange(significance)
        } else if(input$experiment == "Nanopore: WT vs C34F6.10 (RNAi)") {
          data_for_MA <- final_rnai_fndc3 %>% mutate(log2FoldChange = log2(fold_change)) %>% dplyr::arrange(significance)
        } else if(input$experiment == "Nanopore: WT vs C34F6.10 (KO)") {
          data_for_MA <- final_fndc3 %>% mutate(log2FoldChange = log2(fold_change)) %>% dplyr::arrange(significance)
        } else if(input$experiment == "Nanopore: WT vs nspc") {
          data_for_MA <- final_nspc %>% mutate(log2FoldChange = log2(fold_change)) %>% dplyr::arrange(significance)
        } else if(input$experiment == "RNA-seq: WT vs nspc") {
          data_for_MA <- final_N2_NSPC %>% mutate(log2FoldChange = log2(fold_change)) %>% dplyr::arrange(significance)
        } else if(input$experiment == "RNA-seq: tent-5 vs nspc/tent-5") {
          data_for_MA <- final_tm_NSPCtm %>% mutate(log2FoldChange = log2(fold_change)) %>% dplyr::arrange(significance)
        } else if(input$experiment == "RNA-seq: daf-16 vs nspc/daf-16") {
          data_for_MA <- final_daf_NSPCdaf %>% mutate(log2FoldChange = log2(fold_change)) %>% dplyr::arrange(significance)
        } else if(input$experiment == "RNA-seq: WT vs tent-5") {
          data_for_MA <- final_N2_tm %>% mutate(log2FoldChange = log2(fold_change)) %>% dplyr::arrange(significance)
        } else if(input$experiment == "RNA-seq: WT vs daf-16") {
          data_for_MA <- final_N2_daf %>% mutate(log2FoldChange = log2(fold_change)) %>% dplyr::arrange(significance)
        } else if(input$experiment == "RNA-seq: nspc vs nspc/tent-5") {
          data_for_MA <- final_NSPC_NSPCtm %>% mutate(log2FoldChange = log2(fold_change)) %>% dplyr::arrange(significance)
        } else if(input$experiment == "RNA-seq: nspc vs nspc/daf-16") {
          data_for_MA <- final_NSPC_NSPCdaf %>% mutate(log2FoldChange = log2(fold_change)) %>% dplyr::arrange(significance)
        } else if(input$experiment == "RNA-seq: WT vs nspc/tent-5") {
          data_for_MA <- final_N2_NSPCtm %>% mutate(log2FoldChange = log2(fold_change)) %>% dplyr::arrange(significance)
        } else if(input$experiment == "RNA-seq: WT vs nspc/daf-16") {
          data_for_MA <- final_N2_NSPCdaf %>% mutate(log2FoldChange = log2(fold_change)) %>% dplyr::arrange(significance)
        } else if(input$experiment == "RNA-seq: WT vs nspc (RNAi)") {
          data_for_MA <- final_rnai %>% mutate(log2FoldChange = log2(fold_change)) %>% dplyr::arrange(significance)
        } else if(input$experiment == "RNA-seq: WT hermaphrodites vs males") {
          data_for_MA <- final_male %>% mutate(log2FoldChange = log2(fold_change)) %>% dplyr::arrange(significance)
        } else if(input$experiment == "RNA-seq: inactive miniSOG vs active miniSOG (1)") {
          data_for_MA <- final_74_83_1 %>% mutate(log2FoldChange = log2(fold_change)) %>% dplyr::arrange(significance)
        } else if(input$experiment == "RNA-seq: inactive miniSOG vs active miniSOG (2)") {
          data_for_MA <- final_74_83_2 %>% mutate(log2FoldChange = log2(fold_change)) %>% dplyr::arrange(significance)
        } else if(input$experiment == "RNA-seq: inactive miniSOG vs active miniSOG (3)") {
          data_for_MA <- final_74_83_3 %>% mutate(log2FoldChange = log2(fold_change)) %>% dplyr::arrange(significance)
        } else if(input$experiment == "RNA-seq: inactive miniSOG vs active miniSOG (1+2+3)") {
          data_for_MA <- final_74_83_1_2_3 %>% mutate(log2FoldChange = log2(fold_change)) %>% dplyr::arrange(significance)
        } else if(input$experiment == "RNA-seq: WT vs active miniSOG (2)") {
          data_for_MA <- final_N2_83_2 %>% mutate(log2FoldChange = log2(fold_change)) %>% dplyr::arrange(significance)
        } else if(input$experiment == "RNA-seq: WT vs active miniSOG (3)") {
          data_for_MA <- final_N2_83_3 %>% mutate(log2FoldChange = log2(fold_change)) %>% dplyr::arrange(significance)
        } else if(input$experiment == "RNA-seq: WT vs inactive miniSOG (2)") {
          data_for_MA <- final_N2_74_2 %>% mutate(log2FoldChange = log2(fold_change)) %>% dplyr::arrange(significance)
        } else if(input$experiment == "RNA-seq: WT vs inactive miniSOG (3)") {
          data_for_MA <- final_N2_74_3 %>% mutate(log2FoldChange = log2(fold_change)) %>% dplyr::arrange(significance)
        }
        
        group <- data_for_MA %>% filter(symbol %in% group_of_genes)
        
        plot <- ggplot(data_for_MA, aes(x=baseMean, y=log2FoldChange, color=significance)) + ggtitle(colnames(df_of_groups)[i]) + geom_point(size = 3, alpha = 0.3) + scale_color_manual(values=c("gray70", "red")) + theme_light() + geom_point(data=group, size =3, color="black", shape = 21) + ylim(-3,3)+ scale_x_log10()
        
        rv$list[[i]] <- plot
      rv$list_of_names[[i]] <- colnames(df_of_groups)[i]
    }
    updateSelectInput(session, "selection", "Select group from the list:", choices = rv$list_of_names)
  })
  
  # Main plot  
  output$MAplot <- renderPlot({
    
    j <- rv$j
    list <- rv$list
    print(list[[j]])
    
  })
  
  observeEvent(input$selection, {
    selected_group <- input$selection
    index <- which(rv$list_of_names == selected_group)
    rv$j = index
  })
  
  observeEvent(input$next_button, {
    rv$j=rv$j+1
  })
  observeEvent(input$previous_button, {
    rv$j=rv$j-1
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
