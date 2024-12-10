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


final_eggs <- read.csv("final_eggs.csv")
final_tent <- read.csv("final_tent.csv")
final_gld2 <- read.csv("final_gld2.csv")
final_gld4 <- read.csv("final_gld4.csv")
final_gld4_tm <- read.csv("final_gld4_tm.csv")
final_males_herm <- read.csv("final_males_herm.csv")
final_males <- read.csv("final_males.csv")
final_rnai_fndc3 <- read.csv("final_rnai_fndc3.csv")
final_fndc3 <- read.csv("final_fndc3.csv")
final_atx <- read.csv("final_atx.csv")
final_larp <- read.csv("final_larp.csv")
final_nspc <- read.csv("final_nspc.csv")


# Define UI for application that draws a histogram
ui <- fluidPage(
  # Title
  titlePanel("Nanopore data analysis"),
  # Selection of experiments and genes
  fluidRow(
    column(3,
           selectInput("experiment", "Experiment:",
                       choices=c("WT vs tent-5",
                                 "WT vs tent-5 (males)",
                                 "WT vs tent-5 (eggs)",
                                 "WT vs gld-2", 
                                 "WT vs gld-4", 
                                 "WT vs gld-4/tent-5",
                                 "WT hermaphrodites vs males",
                                 "WT vs larp-5 (RNAi)",
                                 "WT vs atx-2 (RNAi)",
                                 "WT vs C34F6.10 (RNAi)",
                                 "WT vs C34F6.10 (KO)",
                                 "WT vs nspc"), 
                       width = '400px'),
    ),
    column(2,
           textInput("gene", "Gene:", width = '250px')     
    ),
    column(5,
           selectizeInput("group", 
                          label = "Group of genes:",
                          choices = c(""), 
                          multiple = T,
                          width = '900px',
                          options = list(delimiter = " ", create = T))
    ),
    column(2, style = "margin-top: 25px;",
           actionButton("clear","Clear selection")
    )
  ),
  # White space   
  fluidRow(
    column(12, style="height:30px")
  ),
  # Main table and plot
  fluidRow(
    column(6,
           DT::dataTableOutput("table")
    ),
    column(6,
           plotOutput("MAplot"))
  ),    
  # White space
  fluidRow(
    column(12,style="height:10px")
  ),
  # Summary of selected gene
  fluidRow(
    column(12,
           titlePanel("Selected gene summary"),
           DT::dataTableOutput("comparison"))
  ),
  # White space
  fluidRow(
    column(12,style="height:10px")
  ),
  fluidRow(
    column(12,style="height:10px")
  )
)
#_________________________________________________
# Define server logic required to draw a histogram
server <- function(input, output,session) {
  
  rv <- reactiveValues(group_of_genes = NULL,list = NULL)
  
  # Main table  
  output$table <- DT::renderDataTable({
    
    if(input$experiment == "WT vs tent-5")
    {table <- final_tent[,c(2:11,14,13,1)] %>% 
      dplyr::mutate(across(c(fold_change,cohen_d), round, digits = 3)) %>% 
      dplyr::mutate(across(c(p.value:padj), signif, digits = 2))}
    else if(input$experiment == "WT vs tent-5 (males)")
    {table <- final_males[,c(2:11,14,13,1)] %>% 
      dplyr::mutate(across(c(fold_change,cohen_d), round, digits = 3)) %>% 
      dplyr::mutate(across(c(p.value:padj), signif, digits = 2))}
    else if(input$experiment == "WT vs tent-5 (eggs)")
    {table <- final_eggs[,c(2:11,14,13,1)] %>% 
      dplyr::mutate(across(c(fold_change,cohen_d), round, digits = 3)) %>% 
      dplyr::mutate(across(c(p.value:padj), signif, digits = 2))}
    else if(input$experiment == "WT vs gld-2")
    {table <- final_gld2[,c(2:11,14,13,1)] %>% 
      dplyr::mutate(across(c(fold_change,cohen_d), round, digits = 3)) %>% 
      dplyr::mutate(across(c(p.value:padj), signif, digits = 2))}
    else if(input$experiment == "WT vs gld-4")
    {table <- final_gld4[,c(2:11,14,13,1)] %>% 
      dplyr::mutate(across(c(fold_change,cohen_d), round, digits = 3)) %>% 
      dplyr::mutate(across(c(p.value:padj), signif, digits = 2))}
    else if(input$experiment == "WT vs gld-4/tent-5")
    {table <- final_gld4_tm[,c(2:11,14,13,1)] %>% 
      dplyr::mutate(across(c(fold_change,cohen_d), round, digits = 3)) %>% 
      dplyr::mutate(across(c(p.value:padj), signif, digits = 2))}
    else if(input$experiment == "WT hermaphrodites vs males")
    {table <- final_males_herm[,c(2:11,14,13,1)] %>% 
      dplyr::mutate(across(c(fold_change,cohen_d), round, digits = 3)) %>% 
      dplyr::mutate(across(c(p.value:padj), signif, digits = 2))}
    else if(input$experiment == "WT vs larp-5 (RNAi)")
    {table <- final_larp[,c(2:11,14,13,1)] %>% 
      dplyr::mutate(across(c(fold_change,cohen_d), round, digits = 3)) %>% 
      dplyr::mutate(across(c(p.value:padj), signif, digits = 2))}
    else if(input$experiment == "WT vs atx-2 (RNAi)")
    {table <- final_atx[,c(2:11,14,13,1)] %>% 
      dplyr::mutate(across(c(fold_change,cohen_d), round, digits = 3)) %>% 
      dplyr::mutate(across(c(p.value:padj), signif, digits = 2))}
    else if(input$experiment == "WT vs C34F6.10 (RNAi)")
    {table <- final_rnai_fndc3[,c(2:11,14,13,1)] %>% 
      dplyr::mutate(across(c(fold_change,cohen_d), round, digits = 3)) %>% 
      dplyr::mutate(across(c(p.value:padj), signif, digits = 2))}
    else if(input$experiment == "WT vs C34F6.10 (KO)")
    {table <- final_fndc3[,c(2:11,14,13,1)] %>% 
      dplyr::mutate(across(c(fold_change,cohen_d), round, digits = 3)) %>% 
      dplyr::mutate(across(c(p.value:padj), signif, digits = 2))}
    else if(input$experiment == "WT vs nspc")
    {table <- final_nspc[,c(2:11,14,13,1)] %>% 
      dplyr::mutate(across(c(fold_change,cohen_d), round, digits = 3)) %>% 
      dplyr::mutate(across(c(p.value:padj), signif, digits = 2))}
    
    
    return(table %>% 
             DT::datatable(rownames = FALSE,
                           selection = list(mode = "single",target = 'cell'),
                           extensions = c('Scroller'), #"Buttons','Select'
                           options = list(deferRender = TRUE,
                                          scrollY = 350,
                                          scroller = TRUE,
                                          dom = 'Blfrtip',
                                          rowId = 0,
                                          # select = list(style = 'single', items = 'cell'),
                                          # buttons = c('selectNone'), 
                                          scrollX = TRUE,
                                          columnDefs = list(list(className = 'dt-center', targets = "_all")))) %>%
             DT::formatStyle(columns = c(1:13), fontSize = '80%'))
  })
  
  # Main plot  
  output$MAplot <- renderPlot({
    
    gene_of_choice <- input$gene
    
    if(input$experiment == "WT vs tent-5")
    {plot <- final_tent %>% dplyr::arrange(significance)}
    else if(input$experiment == "WT vs tent-5 (males)")
    {plot <- final_males %>% dplyr::arrange(significance)} 
    else if(input$experiment == "WT vs tent-5 (eggs)")
    {plot <- final_eggs %>% dplyr::arrange(significance)}
    else if(input$experiment == "WT vs gld-2")
    {plot <- final_gld2 %>% dplyr::arrange(significance)}
    else if(input$experiment == "WT vs gld-4")
    {plot <- final_gld4 %>% dplyr::arrange(significance)}
    else if(input$experiment == "WT vs gld-4/tent-5")
    {plot <- final_gld4_tm %>% dplyr::arrange(significance)}
    else if(input$experiment == "WT hermaphrodites vs males")
    {plot <- final_males_herm %>% dplyr::arrange(significance)}
    else if(input$experiment == "WT vs larp-5 (RNAi)")
    {plot <- final_larp %>% dplyr::arrange(significance)}
    else if(input$experiment == "WT vs atx-2 (RNAi)")
    {plot <- final_atx %>% dplyr::arrange(significance)}
    else if(input$experiment == "WT vs C34F6.10 (RNAi)")
    {plot <- final_rnai_fndc3 %>% dplyr::arrange(significance)}
    else if(input$experiment == "WT vs C34F6.10 (KO)")
    {plot <- final_fndc3 %>% dplyr::arrange(significance)}
    else if(input$experiment == "WT vs nspc")
    {plot <- final_nspc %>% dplyr::arrange(significance)}  
    
    if(is.null(input$group))
    {rv$group_of_genes <- rv$list}
    else
      rv$group_of_genes <- input$group
    
    group <- plot %>% 
      filter(symbol %in% rv$group_of_genes)
    
    if(is.null(rv$group_of_genes) & input$gene == "")
    {return(ggplot(plot, aes(x=baseMean, y=log2(fold_change) , color=significance)) 
            + ggtitle(input$experiment) 
            + geom_point(size = 4 ,alpha = 0.3) 
            + scale_color_manual(values=c("gray70", "red3")) 
            + theme_light()
            + ylim(-3,3)
            + scale_x_log10())}
    else if(is.null(rv$group_of_genes) & input$gene != "")
    {return(ggplot(plot, aes(x=baseMean, y=log2(fold_change), color=significance)) 
            + ggtitle(input$experiment) 
            + geom_point(size = 4, alpha = 0.3) 
            + scale_color_manual(values=c("gray70", "red3")) 
            + theme_light() 
            + geom_point(data=plot %>% filter(grepl(gene_of_choice,symbol)), size =3, color="black", shape = 21) 
            + geom_label_repel(data=plot %>% filter(grepl(gene_of_choice,symbol)), aes(label=symbol, size=NULL, color=NULL), size=5, nudge_y = 0, segment.size  = 0.2, show.legend=FALSE) 
            + ylim(-3,3)
            + scale_x_log10())}
    else 
    {return(ggplot(plot, aes(x=baseMean, y= log2(fold_change) , color=significance)) 
            + ggtitle(input$experiment) + geom_point(size = 4 ,alpha = 0.3) 
            + scale_color_manual(values=c("gray70", "red3")) + theme_light() 
            + geom_point(data=group, size =4, color="black", shape = 21) 
            + geom_label_repel(data=group, aes(label=symbol, size=NULL, color=NULL), size=5, nudge_y = 0, segment.size  = 0.2, show.legend=FALSE) 
            + ylim(-3,3)
            + scale_x_log10())}
  })
  
  # Table with summary  
  output$comparison <- DT::renderDataTable({
    
    columns <- c("symbol",
                 "control_counts",
                 "control_polya", 
                 "test_counts", 
                 "test_polya",
                 "length_diff", 
                 "fold_change", 
                 "padj", 
                 "sig")
    
    if(input$gene != "")
    {group_of_genes <- input$gene}
    else
    {group_of_genes <- input$table_cell_clicked$value}
    
    summary <- data.frame(matrix(ncol = 0, nrow = 13))
    
    filtered_tent <- final_tent %>% filter(symbol %in% group_of_genes)
    filtered_tent <- filtered_tent[,c(2:8,10,11)]
    colnames(filtered_tent) <- columns
    
    filtered_males <- final_males %>% filter(symbol %in% group_of_genes)
    filtered_males <- filtered_males[,c(2:8,10,11)]
    colnames(filtered_males) <- columns
    
    filtered_eggs <- final_eggs %>% filter(symbol %in% group_of_genes)
    filtered_eggs <- filtered_eggs[,c(2:8,10,11)]
    colnames(filtered_eggs) <- columns
    
    filtered_gld2 <- final_gld2 %>% filter(symbol %in% group_of_genes)
    filtered_gld2 <- filtered_gld2[,c(2:8,10,11)]
    colnames(filtered_gld2) <- columns
    
    filtered_gld4 <- final_gld4 %>% filter(symbol %in% group_of_genes)
    filtered_gld4 <- filtered_gld4[,c(2:8,10,11)]
    colnames(filtered_gld4) <- columns
    
    filtered_gld4_tm <- final_gld4_tm %>% filter(symbol %in% group_of_genes)
    filtered_gld4_tm <- filtered_gld4_tm[,c(2:8,10,11)]
    colnames(filtered_gld4_tm) <- columns
    
    filtered_males_herm <- final_males_herm %>% filter(symbol %in% group_of_genes)
    filtered_males_herm <- filtered_males_herm[,c(2:8,10,11)]
    colnames(filtered_males_herm) <- columns
    
    filtered_larp <- final_larp %>% filter(symbol %in% group_of_genes)
    filtered_larp <- filtered_larp[,c(2:8,10,11)]
    colnames(filtered_larp) <- columns
    
    filtered_atx <- final_atx %>% filter(symbol %in% group_of_genes)
    filtered_atx <- filtered_atx[,c(2:8,10,11)]
    colnames(filtered_atx) <- columns
    
    filtered_rnai_fndc3 <- final_rnai_fndc3 %>% filter(symbol %in% group_of_genes)
    filtered_rnai_fndc3 <- filtered_rnai_fndc3[,c(2:8,10,11)]
    colnames(filtered_rnai_fndc3) <- columns
    
    filtered_fndc3 <- final_fndc3 %>% filter(symbol %in% group_of_genes)
    filtered_fndc3 <- filtered_fndc3[,c(2:8,10,11)]
    colnames(filtered_fndc3) <- columns
    
    filtered_nspc <- final_nspc %>% filter(symbol %in% group_of_genes)
    filtered_nspc <- filtered_nspc[,c(2:8,10,11)]
    colnames(filtered_nspc) <- columns
    
    
    summary <- rbind(summary, filtered_tent[1,])
    summary <- rbind(summary, filtered_males[1,])
    summary <- rbind(summary, filtered_eggs[1,])
    summary <- rbind(summary, filtered_gld2[1,])
    summary <- rbind(summary, filtered_gld4[1,])
    summary <- rbind(summary, filtered_gld4_tm[1,])
    summary <- rbind(summary, filtered_males_herm[1,])
    summary <- rbind(summary, filtered_larp[1,])
    summary <- rbind(summary, filtered_atx[1,])
    summary <- rbind(summary, filtered_rnai_fndc3[1,])
    summary <- rbind(summary, filtered_fndc3[1,])
    summary <- rbind(summary, filtered_nspc[1,])
  
    
    summary$experiment <- c("WT vs tent-5",
                            "WT vs tent-5 (males)",
                            "WT vs tent-5 (eggs)",
                            "WT vs gld-2", 
                            "WT vs gld-4", 
                            "WT vs gld-4/tent-5",
                            "WT hermaphrodites vs males",
                            "WT vs larp-5 (RNAi)",
                            "WT vs atx-2 (RNAi)",
                            "WT vs C34F6.10 (RNAi)",
                            "WT vs C34F6.10 (KO)",
                            "WT vs nspc")    
    summary <- summary[c(10,0:9)]
    
    final_summary <- summary %>% 
      dplyr::mutate(across(c(fold_change), round, digits = 3)) %>% 
      dplyr::mutate(across(c(padj), signif, digits = 2)) %>%
      DT::datatable(rownames = FALSE, selection="none",options = list(pageLength = 24,dom = 't', columnDefs = list(list(width = "100px", targets = c(1:9)),list(className = 'dt-center', targets = c(1:9))))) 
    return(final_summary)
  })
  
  # Clear selection button
  observeEvent(input$table_cell_clicked,{
    selected <- input$table_cell_clicked$value
    rv$dList <- c(isolate(rv$dList), isolate(selected))
    rv$list <- rv$dList
  })
  observeEvent(input$clear, {
    updateTextInput(session, "gene", value = "")
    updateTextInput(session, "group", value = "")
    rv$list <- NULL
    rv$dList <- NULL
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
