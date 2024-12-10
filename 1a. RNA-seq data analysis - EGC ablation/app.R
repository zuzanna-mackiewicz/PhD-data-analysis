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


final_74_83_1 <- read.csv("~/datasets/final_74_83_1.csv")
final_74_83_2 <- read.csv("~/datasets/final_74_83_2.csv")
final_74_83_3 <- read.csv("~/datasets/final_74_83_3.csv")
final_74_83_2_3 <- read.csv("~/datasets/final_83_74_2_3.csv")
final_74_83_1_2_3 <- read.csv("~/datasets/final_74_83_1_2_3.csv")
final_N2_83_2 <- read.csv("~/datasets/final_83_N2_2.csv")
final_N2_83_3 <- read.csv("~/datasets/final_83_N2_3.csv")
final_N2_83_2_3 <- read.csv("~/datasets/final_N2_83_2_3.csv")
final_N2_74_2 <- read.csv("~/datasets/final_74_N2_2.csv")
final_N2_74_3 <- read.csv("~/datasets/final_74_N2_3.csv")
final_N2_74_2_3 <- read.csv("~/datasets/final_N2_74_2_3.csv")

final_male <- read.csv("~/datasets/final_male_herm.csv")
final_N2_NSPC <- read.csv("~/datasets/final_N2_NSPC.csv")
final_N2_tm <- read.csv("~/datasets/final_N2_tm.csv")


# Define UI for application that draws a histogram
ui <- fluidPage(
# Title
  titlePanel("RNA-seq data analysis"),
# Selection of experiments and genes
  fluidRow(
    column(3,
      selectInput("experiment", "Experiment:",
                   choices=c("1: inactive miniSOG vs active miniSOG",
                             "2: inactive miniSOG vs active miniSOG",
                             "3: inactive miniSOG vs active miniSOG", 
                             "2+3: inactive miniSOG vs active miniSOG", 
                             "1+2+3: inactive miniSOG vs active miniSOG", 
                             "2: WT vs active miniSOG",
                             "3: WT vs active miniSOG",
                             "2+3: WT vs active miniSOG",
                             "2: WT vs inactive miniSOG",
                             "3: WT vs inactive miniSOG",
                             "2+3: WT vs inactive miniSOG",
                             "WT hermaphrodites vs males",
                             "WT vs nspc",
                             "WT vs tent-5"), 
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

    if(input$experiment == "1: inactive miniSOG vs active miniSOG")
    {table <- final_74_83_1[,c(9,2:3,6:8,1,13:15,10:12)] %>% 
      dplyr::mutate(across(c("ADZ74.1.1":"ADZ83.3.1"), round, digits = 0)) %>% 
      dplyr::mutate(across(c(baseMean:log2FoldChange), round, digits = 3)) %>% 
      dplyr::mutate(across(c(pvalue:padj), signif, digits = 2))}
    else if(input$experiment == "2: inactive miniSOG vs active miniSOG")
    {table <- final_74_83_2[,c(9,2:3,6:8,1,10:15)] %>% 
      dplyr::mutate(across(c("ADZ74.1.2":"ADZ83.3.2"), round, digits = 0)) %>% 
      dplyr::mutate(across(c(baseMean:log2FoldChange), round, digits = 3)) %>% 
      dplyr::mutate(across(c(pvalue:padj), signif, digits = 2))}
    else if(input$experiment == "3: inactive miniSOG vs active miniSOG")
    {table <- final_74_83_3[,c(9,2:3,6:8,1,10:15)] %>% 
      dplyr::mutate(across(c("ADZ74.1.3":"ADZ83.3.3"), round, digits = 0)) %>% 
      dplyr::mutate(across(c(baseMean:log2FoldChange), round, digits = 3)) %>% 
      dplyr::mutate(across(c(pvalue:padj), signif, digits = 2))}
    else if(input$experiment == "2+3: inactive miniSOG vs active miniSOG")
    {table <- final_74_83_2_3[,c(9,2:3,6:8,1,10:12,16:18,13:15,19:21)] %>% 
      dplyr::mutate(across(c("ADZ74.1.3":"ADZ83.3.2"), round, digits = 0)) %>% 
      dplyr::mutate(across(c(baseMean:log2FoldChange), round, digits = 3)) %>% 
      dplyr::mutate(across(c(pvalue:padj), signif, digits = 2))}
    else if(input$experiment == "1+2+3: inactive miniSOG vs active miniSOG")
    {table <- final_74_83_1_2_3[,c(9,2:3,6:8,1,10:12,16:18,25:27,13:15,19:21,22:24)] %>% 
      dplyr::mutate(across(c("ADZ74.1.3":"ADZ83.3.1"), round, digits = 0)) %>% 
      dplyr::mutate(across(c(baseMean:log2FoldChange), round, digits = 3)) %>% 
      dplyr::mutate(across(c(pvalue:padj), signif, digits = 2))}
    else if(input$experiment == "2: WT vs active miniSOG")
    {table <- final_N2_83_2[,c(9,2:3,6:8,1,13:15,10:12)] %>% 
      dplyr::mutate(across(c("N2.1.2":"ADZ83.3.2"), round, digits = 0)) %>% 
      dplyr::mutate(across(c(baseMean:log2FoldChange), round, digits = 3)) %>% 
      dplyr::mutate(across(c(pvalue:padj), signif, digits = 2))}
    else if(input$experiment == "3: WT vs active miniSOG")
    {table <- final_N2_83_3[,c(9,2:3,6:8,1,13:15,10:12)] %>% 
      dplyr::mutate(across(c("N2.1.3":"ADZ83.3.3"), round, digits = 0)) %>% 
      dplyr::mutate(across(c(baseMean:log2FoldChange), round, digits = 3)) %>% 
      dplyr::mutate(across(c(pvalue:padj), signif, digits = 2))}
    else if(input$experiment == "2+3: WT vs active miniSOG")
    {table <- final_N2_83_2_3[,c(9,2:3,6:8,1,13:15,19:21,10:12,16:18)] %>% 
      dplyr::mutate(across(c("N2.1.3":"ADZ83.3.2"), round, digits = 0)) %>% 
      dplyr::mutate(across(c(baseMean:log2FoldChange), round, digits = 3)) %>% 
      dplyr::mutate(across(c(pvalue:padj), signif, digits = 2))}
    else if(input$experiment == "2: WT vs inactive miniSOG")
    {table <- final_N2_74_2[,c(9,2:3,6:8,1,13:15,10:12)] %>% 
      dplyr::mutate(across(c("N2.1.2":"ADZ74.3.2"), round, digits = 0)) %>% 
      dplyr::mutate(across(c(baseMean:log2FoldChange), round, digits = 3)) %>% 
      dplyr::mutate(across(c(pvalue:padj), signif, digits = 2))}
    else if(input$experiment == "3: WT vs inactive miniSOG")
    {table <- final_N2_74_3[,c(9,2:3,6:8,1,13:15,10:12)] %>% 
      dplyr::mutate(across(c("N2.1.3":"ADZ74.3.3"), round, digits = 0)) %>% 
      dplyr::mutate(across(c(baseMean:log2FoldChange), round, digits = 3)) %>% 
      dplyr::mutate(across(c(pvalue:padj), signif, digits = 2))}
     else if(input$experiment == "2+3: WT vs inactive miniSOG")
    {table <- final_N2_74_2_3[,c(9,2:3,6:8,1,13:15,19:21,10:12,16:18)] %>% 
      dplyr::mutate(across(c("N2.1.3":"ADZ74.3.2"), round, digits = 0)) %>% 
      dplyr::mutate(across(c(baseMean:log2FoldChange), round, digits = 3)) %>% 
      dplyr::mutate(across(c(pvalue:padj), signif, digits = 2))}
    else if(input$experiment == "WT hermaphrodites vs males")
    {table <- final_male[,c(9,2,3,6:8,10,13,14,11,12)] %>% 
      dplyr::mutate(across(c(C_elegans_N2_biorep_1:C_elegans_males_N2_biorep_3), round, digits = 0)) %>% 
      dplyr::mutate(across(c(baseMean:log2FoldChange), round, digits = 3)) %>% 
      dplyr::mutate(across(c(pvalue:padj), signif, digits = 2))}
    else if(input$experiment == "WT vs nspc")
    {table <- final_N2_NSPC[,c(9,2:3,6:8,1,10:15)] %>% 
      dplyr::mutate(across(c("N2.1":"NSPC.3"), round, digits = 0)) %>% 
      dplyr::mutate(across(c(baseMean:log2FoldChange), round, digits = 3)) %>% 
      dplyr::mutate(across(c(pvalue:padj), signif, digits = 2))}
    else if(input$experiment == "WT vs tent-5")
    {table <- final_N2_tm[,c(9,2:3,6:8,1,10:15)] %>% 
      dplyr::mutate(across(c("N2.1":"tm.3"), round, digits = 0)) %>% 
      dplyr::mutate(across(c(baseMean:log2FoldChange), round, digits = 3)) %>% 
      dplyr::mutate(across(c(pvalue:padj), signif, digits = 2))}
  
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
    
    if(input$experiment ==  "1: inactive miniSOG vs active miniSOG")
    {plot <- final_74_83_1 %>% dplyr::arrange(sig)}
    else if(input$experiment == "2: inactive miniSOG vs active miniSOG")
    {plot <- final_74_83_2 %>% dplyr::arrange(sig)} 
    else if(input$experiment == "3: inactive miniSOG vs active miniSOG")
    {plot <- final_74_83_3 %>% dplyr::arrange(sig)}
    else if(input$experiment ==  "2+3: inactive miniSOG vs active miniSOG")
    {plot <- final_74_83_2_3 %>% dplyr::arrange(sig)}
    else if(input$experiment ==  "1+2+3: inactive miniSOG vs active miniSOG")
    {plot <- final_74_83_1_2_3 %>% dplyr::arrange(sig)}
    else if(input$experiment == "2: WT vs active miniSOG")
    {plot <- final_N2_83_2 %>% dplyr::arrange(sig)}
    else if(input$experiment ==  "3: WT vs active miniSOG")
    {plot <- final_N2_83_3 %>% dplyr::arrange(sig)}
    else if(input$experiment == "2+3: WT vs active miniSOG")
    {plot <- final_N2_83_2_3 %>% dplyr::arrange(sig)}
    else if(input$experiment == "2: WT vs inactive miniSOG")
    {plot <- final_N2_74_2 %>% dplyr::arrange(sig)}
    else if(input$experiment == "3: WT vs inactive miniSOG")
    {plot <- final_N2_74_3 %>% dplyr::arrange(sig)}
    else if(input$experiment == "2+3: WT vs inactive miniSOG")
    {plot <- final_N2_74_2_3 %>% dplyr::arrange(sig)}
    else if(input$experiment == "WT hermaphrodites vs males")
    {plot <- final_male %>% dplyr::arrange(sig)}
    else if(input$experiment == "WT vs nspc")
    {plot <- final_N2_NSPC %>% dplyr::arrange(sig)}
    else if(input$experiment == "WT vs tent-5")
    {plot <- final_N2_tm %>% dplyr::arrange(sig)}
  
    if(is.null(input$group))
    {rv$group_of_genes <- rv$list}
    else
    rv$group_of_genes <- input$group

    group <- plot %>% 
      filter(symbol %in% rv$group_of_genes)

    if(is.null(rv$group_of_genes) & input$gene == "")
    {return(ggplot(plot, aes(x=baseMean, y= log2FoldChange, color=sig)) 
            + ggtitle(input$experiment) 
            + geom_point(size = 4 ,alpha = 0.3) 
            + scale_color_manual(values=c("gray80", "red3")) 
            + theme_light() 
            + scale_x_log10())}
    else if(is.null(rv$group_of_genes) & input$gene != "")
    {return(ggplot(plot, aes(x=baseMean, y= log2FoldChange, color=sig)) 
            + ggtitle(input$experiment) 
            + geom_point(size = 4, alpha = 0.3) 
            + scale_color_manual(values=c("gray80", "red3")) 
            + theme_light() 
            + geom_point(data=plot %>% filter(grepl(gene_of_choice,symbol)), size =4, color="black", shape = 21) 
            + geom_label_repel(data=plot %>% filter(grepl(gene_of_choice,symbol)), aes(label=symbol, size=NULL, color=NULL), size=5, nudge_y = 0, segment.size  = 0.2, show.legend=FALSE) 
            + scale_x_log10())}
    else 
    {return(ggplot(plot, aes(x=baseMean, y= log2FoldChange, color=sig)) 
           + ggtitle(input$experiment) + geom_point(size = 4 ,alpha = 0.3) 
           + scale_color_manual(values=c("gray80", "red3")) + theme_light() 
           + geom_point(data=group, size =4, color="black", shape = 21) 
           + geom_label_repel(data=group, aes(label=symbol, size=NULL, color=NULL), size=5, nudge_y = 0, segment.size  = 0.2, show.legend=FALSE) 
           + scale_x_log10())}
  })
  
# Table with summary  
  output$comparison <- DT::renderDataTable({
      
    columns <- c("symbol",
                 "baseMean",
                 "log2FoldChange", 
                 "padj", 
                 "sig",
                 "control_1", 
                 "control_2", 
                 "control_3", 
                 "test_1", 
                 "test_2", 
                 "test_3")
      
    if(input$gene != "")
    {group_of_genes <- input$gene}
    else
    {group_of_genes <- input$table_cell_clicked$value}

    summary <- data.frame(matrix(ncol = 0, nrow = 23))
        
    filtered_74_83_1 <- final_74_83_1 %>% filter(symbol %in% group_of_genes)
    filtered_74_83_1 <- filtered_74_83_1[,c(9,2,3,7,8,13:15,10:12)]
    colnames(filtered_74_83_1) <- columns
    
    filtered_74_83_2 <- final_74_83_2 %>% filter(symbol %in% group_of_genes)
    filtered_74_83_2 <- filtered_74_83_2[,c(9,2,3,7,8,10:15)]
    colnames(filtered_74_83_2) <- columns
    
    filtered_74_83_3 <- final_74_83_3 %>% filter(symbol %in% group_of_genes)
    filtered_74_83_3 <- filtered_74_83_3[,c(9,2,3,7,8,10:15)]
    colnames(filtered_74_83_3) <- columns
    
    filtered_N2_83_2 <- final_N2_83_2 %>% filter(symbol %in% group_of_genes)
    filtered_N2_83_2 <- filtered_N2_83_2[,c(9,2,3,7,8,13:15,10:12)]
    colnames(filtered_N2_83_2) <- columns
    
    filtered_N2_83_3 <- final_N2_83_3 %>% filter(symbol %in% group_of_genes)
    filtered_N2_83_3 <- filtered_N2_83_3[,c(9,2,3,7,8,13:15,10:12)]
    colnames(filtered_N2_83_3) <- columns
    
    filtered_N2_74_2 <- final_N2_74_2 %>% filter(symbol %in% group_of_genes)
    filtered_N2_74_2 <- filtered_N2_74_2[,c(9,2,3,7,8,13:15,10:12)]
    colnames(filtered_N2_74_2) <- columns
    
    filtered_N2_74_3 <- final_N2_74_3 %>% filter(symbol %in% group_of_genes)
    filtered_N2_74_3 <- filtered_N2_74_3[,c(9,2,3,7,8,13:15,10:12)]
    colnames(filtered_N2_74_3) <- columns
    
    filtered_male <- final_male %>% filter(symbol %in% group_of_genes)
    filtered_male <- filtered_male[,c(9,2,3,7,8,13,14,17,11,12,16)]
    colnames(filtered_male) <- columns
    
    filtered_N2_NSPC <- final_N2_NSPC %>% filter(symbol %in% group_of_genes)
    filtered_N2_NSPC <- filtered_N2_NSPC[,c(9,2,3,7,8,10:15)]
    colnames(filtered_N2_NSPC) <- columns
    
    filtered_N2_tm <- final_N2_tm %>% filter(symbol %in% group_of_genes)
    filtered_N2_tm <- filtered_N2_tm[,c(9,2,3,7,8,10:15)]
    colnames(filtered_N2_tm) <- columns
        
    summary <- rbind(summary, filtered_74_83_1[1,])
    summary <- rbind(summary, filtered_74_83_2[1,])
    summary <- rbind(summary, filtered_74_83_3[1,])
    summary <- rbind(summary, filtered_N2_83_2[1,])
    summary <- rbind(summary, filtered_N2_83_3[1,])
    summary <- rbind(summary, filtered_N2_74_2[1,])
    summary <- rbind(summary, filtered_N2_74_3[1,])
    summary <- rbind(summary, filtered_male[1,])
    summary <- rbind(summary, filtered_N2_NSPC[1,])
    summary <- rbind(summary, filtered_N2_tm[1,])

    
    summary$experiment <- c("1: inactive miniSOG vs active miniSOG",
                            "2: inactive miniSOG vs active miniSOG",
                            "3: inactive miniSOG vs active miniSOG", 
                            "2: WT vs active miniSOG",
                            "3: WT vs active miniSOG",
                            "2: WT vs inactive miniSOG",
                            "3: WT vs inactive miniSOG",
                            "WT hermaphrodites vs males",
                            "WT vs nspc",
                            "WT vs tent-5")    
    summary <- summary[,c(12,1:11)]
        
    final_summary <- summary %>% 
      dplyr::mutate(across(c(control_1:test_3), round, digits = 0)) %>% 
      dplyr::mutate(across(c(baseMean:log2FoldChange), round, digits = 3)) %>% 
      dplyr::mutate(across(c(padj), signif, digits = 2)) %>% 
      DT::datatable(rownames = FALSE, selection="none",options = list(pageLength = 24, dom = 't', columnDefs = list(list(width = "100px", targets = c(1:11)),list(className = 'dt-center', targets = c(1:11))))) 
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
