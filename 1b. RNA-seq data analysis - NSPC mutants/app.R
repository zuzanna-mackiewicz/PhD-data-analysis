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

final_N2_NSPC <- read.csv("../datasets/final_N2_NSPC.csv")
final_tm_NSPCtm <- read.csv("../datasets/final_tm_NSPCtm.csv")
final_daf_NSPCdaf <- read.csv("../datasets/final_daf_NSPCdaf.csv")
final_N2_tm <- read.csv("../datasets/final_N2_tm.csv")
final_N2_daf <- read.csv("../datasets/final_N2_daf.csv")
final_NSPC_NSPCtm <- read.csv("../datasets/final_NSPC_NSPCtm.csv")
final_NSPC_NSPCdaf <- read.csv("../datasets/final_NSPC_NSPCdaf.csv")
final_N2_NSPCtm <- read.csv("../datasets/final_N2_NSPCtm.csv")
final_N2_NSPCdaf <- read.csv("../datasets/final_N2_NSPCdaf.csv")

final_rnai <- read.csv("../datasets/final_rnai.csv")
final_male <- read.csv("../datasets/final_male_herm.csv")


# Define UI for application that draws a histogram
ui <- fluidPage(
# Title
  titlePanel("RNA-seq data analysis"),
# Selection of experiments and genes
  fluidRow(
    column(3,
      selectInput("experiment", "Experiment:",
                   choices=c("WT vs nspc",
                             "tent-5 vs nspc/tent-5",
                             "daf-16 vs nspc/daf-16",
                             "WT vs tent-5",
                             "WT vs daf-16",
                             "nspc vs nspc/tent-5",
                             "nspc vs nspc/daf-16",
                             "WT vs nspc/tent-5",
                             "WT vs nspc/daf-16",
                             "RNAi: WT vs nspc",
                             "WT hermaphrodites vs males"), 
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
    
    if(input$experiment == "WT vs nspc")
    {table <- final_N2_NSPC[,c(9,2:3,6:8,1,10:15)] %>% 
      dplyr::mutate(across(c("N2.1":"NSPC.3"), round, digits = 0)) %>% 
      dplyr::mutate(across(c(baseMean:log2FoldChange), round, digits = 3)) %>% 
      dplyr::mutate(across(c(pvalue:padj), ~ format(signif(., digits = 2))))}
    else if(input$experiment == "tent-5 vs nspc/tent-5")
    {table <- final_tm_NSPCtm[,c(9,2:3,6:8,1,10:15)] %>% 
      dplyr::mutate(across(c("tm.1":"NSPC.tm.3"), round, digits = 0)) %>% 
      dplyr::mutate(across(c(baseMean:log2FoldChange), round, digits = 3)) %>% 
      dplyr::mutate(across(c(pvalue:padj), ~ format(signif(., digits = 2))))}
    else if(input$experiment == "daf-16 vs nspc/daf-16")
    {table <- final_daf_NSPCdaf[,c(9,2:3,6:8,1,10:15)] %>% 
      dplyr::mutate(across(c("daf.1":"NSPC.daf.3"), round, digits = 0)) %>% 
      dplyr::mutate(across(c(baseMean:log2FoldChange), round, digits = 3)) %>% 
      dplyr::mutate(across(c(pvalue:padj), ~ format(signif(., digits = 2))))}
    else if(input$experiment == "WT vs tent-5")
    {table <- final_N2_tm[,c(9,2:3,6:8,1,10:15)] %>% 
      dplyr::mutate(across(c("N2_1":"tm_3"), round, digits = 0)) %>% 
      dplyr::mutate(across(c(baseMean:log2FoldChange), round, digits = 3)) %>% 
      dplyr::mutate(across(c(pvalue:padj), ~ format(signif(., digits = 2))))}
    else if(input$experiment == "WT vs daf-16")
    {table <- final_N2_daf[,c(9,2:3,6:8,1,10:15)] %>% 
      dplyr::mutate(across(c("N2_1":"daf_3"), round, digits = 0)) %>% 
      dplyr::mutate(across(c(baseMean:log2FoldChange), round, digits = 3)) %>% 
      dplyr::mutate(across(c(pvalue:padj), ~ format(signif(., digits = 2))))}
    else if(input$experiment == "nspc vs nspc/tent-5")
    {table <- final_NSPC_NSPCtm[,c(9,2:3,6:8,1,10:15)] %>% 
      dplyr::mutate(across(c("NSPC_1":"NSPC_tm_3"), round, digits = 0)) %>% 
      dplyr::mutate(across(c(baseMean:log2FoldChange), round, digits = 3)) %>% 
      dplyr::mutate(across(c(pvalue:padj), ~ format(signif(., digits = 2))))}
    else if(input$experiment == "nspc vs nspc/daf-16")
    {table <- final_NSPC_NSPCdaf[,c(9,2:3,6:8,1,10:15)] %>% 
      dplyr::mutate(across(c("NSPC.1":"NSPC.daf.3"), round, digits = 0)) %>% 
      dplyr::mutate(across(c(baseMean:log2FoldChange), round, digits = 3)) %>% 
      dplyr::mutate(across(c(pvalue:padj), ~ format(signif(., digits = 2))))}
    else if(input$experiment == "WT vs nspc/tent-5")
    {table <- final_N2_NSPCtm[,c(9,2:3,6:8,1,10:15)] %>% 
      dplyr::mutate(across(c("N2.1":"NSPC.tm.3"), round, digits = 0)) %>% 
      dplyr::mutate(across(c(baseMean:log2FoldChange), round, digits = 3)) %>% 
      dplyr::mutate(across(c(pvalue:padj), ~ format(signif(., digits = 2))))}
    else if(input$experiment == "WT vs nspc/daf-16")
    {table <- final_N2_NSPCdaf[,c(9,2:3,6:8,1,10:15)] %>% 
      dplyr::mutate(across(c("N2.1":"NSPC.daf.3"), round, digits = 0)) %>% 
      dplyr::mutate(across(c(baseMean:log2FoldChange), round, digits = 3)) %>% 
      dplyr::mutate(across(c(pvalue:padj), ~ format(signif(., digits = 2))))}
    else if(input$experiment == "RNAi: WT vs nspc")
    {table <- final_rnai[,c(9,2:3,6:8,1,10:15)] %>% 
      dplyr::mutate(across(c(empty_1:NSPC_3), round, digits = 0)) %>% 
      dplyr::mutate(across(c(baseMean:log2FoldChange), round, digits = 3)) %>% 
      dplyr::mutate(across(c(pvalue:padj), ~ format(signif(., digits = 2))))}
    else if(input$experiment == "WT hermaphrodites vs males")
    {table <- final_male[,c(9,2,3,6:8,10,13,14,11,12)] %>% 
      dplyr::mutate(across(c(C_elegans_N2_biorep_1:C_elegans_males_N2_biorep_2), round, digits = 0)) %>% 
      dplyr::mutate(across(c(baseMean:log2FoldChange), round, digits = 3)) %>% 
      dplyr::mutate(across(c(pvalue:padj), ~ format(signif(., digits = 2))))}
    
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
                                   # buttons = c('excel'),
                                   scrollX = TRUE,
                                   columnDefs = list(list(className = 'dt-center', targets = "_all")))) %>%
      DT::formatStyle(columns = c(1:13), fontSize = '80%'))
  })
  
# Main plot  
  output$MAplot <- renderPlot({
    
    gene_of_choice <- input$gene
    
    if(input$experiment == "WT vs nspc")
    {plot <- final_N2_NSPC %>% dplyr::arrange(sig)}
    else if(input$experiment == "tent-5 vs nspc/tent-5")
    {plot <- final_tm_NSPCtm %>% dplyr::arrange(sig)}
    else if(input$experiment == "daf-16 vs nspc/daf-16")
    {plot <- final_daf_NSPCdaf %>% dplyr::arrange(sig)}
    else if(input$experiment == "WT vs tent-5")
    {plot <- final_N2_tm %>% dplyr::arrange(sig)}
    else if(input$experiment == "WT vs daf-16")
    {plot <- final_N2_daf %>% dplyr::arrange(sig)}
    else if(input$experiment == "nspc vs nspc/tent-5")
    {plot <- final_NSPC_NSPCtm %>% dplyr::arrange(sig)}
    else if(input$experiment == "nspc vs nspc/daf-16")
    {plot <- final_NSPC_NSPCdaf %>% dplyr::arrange(sig)}
    else if(input$experiment == "WT vs nspc/tent-5")
    {plot <- final_N2_NSPCtm %>% dplyr::arrange(sig)}
    else if(input$experiment == "WT vs nspc/daf-16")
    {plot <- final_N2_NSPCdaf %>% dplyr::arrange(sig)} 
    else if(input$experiment == "RNAi: WT vs nspc")
    {plot <- final_rnai %>% dplyr::arrange(sig)} 
    else if(input$experiment == "WT hermaphrodites vs males")
    {plot <- final_male %>% dplyr::arrange(sig)}


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
            + scale_color_manual(values=c("gray70", "red3")) 
            + theme_light() 
            + scale_x_log10())}
    else if(is.null(rv$group_of_genes) & input$gene != "")
    {return(ggplot(plot, aes(x=baseMean, y= log2FoldChange, color=sig)) 
            + ggtitle(input$experiment) 
            + geom_point(size = 4, alpha = 0.3) 
            + scale_color_manual(values=c("gray70", "red3")) 
            + theme_light() 
            + geom_point(data=plot %>% filter(grepl(gene_of_choice,symbol)), size =3, color="black", shape = 21) 
            + geom_label_repel(data=plot %>% filter(grepl(gene_of_choice,symbol)), aes(label=symbol, size=NULL, color=NULL), size=5, nudge_y = 0, segment.size  = 0.2, show.legend=FALSE) 
            + scale_x_log10())}
    else 
    {return(ggplot(plot, aes(x=baseMean, y= log2FoldChange, color=sig)) 
           + ggtitle(input$experiment) + geom_point(size = 4 ,alpha = 0.3) 
           + scale_color_manual(values=c("gray70", "red3")) + theme_light() 
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
    
    filtered_N2_NSPC <- final_N2_NSPC %>% filter(symbol %in% group_of_genes)
    filtered_N2_NSPC <- filtered_N2_NSPC[,c(9,2,3,7,8,10:15)]
    colnames(filtered_N2_NSPC) <- columns
    
    filtered_tm_NSPCtm <- final_tm_NSPCtm %>% filter(symbol %in% group_of_genes)
    filtered_tm_NSPCtm <- filtered_tm_NSPCtm[,c(9,2,3,7,8,10:15)]
    colnames(filtered_tm_NSPCtm) <- columns
    
    filtered_daf_NSPCdaf <- final_daf_NSPCdaf %>% filter(symbol %in% group_of_genes)
    filtered_daf_NSPCdaf <- filtered_daf_NSPCdaf[,c(9,2,3,7,8,10:15)]
    colnames(filtered_daf_NSPCdaf) <- columns
    
    filtered_N2_tm <- final_N2_tm %>% filter(symbol %in% group_of_genes)
    filtered_N2_tm <- filtered_N2_tm[,c(9,2,3,7,8,10:15)]
    colnames(filtered_N2_tm) <- columns
    
    filtered_N2_daf <- final_N2_daf %>% filter(symbol %in% group_of_genes)
    filtered_N2_daf <- filtered_N2_daf[,c(9,2:3,7:8,10:15)]
    colnames(filtered_N2_daf) <- columns
    
    filtered_NSPC_NSPCtm <- final_NSPC_NSPCtm %>% filter(symbol %in% group_of_genes)
    filtered_NSPC_NSPCtm <- filtered_NSPC_NSPCtm[,c(9,2:3,7:8,10:15)]
    colnames(filtered_NSPC_NSPCtm) <- columns
    
    filtered_NSPC_NSPCdaf <- final_NSPC_NSPCdaf %>% filter(symbol %in% group_of_genes)
    filtered_NSPC_NSPCdaf <- filtered_NSPC_NSPCdaf[,c(9,2:3,7:8,10:15)]
    colnames(filtered_NSPC_NSPCdaf) <- columns
    
    filtered_N2_NSPCtm <- final_N2_NSPCtm %>% filter(symbol %in% group_of_genes)
    filtered_N2_NSPCtm <- filtered_N2_NSPCtm[,c(9,2:3,7:8,10:15)]
    colnames(filtered_N2_NSPCtm) <- columns
    
    filtered_N2_NSPCdaf <- final_N2_NSPCdaf %>% filter(symbol %in% group_of_genes)
    filtered_N2_NSPCdaf <- filtered_N2_NSPCdaf[,c(9,2:3,7:8,10:15)]
    colnames(filtered_N2_NSPCdaf) <- columns
    
    filtered_rnai <- final_rnai %>% filter(symbol %in% group_of_genes)
    filtered_rnai <- filtered_rnai[,c(9,2:3,7:8,10:15)]
    colnames(filtered_rnai) <- columns
    
    filtered_male <- final_male %>% filter(symbol %in% group_of_genes)
    filtered_male <- filtered_male[,c(9,2,3,7,8,13,14,17,11,12,16)]
    colnames(filtered_male) <- columns
        
    summary <- rbind(summary, filtered_N2_NSPC[1,])
    summary <- rbind(summary, filtered_tm_NSPCtm[1,])
    summary <- rbind(summary, filtered_daf_NSPCdaf[1,])
    summary <- rbind(summary, filtered_N2_tm[1,])
    summary <- rbind(summary, filtered_N2_daf[1,])
    summary <- rbind(summary, filtered_NSPC_NSPCtm[1,])
    summary <- rbind(summary, filtered_NSPC_NSPCdaf[1,])
    summary <- rbind(summary, filtered_N2_NSPCtm[1,])
    summary <- rbind(summary, filtered_N2_NSPCdaf[1,])
    summary <- rbind(summary, filtered_rnai[1,])
    summary <- rbind(summary, filtered_male[1,])
    
    
    summary$experiment <- c("WT vs nspc",
                            "tent-5 vs nspc/tent-5",
                            "daf-16 vs nspc/daf-16",
                            "WT vs tent-5",
                            "WT vs daf-16",
                            "nspc vs nscp/tent-5",
                            "nspc vs nspc/daf-16",
                            "WT vs nspc/tent-5",
                            "WT vs nspc/daf-16",
                            "RNAi: WT vs nspc",
                            "WT hermaphrodites vs males")    
    summary <- summary[,c(12,1:11)]
        
    final_summary <- summary %>% 
      dplyr::mutate(across(c(control_1:test_3), round, digits = 0)) %>% 
      dplyr::mutate(across(c(baseMean:log2FoldChange), round, digits = 3)) %>% 
      dplyr::mutate(across(c(padj), ~ format(signif(., digits = 2)))) %>% 
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
