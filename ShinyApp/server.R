function(input, output, session) {
  
  counts_data <- read.csv("countsForShiny.csv")
  gene_options <- counts_data$Protein[1:4919]
  compare_options <- c("TransgenicBioIDvsALL","WTBioIDvsALL","TransgenictdTomvsALL",
                       "WTtdTomvsALL")
  
  # Update selected Comparison ---------
  updateSelectizeInput(session, "Comparison", selected = "TransgenicBioIDvsALL",
                       choices = compare_options, server = TRUE)
  
  fullDE <- reactive({
    selected_compare <- input$Comparison
    fullDE <- read.csv(paste0(selected_compare,".csv"))
  })

  # Update selected gene ------- 
  updateSelectizeInput(session, "gene", choices = gene_options, server = TRUE)
  
  # Select data for comparison of interest to plot count data
  selected_data <- reactive({
    selected_gene <- input$gene
    isolated_gene <- paste("^",selected_gene,"$", sep="") #Specifies the beginning and end of the string so only exact gene matches get plotted
    selected_data <- counts_data[grepl(isolated_gene, counts_data$Protein, ignore.case=TRUE),]
  })
  
  ## Plot Bar plot of Counts. --------
  output$Bar_plot <- renderPlot({
    #Plot with  selected gene
    plot <- selected_data() %>%
      ggplot(aes(x=Group, y=Count, fill=Group)) +
      scale_fill_manual(values=as.character(unique(counts_data$Color))) +
      geom_boxplot(outlier.size=2) +
      geom_jitter(color="black", size=2, position=position_jitter(0.2)) +
      theme(
        legend.position="none",
        plot.title = element_text(size=14, face="bold"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")
      ) +
      ggtitle(input$gene) +
      xlab("") +
      ylab("Counts")
    
    if (nrow(selected_data()) == 0) {
      plot <- NULL
    } else { return(plot) }
    
  })  
  
  ## Create Data tables ----------
  selectedDE <- reactive({
    de_data <- fullDE()
    selected_gene <- input$gene
    selectedDE <- de_data[grepl( selected_gene , de_data$Protein, ignore.case=TRUE),]  
  })
  
  output$DE_data <- DT::renderDataTable({
    DT::datatable(selectedDE(), options= list(paging=FALSE, dom= 't'),rownames = FALSE)
  })
  
  
  
  output$full_data <- DT::renderDataTable({
    DT::datatable(fullDE(), 
                  options = list(lengthMenu = c(5, 30, 50), pageLength = 5),rownames = FALSE)
  })
  ## Downoad Data------------
  output$downloadDE <- downloadHandler(
    filename = function() {
      paste(input$Comparison, Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(fullDE(), file)
    })
  
  counts <- read.csv("counts.csv")
  output$downloadCounts   <- downloadHandler(
    filename = function() {
      paste("AD_BioID_Counts", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(counts, file)
    })
}
