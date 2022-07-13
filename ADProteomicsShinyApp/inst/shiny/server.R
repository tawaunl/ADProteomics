library(dplyr)
library(ggplot2)
library(readr)
library(DT)
library(clusterProfiler)
library(multienrichjam)
library(plotly)
library(processx)
function(input, output, session) {
  
  counts_data <- read.csv("countsForShiny.csv")
  counts <- read.csv("counts.csv") # This is a more user reader format for download
  gene_options <- counts_data$Protein[1:dim(counts)[1]]
  compare_options <- c("Transgenic_BioIDvsALL","WT_BioIDvsALL","Transgenic_tdTOMvsALL",
                       "WT_tdTOMvsALL","Transgenic_BioIDvsTransgenic_tdTOM","Transgenic_BioIDvsWT_tdTOM",
                       "WT_BioIDvsTransgenic_BioID","WT_BioIDvsWT_tdTOM")
  
  # Update selected Comparison ---------
  updateSelectizeInput(session, "Comparison", selected = "Transgenic_BioIDvsALL",
                       choices = compare_options, server = TRUE)
  
  fullDE <- reactive({
    selected_compare <- input$Comparison
    if(input$interleaved==TRUE){
      fullDE <- read.csv(paste0(selected_compare,".csv"))
      x <- fullDE$Log2FC*-1
      fullDE$Log2FC <- x
      fullDE <- data.frame(fullDE)
    }
    
    else{fullDE <- read.csv(paste0(selected_compare,".csv")) }
    
  })
  
  ### Choose correct Pathways table for comparison ----------
  
  output$de_enriched <- renderUI({
    list(
      checkboxInput("de_enriched",
                    "Check for de-enriched pathways", value = FALSE)
    )
  })
  output$barplot <- renderUI({
    list(
      checkboxInput("barplot",
                    "Check to plot pathways", value = FALSE)
    )
  })
  
  ### PAthway Analysis tables and plot  ---------
  
  pathways <- reactive({
    selected_compare <- input$Comparison
    if(input$de_enriched==FALSE){
      pathways <- read.csv(paste0(selected_compare,"_enrichedPathways.csv"))}
    
    if(input$de_enriched==TRUE){
      pathways <- read.csv(paste0(selected_compare,"_de-enrichedPathways.csv"))}

     pathways <- pathways[,c(2:8,10)]
  })
  
  output$pathways_data <- DT::renderDataTable({
    DT::datatable(pathways(),
                  options = list(lengthMenu = c(10,15,20, 50), pageLength =15 ,rowCallback = JS(
                    "function(row, data) {",
                    "for (i = 1; i < data.length; i++) {",
                    "if (data[i]>1000 | data[i]<1){",
                    "$('td:eq('+i+')', row).html(data[i].toExponential(1));",
                    "}",
                    "}",
                    "}")),rownames = FALSE)
  })
  
  output$pathwayplot <- renderPlot({
    selected_compare <- input$Comparison
    if(input$de_enriched==FALSE){
      data <- read.csv(paste0(selected_compare,"_enrichedPathways.csv"))}
    
    if(input$de_enriched==TRUE){
      data <- read.csv(paste0(selected_compare,"_de-enrichedPathways.csv"))}  
    
    data<- enrichDF2enrichResult(data)
    barplot(data,  showCategory=20)
  })
  
  # Update selected gene ------- 
  updateSelectizeInput(session, "gene",selected = "GFAP", choices = gene_options, server = TRUE)
  
  # Select data for comparison of interest to plot count data---------
  selected_data <- reactive({
    selected_gene <- input$gene
    selected_compare <- input$Comparison
    
    if(grepl( "ALL", selected_compare, fixed = TRUE)){ isolated_comp <- c("Transgenic_BioID2P2AtdTOM","Transgenic_tdTOM",
                                                                          "WT_BioID2P2AtdTOM","WT_tdTOM")}
    if(selected_compare == "Transgenic_BioIDvsTransgenic_tdTOM"){
      isolated_comp <- c("Transgenic_BioID2P2AtdTOM","Transgenic_tdTOM")}
    if(selected_compare == "Transgenic_BioIDvsWT_tdTOM"){
      isolated_comp <- c("Transgenic_BioID2P2AtdTOM","WT_tdTOM")}
    if(selected_compare == "WT_BioIDvsTransgenic_BioID"){
      isolated_comp <- c("WT_BioID2P2AtdTOM","Transgenic_BioID2P2AtdTOM")}
    if(selected_compare == "WT_BioIDvsWT_tdTOM"){
      isolated_comp <- c("WT_BioID2P2AtdTOM","WT_tdTOM")}
    
    isolated_gene <- paste("^",selected_gene,"$", sep="") #Specifies the beginning and end of the string so only exact gene matches get plotted
    selected_data <- counts_data[grepl(isolated_gene, counts_data$Protein, ignore.case=TRUE),]
    selected_data %>% filter(Group %in% isolated_comp)
    
    
  })
  
  ## Plot Bar plot of Counts. --------
  output$Bar_plot <- renderPlotly({
    
    if (startsWith(input$Comparison,"Transgenic_B")) {
      group_order <- factor(levels=c("Transgenic_BioID2P2AtdTOM","Transgenic_tdTOM",
                                     "WT_BioID2P2AtdTOM","WT_tdTOM"))
      if(input$interleaved==TRUE){
        group_order <- factor(levels=c("Transgenic_tdTOM",
                                       "WT_BioID2P2AtdTOM","WT_tdTOM","Transgenic_BioID2P2AtdTOM"))
      }}
    if (startsWith(input$Comparison,"Transgenic_t")) {
      group_order <- factor(levels=c("Transgenic_tdTOM","Transgenic_BioID2P2AtdTOM",
                                     "WT_BioID2P2AtdTOM","WT_tdTOM"))
      if(input$interleaved==TRUE){
        group_order <- factor(levels=c("WT_BioID2P2AtdTOM","WT_tdTOM",
                                       "Transgenic_BioID2P2AtdTOM","Transgenic_tdTOM"))
      }}
    if (startsWith(input$Comparison,"WT_B") ) {
      group_order <- factor(levels=c("WT_BioID2P2AtdTOM","WT_tdTOM",
                                     "Transgenic_BioID2P2AtdTOM","Transgenic_tdTOM"))
      if(input$interleaved==TRUE){
        group_order <- factor(levels=c("WT_tdTOM","Transgenic_BioID2P2AtdTOM",
                                       "Transgenic_tdTOM","WT_BioID2P2AtdTOM"))
      }}
    if (startsWith(input$Comparison,"WT_t") ) {
      group_order <- factor(levels=c("WT_tdTOM","WT_BioID2P2AtdTOM",
                                     "Transgenic_BioID2P2AtdTOM","Transgenic_tdTOM")) 
      if(input$interleaved==TRUE){
        group_order <- factor(levels=c("Transgenic_BioID2P2AtdTOM","Transgenic_tdTOM",
                                       "WT_BioID2P2AtdTOM","WT_tdTOM"))
      }}
    
    plotData <- selected_data() %>%
      arrange(Count) %>%
      mutate(Group = factor(Group, levels=levels(group_order))) %>%
      arrange(Group)
    
    #Plot with  selected gene
    plot <- plot_ly(plotData,type = "box",y=plotData$Count,x=plotData$Group,
                    boxpoints="all",jitter = 0.3,
                    pointpos = 0,text = c(paste0(plotData$Sample,"\n Count =",plotData$Count)),
                    hoverinfo = 'text',marker = list(plotData$Color),color=plotData$Color,
                    showlegend=FALSE) 
    
      
    
    if (nrow(selected_data()) == 0) {
      plot <- NULL
    } else { return(plot) }
    
    print(plot%>%
            layout(title = input$gene,
                   xaxis = list(autotypenumbers = 'strict', title = 'Group'),
                   yaxis = list(title = 'log(Normalized Protein Intensity)'),
                   plot_bgcolor='#ffff',
                   xaxis = list(
                     zerolinecolor = '#ffff',
                     zerolinewidth = 2,
                     gridcolor = 'ffff'),
                   yaxis = list(
                     zerolinecolor = '#ffff',
                     zerolinewidth = 2,
                     gridcolor = 'ffff')))
  })  
  
  ## Create Data tables ----------
  selectedDE <- reactive({
    de_data <- fullDE()
    selected_gene <- input$gene
    selectedDE <- de_data[grepl(selected_gene, de_data$Protein, ignore.case=TRUE),]
  })
  
  
  output$DE_data <- DT::renderDataTable({
    DT::datatable(selectedDE(), options= list(paging=FALSE, dom= 't'),rownames = FALSE)
  })
  
  
  
  output$full_data <- DT::renderDataTable({
    DT::datatable(fullDE(), 
                  options = list(lengthMenu = c(5, 30, 50), pageLength = 5),rownames = FALSE)
  })
  
  
  ### Create 4-way plots ---------------------
  prepare4way <- function(markers.1,markers.2){
    x <- markers.1$Protein
    y <- markers.2$Protein
    rownames(markers.1)<-markers.1$Protein
    rownames(markers.2) <- markers.2$Protein
    commonGenes <- intersect(x,y)
    sub1 <- markers.1[commonGenes,]
    sub2 <- markers.2[commonGenes,]
    df <- data.frame(markers.1,markers.2)
    return(df)
  }
  
  plot4way <- function(df, comparisonName,xlab, ylab){
    rsq <- round(cor(x=df$Log2FC,y=df$Log2FC.1 )^2,3)
    
    plot <- plot_ly(df,type = "scatter",mode='markers', y=df$Log2FC.1,x=df$Log2FC,
                    text = c(paste0("Protein: ",df$Protein,
                                    "\n x-FC = ",df$Log2FC,
                                    "\n y-FC = ",df$Log2FC.1,
                                    "\n Gene =",df$Gene,
                                    "\n Description=",df$Description)),
                    hoverinfo = 'text',marker = list(color="black"),
                    showlegend=FALSE) 
     
     plot <- plot %>% layout(title= paste0(comparisonName,' comparison',
                                   '<br>','<sup>','R2 = ',rsq),
                     xaxis = list(autotypenumbers = 'strict', title = paste("log2(FoldChange)",xlab)),
                     yaxis = list(title = paste("log2(FoldChange)",ylab)),
                     plot_bgcolor='#ffff',
                     xaxis = list(
                       zerolinecolor = '#ffff',
                       zerolinewidth = 2,
                       gridcolor = 'ffff'),
                     yaxis = list(
                       zerolinecolor = '#ffff',
                       zerolinewidth = 2,
                       gridcolor = 'ffff'))
     return(plot)
      
  }
  
  updateSelectizeInput(session, "x", selected = "Transgenic_BioIDvsALL",
                       choices = compare_options, server = TRUE)
  
  updateSelectizeInput(session, "y", selected = "WT_BioIDvsALL",
                       choices = compare_options, server = TRUE)
  xDE <- reactive({
    selected_compare <- input$x
    xDE <- read.csv(paste0(selected_compare,".csv"))
    
  })
  
  yDE <- reactive({
    selected_compare <- input$y
    yDE <- read.csv(paste0(selected_compare,".csv"))
    
  })
  
  fourWayinput <- reactive({
    df <- prepare4way(xDE(),yDE())
    comparisonName <- paste(input$x,"vs.",input$y)
    fourWayinput <- plot4way(df,comparisonName,xlab=input$x,ylab=input$y)
    
  })
  
  Quadrant <- reactive({
    df <- prepare4way(xDE(),yDE())
    
  })
  output$four_way <- renderPlotly({
    
    qdf <- data.frame(Quadrant())
    q1 <- qdf[which(qdf$Log2FC<0 & qdf$Log2FC.1>0),]
    q2 <- qdf[which(qdf$Log2FC>0 & qdf$Log2FC.1>0),]
    q3 <- qdf[which(qdf$Log2FC<0 & qdf$Log2FC.1<0),]
    q4 <- qdf[which(qdf$Log2FC>0 & qdf$Log2FC.1<0),]
    
    if(input$quad=="Q1"){
      plot <- fourWayinput() %>% add_trace(q1,name="Q1", mode='markers', y=q1$Log2FC.1,x=q1$Log2FC,
                                         text = c(paste0("Protein: ",q1$Protein,
                                                         "\n x-FC = ",q1$Log2FC,
                                                         "\n y-FC = ",q1$Log2FC.1,
                                                         "\n Gene =",q1$Gene,
                                                         "\n Description=",q1$Description)),
                                         hoverinfo = 'text',marker = list(color="red"),
                                         showlegend=FALSE)  
    }
    if(input$quad=="Q2"){
      plot <- fourWayinput() %>% add_trace(q2,name="Q2", mode='markers', y=q2$Log2FC.1,x=q2$Log2FC,
                                        text = c(paste0("Protein: ",q2$Protein,
                                                        "\n x-FC = ",q2$Log2FC,
                                                        "\n y-FC = ",q2$Log2FC.1,
                                                        "\n Gene =",q2$Gene,
                                                        "\n Description=",q2$Description)),
                                        hoverinfo = 'text',marker = list(color="red"),
                                        showlegend=FALSE) 
    }
    if(input$quad=="Q3"){
      plot <- fourWayinput() %>% add_trace(q3,name="Q3", mode='markers', y=q3$Log2FC.1,x=q3$Log2FC,
                                        text = c(paste0("Protein: ",q3$Protein,
                                                        "\n x-FC = ",q3$Log2FC,
                                                        "\n y-FC = ",q3$Log2FC.1,
                                                        "\n Gene =",q3$Gene,
                                                        "\n Description=",q3$Description)),
                                        hoverinfo = 'text',marker = list(color="red"),
                                        showlegend=FALSE) 
    }
    if(input$quad=="Q4"){
      plot <- fourWayinput() %>% add_trace(q4,name="Q4", mode='markers', y=q4$Log2FC.1,x=q4$Log2FC,
                                        text = c(paste0("Protein: ",q4$Protein,
                                                        "\n x-FC = ",q4$Log2FC,
                                                        "\n y-FC = ",q4$Log2FC.1,
                                                        "\n Gene =",q4$Gene,
                                                        "\n Description=",q4$Description)),
                                        hoverinfo = 'text',marker = list(color="red"),
                                        showlegend=FALSE) 
    }
    print(plot)
    
  })
  
  
  ## Download Data------------
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
  
  ##### Download Pathway table ------------
  output$downloadPathways <- downloadHandler(
    filename = function() {
      if(input$de_enriched==TRUE){direction="_de-enriched"}
      else{direction="_enriched" }
      paste(input$Comparison, Sys.Date(), direction,"Pathways.csv", sep="")
    },
    content = function(file) {
      write.csv(pathways(), file)
    })
  
  ## Download 4 way plot -------------
  # output$download4way <- downloadHandler(
  #   filename = function() {
  #     comparisonName <- paste(input$x,"vs.",input$y)
  # 
  #     paste(comparisonName, Sys.Date(),"4wayPlot.png", sep="")
  #   },
  #   content = function(file) {
  #     orca(fourWayinput(),file)
  #   })
  ### Download 4 way quadrants ------------
  quad_options <- c("Q1","Q2","Q3","Q4")
  updateSelectizeInput(session, "quad", selected = "Q2",
                       choices = quad_options, server = TRUE)
  
  getquadrantGenes <- function(markers.1,markers.2,quadrant){
    if(quadrant == "Q1"){
      x <- which(markers.1$Log2FC < 0)
      y <- which(markers.2$Log2FC > 0 )
    }
    if(quadrant == "Q2"){
      x <- which(markers.1$Log2FC > 0)
      y <- which(markers.2$Log2FC > 0 )
    }
    if(quadrant =="Q3"){
      x <- which(markers.1$Log2FC < 0)
      y <- which(markers.2$Log2FC < 0 )
    }
    if(quadrant=="Q4"){
      x <- which(markers.1$Log2FC > 0)
      y <- which(markers.2$Log2FC < 0 )
    }
    markers.1 <- markers.1[x,]
    markers.2 <- markers.2[y,]
    rownames(markers.1) <- markers.1$Protein
    rownames(markers.2) <- markers.2$Protein
    x <- rownames(markers.1)
    y <- rownames(markers.2)
    
    commonGenes <- intersect(x,y)
    sub1 <- markers.1[commonGenes,]
    sub2 <- markers.2[commonGenes,]
    df <- data.frame(sub1,sub2[,c(4,5)])
    colnames(df)<-c("Gene","Protein","Description",paste0("p_Value_",input$x),
                   paste0("Log2FC_",input$x),
                   paste0("p_Value_",input$y),
                   paste0("Log2FC_",input$y))
    return(df)
  }
  quadInput <- reactive({
    df <- getquadrantGenes(xDE(),yDE(),input$quad)
  })
  output$downloadQuadrant <- downloadHandler(
    filename = function() {
      comparisonName <- paste(input$x,"vs.",input$y)
      paste(comparisonName,input$quad,"_genes.csv", sep="")
    },
    content = function(file) {
      write.csv(quadInput(), file,row.names = FALSE)
    })
  
  
}
