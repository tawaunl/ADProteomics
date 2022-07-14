#
# Shiny App to interactively display AD proteomics Data
# Written by: Tawaun Lucas
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
# 

library(shiny)
library(dplyr)
library(ggplot2)
library(readr)
library(DT)
library(shinyjs)
library(gghighlight)
library(plotly)
library(clusterProfiler)

ui <- navbarPage("AD Proteomics Comparisons",
                 
                 tabPanel("Single Comparison",
                          selectizeInput("Comparison", "Select comparison:",
                                         selected = "TransgenicBioIDvsALL", choices = NULL), 
                          checkboxInput("interleaved",
                                        "Check to reverse comparison", value = FALSE),
                          tabsetPanel(
                            tabPanel(title = "Bar Plot", tabName = "bar_plots",
                                     br(),
                                     print(""),
                                     downloadButton("downloadCounts", "Download Full Counts Table"),
                                     br(),br(),br(),
                                     # Gene selector
                                     selectizeInput("gene", "Enter a Protein name:", selected = "GFAP", choices = NULL), 
                                     br(),
                                     plotlyOutput("Bar_plot",height = 650,width = 850),
                                     hr(),
                                     fluidRow(DT::dataTableOutput("DE_data")),
                                     br(),br(),br(),
                                     br()),
                            
                            tabPanel(title = "DE Table", tabName = "data_table",
                                     br(),downloadButton("downloadDE", "Download Full DE Dataset"),
                                     print(""),
                                     fluidRow(DT::dataTableOutput("full_data")),
                                     br(),br(),
                                     br(),br(),br()
                            ),
                            
                            tabPanel(title = "Pathways", tabName = "pathways",
                                     br(),
                                     downloadButton("downloadPathways", "Download Full Pathways Table"),
                                     br(),br(),
                                     fluidRow(
                                       column(width=4,uiOutput("de_enriched"))),
                                     print(""),
                                     useShinyjs(),
            
                                       radioButtons("visuBtn", NULL,
                                                    choices = c(Table = "table", Plot = "plot")),
                                     conditionalPanel(
                                       condition = "input.visuBtn == 'table'",
                                       DT::dataTableOutput("pathways_data")),
                                     conditionalPanel(
                                       condition = "input.visuBtn == 'plot'",
                                       plotOutput('pathwayplot', height = 700,width = 900)),
                                     br(),br(),br()
                            ))),
                 
                 tabPanel("4-Way Plot",
                          fixedRow( 
                            column(4,
                                   selectizeInput("x", "Select x- axis comparison:",
                                                  selected = "TransgenicBioIDvsALL", choices = NULL)),
                            column(4,
                                   selectizeInput("y", "Select y- axis comparison:",
                                                  selected = "TransgenicBioIDvsALL", choices = NULL))),
                          
                          #Downlaoad buttons
                          # downloadButton("download4way", "Download Plot"),
                          fixedRow( 
                            column(3,
                                   selectizeInput("quad", "Select Quadrant",
                                                  selected = "Q2", choices = NULL)),
                            column(3,
                                   downloadButton("downloadQuadrant", "Download Quadrant Genes"))),
                          
                          br(),br(),br(),
                          #Plot
                          br(),
                          plotlyOutput("four_way",height = 900,width = 850)
                          )
   
               
)


