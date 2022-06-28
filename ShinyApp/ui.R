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

ui <- fluidPage(
  
  titlePanel("AD Proteomics Comparisons"),
  selectizeInput("Comparison", "Select comparison:",
                 selected = "TransgenicBioIDvsALL", choices = NULL), 
  
  mainPanel(
    tabsetPanel(
      tabPanel(title = "Bar Plot", tabName = "bar_plots",
               br(),
               print(""),
               downloadButton("downloadCounts", "Download Full Counts Table"),
               br(),br(),br(),
               # Gene selector
               selectizeInput("gene", "Enter a Protein name:", selected = "GFAP", choices = NULL), 
               br(),
               plotOutput("Bar_plot"), hr(),
               fluidRow(DT::dataTableOutput("DE_data")),
               br(),br(),br(),
               br(),
      ),
      tabPanel(title = "Table", tabName = "data_table",
               br(),downloadButton("downloadDE", "Download Full DE Dataset"),
               print(""),
               fluidRow(DT::dataTableOutput("full_data")),
               br(),br(),
               br(),br(),br(),
               
      )
    )
  )
  )

