# Creating Volcano plots 
library(EnhancedVolcano)
library(MSstatsTMT)
library(UniProt.ws)
library(MSstats)
library(readr)

setwd("/gstore/scratch/u/lucast3/Alzheimers_Proteomics/")

## load data
### tgBioId vs all
data <- read_delim("Tg_BioIDvsAll.csv", 
                       delim = ";", escape_double = FALSE,
                       locale = locale(decimal_mark = ",", grouping_mark = "."),
                       trim_ws = TRUE)

### Plot Volcanco
x <- data$GeneName
data <- data[-which(startsWith(as.character(x),"Krt")==TRUE),]
y1 <- ceiling(max(-log10(data$adj.pvalue))*2)/2


plot <- EnhancedVolcano(data,
                lab = data$GeneName,
                x = 'log2FC',
                y = 'adj.pvalue',
                title = 'TgBioID vs. All', caption='',
                subtitle = paste("FDR cutoff = 0.05 \nFC cutoff = 1"),
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 1.0,
                labSize = 4.0,
                col=c('black', 'black', 'black', 'red3'),
                ylim = c(0,y1),
                ylab = bquote(~-Log[10]~ 'FDR'),
                drawConnectors = TRUE,border = 'full',colConnectors = 'grey30',
                widthConnectors = .5, borderWidth = 1.0, borderColour = 'black') + theme(legend.position="none")

pdf(file = "/gstore/scratch/u/lucast3/Alzheimers_Proteomics/Plots/Tg_BioIDvsALL_Volcano.pdf",width = 11.75,height=10.04)
print(plot)
dev.off()


### TG_tdTomvsALL vs all -------------------------------------------------
data <- read_delim("TG_tdTomvsALL.csv", 
                       delim = ";", escape_double = FALSE,
                       locale = locale(decimal_mark = ",", grouping_mark = "."),
                       trim_ws = TRUE)

### Plot Volcanco
x <- data$GeneName
data <- data[-which(startsWith(as.character(x),"Krt")==TRUE),]
y1 <- ceiling(max(-log10(data$adj.pvalue))*2)/2
EnhancedVolcano(data,
                lab = data$GeneName,
                x = 'log2FC',
                y = 'adj.pvalue',
                title = 'TG_tdTom vs. All', caption='',
                subtitle = paste("FDR cutoff = 0.05 \nFC cutoff = 1"),
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 1.0,
                labSize = 4.0,
                col=c('black', 'black', 'black', 'red3'),
                ylim = c(0,y1),
                ylab = bquote(~-Log[10]~ 'FDR'),
                drawConnectors = TRUE,border = 'full',colConnectors = 'grey30',
                widthConnectors = .5, borderWidth = 1.0, borderColour = 'black') + theme(legend.position="none")

pdf(file = "/gstore/scratch/u/lucast3/Alzheimers_Proteomics/Plots/Tg_tdTomvsAll_Volcano.pdf",width = 11.75,height=10.04)
print(plot)
dev.off()



### WtBioID vs all -----------------------------------------------------------
data <- read_delim("WtBvsALL.csv", 
                       delim = ";", escape_double = FALSE,
                       locale = locale(decimal_mark = ",", grouping_mark = "."),
                       trim_ws = TRUE)

### Plot Volcanco
x <- data$GeneName
data <- data[-which(startsWith(as.character(x),"Krt")==TRUE),]
y1 <- ceiling(max(-log10(data$adj.pvalue))*2)/2
EnhancedVolcano(data,
                lab = data$GeneName,
                x = 'log2FC',
                y = 'adj.pvalue',
                title = 'WtBioID vs. All', caption='',
                subtitle = paste("FDR cutoff = 0.05 \nFC cutoff = 1"),
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 1.0,
                labSize = 4.0,
                col=c('black', 'black', 'black', 'red3'),
                ylim = c(0,y1),
                ylab = bquote(~-Log[10]~ 'FDR'),
                drawConnectors = TRUE,border = 'full',colConnectors = 'grey30',
                widthConnectors = .5, borderWidth = 1.0, borderColour = 'black') + theme(legend.position="none")

pdf(file = "/gstore/scratch/u/lucast3/Alzheimers_Proteomics/Plots/WtBioIDvsAll_Volcano.pdf",width = 11.75,height=10.04)
print(plot)
dev.off()


### WttdTomato vs all ---------------------------------------------------
data <- read_delim("WtTdvsALL.csv", 
                       delim = ";", escape_double = FALSE,
                       locale = locale(decimal_mark = ",", grouping_mark = "."),
                       trim_ws = TRUE)

### Plot Volcanco
x <- data$GeneName
data <- data[-which(startsWith(as.character(x),"Krt")==TRUE),]
y1 <- ceiling(max(-log10(data$adj.pvalue))*2)/2
EnhancedVolcano(data,
                lab = data$GeneName,
                x = 'log2FC',
                y = 'adj.pvalue',
                title = 'WTtdTomato vs. All', caption='',
                subtitle = paste("FDR cutoff = 0.05 \nFC cutoff = 1"),
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 1.0,
                labSize = 4.0,
                col=c('black', 'black', 'black', 'red3'),
                ylim = c(0,y1),
                ylab = bquote(~-Log[10]~ 'FDR'),
                drawConnectors = TRUE,border = 'full',colConnectors = 'grey30',
                widthConnectors = .5, borderWidth = 1.0, borderColour = 'black') + theme(legend.position="none")

pdf(file = "/gstore/scratch/u/lucast3/Alzheimers_Proteomics/Plots/WtTdvsAll_Volcano.pdf",width = 11.75,height=10.04)
print(plot)
dev.off()