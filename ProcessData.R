library(MSstatsTMT)
library(EnhancedVolcano)
library(BlueCopper2GP)
library(MultiAssayExperiment)
library(MSstats)
library(lme4)
library(UniProt.ws)
mouseUp <- UniProt.ws(10090)
setwd("/gstore/scratch/u/lucast3/Alzheimers_Proteomics/")
### ------ Read in MSstats input from Blue Copper -------------

# IDs 2129-full dataset, 2158- No129C134C135N
input <- getBCResource(id="2158", tag="MSstats Input")

quant.msstats <- proteinSummarization(input,
                                      method="msstats",
                                      global_norm=TRUE,
                                      reference_norm=TRUE,
                                      remove_norm_channel = TRUE,
                                      remove_empty_channel = TRUE)
x <- quant.msstats$FeatureLevelData
x <- x[-which(startsWith(as.character(x[,1]),"KRT")==TRUE),]
quant.msstats$FeatureLevelData <- x

x <- quant.msstats$ProteinLevelData
x <- x[-which(startsWith(as.character(x[,5]),"KRT")==TRUE),]
quant.msstats$ProteinLevelData <- x

mouseUp <- UniProt.ws(10090)

save(quant.msstats, file = "/gstore/scratch/u/lucast3/Alzheimers_Proteomics/ALZdata.Rdata")


## Plot QC for all proteins summarized

dataProcessPlotsTMT(data=quant.msstats, 
                     type='QCPlot',
                     width = 21, 
                     height = 7,
                    which.Protein = "allonly")


### ------ Group comparisons ---------
test.pairwise <- groupComparisonTMT(quant.msstats)
ProteinID <- vector()
Entry <- vector()
for(p in 1:length(levels(test.pairwise$ComparisonResult$Protein))){
  p.name <- levels(test.pairwise$ComparisonResult$Protein)[p]
  Entry[p] <- strsplit(p.name,"_")[[1]][1]
  id <- unlist(strsplit(strsplit(p.name,"_")[[1]][2],"|"))
  ProteinID[p] <- paste(id[7:length(id)],collapse = "")
}

labels <- select(mouseUp,
                 keys = ProteinID,
                 columns = c("GENENAME","PROTEIN-NAMES"),
                 keytype = "UNIPROTKB_ID" )
genes <- vector()
for(p in 1:length(ProteinID)){
  if(length(which(ProteinID[p]==labels[,1])) > 0){
    genes[p] <- labels[which(ProteinID[p]==labels[,1]),2]
  }
  else{
    genes[p] <- NA
  }
}

test.pairwise[["ComparisonResult"]]$GeneName <- test.pairwise$ComparisonResult$Protein
levels(test.pairwise[["ComparisonResult"]]$GeneName) <- genes


# Comarison level:"Transgenic_BioID2P2AtdTOM" "Transgenic_tdTOM" "WT_BioID2P2AtdTOM" "WT_tdTOM"
c.names <- c("Transgenic_BioID2P2AtdTOM","Transgenic_tdTOM", "WT_BioID2P2AtdTOM","WT_tdTOM" )

### Tg_bioID vs all --------------------------------------
comparisonTGB <- matrix(c(1,-1,-1,-1),nrow=1)
#Set names of each row and columns
row.names(comparisonTGB) <- "Tg_BioIDvsAll"
colnames(comparisonTGB) <- c.names

test.TGB <- groupComparisonTMT(data = quant.msstats, contrast.matrix = comparisonTGB, moderated = TRUE)

ProteinID <- vector()
Entry <- vector()
for(p in 1:length(levels(test.TGB$ComparisonResult$Protein))){
  p.name <- levels(test.TGB$ComparisonResult$Protein)[p]
  Entry[p] <- strsplit(p.name,"_")[[1]][1]
  id <- unlist(strsplit(strsplit(p.name,"_")[[1]][2],"|"))
  ProteinID[p] <- paste(id[7:length(id)],collapse = "")
}

labels <- select(mouseUp,
                 keys = ProteinID,
                 columns = c("GENENAME","PROTEIN-NAMES"),
                 keytype = "UNIPROTKB_ID" )
genes <- vector()
for(p in 1:length(ProteinID)){
  if(length(which(ProteinID[p]==labels[,1])) > 0){
    genes[p] <- labels[which(ProteinID[p]==labels[,1]),2]
  }
  else{
    genes[p] <- NA
  }
}


test.TGB[["ComparisonResult"]]$GeneName <- test.TGB$ComparisonResult$Protein
levels(test.TGB[["ComparisonResult"]]$GeneName) <- genes
#plot 
MSstatsConvert::MSstatsLogsSettings(FALSE, FALSE, FALSE)
groupComparisonPlots(test.TGB$ComparisonResult,type="VolcanoPlot",
                     address = "Tg_BioIDvsAll_")
groupComparisonPlots(test.TGB$ComparisonResult,type="ComparisonPlot",
                     address = "Tg_BioIDvsAll_")
#save table
write.csv2(test.TGB$ComparisonResult,file = "Tg_BioIDvsAll.csv")

###  ""Tg_tdTom""vs all --------------------------------------
comparisonTGT <- matrix(c(-1,1,-1,-1),nrow=1)
#Set names of each row and columns
row.names(comparisonTGT) <- "Tg_tdTomvsAll"
colnames(comparisonTGT) <- c.names

test.TGT <- groupComparisonTMT(data = quant.msstats, contrast.matrix = comparisonTGT, moderated = TRUE)

ProteinID <- vector()
Entry <- vector()
for(p in 1:length(levels(test.TGT$ComparisonResult$Protein))){
  p.name <- levels(test.TGT$ComparisonResult$Protein)[p]
  Entry[p] <- strsplit(p.name,"_")[[1]][1]
  id <- unlist(strsplit(strsplit(p.name,"_")[[1]][2],"|"))
  ProteinID[p] <- paste(id[7:length(id)],collapse = "")
}

labels <- select(mouseUp,
                 keys = ProteinID,
                 columns = c("GENENAME","PROTEIN-NAMES"),
                 keytype = "UNIPROTKB_ID" )
genes <- vector()
for(p in 1:length(ProteinID)){
  if(length(which(ProteinID[p]==labels[,1])) > 0){
    genes[p] <- labels[which(ProteinID[p]==labels[,1]),2]
  }
  else{
    genes[p] <- NA
  }
}

test.TGT[["ComparisonResult"]]$GeneName <- test.TGT$ComparisonResult$Protein
levels(test.TGT[["ComparisonResult"]]$GeneName) <- genes
#plot
MSstatsConvert::MSstatsLogsSettings(FALSE, FALSE, FALSE)
groupComparisonPlots(test.TGT$ComparisonResult,type="VolcanoPlot",
                     address = "TG_tdTomvsALL_")
groupComparisonPlots(test.TGT$ComparisonResult,type="ComparisonPlot",
                     address = "TG_tdTomvsALL_")
#save table
write.csv2(test.TGT$ComparisonResult,file = "TG_tdTomvsALL.csv")

### WT_BIOID vs all -----------------------------------------
comparisonWtB <- matrix(c(-1,-1,1,-1),nrow=1)
#Set names of each row and columns
row.names(comparisonWtB) <- "WtBvsALL"
colnames(comparisonWtB) <- c.names

test.WtB <- groupComparisonTMT(data = quant.msstats, contrast.matrix = comparisonWtB, moderated = TRUE)

ProteinID <- vector()
Entry <- vector()
for(p in 1:length(levels(test.WtB$ComparisonResult$Protein))){
  p.name <- levels(test.WtB$ComparisonResult$Protein)[p]
  Entry[p] <- strsplit(p.name,"_")[[1]][1]
  id <- unlist(strsplit(strsplit(p.name,"_")[[1]][2],"|"))
  ProteinID[p] <- paste(id[7:length(id)],collapse = "")
}

labels <- select(mouseUp,
                 keys = ProteinID,
                 columns = c("GENENAME","PROTEIN-NAMES"),
                 keytype = "UNIPROTKB_ID" )
genes <- vector()
for(p in 1:length(ProteinID)){
  if(length(which(ProteinID[p]==labels[,1])) > 0){
    genes[p] <- labels[which(ProteinID[p]==labels[,1]),2]
  }
  else{
    genes[p] <- NA
  }
}

test.WtB[["ComparisonResult"]]$GeneName <- test.WtB$ComparisonResult$Protein
levels(test.WtB[["ComparisonResult"]]$GeneName) <- genes
#plot 
MSstatsConvert::MSstatsLogsSettings(FALSE, FALSE, FALSE)
groupComparisonPlots(test.WtB$ComparisonResult,type="VolcanoPlot",
                     address = "WtBvsALL_")
groupComparisonPlots(test.WtB$ComparisonResult,type="ComparisonPlot",
                     address = "WtBvsALL_")
#save table
write.csv2(test.WtB$ComparisonResult,file = "WtBvsALL.csv")

### wtTdTom vs All---------------------------------------------------
comparisonWtTd <- matrix(c(-1,-1,-1,1),nrow=1)
#Set names of each row and columns
row.names(comparisonWtTd) <- "WtTdvsALL"
colnames(comparisonWtTd) <- c.names

test.WtTd <- groupComparisonTMT(data = quant.msstats, contrast.matrix = comparisonWtTd, moderated = TRUE)

ProteinID <- vector()
Entry <- vector()
for(p in 1:length(levels(test.WtTD$ComparisonResult$Protein))){
  p.name <- levels(test.WtTD$ComparisonResult$Protein)[p]
  Entry[p] <- strsplit(p.name,"_")[[1]][1]
  id <- unlist(strsplit(strsplit(p.name,"_")[[1]][2],"|"))
  ProteinID[p] <- paste(id[7:length(id)],collapse = "")
}

labels <- select(mouseUp,
                 keys = ProteinID,
                 columns = c("GENENAME","PROTEIN-NAMES"),
                 keytype = "UNIPROTKB_ID" )
genes <- vector()
for(p in 1:length(ProteinID)){
  if(length(which(ProteinID[p]==labels[,1])) > 0){
    genes[p] <- labels[which(ProteinID[p]==labels[,1]),2]
  }
  else{
    genes[p] <- NA
  }
}

test.WtTd[["ComparisonResult"]]$GeneName <- test.WtTd$ComparisonResult$Protein
levels(test.WtTd[["ComparisonResult"]]$GeneName) <- genes
#plot 
MSstatsConvert::MSstatsLogsSettings(FALSE, FALSE, FALSE)
groupComparisonPlots(test.WtTd$ComparisonResult,type="VolcanoPlot",
                     address = "WtTdvsALL_")
groupComparisonPlots(test.WtTd$ComparisonResult,type="ComparisonPlot",
                     address = "WtTdvsALL_")
#save table
write.csv2(test.WtTd$ComparisonResult,file = "WtTdvsALL.csv")


### ----- plot Group Comparisons -------------------

#Volcano plot - Creates pdf with all comparisons
groupComparisonPlots(test.pairwise$ComparisonResult,type="VolcanoPlot",
                     address = "AllPariwise_")

#Heatmap plot - Creates pdf with all comparisons
groupComparisonPlots(test.pairwise$ComparisonResult,type="Heatmap",
                     address = "AllPariwise_")

#Comparison Plots plot - Creates pdf with all comparisons
groupComparisonPlots(test.pairwise$ComparisonResult,type="ComparisonPlot",
                     address = "AllPariwise_")

#save table
write.csv2(test.pairwise$ComparisonResult,file = "AllPairwise.csv")


