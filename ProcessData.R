library(MSstatsTMT)
library(EnhancedVolcano)
library(BlueCopper2GP)
library(MultiAssayExperiment)
library(MSstats)
library(lme4)
library(UniProt.ws)
mouseUp <- UniProt.ws(10090)

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

ProteinIDFeat <- vector()
for(p in 1:length(levels(quant.msstats$FeatureLevelData$ProteinName))){
  p.name <- levels(quant.msstats$FeatureLevelData$ProteinName)[p]
  levels(quant.msstats$FeatureLevelData$ProteinName)[p] <- strsplit(p.name,"_")[[1]][1]
  id <- unlist(strsplit(strsplit(p.name,"_")[[1]][2],"|"))
  ProteinIDFeat[p] <- paste(id[7:length(id)],collapse = "")
}

ProteinIDPro <- vector()
for(p in 1:length(levels(quant.msstats$ProteinLevelData$Protein))){
  p.name <- levels(quant.msstats$ProteinLevelData$Protein)[p]
  levels(quant.msstats$ProteinLevelData$Protein)[p] <- strsplit(p.name,"_")[[1]][1]
  id <- unlist(strsplit(strsplit(p.name,"_")[[1]][2],"|"))
  ProteinIDPro[p] <- paste(id[7:length(id)],collapse = "")
}

mouseUp <- UniProt.ws(10090)

labels <- select(mouseUp,
                 keys = ProteinIDFeat,
                 columns = c("GENENAME","PROTEIN-NAMES"),
                 keytype = "UNIPROTKB_ID" )


quant.msstats$FeatureLevelData$Gene <- quant.msstats$FeatureLevelData$ProteinName
quant.msstats$FeatureLevelData$ID <- quant.msstats$FeatureLevelData$ProteinName
levels(quant.msstats$FeatureLevelData$Gene) <- labels$GENENAME
levels(quant.msstats$FeatureLevelData$ID) <- ProteinIDPro

labels <- select(mouseUp,
                 keys = ProteinIDPro,
                 columns = c("GENENAME","PROTEIN-NAMES"),
                 keytype = "UNIPROTKB_ID" )

quant.msstats$ProteinLevelData$Gene <- quant.msstats$ProteinLevelData$Protein
quant.msstats$ProteinLevelData$ID <- quant.msstats$ProteinLevelData$Protein
levels(quant.msstats$ProteinLevelData$Gene) <- labels$GENENAME
levels(quant.msstats$ProteinLevelData$ID) <- ProteinIDPro

save(quant.msstats, file = "/gstore/scratch/u/lucast3/Alzheimers_Proteomics/ALZdata.Rdata")



dataProcessPlotsTMT(data=quant.msstats,
                    type = 'ProfilePlot',
                    width = 21, 
                    height = 7)

## Plot QC for all proteins summarized

dataProcessPlotsTMT(data=quant.msstats, 
                     type='QCPlot',
                     width = 21, 
                     height = 7,
                    which.Protein = "allonly")


### ------ Group comparisons ---------
test.pairwise <- groupComparisonTMT(quant.msstats)

# Comarison level:"Transgenic_BioID2P2AtdTOM" "Transgenic_tdTOM" "WT_BioID2P2AtdTOM" "WT_tdTOM"
c.names <- c("Tg_BioID","Tg_tdTom", "WT_BioID","WT_tdTom")

### Tg_bioID vs all --------------------------------------
comparisonTGB <- matrix(c(1,-1,-1,-1),nrow=1)
#Set names of each row and columns
row.names(comparisonTGB) <- "Tg_BioIDvsAll"
colnames(comparisonTGB) <- c.names

test.TGB <- groupComparisonTMT(data = quant.msstats, contrast.matrix = comparisonTGB, moderated = TRUE)
#plot
groupComparisonPlots(test.TGB$ComparisonResult,type="VolcanoPlot",
                     address = "Tg_BioIDvsAll_")
groupComparisonPlots(test.TGB$ComparisonResult,type="ComparisonPlot",
                     address = "Tg_BioIDvsAll_")
#save table
write.csv2(test.TGB$ComparisonResult,file = "Tg_BioIDvsAl.csv")

###  ""Tg_tdTom""vs all --------------------------------------
comparisonTGT <- matrix(c(-1,1,-1,-1),nrow=1)
#Set names of each row and columns
row.names(comparisonTGT) <- "Tg_tdTomvsAll"
colnames(comparisonTGT) <- c.names

test.TGT <- groupComparisonTMT(data = quant.msstats, contrast.matrix = comparisonTGT, moderated = TRUE)
#plot
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
#plot
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
#plot
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


