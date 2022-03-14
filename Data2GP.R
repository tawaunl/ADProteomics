# Script to get data from blue copper to GPSA in MultiAssay experiment.

### ------- Load Librarires --------------
library(BlueCopper2GP)
library(MultiAssayExperiment)


### ------ Load data -----------

proj <- BlueCopperProject("2129")

pro.url <- "http://bluecopper.gene.com/bluecopper/permalink/0790e005-b44d-459a-afa6-591b7ad98e7d"
pep.url <- "http://bluecopper.gene.com/bluecopper/permalink/a79d1eb2-479a-41d7-854f-dd88b277c67b"
name.s <- list(list(type= "genome",id="GRCm38"))
### ------ Prepare Data for GPSA

data <- prepareGipsyDataset(id="2091",
                    protein.url = pro.url,
                    peptide.url = pep.url,
                    authors = "lucast3",
                    technology = "TandemMassTagSpectrometry",
                    platform = "TandemMassTag",
                    organism = "Mus musculus",
                    namespace = name.s,
                    manual = TRUE,
                    bfc = FALSE,
                    update = FALSE)

## convert DFrames to Matrices
data@ExperimentList@listData$feature_level@assays@data@listData$raw_feature_intensity <- as(data@ExperimentList@listData$feature_level@assays@data@listData$raw_feature_intensity,"Matrix")
data@ExperimentList@listData[["protein_level"]]@assays@data@listData[["normalized_protein_intensity"]] <- as(data@ExperimentList@listData[["protein_level"]]@assays@data@listData[["normalized_protein_intensity"]],"Matrix")
data@ExperimentList@listData[["protein_level"]]@assays@data@listData[["log_normalized_protein_intensity"]] <- as(data@ExperimentList@listData[["protein_level"]]@assays@data@listData[["log_normalized_protein_intensity"]],"Matrix")

### Upload Dataset to DatasetDB ----------

library(dsassembly)

dsassembly::saveDataset(data, testing = TRUE)
