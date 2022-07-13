## Lets do the Pathway analysis for the shiny APP

### load libraries

library(clusterProfiler)


compare_options <- c("Transgenic_BioIDvsALL","WT_BioIDvsALL","Transgenic_tdTOMvsALL",
                     "WT_tdTOMvsALL","Transgenic_BioIDvsTransgenic_tdTOM","Transgenic_BioIDvsWT_tdTOM",
                     "WT_BioIDvsTransgenic_BioID","WT_BioIDvsWT_tdTOM")

mm <- search_kegg_organism('musculus', by='scientific_name')
wd <- "/gstore/scratch/u/lucast3/Alzheimers_Proteomics/AD_Proteomics/ADProteomicsShinyApp/inst/shiny/"

for (comparison in compare_options) {
  
  ## Read in Data
  data <- read.csv(file= paste0(wd, comparison,".csv"))
  
  ## Find instances of genes enriched and de- enriched
  enriched  <- which(data$Log2FC > .5)
  deenriched <- which(data$Log2FC < 0)
  
  # subset genes to comput pathway analysis
  enrichedGenes    <- data$Gene[enriched]
  de_enrichedGenes <- data$Gene[deenriched]
  
  
  #convert gene symbols to Entrez Ids
  
  enrichedGenes <- bitr(enrichedGenes,
                        fromType = "SYMBOL",
                        toType   = c("ENSEMBL","ENTREZID"),
                        OrgDb    = "org.Mm.eg.db")
  de_enrichedGenes <- bitr(de_enrichedGenes,
                        fromType = "SYMBOL",
                        toType   = c("ENSEMBL","ENTREZID"),
                        OrgDb    = "org.Mm.eg.db")
  
  
  ## compute Pathway analysis
  if(dim(enrichedGenes)[1] > 0){
  enrichedPathways <- enrichKEGG(gene         = enrichedGenes$ENTREZID,
                                 organism     = 'mmu',
                                 pvalueCutoff = 0.05)
  write.csv(enrichedPathways, file=paste0(wd,comparison,"_enrichedPathways.csv"))
  save(enrichedPathways,file=paste0(wd,comparison,"_enrichedPathways.RData"))
  }
  
  if(dim(de_enrichedGenes)[1] > 0){
  de_enrichedPathways <- enrichKEGG(gene         = de_enrichedGenes$ENTREZID,
                                 organism        = 'mmu',
                                 pvalueCutoff    = 0.05)
  
  write.csv(de_enrichedPathways, file=paste0(wd,comparison,"_de-enrichedPathways.csv"))
  save(de_enrichedPathways,file=paste0(wd,comparison,"_de-enrichedPathways.RData"))
  }

}

