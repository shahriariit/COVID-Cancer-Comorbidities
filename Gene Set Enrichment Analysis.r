
setwd("D:/corona")

library(RCurl)
library(GEOquery)
library(limma)
library(topGO)
library(genefilter)

gseRes <- read.csv("COVID_PBMC_File1_3.csv")
subtab_result <- subset(gseRes, select=c("GENE_SYMBOL","padj","pvalue","Tag","baseMean","log2FoldChange"))
genList <- subtab_result$log2FoldChange
names(genList) <- subtab_result$GENE_SYMBOL

topDifGenes <- function(value){
  return(abs(value)>1.00)
}

COVID_GOdata <- new("topGOdata",
                 description = "COVID study",
                 ontology = "BP", 
                 allGenes = genList,
                 geneSel = topDifGenes,
                 annot = annFUN.org, 
                 ID = "symbol", 
                 mapping = "org.Hs.eg.db",
                 nodeSize = 10)

n_sig <- sum(topDifGenes(genList))
sig <- sigGenes(COVID_GOdata)
uig <- usedGO(COVID_GOdata)

fisherRes <- runTest(COVID_GOdata, algorithm = "classic", statistic = "fisher")
ksRes <- runTest(COVID_GOdata, algorithm = "classic", statistic = "ks")

allRes <- GenTable(COVID_GOdata, classic = fisherRes, KS = ksRes, orderBy = "classic", topNodes = 30)
showSigOfNodes(COVID_GOdata, score(fisherRes), firstSigNodes = 5, useInfo = "all")
printGraph(COVID_GOdata, fisherRes, firstSigNodes = 5, fn.prefix = "COVID_GSE147507", useInfo = "all", pdfSW = TRUE)

allterms <- allRes$GO.ID
genes <- genesInTerm(COVID_GOdata,allterms)

for (i in 1:length(allterms))
{
  term <- allterms[i]
  term_genes <- genes[term][[1]]
  fact <- term_genes %in% sig
  term_genes_2 <- term_genes[fact == TRUE]
  term_genes_2 <- paste(term_genes_2, collapse=',')
  cat(term,"genes:",term_genes_2,"\n", append = TRUE, file = "COVID_GSE147507_correspondence.txt" )
}