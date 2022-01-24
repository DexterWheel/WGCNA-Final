
library(BiocManager)
library(topGO)
library(ALL)



#first the topGOdata object must be created

#This object will contain:

#gene identifiers and their scores (p-values)

#GO annotations

#GO hierarchical structure

#all other information needed to perform the desired enrichment analysis

data(ALL)
data(geneList)

#geneList: 323, of genes and their corresponding p-values
View(geneList)

#The next data one needs are the GO terms, and the mapping that associate each gene with one or more GO term(s)

affyLib <- paste(annotation(ALL), "db", sep = ".")
library(package = affyLib, character.only = TRUE)



sampleGOdata <- new("topGOdata",
                    + description = "Simple session", ontology = "BP",
                    + allGenes = geneList, geneSel = topDiffGenes,
                    + nodeSize = 10,
                    + annot = annFUN.db, affyLib = affyLib)






