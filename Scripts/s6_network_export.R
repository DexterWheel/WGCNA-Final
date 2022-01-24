#==== 
#6.  Exporting a gene network to external visualization software
#====
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Load the expression and trait data saved in the first part
lnames = load(file = "data-processed/wheat-dataInput.RData")
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = "data-processed/wheat-networkConstruction-auto.RData")
lnames

#==== 
#6.a   Exporting to VisANT
#====

### The package provides a convenient function for exporting the network to VisANT [1]. 
###We illustrate a simple export of the full weighted network of a single module.

# We need the Topological overlap matrix (TOM) to export the data, however due to the size of the dataset it is not possible to create a single TOM file with my current machine. I would need at least 128gb of ram, potentially less on linux, to complete the analysis in one block.

#thankfully the function Vector TOM allows us to create a TOM using just the genes from our chosen modules while taking into account the entire dataset.

#This allows us to export this data to external software

# Select module
module = "brown"
# Select module probes
probes = names(data_processed)
inModule = (moduleColors==module)
modProbes = probes[inModule]

#modProbes provides us with the list of genes we need to select from the original dataset 

#subsetting the original dataset for just the brown module genes
brown = data_processed[modProbes]

TOM = vectorTOM(data_processed,
                brown, 
                subtract1 = FALSE, 
                blockSize = 5000, 
                corFnc = "cor", corOptions = "use = 'p'", 
                networkType = "signed", 
                power = 18, 
                verbose = 1, indent = 0)

save(TOM, file = "data-processed/TOM-brown.RData")

# Read in the annotation file

lnames = load(file = "data-project/annotation-final.RData")



# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into an edge list file VisANT can read
vis = exportNetworkToVisANT(modTOM,
                            file = paste("export/VisANTInput-", module, ".txt", sep=""),
                            weighted = TRUE,
                            threshold = 0,
                            probeToGene = data.frame(annotation_final$ensembl_gene_id, annotation_final$gene_symbol) )

###Because the brown module is rather large, 
#we can restrict the genes in the output to say the 30 top hub genes in the module:

nTop = 30;
IMConn = softConnectivity(data_processed[, modProbes]);
top = (rank(-IMConn) <= nTop)
vis = exportNetworkToVisANT(modTOM[top, top],
                            file = paste("export/VisANTInput-", module, "-top30.txt", sep=""),
                            weighted = TRUE,
                            threshold = 0,
                            probeToGene = data.frame(annot$ensembl_gene_id, annot$gene_symbol) )
#==== 
#6.b   Exporting to Cytoscape
#====    


### Cytoscape [2] allows the user to input an edge file and a node file, allowing the user to specify for example the link weights and the node colors.

### Here we demonstrate the output of two modules, the red and brown ones, to Cytoscape.

# Recalculate topological overlap if needed
TOM = TOMsimilarityFromExpr(data_processed, power = 18)


lnames = load(file = "data-processed/wheatTOM-block.1.RData")

TOM1 = TOM

lnames = load(file = "data-processed/wheatTOM-block.2.RData")

TOM2 = TOM

lnames = load(file = "data-processed/wheatTOM-block.3.RData")

TOM3 = TOM

lnames = load(file = "data-processed/wheatTOM-block.4.RData")

TOM4 = TOM

TOMs = c(TOM1, TOM2, TOM3, TOM4)


# Read in the annotation file
lnames = load(file = "data-project/annotation-final.RData")

annot = annotation_final
# Select modules
modules = c("brown", "red")

# Select module probes
probes = names(data_processed)
inModule = is.finite(match(moduleColors, modules))
modProbes = probes[inModule]
modGenes = annot$gene_symbol[match(modProbes, annot$ensembl_gene_id)]

# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule]

dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("export/CytoscapeInput-edges-", 
                                                paste(modules, collapse="-"),
                                                ".txt", sep=""),
                               nodeFile = paste("export/CytoscapeInput-nodes-", 
                                                paste(modules, collapse="-"),
                                                ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule])


#Note that network input to Cytoscape is a bit more involved and the user should take care to select all necessary
#options for the edge and node files to be interpreted correctly. We refer the reader to Cytoscape documentation for
#all the necessary details.