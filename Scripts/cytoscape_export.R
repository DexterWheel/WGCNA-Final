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

# Recalculate topological overlap
TOM = TOMsimilarityFromExpr(data_processed, power = 18, networkType = "signed", TOMType = "signed");
# Read in the annotation file

lnames = load(file = "data-project/annotation-final.RData")

# Select module
module = "brown";
# Select module probes
probes = names(data_processed)
inModule = (moduleColors==module);
modProbes = probes[inModule];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
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
TOM = TOMsimilarityFromExpr(data_processed, power = 6)

# Read in the annotation file
lnames = load(file = "data-project/annotation-final.RData")

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