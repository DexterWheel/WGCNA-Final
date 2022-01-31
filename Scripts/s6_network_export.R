#==== 
#6.  Exporting a gene network to external visualization software
#====
# Load the WGCNA package
library(WGCNA)
library(dplyr)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Load the expression and trait data saved in the first part
lnames = load(file = "data-processed/wheat-dataInput.RData")
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = "data-processed/wheat-networkConstruction-auto.RData")
lnames

# Read in the annotation file
lnames = load(file = "data-project/annotation-final.RData")
#==== 
#6.b   Exporting to Cytoscape
#====

### Cytoscape [2] allows the user to input an edge file and a node file, allowing the user to specify for example the link weights and the node colors.

### Here we demonstrate the output of two modules, the red and brown ones, to Cytoscape.

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

#subsetting the original datasets for just the brown module genes
brown = data_processed[modProbes]
brown_annot= annotation_final[match(modProbes, annotation_final$ensembl_gene_id, nomatch=0),]
mod_go = brown_annot$go_id[match(modProbes, brown_annot$ensembl_gene_id)]

#####
#recreating the TOM with just the module genes
#####

powers = c(c(1:10), seq(from = 12, to=20, by=2))

sft = pickSoftThreshold(brown, powerVector = powers, verbose = 5, networkType = "signed")

sizeGrWindow(9, 5)
par(mfrow = c(1, 1))
cex1 = 0.9

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], 
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab = "Soft Threshold (power)", 
     ylab = "Scale Free Topology Model Fit,signed R^2",type = "n",
     main = paste("Scale independence"));

text(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], 
     sft$fitIndices[,5],
     xlab = "Soft Threshold (power)", 
     ylab = "Mean Connectivity", type = "n",
     main = paste("Mean connectivity"))

text(sft$fitIndices[,1], 
     sft$fitIndices[,5], 
     labels=powers, 
     cex=cex1,col="red")
#####
#power still 18

# Recalculate topological overlap
TOM = TOMsimilarityFromExpr(brown, power = 18)

######
#TOM = vectorTOM(data_processed,
 #               brown, 
  #              subtract1 = FALSE, 
   #             blockSize = 5000, 
    #            corFnc = "cor", corOptions = "use = 'p'", 
     #           networkType = "signed", 
      #          power = 18, 
       #         verbose = 1, indent = 0)

#save(TOM, file = "data-processed/TOM-brown.RData")
######

dimnames(TOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(TOM,
                               edgeFile = paste("export/CytoscapeInput-edges-",
                                                paste(module, collapse="-"),
                                                ".txt", sep=""),
                               nodeFile = paste("export/CytoscapeInput-nodes-", 
                                                paste(module, collapse="-"),
                                                ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.5,
                               nodeNames = modProbes,
                               altNodeNames = mod_go,
                               nodeAttr = moduleColors[inModule])


#Note that network input to Cytoscape is a bit more involved and the user should take care to select all necessary
#options for the edge and node files to be interpreted correctly. We refer the reader to Cytoscape documentation for
#all the necessary details.