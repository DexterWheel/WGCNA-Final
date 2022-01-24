lnames = load(file = "data-processed/wheat-dataInput.RData")

#====
# S2.1 : soft threshold selection
#====

powers = c(c(1:10), seq(from = 12, to=20, by=2))

powers

# Call the network topology analysis function
sft = pickSoftThreshold(data_processed, powerVector = powers, verbose = 5, networkType = "signed")

# Plot the results:
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

# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")

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

#====
# S2.2 : One-step network construction and module detection
#====

#Constructing the gene network and identifying modules is now a simple function call:

net = blockwiseModules(data_processed, power = 18, maxBlockSize = 20000,
                       TOMType = "signed", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "data-processed/wheatTOM",
                       verbose = 3)
#We now return to the network analysis. 
#To see how many modules were identified and what the module sizes are, one can use table(net$colors).

table(net$colors)

#This indicates that there are 68 modules, 
#labeled 1 through 68 in order of descending size, 
#with sizes ranging from 6754 to 38 genes. 
#The label 0 is reserved for genes outside of all modules.


#hierarchical clustering dendrogram (tree) used for the module identification 
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath

plotDendroAndColors(net$dendrograms[[1]], 
                    mergedColors[net$blockGenes[[1]]],
                    "Module colors", main = "Gene dendrogram and module colors in block 1",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

plotDendroAndColors(net$dendrograms[[2]], 
                    mergedColors[net$blockGenes[[2]]],
                    "Module colors", main = "Gene dendrogram and module colors in block 2",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

plotDendroAndColors(net$dendrograms[[3]], 
                    mergedColors[net$blockGenes[[3]]],
                    "Module colors", main = "Gene dendrogram and module colors in block 3",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

plotDendroAndColors(net$dendrograms[[4]], 
                    mergedColors[net$blockGenes[[4]]],
                    "Module colors", main = "Gene dendrogram and module colors in block 4",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

plotDendroAndColors(net$dendrograms[[0]], 
                    mergedColors[net$blockGenes[[0]]],
                    "Module colors", main = "Gene dendrogram and module colors in block 4",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#Important     #We note that if the user would like to change some of the tree cut, module membership, and module merging criteria, the package provides the function recutBlockwiseTrees that can apply modified criteria without having to recompute the network and the clustering dendrogram. This may save a substantial amount of time.


#We now save the module assignment and module eigengene information necessary for subsequent analysis.

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, 
     moduleColors, 
     geneTree, net,
     file = "data-processed/wheat-networkConstruction-auto.RData")
