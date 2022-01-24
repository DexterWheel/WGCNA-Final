#==== 
#5.  Network visualization using WGCNA functions
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
nGenes = ncol(data_processed)
nSamples = nrow(data_processed)

#==== 
#5.a   Visualizing the gene network
#====

###One way to visualize a weighted network is to plot its heatmap, Fig. 1.

###Each row and column of the heatmap correspond to a single gene.

###The heatmap can depict adjacencies or topological overlaps, with light colors denoting low adjacency (overlap) and darker colors higher adjacency (overlap).

###In addition, the gene dendrograms and module colors are plotted along the top and left side of the heatmap.

###The package provides a convenient function to create such network plots;

###Fig. 1 was created using the following code.

###This code can be executed only if the network was calculated using a single-block approach (that is, using the 1-step automatic or the step-by-step tutorials). 

###If the networks were calculated using the block-wise approach, the user will need to modify this code to perform the visualization in each block separately.

###The modification is simple and we leave it as an exercise for the interested reader.

###Calculate topological overlap anew: this could be done more efficiently by saving the TOM calculated during module detection, but let us do it again here.

lnames = load(file = "data-processed/wheatTOM-block.1.RData")
lnames = load(file = "data-processed/wheatTOM-block.2.RData")
lnames = load(file = "data-processed/wheatTOM-block.3.RData")
lnames = load(file = "data-processed/wheatTOM-block.4.RData")

dissTOM = 1-TOM
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^18
# Set diagonal to NA for a nicer plot
diag(TOM) = NA
# Call the plot function
sizeGrWindow(9,9)
TOMplot(TOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")

nSelect = 2000
# For reproducibility, we set the random seed
set.seed(10);
select = sample(nGenes, size = nSelect);
selectTOM = dissTOM[select, select];
# There’s no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select];
# Open a graphical window
sizeGrWindow(9,9)
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing
# the color palette; setting the diagonal to NA also improves the clarity of the plot
plotDiss = selectTOM^18;
diag(plotDiss) = NA;
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")

#==== 
#5.b   Visualizing the network of eigengenes
#====  

### It is often interesting to study the relationships among the found modules. 
###One can use the eigengenes as representative profiles and quantify module similarity by eigengene correlation.

###The package contains a convenient function plotEigengeneNetworks that generates a summary plot of the eigengene network. 

###It is usually informative to add a clinical trait (or multiple traits) to the eigengenes to see how the traits fit into the eigengene network:


# Recalculate module eigengenes
MEs = moduleEigengenes(data_processed, moduleColors)$eigengenes

# Isolate dev from the clinical traits
dev = as.data.frame(datTraits$dev);
names(dev) = "dev"
# Add the dev to existing module eigengenes
MET = orderMEs(cbind(MEs, dev))

# Plot the relationships among the eigengenes and the trait
sizeGrWindow(5,7.5);
par(cex = 0.9)
plotEigengeneNetworks(MET, "", 
                      marDendro = c(0,4,1,2), 
                      marHeatmap = c(3,4,1,2), 
                      cex.lab = 0.8, 
                      xLabelsAngle = 90)


### The function produces a dendrogram of the eigengenes and trait(s), and a heatmap of their relationships.

### To split the dendrogram and heatmap plots, we can use the following code

# Plot the dendrogram
sizeGrWindow(6,6)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", 
                      marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)

# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", 
                      marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, 
                      xLabelsAngle = 90)


#The eigengene dendrogram and heatmap identify groups of correlated eigengenes termed meta-modules. 
#For example, the dendrogram indicates that red, brown and bluw modules are highly related; 
#their mutual correlations are stronger than their correlations with weight. 
#On the other hand, the salmon module, which is also significantly correlated with weight, is not part of the same meta-module as the red, brown and blue modules, at least if meta-modules are defined as tight custers of modules (for example, modules with a correlation of eigengenes of at least 0.5).





