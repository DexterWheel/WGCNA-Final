#==== 
# S0 : Sources
#====


# Paper : https://europepmc.org/article/MED/25853487 

  # Title : Deep transcriptome sequencing provides new insights into the structural and functional organization of the wheat genome.

    # ref : Pingault L, Choulet F, Alberti A, et al. 
        #   Deep transcriptome sequencing provides new insights into the structural and functional organization of the wheat genome.
        #   Genome Biology. 2015 Feb;16:29.
        #   DOI: 10.1186/s13059-015-0601-9. 
        #   PMID: 25853487; PMCID: PMC4355351.


#Database : Expression atlas: https://www.ebi.ac.uk/gxa/experiments/E-MTAB-4484/Results

  # Title : Transcription profiling by high throughput sequencing of five different organs of wheat 
        #   (leaf, shoot, root, spike, and grain) at three developmental stages each

    # Method : RNA-Seq mRNA baseline

      # Organism : Triticum aestivum

#The hope is to specify tissue defining pathways and genes

#==== 
# S1.1 : Data Import and reformatting
#====

library(WGCNA)
options(stringsAsFactors = FALSE)

  #Import raw fpkm data
  #fpkm: Fragments Per Kilobase of transcript per Million mapped reads
  #"In RNA-Seq, the relative expression of a transcript is proportional to the number of cDNA fragments that originate from it. 
  #Paired-end RNA-Seq experiments produce two reads per fragment, but that doesn't necessarily mean that both reads will be mappable. 
  #For example, the second read is of poor quality.
  #If we were to count reads rather than fragments, we might double-count some fragments but not others, leading to a skewed     expression value. 
  #Thus, FPKM is calculated by counting fragments, not reads"
data_raw0 = read.table("data-project/E-MTAB-4484-query-results.fpkms.tsv", header = T, fill = T)
data_raw0[1:10, 1:15]


  #Gene_name needs to be removed as it is unnecessary. Gene ids and tissue and development data need to transposed
data_raw1 = as.data.frame( t(data_raw0[, -c(1:2) ]) )
  #naming columns after the genes they represent
names(data_raw1) = data_raw0$Gene_ID
  #double checking that the names have been properly changed
data_raw1[1:15, 1:2]


#==== 
# S1.2 : Trimming data
#====

#####
# Need to convert some of the Ids so that they are included in annotations
#####

  # we need to remove genes and samples with too many missing values

gsg = goodSamplesGenes(data_raw1, verbose = 3)
gsg$allOK

  #allok = False and so samples need to be removed

if (!gsg$allOK){
  #Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes)>0)
      printFlush(paste("Removing genes:", paste(names(data_raw1)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(data_raw1)[!gsg$goodSamples], collapse = ", ")))
  #Remove the offending genes and samples from the data:
  data_raw1 = data_raw1[gsg$goodSamples, gsg$goodGenes]
    }

  #running the code again shows that all undesirable values have been removed
gsg = goodSamplesGenes(data_raw1, verbose = 3)
gsg$allOK


### Next we cluster the samples to see if there are any obvious outliers. 
  #in contrast to clustering genes that will come later)

sampleTree = hclust(dist(data_raw1), method = "average")

plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

  #The plot suggests that there are two major clusters. Apart from Stem BBCH 55 all tissues appear to cluster together. 
  #It has been decided to retain Stem BBCH 55 as there may be a developmental reason for the unexpected difference in clustering

#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9)

#====
# S1.3 : Trait data
#====

traitData = read.table("data-project/trait-data-tissue-split.txt", header = T)

dim(traitData); names(traitData)

  #
samples = rownames(data_raw1)
traitRows = match(samples, traitData$ref)

datTraits = traitData[traitRows, -7]

#renaming the rows 
rownames(datTraits) = traitData[traitRows, 7]

  # Re-cluster samples
sampleTree2 = hclust(dist(data_raw1), method = "average")

  # Convert traits to a color representation
traitColors = numbers2colors(datTraits, signed = T)

  # Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")

  # trait data visualises the findings mentioned above

###save expression and trait data for use in the next steps

save(data_raw1, datTraits, file = "data-processed/wheat-dataInput.RData")

#====
# S2.1 : soft threshold selection
#====

powers = c(c(1:10), seq(from = 12, to=20, by=2))

powers

# Call the network topology analysis function
sft = pickSoftThreshold(data_raw1, powerVector = powers, verbose = 5, networkType = "signed")

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

net = blockwiseModules(data_raw1, power = 18, maxBlockSize = 20000,
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
                    "Module colors",
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
     geneTree,
     file = "data-processed/wheat-networkConstruction-auto.RData")


#====
# S3.1 :  Quantifying module–trait associations
#====

#In this analysis we would like to identify modules that are significantly associated with the measured traits.

#Since we already have a summary profile (eigengene) for each module, 
#we simply correlate eigengenes with external traits and look for the most significant associations:

#It is not likely that we will find too much of biological significance here as we are more interested in the ontological data


lnames = load(file = "data-processed/wheat-networkConstruction-auto.RData")
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = "data-processed/wheat-dataInput.RData")
lnames

# Define numbers of genes and samples
nGenes = ncol(data_raw1)
nSamples = nrow(data_raw1)

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(data_raw1, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

#Since we have a moderately large number of modules and traits, 
#a suitable graphical representation will help in reading the table. 
#We color code each association by the correlation value:

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")

dim(textMatrix) = dim(moduleTraitCor)

par(mar = c(8, 18, 5, 3))

# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1,
               cex.legendLabel = 2,
               cex.lab = 2,
               cex.main = 2,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

#====
# S3.2 : Gene Significance and Module Membership
#====    

###We quantify associations of individual genes with our trait of interest (weight) by defining Gene Significance GS as (the absolute value of) the correlation between the gene and the trait. 

###For each module, we also define a quantitative measure of module membership MM as the correlation of the module eigengene and the gene expression profile.

###This allows us to quantify the similarity of all genes on the array to every module.
#####
# tissue
#####
# Define variable weight containing the weight column of datTrait
leaf = as.data.frame(datTraits$leaf)
#rename column
names(leaf) = "leaf"

# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(data_raw1, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

geneTraitSignificance = as.data.frame(cor(data_raw1, leaf, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) = paste("GS.", names(leaf), sep="")
names(GSPvalue) = paste("p.GS.", names(leaf), sep="")


#  Intramodular analysis: identifying genes with high GS and MM
      

###Using the GS and MM measures, we can identify genes that have a high significance for weight as well as high module membership in interesting modules.

###As an example, we look at the brown module that has the highest association with weight.

###We plot a scatterplot of Gene Significance vs. Module Membership in the brown module:


module = "brown"
column = match(module, modNames)

moduleGenes = moduleColors==module

sizeGrWindow(7, 7)

par(mfrow = c(1,1))

verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for tissue",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

#This plot. Clearly, GS and MM are highly correlated, illustrating that genes highly significantly associated with a trait are often also the most important (central) elements of modules associated with the trait.

#The reader is encouraged to try this code with other significance trait/module correlation (for example, the magenta, midnightblue, and red modules with weight).

#####
# development
#####

# Define variable weight containing the weight column of datTrait
dev = as.data.frame(datTraits$dev)
#rename column
names(dev) = "dev"

# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(data_raw1, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

geneTraitSignificance = as.data.frame(cor(data_raw1, dev, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) = paste("GS.", names(dev), sep="")
names(GSPvalue) = paste("p.GS.", names(dev), sep="")


module = "brown"
column = match(module, modNames);

moduleGenes = moduleColors==module

par(mfrow = c(1,1))

verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for development",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

#====
# S3.3    Summary output of network analysis results
#====      

### We have found modules with high association with our trait of interest, 
#and have identified their central players by the Module Membership measure.

###We now merge this statistical information with gene annotation and write out a file that
#summarizes the most important results and can be inspected in standard spreadsheet software 
#such as MS Excel or Open Office Calc.

###Our expression data are only annotated by probe ID names: the command
names(data_raw1)
###will return all probe IDs included in the analysis. Similarly,

names(data_raw1)[moduleColors=="brown"]
###will return probe IDs belonging to the brown module. 

###To facilitate interpretation of the results, 
#we use a probe annotation file

lnames = load(file = "data-project/annotation-final.RData")

#renaming the rows 
rownames(datTraits) = traitData[traitRows, 7]

dim(annotation_final)
names(annotation_final)
probes = names(data_raw1)
probes2annot = match(probes, annotation_final$ensembl_gene_id)


# The following is the number of probes without annotation:
sum(is.na(probes2annot))
# Should return 0.

###We now create a data frame holding the following information for all probes: 
#probe ID, gene symbol, Locus Link ID (Entrez code), module color, gene significance for dev, and module membership and p-values in all modules. 

###The modules will be ordered by their significance for dev, 
#with the most significant ones to the left.

# Create the starting data frame
geneInfo0 = data.frame(ensembl_gene_id = probes,
                       wheatexp_gene = annotation_final$wheatexp_gene[probes2annot],
                       description = annotation_final$description[probes2annot],
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)

# Order modules by their significance for dev
modOrder = order(-abs(cor(MEs, dev, use = "p")))

###-abs?#
##

# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership)){
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]])
  
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.dev));
geneInfo = geneInfo0[geneOrder, ]

#This data frame can be written into a text-format spreadsheet, for example by
write.csv(geneInfo, file = "data-processed/geneInfo_dev.csv")

#The reader is encouraged to open and view the file in a spreadsheet software, or inspect it directly within R using the command 
fix(geneInfo)

#==== 
# S4.0  Interfacing network analysis with other data such as functional annotation and gene ontology
#====      

library(WGCNA)
library(org.Mm.eg.db)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

lnames = load(file = "data-processed/wheat-networkConstruction-auto.RData")
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = "data-processed/wheat-dataInput.RData")

###Our previous analysis has identified several modules (labeled brown, red, and salmon) that are highly associated with weight.

###To facilitate a biological interpretation, we would like to know the gene ontologies of the genes in the modules, whether they are significantly enriched in certain functional categories etc

#==== 
# S4.1   Output gene lists for use with online software and services
#====     

###One option is to simply export a list of gene identifiers 
#that can be used as input for several popular gene ontology
#and functional enrichment analysis suites such as David or AmiGO. 

###For example, we write out the LocusLinkID (entrez) codes for the brown module into a file:
getwd()
# Read in the probe annotation
lnames = load(file = "data-project/annotation-final.RData")

# Match probes in the data set to the probe IDs in the annotation file
probes = names(data_raw1)
probes2annot = match(probes, annotation_final$ensembl_gene_id)

# Get the corresponding Locuis Link IDs
allLLIDs = annotation_final$ensembl_gene_id[probes2annot]


#########
# $ Choose interesting modules
########

intModules = c("brown", "red", "salmon")

for (module in intModules){
  # Select module probes
  modGenes = (moduleColors==module)
  # Get their entrez ID codes
  modLLIDs = allLLIDs[modGenes]
  # Write them into a file
  fileName = paste("data-processed/LocusLinkIDs-", module, ".txt", sep="")
  
  write.table(as.data.frame(modLLIDs), file = fileName,
              row.names = FALSE, col.names = FALSE)
}

# As background in the enrichment analysis, we will use all probes in the analysis.
fileName = paste("data-processed/LocusLinkIDs-all.txt", sep="");
write.table(as.data.frame(allLLIDs), file = fileName,
            row.names = FALSE, col.names = FALSE)

#==== 
#4.b   Enrichment analysis directly within R
#====

###The WGCNA package now contains a function to perform GO enrichment analysis using a simple, single step.

###To run the function, Biconductor packages GO.db, AnnotationDBI, and the appropriate organism-specific annotation package(s) need to be installed before running this code.

###The organism-specific packages have names of the form org.Xx.eg.db, where Xx stands for organism code, for example, Mm for mouse, Hs for human, etc.

###The only exception is yeast, for which no org.Xx.eg.db package is available; instead, the package carries the name org.Sc.sgd.db.

###Please visit the Bioconductor main page at http://www.bioconductor.org to download and install the required packages.

###In our case we are studying gene expressions from wheat, so this code needs the package ....

###Calling the GO enrichment analysis function GOenrichmentAnalysis is very simple. 

###The function takes a vector of module labels, and the Entrez (a.k.a. Locus Link) codes for the genes whose labels are given.


GOenr = GOenrichmentAnalysis(moduleColors, allLLIDs, organism = "wheat", nBestP = 10)

###   ###   ###   ###  Q  ###   ###   ###   ###

### enrichmentAnalysis is apparently the new verion on GOenrichment analysis
#The function runs for awhile and returns a long list, the most interesting component of which is:

tab = GOenr$bestPTerms[[4]]$enrichment

#This is an enrichment table containing the 10 best terms for each module present in moduleColors.

#Names of the columns within the table can be accessed by:
names(tab)

###We refer the reader to the help page of the function within R (available using ?GOenrichmentAnalysis at the R prompt) for details of what each column means. Because the term definitions can be quite long, the table is a bit difficult to display on the screen. 

#For readers who prefer to look at tables in Excel or similar spreadsheet software, it is best to
#save the table into a file and open it using their favorite tool:

write.table(tab, file = "data-processed/GOEnrichmentTable.csv", sep = ",", quote = TRUE, row.names = FALSE)


#On the other hand, to quickly take a look at the results, one can also abridge the table a bit and display it directly on screen:

keepCols = c(1, 2, 5, 6, 7, 12, 13);
screenTab = tab[, keepCols];
# Round the numeric columns to 2 decimal places:
numCols = c(3, 4);
screenTab[, numCols] = signif(apply(screenTab[, numCols], 2, as.numeric), 2)
# Truncate the the term name to at most 40 characters
screenTab[, 7] = substring(screenTab[, 7], 1, 40)
# Shorten the column names:
colnames(screenTab) = c("module", "size", "p-val", "Bonf", "nInTerm", "ont", "term name");
rownames(screenTab) = NULL;
# Set the width of R’s output. The reader should play with this number to obtain satisfactory output.
options(width=95)
# Finally, display the enrichment table:
screenTab
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
nGenes = ncol(data_raw1)
nSamples = nrow(data_raw1)

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
MEs = moduleEigengenes(data_raw1, moduleColors)$eigengenes

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
TOM = TOMsimilarityFromExpr(data_raw1, power = 6);
# Read in the annotation file

lnames = load(file = "data-project/annotation-final.RData")

# Select module
module = "brown";
# Select module probes
probes = names(data_raw1)
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
IMConn = softConnectivity(data_raw1[, modProbes]);
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
TOM = TOMsimilarityFromExpr(data_raw1, power = 6)

# Read in the annotation file
lnames = load(file = "data-project/annotation-final.RData")

# Select modules
modules = c("brown", "red")

# Select module probes
probes = names(data_raw1)
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