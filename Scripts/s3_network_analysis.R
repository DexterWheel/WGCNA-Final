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
nGenes = ncol(data_processed)
nSamples = nrow(data_processed)

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(data_processed, moduleColors)$eigengenes
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
# Define variable leaf containing the leaf column of datTrait
leaf = as.data.frame(datTraits$leaf)
#rename column
names(leaf) = "leaf"

# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(data_processed, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

geneTraitSignificance = as.data.frame(cor(data_processed, leaf, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) = paste("GS.", names(leaf), sep="")
names(GSPvalue) = paste("p.GS.", names(leaf), sep="")


#  Intramodular analysis: identifying genes with high GS and MM


###Using the GS and MM measures, we can identify genes that have a high significance for weight as well as high module membership in interesting modules.

###As an example, we look at the brown module that has the highest association with leaf.

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

#####
# development: without split
#####

# Define variable weight containing the weight column of datTrait
dev = as.data.frame(datTraits$dev)
#rename column
names(dev) = "dev"

# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(data_processed, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

geneTraitSignificance = as.data.frame(cor(data_processed, dev, use = "p"))
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
#####
# development: with split
#####

# Define variable weight containing the weight column of datTrait
frt_dev = as.data.frame(datTraits$frt_dev)
#rename column
names(frt_dev) = "frt_dev"

# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(data_processed, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

geneTraitSignificance = as.data.frame(cor(data_processed, dev, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) = paste("GS.", names(dev), sep="")
names(GSPvalue) = paste("p.GS.", names(dev), sep="")


module = "coral1"
column = match(module, modNames);

moduleGenes = moduleColors==module

par(mfrow = c(1,1))

verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for development",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

#results from this suggest that there may not be enough data for this kind of analysis for some developmental stages

#====
# S3.3    Summary output of network analysis results
#====

### We have found modules with high association with our trait of interest, 
#and have identified their central players by the Module Membership measure.

###We now merge this statistical information with gene annotation and write out a file that
#summarizes the most important results and can be inspected in standard spreadsheet software 
#such as MS Excel or Open Office Calc.

###Our expression data are only annotated by probe ID names: the command
names(data_processed)
###will return all probe IDs included in the analysis. Similarly,

names(data_processed)[moduleColors=="brown"]
###will return probe IDs belonging to the brown module. 

###To facilitate interpretation of the results, 
#we use a probe annotation file

lnames = load(file = "data-project/annotation-final.RData")

lnames = load(file = "annotation.RData")

dim(annotation_final)
names(annotation_final)
probes = names(data_processed)
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
                       go_id = annotation_final$go_id[probes2annot],
                       description = annotation_final$description[probes2annot],
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)

# Order modules by their significance for dev
modOrder = order(-abs(cor(MEs, leaf, use = "p")))

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
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.leaf));
geneInfo = geneInfo0[geneOrder, ]

#This data frame can be written into a text-format spreadsheet, for example by
write.csv(geneInfo, file = "data-processed/geneInfo_dev.csv")

#The reader is encouraged to open and view the file in a spreadsheet software, or inspect it directly within R using the command 
fix(geneInfo)
