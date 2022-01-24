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

###To facilitate a biological interpretation, we would like to know the gene ontologies of the genes in our modules, whether they are significantly enriched in certain functional categories etc

#==== 
# S4.1   Output gene lists for use with online software and services
#====     

###One option is to simply export a list of gene identifiers 
#that can be used as input for several popular gene ontology
#and functional enrichment analysis suites such as David or AmiGO. 

###For example, we write out the ensembl_gene_id codes for the modules into a file:
getwd()
# Read in the probe annotation created earlier via the biomart database
lnames = load(file = "data-project/annotation-final.RData")

# Match probes in the data set to the probe IDs in the annotation file
probes = names(data_processed)
probes2annot = match(probes, annotation_final$ensembl_gene_id)

# Get the corresponding ensembl_gene_id
allLLIDs = annotation_final$ensembl_gene_id[probes2annot]


#########
# $ Choose interesting modules
########

intModules = c("brown", "red", "salmon", "coral1")

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

#the above code just selects and makes a list of genes within the brown red and salmon modules via their ensmbl gene ids

# As background in the enrichment analysis, we will use all probes in the analysis.
fileName = paste("data-processed/ensemblgeneid-all.txt", sep="");
write.table(as.data.frame(allLLIDs), file = fileName,
            row.names = FALSE, col.names = FALSE)

#==== 
#4.b   Enrichment analysis 
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

### enrichmentAnalysis is apparently the new verion of GOenrichment analysis
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