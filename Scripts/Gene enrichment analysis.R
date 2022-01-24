library(topGO)
library(ALL)
library(WGCNA)
#####
# Corresponding section in master script is s4
#####

  #Analysis of the functions of genes within modules will provide an insight into the biological pathways that these modules are involved in.
  #This could tell us which genes are most integral to these pathways as well as how they could be manipulated directly or indirectly. Modules that appear to be closely co-regulated may be discovered.

  #To do this analysis we need to map the functions of the genes in the network to their gene ontology (GO) terms
  #We need to do an Enrichment analysis

  #####
  #A typical session can be divided into three steps:

  #1. Data preparation: List of genes identifiers, gene scores, list of differentially expressed genes, as well as gene-to-GO annotations are all collected and stored in a single R object.(topGOdata)

  #2. Running the enrichment tests: Using the object created in the first step the user can perform enrichment
  #analysis using any feasible mixture of statistical tests and methods that deal with the GO topology.

  #3. Analysis of the results: The results obtained in the second step are analysed using summary functions
  #and visualisation tools

#####
#1. data prep
#####

#first the topGOdata object must be created

#This object will contain:

    #gene identifiers (ensmbl gene id) and their scores (p-values)

    #GO annotations (obtained via biomart)

    #GO hierarchical structure ()

    #all other information needed to perform the desired enrichment analysis

#####
# loading data input/network
#####
options(stringsAsFactors = FALSE)

#load in data input which contains just the genes, the pvalues and which individuals they were taken from
lnames = load(file = "data-processed/wheat-dataInput.RData")
#load in network created in s2
lnames = load(file = "data-processed/wheat-networkConstruction-auto.RData")

#the network data will tell us which genes are found in which modules

data_processed2 <- t(data_processed)

#####
#annotations/GO terms
#####

lnames = load(file = "annotation.RData")

#annot_list is a list data frames. These dataframes each contain ensembl gene ids and data for a particular attribute from the Biomart database

#these dataframes have been kept separate due to differing sizes. This is because some attributes such as GO ids have multiple matches corresponding to a single gene

#####

#next we need to create topGOdata. This object will "contain all gene identifiers and their scores, the GO annotations, the GO hierarchical structure and all other information needed to perform the desired enrichment analysis."

sampleGOdata <- new("topGOdata",
                    + description = "Simple session", ontology = "BP",
                    + allGenes = data_processed, geneSel = topDiffGenes,
                    + nodeSize = 10,
                    + annot = annFUN.db, affyLib = affyLib)

#####
#2. Enrichment test
#####

#the enrichment test must be carried out on each module individually
#there are a total of ... modules
#this will produce a lot of data so we should automate the testing process to reduce the amount of code. A selection of "interesting" modules could also be selected for analysis to further reduce the amount of information to interpret

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
probes = names(data_raw1)
probes2annot = match(probes, annotation_final$ensembl_gene_id)

# Get the corresponding ensembl_gene_id
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

#the above code just selects and makes a list of genes within the brown red and salmon modules via their ensmbl gene ids

# As background in the enrichment analysis, we will use all probes in the analysis.
fileName = paste("data-processed/ensemblgeneid-all.txt", sep="");
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

#####

#The user needs to provide the gene universe, GO annotations and either a
# criteria for selecting interesting genes (e.g. differentially expressed genes) from the gene universe or a score
# associated with each gene.


data(ALL)


# first we need the GO terms and their respective gene Ids



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