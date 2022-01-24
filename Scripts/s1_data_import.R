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
data_processed = as.data.frame( t(data_raw0[, -c(1:2) ]) )
#naming columns after the genes they represent
names(data_processed) = data_raw0$Gene_ID
#double checking that the names have been properly changed
data_processed[1:15, 1:2]


#==== 
# S1.2 : Trimming data
#====

#####
# Need to convert some of the Ids so that they are included in annotations
#####

# we need to remove genes and samples with too many missing values

gsg = goodSamplesGenes(data_processed, verbose = 3)
gsg$allOK

#allok = False and so samples need to be removed

if (!gsg$allOK){
  #Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(data_processed)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(data_processed)[!gsg$goodSamples], collapse = ", ")))
  #Remove the offending genes and samples from the data:
  data_processed = data_processed[gsg$goodSamples, gsg$goodGenes]
}

#running the code again shows that all undesirable values have been removed
gsg = goodSamplesGenes(data_processed, verbose = 3)
gsg$allOK


### Next we cluster the samples to see if there are any obvious outliers. 
#in contrast to clustering genes that will come later)

sampleTree = hclust(dist(data_processed), method = "average")

plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

#The plot suggests that there are two major clusters. Apart from Stem BBCH 55 all tissues appear to cluster together. 
#It has been decided to retain Stem BBCH 55 as there may be a developmental reason for the unexpected difference in clustering

#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9)

#====
# S1.3 : Trait data
#====

library(stringr)# allows characters to be extracted based on user set conditions

#two versions of the trait data are to be used.
#One which separated by just tissue with dev stage as one column
#and one which separates by both tissue and developmental stage

trait_list = c("data-project/trait-data-tissue-split.txt", "data-project/trait-data-tissue-dev-split-.txt")

for (i in trait_list){
  
  temp = print(i)
  
  newname = str_extract(temp, "tissue\\s*(.*?)\\s*split")
  
  traitData = read.table(i, header = T)
  
  dim(traitData); names(traitData)
  
  samples = rownames(data_processed)
  traitRows = match(samples, traitData$ref)
  datTraits = traitData[traitRows,]
  datTraits = subset(datTraits, select=-(ref))
  
  #renaming the rows 
  rownames(datTraits) = traitData[traitRows, "ref"]
  
  # Re-cluster samples
  sampleTree2 = hclust(dist(data_processed), method = "average")
  
  # Convert traits to a color representation
  traitColors = numbers2colors(datTraits, signed = T)
  
  # Plot the sample dendrogram and the colors underneath.
  plotDendroAndColors(sampleTree2, traitColors,
                      groupLabels = names(datTraits),
                      main = "Sample dendrogram and trait heatmap")
  
  assign(newname, datTraits)
  }


datTraits_dev = `tissue-dev-split`

datTraits_tissue = `tissue-split`
# trait data visualises the findings mentioned above

###save expression and trait data for use in the next steps

#datTraits is also saved in case I decide to just use one of the trait files

save(data_processed, datTraits, datTraits_tissue, datTraits_dev, file = "data-processed/wheat-dataInput.RData")

lnames = load(file = "data-processed/wheat-dataInput.RData")
