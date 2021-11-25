library(biomaRt)
library(rlist)
library(tidyverse)

#====
# Data extraction
#====
marts = listMarts(host="https://plants.ensembl.org")
#found mart: plants_mart

datasets = listDatasets(useMart("plants_mart", host="https://plants.ensembl.org"))
#found dataset: taestivum_eg_gene

mart.hs = useMart("plants_mart", "taestivum_eg_gene", host="https://plants.ensembl.org")

filters = listFilters(mart.hs)
#wheatexp_gene ensembl_gene_id

attributes = listAttributes(mart.hs)
#####
att = c("description", "source", "source_description",
         "name_1006", "definition_1006", "namespace_1003",
         "interpro", "interpro_description", "interpro_short_description",
         "external_gene_name", "chromosome_name", "go_id",
         "pdb", "uniparc", "uniprotsptrembl", "uniprotswissprot",
         "wheatexp_gene", "wheatexp_trans", "rfam", "sfld",
         "embl", "ensembl_peptide_id", "ensembl_exon_id",
         "ensembl_transcript_id", "entrezgene_accession",
         "entrezgene_description", "entrezgene_id",
         "external_gene_source", "external_synonym",
         "gene_biotype", "kegg_enzyme", "peptide",
         "peptide_location", "pfam", "pirsf", "refseq_dna",
         "refseq_peptide", "plant_reactome_pathway",
         "plant_reactome_reaction", "sift_prediction_2076",
         "sift_score_2076", "superfamily", "tigrfam", "transcript_biotype")
##### 

mart.hs = useMart("plants_mart", "taestivum_eg_gene", host="https://plants.ensembl.org")

attributes = listAttributes(mart.hs)
#making a list of values for the search
lnames = load(file = "data-project/annotation.RData")
lnames = load(file = "data-processed/wheat-dataInput.RData")

data_raw2 = as.data.frame(t(data_raw1))

IDs = rownames(data_raw2)

annot_list = list()

for (i in att){
  
  #searching for data for each attribute sequentially
  test = getBM(attributes = c("ensembl_gene_id", i),
               filters = "ensembl_gene_id",
               values = IDs,
               mart = mart.hs)
  
  #traitRows = match(IDs, test$ensembl_gene_id)
  
  #test = test[traitRows,]
  
  list = list(test)
  
  #appending the dataframe to an existing list of dataframes
  annot_list = c(annot_list, list)
  }

#naming the data frames after their attributes
names(annot_list) = att

save(annot_list, file = "data-project/annotation2.RData")

#====
# tidying data
#====

lnames = load(file = "data-project/annotation.RData")
lnames = load(file = "data-processed/wheat-dataInput.RData")

data_raw2 = as.data.frame(t(data_raw1))

samples = rownames(data_raw2)

names(annot_list)

annot_list2 = list()

for (names in annot_list){
  
  i = print(names)
  
  traitRows = match(samples, i$ensembl_gene_id)
  
  a = names[traitRows,]
  
  list = list(a)
  
  #appending the dataframe to an existing list of dataframes
  annot_list2 = c(annot_list2, list)
  
  }


annotation_final <- annot_list2 %>% reduce(left_join, by = "ensembl_gene_id")

save(annotation_final, file = "data-project/annotation-final.RData")



