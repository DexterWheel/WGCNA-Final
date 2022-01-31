#### Important ####
#this code does not account for the fact that some attributes such as GO ids can have multiple matches to the same gene

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

data_raw2 = as.data.frame(t(data_processed))

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

save(annot_list, file = "data-project/annotation.RData")

#====
# tidying data
#====

lnames = load(file = "data-processed/annotation.RData")

#first we need to select the attributes we want to keep, due to issues with biomart we can no longer access it, so we must instead trim the data we have already recieved

att = c("go_id", "description", "gene_biotype", "transcript_biotype", "name_1006", "interpro", "interpro_short_description",
        "chromosome_name",  "pdb", "uniparc", "uniprotsptrembl", "uniprotswissprot",
        "external_gene_source",  "superfamily", "source", "source_description")

annot_list2 = list()

for (i in att){
  
  name = print(i)
  temp = annot_list[[i]]
  
  temp_list= list(temp)
  
  annot_list2 = c(annot_list2,temp_list)
}

names(annot_list2) = att
#we then need to consolidate duplicate rows together
dfs = names(annot_list2)

annot_list3 = list()

for (i in dfs){
  
  x = annot_list2[[i]]
  
  y = annot_list2[[i]][[2]]
  
  columns = colnames(x)
  
  attribute = columns[2]
  
  temp = aggregate(y ~ ensembl_gene_id, data = x, paste, collapse = ",")
  
  names(temp)[names(temp) == "y"] = attribute
  
  lists = list(temp)
  
  #appending the dataframe to an existing list of dataframes
  annot_list3 = c(annot_list3, lists)
  }

names(annot_list3) = att
annot_list4 = list()

for (i in dfs){
  
  x = annot_list2[[i]]
  
  y = annot_list2[[i]][[2]]
  
  columns = colnames(x)
  
  attribute = columns[2]
  
  temp = aggregate(y ~ ensembl_gene_id, data = x, paste)
  
  names(temp)[names(temp) == "y"] = attribute
  
  lists = list(temp)
  
  #appending the dataframe to an existing list of dataframes
  annot_list4 = c(annot_list3, lists)
}

names(annot_list4) = att

go = annot_list4[[1]]


#loading in input data to use the ids as reference
lnames = load(file = "data-processed/wheat-dataInput.RData")

data_raw2 = as.data.frame(t(data_processed))

samples = rownames(data_raw2)

#making an empty list

annot_list4 = list()

for (names in annot_list3){
  
  i = names
  
  traitRows = match(samples, names$ensembl_gene_id)
  
  a = names[traitRows,]
  
  list = list(a)
  
  #appending the dataframe to an existing list of dataframes
  annot_list4 = c(annot_list4, list)
  
  }

names(annot_list4) = att

annotation_final = annot_list4 %>% reduce(left_join, by = "ensembl_gene_id")

save(annotation_final, go, file = "data-project/annotation-final.RData")


