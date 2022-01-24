library(biomaRt)

marts <- listMarts(host="https://plants.ensembl.org")
#found mart: plants_mart

datasets <- listDatasets(useMart("plants_mart", host="https://plants.ensembl.org"))
#found dataset: taestivum_eg_gene

mart.hs <- useMart("plants_mart", "taestivum_eg_gene", host="https://plants.ensembl.org")

attributes <- listAttributes(mart.hs)

filters <- listFilters(mart.hs)
#wheatexp_gene ensembl_gene_id

library(biomaRt)
mart.hs <- useMart("plants_mart", "taestivum_eg_gene", host="https://plants.ensembl.org")

#making a list of values for the search
data_raw0 = read.table("data-project/E-MTAB-4484-query-results.fpkms.tsv", header = T, fill = T)
IDs <- data_raw0[,1]

mart_finished <- getBM(attributes = c("ensembl_gene_id",
                                   "description",
                                   "external_gene_name",
                                   "uniprotswissprot",
                                   "pdb",
                                   "chromosome_name"), 
                    filters = "ensembl_gene_id", 
                    values = IDs,
                    mart = mart.hs)

"ensembl_gene_id"

att <- c("description", "source", "source_description",
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


library(biomaRt)
mart.hs <- useMart("plants_mart", "taestivum_eg_gene", host="https://plants.ensembl.org")

attributes <- listAttributes(mart.hs)
#making a list of values for the search
data_raw0 = read.table("data-project/E-MTAB-4484-query-results.fpkms.tsv", header = T, fill = T)
IDs <- data_raw0[,1]

annot <- data.frame(matrix(ncol=0,nrow=0))






for (i in att) {
  a <- get(i)
  for (j in all.files) {
    b <- get(j)
    d <- dataframe(mergefun(a, b))
    newname <- getBM(attributes = c("ensembl_gene_id",i),
                     filters = "ensembl_gene_id",
                     values = IDs,
                     mart = mart.hs)
    names(d) <- c("description", "source", "source_description",
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
    assign(newname,d)
    
  }  
}



library(magicfor)               
magic_for(getBM, silent = TRUE) 

for (i in att){
  
  test <- getBM(attributes = c("ensembl_gene_id",i),
                     filters = "ensembl_gene_id",
                     values = IDs,
                     mart = mart.hs)
  
  test <- annot
  ann_[i] <- annot
  
  if (nrow(test) = 0) {
    print(i)
    }
   }


nrow()







