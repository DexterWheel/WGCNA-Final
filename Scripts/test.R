
library(biomaRt)
library(rlist)

#making a list of values for the search
lnames = load(file = "data-project/annotation.RData")
lnames = load(file = "data-processed/wheat-dataInput.RData")

data_raw2 = as.data.frame(t(data_raw1))
IDs = rownames(data_raw2)

mart.hs = useMart("plants_mart", "taestivum_eg_gene", host="https://plants.ensembl.org")

test = getBM(attributes = c("ensembl_gene_id", "transcript_biotype"),
               filters = "ensembl_gene_id",
               values = IDs,
               mart = mart.hs)
