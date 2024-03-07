library(parallel)
library(iterators)
library(foreach)
library(doParallel)
library(MASS)
library(biomaRt)

# add and set directory of the folder where to save all the update information
# make sure to first copy the genes_and_ids_all.csv file in this directory
directory <- ''
setwd(directory)
df <- read.csv( 'genes_and_ids_all.csv', header=T)
genes_ids <- unique(df$id)
splits <- split(genes_ids, ceiling(seq_along(genes_ids)/100))

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
gene_ids_all <- genes_ids 

goids = getBM(attributes = c('entrezgene_id', 'go_id'), 
              filters = 'entrezgene_id', 
              values = gene_ids_all[1:30], 
              mart = ensembl)
goids

# Parallel implementation
numCores <- detectCores()
registerDoParallel(numCores)

system.time({
  res<- foreach::foreach(i=1:length(splits), .combine=rbind) %dopar%{
    getBM(attributes = c('entrezgene_id', 'go_id'), 
          filters = 'entrezgene_id', 
          values = splits[[i]], 
          mart = ensembl) }
})

length(unique(res$entrezgene_id))
sum(is.na(res$go_id))
# remove NAs
res_new<- res[!is.na(res$entrezgene_id), ]
res_new <- res_new[!is.na(res_new$go_id), ]
res_new <- res_new[res_new$go_id != "", ]
length(unique(res_new$entrezgene_id)) 
write.table(res_new, 'go_annotations_new.csv', row.names = F, col.names = F, sep=",", quote = F)




