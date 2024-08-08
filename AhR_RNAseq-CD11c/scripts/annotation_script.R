####R code for fetching gene names with biomaRt #####
require(biomaRt)
annot.table  <- data.frame() # This is your annotation table, it can't
be NULL
# Collect ensembl IDs from annot.table before converting to normal gene names.
ensembl_ids <- character() # Get this from may be from annot.table
# Prepare gene table with some simple caching to avoid stressing the
Ensembl server by many repeated runs
genes.table = NULL
if (!file.exists("cache.genes.table")) {
  message("Retrieving genes table from Ensembl...")
  mart <- useMart("ensembl")
  #listDatasets(mart=mart)
  mart <- useDataset("hsapiens_gene_ensembl", mart = mart)
  genes.table <- getBM(filters= "ensembl_gene_id",
                       attributes= c("ensembl_gene_id",
                                     "external_gene_id", "description"), values= ensembl_ids, mart= mart)
  save(genes.table, file= "cache.genes.table")
} else {
  load("cache.genes.table")
  message("Reading gene annotation information from cache file:
cache/cache.genes.table
            Remove the file if you want to force retrieving data from
Ensembl")
}

# Merging two tables, syntax showed here are in full forms, by.x and
by.y can be simplified
annot.table <- merge(x = annot.table, y = genes.table,
                     by.x = "ensembl_id", by.y =
                       "ensembl_gene_id", all.x = T, all.y = F )
# Do something with the table, i.e export to TSV