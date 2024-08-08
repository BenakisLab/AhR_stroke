pacman::p_load("tximport", "tximportData", "readr", "DESeq2", "GenomicFeatures", 
               "tximeta", "org.Mm.eg.db", "AnnotationDbi","tidyverse",
               "clusterProfiler", "pheatmap", "SOAR", "ComplexHeatmap")

# functions 
annotate_degs <- function(db, deg, keys, keytype, multivals) {
  deg_annot <- deg
  
  deg_annot[["gene_symbol"]] <- mapIds(db, keys = keys, column = "SYMBOL", keytype = keytype, multiVals = multivals)
  
  deg_annot[["entrez"]] <- mapIds(db, keys = keys, column = "ENTREZID", keytype = keytype, multiVals = multivals)
  
  deg_annot[["gene name"]] <- mapIds(db, keys = keys, column = "GENENAME", keytype = keytype, multiVals = multivals)
  
  return(deg_annot)
}

order_extract_significant <- function(res_deseq2){
  
  df <- as.data.frame(res_deseq2)
  
  df_out <- filter(df, ) # use tidyeval maybe 
}

# return 
subset_genes <- function(df, term, description_col, gene_id_col) {
  
  formatted_df <- df %>% filter(grepl(term, {{ description_col }} ))
  
  formatted_df[[gene_id_col]] <- gsub("/", ",", as.character(formatted_df[[gene_id_col]]))
  
  genes <- filter(formatted_df, Count > 0) %>%  select(gene_id_col)
  
  unique_genes <- unique(unlist(strsplit(genes[[gene_id_col]], ",")))
  
  unique_genes_df <- data.frame("Gene_list"  = unique_genes)
  
  return_list <- list("gene_list" = unique_genes, "gene_df" = unique_genes_df)
  return(return_list)
}

cal_z_score <- function(mat){
  t(scale(t(mat)))
}



format_matrix <- function(df, subset, id_col, col_order, drop_cols) {
  if (missing(subset)){
    mat <- df
  }
  else {
    mat <- df[df[[id_col]] %in% subset, ] 
  }
  
  rownames(mat) <- mat[[id_col]]
  # select numeric columns - has to be a nicer way of doing this 
  mat <- data.matrix(mat)
  
  mat <- mat[ , !apply(is.na(mat), 2, all)]
  
  if (exists("col_order") == T) {
    mat <- mat[, col_order]
  }
  
  if (exists("drop_cols") == T) {
    mat <- mat[,-drop_cols]
  }
  
  return(mat)
}

getGO <- function(gene_col, OrgDb, ont, level, df_out) {
  
  go <- groupGO(gene = gene_col, OrgDb = OrgDb, ont = ont, level = level, 
                readable = T)
  
  go_df <- as.data.frame(go@result)
  
  return(go)
}

create_metadata_df <- function(vars, group_sizes, column_names_df) {
  
  sample_col <- data.frame(sample = rep(vars), group_sizes)
  
  row(sample_col) <- colnames(column_names_df)
  
  return(sample_col)
}

metadata <- read.table("../RNAseq_AV_Mice.txt", header=T, row.names = 1)
# 
quant_path <- "../quants"
files <- file.path(quant_path, rownames(metadata), "quant.sf")
names(files) <- rownames(metadata)

#inspecting files
files
#
#txdb <- makeTxDbFromEnsembl(organism = "Mus musculus")
#saveRDS(txdb, file= "txdb.rds")
txdb <- readRDS("txdb.rds")
#k <- keys(txdb, keytype = "TXNAME")
#saveRDS(k, file = "k.rds")
k - readRDS("k.rds")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME") #try ignore txversion 
saveRDS(tx2gene, file = "tx2gene.rds")
#tx2gene <- readRDS("tx2gene.rds")

#this function will output gene-level matrices. We can avoid gene-level summarization by setting txOut=TRUE
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = T)
saveRDS(txi, file = "txi.rds")
# txi <- readRDS("txi.rds")
# create deseq object
dds <- DESeqDataSetFromTximport(txi, colData = metadata, design = ~ Genotype)

# filter to remove rows with very few reads - not strictly necessary, see here: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#pre-filtering
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# check factor levels are in the correct order - actually not required since genotypes alphabetical but good to have for future use
dds$Genotype 
dds$Genotype <- relevel(dds$Genotype, ref = "fl_fl")

# differential expression
dds <- DESeq(dds)
saveRDS(dds, "dds.RDS")
#dds <- readRDS("dds.RDS")

res_tgtg <- results(dds, contrast = c("Genotype", "tg_tg", "fl_fl"))
res_tgtg.annot <- annotate_degs(org.Mm.eg.db, res_tgtg, keys = row.names(res_tgtg), keytype = "ENSEMBL", multivals = "first")

res_tgwt <- results(dds, contrast = c("Genotype", "tg_wt", "fl_fl"))
res_tgwt.annot <- annotate_degs(org.Mm.eg.db, res_tgwt, keys = row.names(res_tgtg), keytype = "ENSEMBL", multivals = "first")

res_atf6 <- results(dds, contrast = c("Genotype", "tg_tg", "tg_wt"))
res_atf6.annot <- annotate_degs(org.Mm.eg.db, res_atf6, keys = row.names(res_tgtg), keytype = "ENSEMBL", multivals = "first")

# filter by p-val and fold change
padj.cutoff <- 0.05
lfc.cutoff <- 1

sig_tgtg <- as.data.frame(res_tgtg.annot) 
sig_tgtg <- filter(sig_tgtg, padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff) %>% arrange(desc(log2FoldChange))

sig_tgwt <- as.data.frame(res_tgwt.annot)
sig_tgwt <- dplyr::filter(sig_tgwt, padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)

sig_atf6 <- as.data.frame(res_atf6.annot)
sig_atf6 <- dplyr::filter(sig_atf6, padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)

# clusterprofiler - GO classification
tgtg.BP.4 <- getGO(gene_col = sig_tgtg$entrez, OrgDb = org.Mm.eg.db, ont = "BP", level = 4)
tgtg.BP.6 <- getGO(gene_col = sig_tgtg$entrez, OrgDb = org.Mm.eg.db, ont = "BP", level = 6)

tgwt.BP.4 <- getGO(gene_col = sig_tgwt$entrez, OrgDb = org.Mm.eg.db, ont = "BP", level = 4)
tgwt.BP.6 <- getGO(gene_col = sig_tgwt$entrez, OrgDb = org.Mm.eg.db, ont = "BP", level = 6)

atf6.BP.4 <- getGO(gene_col = sig_atf6$entrez, OrgDb = org.Mm.eg.db, ont = "BP", level = 4)
atf6.BP.6 <- getGO(gene_col = sig_atf6$entrez, OrgDb = org.Mm.eg.db, ont = "BP", level = 6)



# this functions needs df so make sure to access result slot in getGO class
aa_genes_tgtg <- subset_genes(tgtg.BP.4@result, term = "amino acid", description_col = Description, gene_id_col = "geneID")
aa_genes_tgwt <- subset_genes(tgwt.BP.4@result, term = "amino acid", description_col = Description, gene_id_col = "geneID")
aa_genes_atf6 <- subset_genes(atf6.BP.4@result, term = "amino acid", description_col = Description, gene_id_col = "geneID")

# fa gene specific 
peptide_genes_tgtg <- subset_genes(tgtg.BP.6@result, term = "peptide", description_col = Description, gene_id_col = "geneID")
peptide_genes_tgwt <- subset_genes(tgwt.BP.6@result, term = "peptide", description_col = Description, gene_id_col = "geneID")
peptide_genes_atf6 <- subset_genes(atf6.BP.6@result, term = "peptide", description_col = Description, gene_id_col = "geneID")

# normalized read counts - abstract some of this into function 

norm_counts <- DESeq2::counts(dds, normalized=T)
norm_counts.annot <- merge(as.data.frame(res_tgtg.annot[7]), as.data.frame(norm_counts), by="row.names", sort=FALSE) 
norm_counts.annot$gene_symbol %>% replace_na("unknown")
rownames(norm_counts.annot) <- norm_counts.annot[,1]
#write.table(norm_counts.annot, "normalised_counts.tab", sep = "\t", quote = F, row.names = F)

# prepare data for plotting 
# drop F12356 from datafram and reorder columns 
col.order <- c("F12361",	"F12362",	"F12363",	"F12364",	"F12365", "F12366", 
               "F12367",	"F12368", "F12369", "F12370", "F12371", "F12372",
               "F12354", "F12355", "F12356", "F12357",	"F12358",	"F12359", "F12360")

# format matrices for heatmap

mat_tgtg.aa <- format_matrix(norm_counts.annot, subset = aa_genes_tgtg$gene_list, id_col = "gene_symbol", col_order = col.order, drop_cols = 15)
mat_tgwt.aa <- format_matrix(norm_counts.annot, subset = aa_genes_tgwt$gene_list, id_col = "gene_symbol", col_order = col.order, drop_cols = 15)
mat_atf6.aa <- format_matrix(norm_counts.annot, subset = aa_genes_atf6$gene_list, id_col = "gene_symbol", col_order = col.order, drop_cols = 15)


mat_tgtg.pep <- format_matrix(norm_counts.annot, subset = peptide_genes_tgtg$gene_list, id_col = "gene_symbol", col_order = col.order, drop_cols = 15)
mat_tgwt.pep <- format_matrix(norm_counts.annot, subset = peptide_genes_tgwt$gene_list, id_col = "gene_symbol", col_order = col.order, drop_cols = 15)
mat_atf6.pep <- format_matrix(norm_counts.annot, subset = peptide_genes_atf6$gene_list, id_col = "gene_symbol", col_order = col.order, drop_cols = 15)
# create names list of matrices and apply z score function returning list  
mat_list <- list("mat_atf6_aa_zscore" = mat_atf6.aa,  "mat_tgtg_aa_zscore" = mat_tgtg.aa, "mat_tgwt_aa_zscore" = mat_tgwt.aa,
                 "mat_atf6_pep_zscore" = mat_atf6.pep,  "mat_tgtg_pep_zscore" = mat_tgtg.pep, "mat_tgwt_pep_zscore" = mat_tgwt.pep)
mat_zscores <- lapply(mat_list, cal_z_score)

# generating metadata could also be moved into function 
sample_col <- data.frame(sample = rep(c("fl/fl", "tg/wt", "tg/tg"), c(6,6, 6)))

# remove outlier sample here
#row.names(sample_col) <- colnames(matrix_in_outlier_rem)



# heatmaps 

colours <- list("sample"=c("fl/fl"="royalblue","tg/wt"="gold", "tg/tg"="red2"))
colAnn <- HeatmapAnnotation(df=sample_col, which="col", col=colours)

# amino acids 
ComplexHeatmap::Heatmap(mat_zscores$mat_atf6_aa_zscore, cluster_columns = F, cluster_rows = T, top_annotation=colAnn)


ComplexHeatmap::Heatmap(mat_zscores$mat_tgtg_aa_zscore, cluster_columns = F, cluster_rows = T, top_annotation=colAnn)
# repeat for tgtg - when you have the time rewrite this as a function  

ComplexHeatmap::Heatmap(mat_zscores$mat_tgwt_aa_zscore, cluster_columns = F, cluster_rows = T, top_annotation=colAnn)

# peptides 
ComplexHeatmap::Heatmap(mat_zscores$mat_atf6_pep_zscore, cluster_columns = F, cluster_rows = T, top_annotation=colAnn)


ComplexHeatmap::Heatmap(mat_zscores$mat_tgtg_pep_zscore, cluster_columns = F, cluster_rows = T, top_annotation=colAnn)
# repeat for tgtg - when you have the time rewrite this as a function  

ComplexHeatmap::Heatmap(mat_zscores$mat_tgwt_pep_zscore, cluster_columns = F, cluster_rows = T, top_annotation=colAnn)

