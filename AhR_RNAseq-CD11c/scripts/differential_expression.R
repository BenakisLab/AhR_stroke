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
  
  go <- groupGO(gene = gene_col, OrgDb = OrgDb, ont = ont, level = level, readable = T )
  
  go_df <- as.data.frame(go@result)
  
  return(go_df)
}

create_metadata_df <- function(vars, group_sizes, column_names_df) {
  
  sample_col <- data.frame(sample = rep(vars), group_sizes)
  
  row(sample_col) <- colnames(column_names_df)
  
  return(sample_col)
}