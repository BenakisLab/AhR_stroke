
# shotgun utils
library(phyloseq)
library(tidyverse)


common_cols <- function(list_df, n=length(list_df)) {
  col_list <- map(list_df,names)
  
  if (n < length(list_df)){
    n_cols <- table(unlist(col_list))
    return(names(n_cols[n_cols >= n]))
  }
  
  return(Reduce(intersect,col_list))
}


combine_counts <- function(counts_list, taxa_are_cols=T, prev_thresh=0.3) {
  if (taxa_are_cols == F){
    counts_list <- map(counts_list, function(x) t(x) %>% as.data.frame)
  }
  # purely presence absence
  shared_taxa <- common_cols(counts_list, n=round(length(counts_list) * prev_thresh))
  
  combined_counts <- bind_rows(counts_list) %>% 
    select(all_of(shared_taxa))
  return(combined_counts)
}

combine_meta <- function(meta_list, col_list) {
  
  combined_meta <- bind_rows(meta_list) %>% 
    select(all_of(col_list))
}

return_level <- function(level) {
  if (level %in% c("Species", "species", "s")){
    return("s")
  } else if (level %in% c("Genus", "genus", "g")){
    return("g")
  } else if (level %in% c("Family", "family", "f")){
    return("f")
  } else if (level %in% c("Order", "order", "o")){
    return("o")
  } else if (level %in% c("Class", "class", "c")){
    return("c")
  } else if (level %in% c("Phylum", "phylum", "p")){
    return("p")
  } else if (level %in% c("Kingdom", "kingdom", "k")){
    return("k")
  } else {
    stop("Not a valid taxonomic rank")
  }
}

select_rank <- function(merged_table, level) {
  
  level <- return_level(level)
  
  tax_table <-  merged_table %>% 
    rownames_to_column("taxa") %>% 
    dplyr::select(taxa)
  
  tax_table_sep <- tax_table %>% 
    rowwise() %>% 
    separate(taxa, into = c("Kingdom", "Phylum", "Class", "Order", 
                            "Family","Genus", "Species", 
                            "Subspecies"), remove = F, sep = "\\|")
  tax_table_fill <- t(zoo::na.locf(t(tax_table_sep))) %>%
    as.data.frame()
  
  tax_table_out <- tax_table_fill %>% 
    rowwise() %>% 
    separate(Subspecies, into = c("Level", "tmp"), sep = "_") %>% 
    select(-tmp) %>% 
    filter(Level == level) %>% 
    column_to_rownames("taxa") %>% 
    select(-Level)
  
  merged_table_out <- merged_table %>% 
    rownames_to_column("taxa") %>% 
    dplyr::filter(taxa %in% rownames(tax_table_out)) %>% 
    column_to_rownames("taxa")
  
  return_list <- list(counts = merged_table_out, tax_table = tax_table_out)
  return(return_list)
}

import_pseq_metag <- function(merged_table_path, metapath, level) {
  # read files
  merged_table <- read_tab_delim_metag(merged_table_path)
  metadata <- read_tab_delim_metag(metapath)
  
  merged_table_rank <- select_rank(merged_table, level=level)
  
  out <- load_phylo(merged_table_rank$counts, merged_table_rank$tax_table, metadata)
  
  return(out)
}


