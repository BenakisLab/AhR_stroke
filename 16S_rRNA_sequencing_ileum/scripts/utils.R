library(rstatix)

pairwise_wilcox_ps <- function(ps, factor, case, control, p_adjust="BH"){
  
  featmat <- ps_to_asvtab(ps) %>% 
    as.matrix()
  meta <- meta_to_df(ps) %>% 
    rownames_to_column("SampleID")
  # initialise pval matrix 
  p.val <- matrix(NA, nrow=nrow(featmat), ncol=1, 
                  dimnames=list(row.names(featmat), "p.val"))
  
  for (row in rownames(featmat)){
    x <- featmat[row, meta %>% 
                   filter(.data[[ factor ]] ==control) %>% pull(SampleID)]
    y <- featmat[row, meta %>% 
                   filter(.data[[ factor ]] ==case) %>% pull(SampleID)]
    
    
    p.val[row, ] <- wilcox.test(x, y, exact=FALSE)$p.value
  }
  
  p.val.adj <- p.val %>% 
    as.data.frame() %>% 
    adjust_pvalue("p.val", method = p_adjust)
  
  return(p.val.adj)
}

filter_rownames <- function(df, filt_vector) {
  # wrapper script around filter to handle rownames
  df_filt <- df %>% 
    rownames_to_column(var="id") %>%
    filter(id %in% filt_vector) %>%
    column_to_rownames(var="id")
  return(df_filt)
}  

filter_na <- function(df, expr){
  df %>% filter({{expr}} %>% replace_na(T))
}

t_df <- function(x) {
  return(as.data.frame(t(x)))
}