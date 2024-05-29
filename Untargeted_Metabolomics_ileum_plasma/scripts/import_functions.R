#' Functions to handle imputation, and scaling of metabolome data, while handling
#' feature metadata

metabo_list <- function(df, meta, sampleid_col){
  
  sample_names <- meta %>% pull(.data[[sampleid_col]])
  
  features <- df %>% 
    select(all_of(sample_names))
  
  feat_meta <- df %>% 
    select(!all_of(sample_names))
  
  return_list <- list("Features" = features, "Feature_Metadata" = feat_meta, "Sample_Metadata" = meta)
  return(return_list)
}

normalise_sample_weight <- function(metabo_list, weight_col, factor=1) {
  
  weights <- metabo_list$Sample_Metadata %>% 
    pull(.data[[ weight_col ]])
  names(weights) <- metabo_list$Sample_Metadata[["SampleID"]]
  
  if (ncol(metabo_list$Features) != length(weights)) {
    stop("Number of weights must be equal to the number of samples")
  } 
  if (sum(names(weights) %in% colnames(metabo_list$Features)) != length(weights)) {
    stop("Sample names in the weight vector must match sample names in the data")
  }
  
  normed <- as.matrix(sapply(names(weights), function(w) metabo_list$Features[,w]/weights[w]*factor))
  colnames(normed) <- names(weights)
  metabo_list$Weight_normed_feat <- normed
  return(metabo_list)
}

half_min_imp <- function(metabo_list) {
  
  # transpose
  df_t <- as.data.frame(t(metabo_list$Features))
  
  # colwise min, fill NA with 0.5 * min 
  half_min <- function(x) replace(x, is.na(x), min(x, na.rm = T) * 0.5)
  df_imp <- replace(df_t, TRUE, lapply(df_t, half_min))
  
  # re-tranpose and add back to metabolist
  metabo_list$Features <- as.data.frame(t(df_imp))
  
  return(metabo_list)
}
#######################################

# adjust so you can specify number of samples as well as % 
filter_prevalence <- function(metabo_list, filt_type="prev_trh", prev_trh=0.5, min_samples=0) {
  
  if (filt_type == "prev_trh"){ 
    feat_filt <- metabo_list$Features[which(rowMeans(is.na(metabo_list$Features)) <= prev_trh), ]
  }
  else if (filt_type == "min_samples"){
    
    if (min_samples == 0){
      stop("No value specified for minimum samples")
    } 
    else {
      feat_filt <- metabo_list$Features[which(rowSums(is.na(metabo_list$Features)) <= min_samples), ]
    }
  }
  
  feat_meta_filt <- metabo_list$Feature_Metadata[rownames(feat_filt), ]
  
  return_list <- list("Features" = feat_filt, 
                      "Feature_Metadata" = feat_meta_filt,
                      "Sample_Metadata" = metabo_list$Sample_Metadata)
  return(return_list)
}