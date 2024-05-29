#' Functions to handle imputation, and scaling of metabolome data, while handling
#' feature metadata
library(tidyverse)
library(POMA)


filter_rownames <- function(df, filt_vector) {
  # wrapper script around filter to handle rownames
  df_filt <- df %>% 
    rownames_to_column(var="id") %>%
    filter(id %in% filt_vector) %>%
    column_to_rownames(var="id")
  return(df_filt)
}

metabo_list <- function(df, meta, sampleid_col) {
  sample_names <- meta %>% pull(.data[[sampleid_col]])
  
  features <- df %>%
    select(all_of(sample_names))
  
  feat_meta <- df %>%
    select(!all_of(sample_names))
  
  return_list <-
    list(
      "features" = features,
      "feature.metadata" = feat_meta,
      "sample.metadata" = meta
    )
  return(return_list)
}

normalise_sample_weight <-
  function(feat, metadata, weight_col, factor = 1) {
    weights <- metadata %>%
      pull(.data[[weight_col]])
    names(weights) <- metadata[["SampleID"]]
    
    if (ncol(feat) != length(weights)) {
      stop("Number of weights must be equal to the number of samples")
    }
    if (sum(names(weights) %in% colnames(feat)) != length(weights)) {
      stop("Sample names in the weight vector must match sample names in the data")
    }
    
    normed <-
      as.matrix(sapply(names(weights), function(w)
        feat[, w] / weights[w] * factor))
    colnames(normed) <- names(weights)
    rownames(normed) <- rownames(feat)
    return(feat)
  }

half_min_imp <- function(metabo_list) {
  # transpose
  df_t <- as.data.frame(t(metabo_list$features))
  
  # colwise min, fill NA with 0.5 * min
  half_min <-
    function(x)
      replace(x, is.na(x), min(x, na.rm = T) * 0.5)
  df_imp <- replace(df_t, TRUE, lapply(df_t, half_min))
  
  # re-tranpose and add back to metabolist
  metabo_list$features <- as.data.frame(t(df_imp))
  
  return(metabo_list)
}
#######################################

# adjust so you can specify number of samples as well as %
filter_prevalence <-
  function(metabo_list,
           filt_type = "prev_trh",
           prev_trh = 0.5,
           max_na = 1) {
    if (filt_type == "prev_trh") {
      feat_filt <-
        metabo_list$features[which(rowMeans(is.na(metabo_list$features)) < prev_trh),]
    }
    else if (filt_type == "min_samples") {
      if (max_na == 0) {
        stop("No value specified for max na")
      }
      else {
        min_samples <- dim(metabo_list$features)[2] - max_na
        feat_filt <-
          metabo_list$features[which(rowSums(is.na(metabo_list$features)) < min_samples),]
      }
    }
    
    feat_meta_filt <-
      metabo_list$feature.metadata[rownames(feat_filt),]
    
    return_list <- list(
      "features" = feat_filt,
      "feature.metadata" = feat_meta_filt,
      "sample.metadata" = metabo_list$sample.metadata
    )
    return(return_list)
  }


detect_outliers <- function(metabo_list, ...) {
  # use POMA for this for now
  SE <- PomaSummarizedExperiment(target = metabo_list$sample.metadata, features = t(metabo_list$norm.features))
  outliers <- PomaOutliers(SE, do = "analyze", ...)
  
  return(outliers$outliers)
}


pqn <- function(metabo_list, method = "median") {
  feat <- metabo_list$features
  
  if (method == "mean") {
    ref <- as.numeric(rowMeans(feat))
  }
  if (method == "median") {
    ref <- as.numeric(apply(feat, 1, median))
  }
  
  # do the actual normalisation
  # replace with apply or map
  feat_normed <-
    apply(feat, 2,  function(x)
      as.numeric(x) / median(as.numeric(x) / ref))
  
  
  metabo_list$norm.features <- feat_normed
  # reset rownames from feature metadata
  rownames(metabo_list$norm.features) <- rownames(metabo_list$features)
  
  return(metabo_list)
}

scale_z_score <- function(feat,
                     samples_are_columns = T,
                     na.rm = TRUE) {
  if (samples_are_columns == T) {
    scaled_feat <-
      t(apply(feat, 1, function(x)
        (x - mean(x, na.rm = na.rm)) / sd(x, na.rm = na.rm)))
  }
  else {
    scaled_feat <-
      apply(data, 2, function(x)
        (x - mean(x, na.rm = na.rm)) / sd(x, na.rm = na.rm))
  }
  return(scaled_feat)
}
