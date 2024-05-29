library(tidyverse)
library(coin)
library(docstring)
library(roxygen2)
library(EnhancedVolcano)

get_padj <- function(data,
                     groups,
                     feature.orient = 1,
                     test.method = t.test,
                     p.adjust.method = "BH") {
  raw_p <- apply(data, feature.orient,
                 function(x)
                   tryCatch(
                     test.method(x ~ groups)$p.value,
                     error = function(x)
                       return(1)
                   ))
  adj_p <- p.adjust(raw_p, p.adjust.method)
  names(adj_p) <- dimnames(data)[[feature.orient]]
  
  return(adj_p)
}

# FIX ISSUE WITH REVERSE GROUPS
fold.change <- function(metabo_list, group_col, group.order = NULL, log.data=TRUE) {
  # get groups
  groups <- metabo_list$sample.metadata[[group_col]]
  names(groups) <- metabo_list$sample.metadata[[1]]
  
  if (is.null(group.order)) {
    uni.groups <- sort(unique(groups))
  } else {
    uni.groups <- group.order
  }
  n.groups <- length(uni.groups)
  
  
  # NOTE: only geometric mean if log-data
  if (log.data == TRUE){
    # extract data
    data <- list(FC = metabo_list$log.features,
                 pval = metabo_list$scaled.features)
    
    # calculate means on log data
    means <- sapply(uni.groups,
                         function(g) {
                           sub <- data$FC[, names(groups[groups == g])]
                           if (is.null(dim(sub))) {
                             warning(paste0("Only one sample in group ", g, "!"))
                             return(sub)
                           }
                           rowMeans(sub, na.rm = TRUE)
                         })
    # calc FC using log2 mean(x) - log2 mean(y)
    comps <- sapply(1:(n.groups - 1),
                    function(i) {
                      sapply((i + 1):n.groups,
                             function(j)
                               return(means[, i] - means[, j]))
                    })
  } else {
    # extract data
    data <- list(FC = metabo_list$norm.features,
                 pval = metabo_list$scaled.features)
    # calculate means on normalised data
    means <- sapply(uni.groups,
                    function(g) {
                      sub <- data$FC[, names(groups[groups == g])]
                      if (is.null(dim(sub))) {
                        warning(paste0("Only one sample in group ", g, "!"))
                        return(sub)
                      }
                      rowMeans(sub, na.rm = TRUE)
                    })
    # calc FC using log2(mean(x) / mean(y))
    comps <- sapply(1:(n.groups - 1),
                    function(i) {
                      sapply((i + 1):n.groups,
                             function(j)
                               return(log2(means[, i] / means[, j])))
                      
                    })
    
   }

  # else if not log transformed do log2(x / y)
  p.adj <- get_padj(data$pval,
                    groups,
                    test.method = t.test,
                    p.adjust.method = "fdr")
  
  res <- cbind(comps, p.adj) %>% 
    as.data.frame()
  colnames(res) <- c("Log2FC", "p.adj")
  
  res_annot <- merge(res, metabo_list$feature.metadata, by=0) %>% 
    column_to_rownames("Row.names")
  
  return(res_annot)
}





plot_volcano <- function(res, groups, names, pthresh=0.05, FCthresh=0.58, ...) {
  
  keyvals <- ifelse(
    res$Log2FC < -FCthresh & res$p.adj < pthresh, "#016FB9",
    ifelse(res$Log2FC > FCthresh & res$p.adj < pthresh, "#FB3640",
           'grey'))
  keyvals[is.na(keyvals)] <- 'grey'
  names(keyvals)[keyvals == '#016FB9'] <- 'Decreased'
  names(keyvals)[keyvals == 'grey'] <- 'Unchanged'
  names(keyvals)[keyvals == '#FB3640'] <- 'Increased'
  
  
  #xlimits <- c(ceiling(min(res$Log2FC)), ceiling(max(res$Log2FC))) 
  #ylimits <- c(ceiling(min(-log10(res$p.adj))), ceiling(max(-log10(res$p.adj))) )
  
  EnhancedVolcano(res, x="Log2FC", y="p.adj", 
                  lab=names,
                  labSize = 4,
                  colCustom = keyvals,
                  colAlpha=0.7, 
                  caption=NULL, 
                  gridlines.major = F,
                  gridlines.minor = F,
                  pCutoff = pthresh, 
                  FCcutoff = FCthresh, 
                  title = paste(groups[1], 
                                "versus", 
                                groups[2], 
                                sep=" "), 
                  subtitle = NULL,...)
}

subpathway_enrichment <- function(fc_res, thresh) {
  subpathways <- unlist(unique(fc_res["SUB PATHWAY"]))
  k <- map(subpathways, function(x) fc_res %>% 
             filter(`SUB PATHWAY` == x) %>% 
             filter(p.adj <=  thresh) %>% 
             nrow())
  names(k) <- subpathways
  
  m <- map(subpathways, function(x) fc_res %>% 
             filter(`SUB PATHWAY` == x) %>% 
             nrow())
  names(m) <- subpathways
  
  n <- fc_res %>% filter(p.adj <= thresh) %>% 
    nrow()
  N <- nrow(fc_res)
  
  enrichment_scores <- list()
  
  for (i in subpathways){
    score <- (k[[i]] / m[[i]]) / ((n - k[[i]]) / (N - m[[i]]))
    enrichment_scores[i] <- score
  }
  
  res <- data.frame(SUB_PATHWAY = names(enrichment_scores),
                    Enrichment = unlist(enrichment_scores))
  return(res)
}
