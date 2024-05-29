if(!require("pacman")) {
  install.packages("pacman", repos = "http://cran.us.r-project.org")
}
pacman::p_load(
  BiocManager,
  factoextra,
  ggfortify,
  tidyverse,
  ggpubr,
  ggsci,
  rstatix,
  phyloseq,
  zoo,
  mixOmics,
  lemon,
  EnhancedVolcano,
  cowplot
)


####################### COLOR PALETTES #########################################
stacked_bar.palette <- c(
  "#E64B35FF",
             "#4DBBD5FF",
             "#00A087FF",
             "#3C5488FF",
             "#F39B7FFF",
             "#8491B4FF",
             "#91D1C2FF",
             "#DC0000FF",
             "#7E6148FF",
             "#B09C85FF",
             "#E4E9B2",
             "#F9A620",
             "#054A29",
             "#52414C"
)

######################### DATA WRANGLING #########################

get_feat_tables <- function(func_profile) {
  feat_tables <- within(func_profile, rm("Metadata"))
  return(feat_tables)
}

filter_rownames <- function(df, filt_vector) {
  # wrapper script around filter to handle rownames
  df_filt <- df %>% 
    rownames_to_column(var="id") %>%
    filter(id %in% filt_vector) %>%
    column_to_rownames(var="id")
  return(df_filt)
}

asin_sqrt <- function(x) {
  mean_sum <- mean(colSums(x))
  if (mean_sum == 1) {
    ps_t <- asin(sqrt(x))
  } else if (mean_sum == 100) {
    x <- x / 100
    ps_t <- asin(sqrt(x))
  } else {
    stop("not proportion data, exiting")
  } 
}

transform_func <- function(func_profile, transform, features_are_rows=TRUE, offset=1e-6) {
  
  if ("Metadata" %in% names(func_profile)){
    feat_tables <- get_feat_tables(func_profile)
  } else {
    feat_tables <- func_profile
  }
  
  if (features_are_rows == FALSE){
    feat_tables <- map(feat_tables, function(x) t_df(x)) 
  }
  if (transform %in% c("relative", "arcsin", "log10")) {
    if (transform == "relative") {
      feat_tables <- map(feat_tables, function(x) t_df(100 * t_df(x) / colSums(x))) 
    }
    else if (transform == "arcsin") {
        feat_tables <- map(feat_tables, function(x) asin_sqrt(x))
    } else if (transform == "log10"){
      feat_tables <- map(feat_tables, function(x) log10(x + offset))
    }
  } else {
      stop("Not a valid transform")
 } 
 
  
  
  if (features_are_rows == TRUE){
    return(feat_tables)
  } else {
    feat_tables <- map(feat_tables, function(x) t_df(x))
    feat_tables$Metadata <- func_profile$Metadata
    return(feat_tables)
  }
}

filt_prev_abund  <- function(x, abund, prev, orientation="cols") {
  if (orientation %in% c("cols", "columns", "c", "r", "rows")){
    if (orientation %in% c("cols", "columns", "c")){
      filt <- x[which(rowSums(x > abund) 
                      >= prev * ncol(x)), ]
      return(filt)
    } else {
      filt <- x[, which(colSums(x > abund) 
                      >= prev * nrow(x))]
      return(filt)
    }
    
  }

  return(filt)
}

filter_func <- function(func_profile, abund, prev, renorm=T) {
  
  feat_tables <- get_feat_tables(func_profile)
  func_profile_filt <- map(feat_tables, function(x) filt_prev_abund(x, abund, prev, orientation="rows"))
  
  func_profile_filt$Metadata <- func_profile$Metadata
  
  if (renorm == T){
   func_profile_filt <- transform(func_profile_filt, "relative", features_are_rows = F)
  }
  return(func_profile_filt)
}

annotate_func <- function(func_profile, feat_type=NULL) {
  
  if (feat_type %in% c("EC", "KO", "pathways")){
    if (feat_type == "KO"){
      mapfile <- read_tsv("D:/Users/adam-/Data/Metagenomics/data/ko_info.tsv", 
                      col_names = c("feat", "name"))
    } else if (feat_type == "EC"){
      mapfile <- read_tsv("D:/Users/adam-/Data/Metagenomics/data/ec_level4_info.tsv", 
                          col_names = c("feat", "name"))
    } else {
      mapfile <- read_tsv("D:/Users/adam-/Data/Metagenomics/data/metacyc_pathways_info.txt", 
                             col_names = c("feat", "name"))
    }
    
  }
  
  lookup <- mapfile$feat
  names(lookup) <- mapfile$name
  
  feat <- func_profile[[feat_type]]
  
  feat$name <- lapply(rownames(feat), function(x) names(lookup)[match(x, lookup)])

  feat <- feat %>% rownames_to_column("id") %>%
    unite(feat_type, all_of(c("id", "name")), sep = "_", remove = T) %>%
    column_to_rownames(var="feat_type") %>% 
    t_df()

  func_profile[feat_type] <- feat
}

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

t_df <- function(x) {
  return(as.data.frame(t(x)))
}

#########################    IMPORTING DATA #########################
read_tab_delim_metag <- function(df, comment_char="#") {
  # read all tab delimited files using these params
  df_out <-
    read.table(
      df,
      row.names = 1,
      header = 1,
      sep = "\t",
      check.names = F,
      quote = "",
      comment.char = comment_char
    )
  return(df_out)
}

import_func_profile <- function(meta, un, ec, ko, py, filter_ungrouped=T, renorm=T, 
                                read=F, fix_names=T){
  
  if (read == TRUE){
    
    un <- read_tab_delim_metag(un, comment_char = "")
    ec <- read_tab_delim_metag(ec, comment_char = "")
    ko <- read_tab_delim_metag(ko, comment_char = "")
    py <- read_tab_delim_metag(py, comment_char = "")
    meta <- read_tab_delim_metag(meta, comment_char = )
  }
  
  feat <- list(un, ec, ko, py)
  names(feat) <- c("UNIREF", "EC", "KO", "pathways")
  
  if (fix_names == T){
   colnames(feat$UNIREF) <- gsub("_merged_clean_Abundance-RPKs", "", colnames(feat$UNIREF))
   colnames(feat$EC) <- gsub("_merged_clean_Abundance-RPKs", "", colnames(feat$EC))
   colnames(feat$KO) <- gsub("_merged_clean_Abundance-RPKs", "", colnames(feat$KO))
   colnames(feat$pathways) <- gsub("_merged_clean_Abundance", "", colnames(feat$pathways))
    
  }
  
  sample_names <- meta %>% rownames()
  
  feat <- map(feat, function(x) x %>% 
                select(all_of(sample_names)) )
  
  if (filter_ungrouped == TRUE){
    feat <- map(feat, function(x) x %>% 
                  rownames_to_column("tmp") %>% 
                  filter(! tmp %in% c("UNGROUPED", "UNMAPPED", "UNINTEGRATED")) %>% 
                  column_to_rownames("tmp"))
    if (renorm == T){
      feat <- map(feat, function(x) t_df(100 * t_df(x) / colSums(x)))
    }
  }
  # annotate
  ## skip UNIREF as these are already annotated 
  feat <- map(c("EC","KO", "pathways"), function(x) annotate_func(func_profile = feat, 
                                                        feat_type = x))
  names(feat) <- c("EC", "KO", "pathways")
  
  # add EC prefix 
  colnames(feat[["EC"]]) <- paste0("EC:", colnames(feat[["EC"]]))
  
  return_list <- list("UNIREF" = feat[["UNIREF"]],
                      "EC" = feat[["EC"]],
                      "KO" = feat[["KO"]],
                      "pathways" = feat[["pathways"]],
                      "Metadata" = meta)
  return(return_list)

}


######################### NORMALISATION & FILTERING ############################

subset_samples_func <- function(func_profile, filter_statement) {
  
  keep <- func_profile$Metadata %>% 
    filter(eval(rlang::parse_expr(filter_statement))) %>% 
    rownames()
  
  func_profile_filt <- map(func_profile, function(x) filter_rownames(x, keep))
  
  return(func_profile_filt)
  
}


prune_samples_func  <- function(func_profile, keep) {
  
  func_profile_filt <- map(func_profile, function(x) filter_rownames(x, keep))
  return(func_profile_filt)
  
}

################# DIFFERENTIAL ABUNDANCE #######################################

func_fold_change <- function(func_profile, feat_type, log.n0=1e-08, mult.corr="fdr",
                             var, case, ctrl, fc_method="fc") {
  
  if (feat_type %in% c("UNIREF","uniref", "EC", "ec", "KO", "ko", "PY", "PWY", "pathways")) {
  
  
  feat <- t_df(func_profile[[feat_type]])
  meta <- func_profile[["Metadata"]] %>% 
    rownames_to_column("id")
  
 
  # infer name of sample column 
  sampleid_col <- colnames(meta)[1]
  
  stopifnot(all(meta[sampleid_col] %>% pull() 
                %in% colnames(feat)))
  
  # create pval matrix with row for each feature
  p.val <- matrix(NA, nrow=nrow(feat), ncol=1, 
                  dimnames=list(row.names(feat), c("pval")))
  # duplicate pval matrix to store gfc 
  fc <- p.val
  colnames(fc) <- c("FC") 
  
  
  cat("Calculating pvalue and FC/GFC for each feature...\n")
  pb <- txtProgressBar(max=nrow(feat), style=3)
  
  # calculate wilcoxon test and effect size for each feature and study
  for (f in row.names(feat)) {
    
    x <- unlist(feat[f, meta %>% 
                         filter(.data[[var]]==case) 
                       %>% pull(sampleid_col)])
    
    y <- unlist(feat[f, meta %>% 
                         filter(.data[[var]]==ctrl) 
                       %>% pull(sampleid_col)])
    
    # Wilcoxon
    p.val[f, ] <- wilcox.test(x, y, exact=FALSE)$p.value
    
    # FC
    if (fc_method %in% c("gfc", "GFC")){
      q.p <- quantile(log2(x+log.n0), probs=seq(.1, .9, .05))
      q.n <- quantile(log2(y+log.n0), probs=seq(.1, .9, .05))
      fc[f,] <- sum(q.p - q.n)/length(q.p)
    } else if (fc_method %in% c("FC", "fc")){
      q.p <- log2(x+log.n0)
      q.n <- log2(y+log.n0)
      fc[f,] <- median(q.p - q.n)
    }
    
    # progressbar
    setTxtProgressBar(pb, (pb$getVal()+1))
  }
  cat('\n')
  
  # multiple hypothesis correction
  p.adj <- data.frame(apply(p.val, MARGIN=2, FUN=p.adjust, method=mult.corr),
                      check.names = FALSE)
  ### Return datafram containing pval and log2FC
  
  df_res <- data.frame(p.adj,fc)
  colnames(df_res) <- c("p.adj", "log2FC")
  
  
  comparison <- sort(c(case, ctrl))
  
  return_list <- list(res = df_res, groups = comparison)
  }
}

# rewrite for func profile
pairwise_wilcox_feattab <- function(ps, factor, case, control, 
                                    sample_col="sample_id", p_adjust="BH"){
  
  feattab <- ps_to_feattab(ps)
  meta <- meta_to_df(ps, rownames = F)
  # initialise pval matrix 
  p.val <- matrix(NA, nrow=nrow(feattab), ncol=1, 
                  dimnames=list(row.names(feattab), "p.val"))
  
  featmat <- as.matrix(feattab)
  
  for (row in rownames(featmat)){
    x <- featmat[row, meta %>% 
                   filter(.data[[ factor ]] ==control) %>% pull(.data[[ sample_col ]])]
    y <- featmat[row, meta %>% 
                   filter(.data[[ factor ]] ==case) %>% pull(.data[[ sample_col ]])]
    
    
    p.val[row, ] <- wilcox.test(x, y, exact=FALSE)$p.value
  }
  
  p.val.adj <- p.val %>% 
    as.data.frame() %>% 
    adjust_pvalue("p.val", method = p_adjust)
  
  return(p.val.adj)
}



######################### PLOTTING AND STATS #########################

# calculate colour palettes
calc_pal <- function(ps, group_variable) {
  # update function to accept nonps objects and add support for continuous palettes
  meta <- meta_to_df(ps)
  groups <- unique(meta[, group_variable])
  
  if (length(groups) < 10) {
    pal <- pal_npg()(length(groups))
  }
  else {
    print("Exceeded colour palette limit, use custom palette")
  }
}


# create statistical formula
xyform <- function (y_var, x_vars) {
  # y_var: a length-one character vector
  # x_vars: a character vector of object names
  as.formula(sprintf("%s ~ %s", y_var, paste(x_vars, collapse = " + ")))
}

round_any <- function(x, accuracy, f=ceiling) {
  f(x/ accuracy) * accuracy
}

# plot boxplot with stats
plot_boxplot <- function(df,
                         variable_col,
                         value_col,
                         comparisons_list,
                         fill_var = variable_col,
                         xlab = variable_col,
                         ylab = value_col,
                         p_title = NULL,
                         multiple_groups = FALSE,
                         cols = NULL,
                         group.order = NULL,
                         paired = FALSE,
                         ...) {
  # extend color palette with transparent value - required due to way we are
  # layering plot
  if (is.null(cols)) {
    cols <- pal_npg()(length(unique(df[, variable_col])))
  }
  cols <- c(cols, "transparent")
  
  if (!is.null(group.order)) {
    df[, variable_col] <-
      factor(df[, variable_col], levels = group.order)
  }
  
  formula <- xyform(value_col, variable_col)
  
  if (multiple_groups == TRUE) {
    if (paired == TRUE) {
      stat_variance <- df %>%
        friedman_test(formula)
      stat_test <- df %>%
        pairwise_wilcox_test(
          formula,
          comparisons = comparisons_list,
          p.adjust.method = "BH",
          paired = TRUE
        ) %>%
        add_significance() %>%
        add_xy_position(x = variable_col) %>%
        filter(p.adj < 0.05)
    }
    else {
      stat_variance <- df %>%
        kruskal_test(formula)
      stat_test <- df %>%
        pairwise_wilcox_test(formula,
                             comparisons = comparisons_list,
                             p.adjust.method = "BH") %>%
        add_significance() %>%
        add_xy_position(x = variable_col) %>%
        filter(p.adj < 0.05)
    }
  }
  else if (multiple_groups == FALSE) {
    if (paired == TRUE) {
      stat_test <- df %>%
        wilcox_test(formula, paired = TRUE) %>%
        add_significance() %>%
        add_xy_position(x = variable_col) %>%
        filter(p < 0.05)
    }
    else {
      stat_test <- df %>%
        wilcox_test(formula) %>%
        add_significance() %>%
        add_xy_position(x = variable_col) %>%
        filter(p < 0.05)
    }
    
  }
  
  # aes string accepts strings as column names, this code plots boxplot and adds error bars
  plot <- ggplot(
    df,
    aes_string(
      x = variable_col,
      y = value_col,
      fill = variable_col,
      color = variable_col
    )
  ) +
    geom_boxplot(
      color = "black",
      alpha = 0.8,
      outlier.shape = 5,
      outlier.size = 1
    ) +
    geom_point(size = 1.5, position = position_jitterdodge()) +
    labs(x = xlab, y = ylab) +
    stat_boxplot(color = "black",
                 geom = "errorbar",
                 width = 0.2)
  # creates new 'finalised plot' and adds statistical significance, labels and adjusts theme and title
  final_plot <- plot +
    theme_classic2() +
    ggtitle(p_title) +
    theme(
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.title.x = element_text(size = 18),
      axis.title.y = element_text(size = 18),
      axis.line = element_line(colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.border = element_blank(),
      plot.title = element_text(hjust = 0.5),
      legend.position = "None"
    ) +
    scale_fill_manual(values = cols) +
    scale_color_manual(values = cols) +
    rotate_x_text(angle = 45)
  
  if (dim(stat_test)[1] == 0) {
    plot_out <- final_plot
  }
  else {
    
    # get ymax for calculating limits and breaks
    datamax <- max(df[, value_col])
    statmax <- max(stat_test["y.position"])
    if (datamax > 100){
      if (datamax < 250) {
        ybreak <- 50
      }
      else {
        ybreak <- 100
      }
      ylimit <- round_any(statmax, accuracy = 100)
      ybreakmax <- round_any(datamax, accuracy = 100)
    } else if (datamax < 100) {
      ybreak <- 10
      ylimit <- round_any(statmax, accuracy = 10)
      ybreakmax <- round_any(datamax, accuracy = 10)
    } else if (datamax < 10) {
      ybreak <- 1
      ylimit <- round_any(statmax, accuracy = 1)
      ybreakmax <- round_any(datamax, accuracy = 1)
    }
    
    if (multiple_groups == T) {
      plot_out <- final_plot +
        stat_pvalue_manual(
          stat_test,
          label = "p.adj.signif",
          hide.ns = T,
          inherit.aes = FALSE,
          ...
        ) +
        scale_y_continuous(breaks = seq(0, ybreakmax,by=ybreak), limits = c(0, ylimit)) +
        coord_capped_cart(left='top', expand = F)
    }
    else {
      plot_out <- final_plot +
        stat_pvalue_manual(
          stat_test,
          label = "p.signif",
          hide.ns = T,
          inherit.aes = FALSE,
          ...
        ) +
        scale_y_continuous(breaks = seq(0, ybreakmax, by=ybreak), limits = c(0, ylimit)) +
        coord_capped_cart(left='top', expand = F)
      
    }
  }
  
  return(plot_out)
}

# plot scatter plot with correlation if desired
plot_scatter <- function(df,
                         x,
                         y,
                         point_color,
                         line_color,
                         fill_color,
                         xlab,
                         ylab,
                         corr.method = NULL,
                         ...) {
  p <-
    ggplot(data = df, mapping = aes(x = .data[[x]], y = .data[[y]])) +
    geom_point(aes(color = point_color), size = 2.5) +
    geom_smooth(method = "lm",
                color = line_color,
                fill = fill_color) +
    theme_bw() +
    theme(
      legend.position = "None",
      axis.title.x = element_text(size = 14),
      axis.text.x = element_text(size = 12),
      axis.title.y = element_text(size = 14),
      axis.text.y = element_text(size = 12)
    ) +
    xlab(xlab) +
    ylab(ylab)
  
  if (!is.null(corr.method)) {
    p <- p + stat_cor(method = corr.method, ...)
    return(p)
  }
  else {
    return(p)
  }
}

plot_taxonomic_comp  <- function(ps, tax_level, var, ord=NULL, n_taxa=10) {
  
  top_taxa <- c(get_top_n(ps, n=n_taxa, level = tax_level), "other")
  taxa_pal <- c(stacked_bar.palette[1:length(top_taxa) -1], "#DCDCDC")
  names(taxa_pal) <- top_taxa
  
  if (! is.null(ord)){
    ps <- ps %>% 
      ps_mutate(plot_var = factor(.data[[ var ]], levels = ord))
  }
  
  comp_fig <- ps %>% 
    tax_fix(unknowns = c("unassigned")) %>%
    comp_barplot(
      tax_level = tax_level, n_taxa = n_taxa,
      bar_outline_colour = "black", bar_width = 0.9,
      palette = taxa_pal) +
    facet_grid(.~factor(plot_var),
               scales = "free", space = "free"
    ) + 
    theme_cowplot() +
    scale_y_continuous(labels = scales::percent, expand = c(0,0)) +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=22),
          axis.text.y = element_text(size=18),
          axis.ticks = element_blank())
  
  return(comp_fig)
}

plot_pca <- function(func_profile,what, var, colours, add_ellipse=F) {
  
  # need to check what variable heres
  
  # calc PCA
  metadata_filt <- func_profile[["Metadata"]] %>% select(all_of(var))
  dat <- merge(func_profile[[what]], metadata_filt, by=0) %>% 
    column_to_rownames("Row.names") 
  # this will select numerical metadata
  pca <- prcomp(dat %>% select(where(is.numeric)), scale. = T)
 
  
  # Plot
  scree_plot <- fviz_eig(pca)
  pca_plot <- autoplot(pca, data = dat,
                       colour = var) + 
    geom_point(aes(color= .data[[ var ]] ), size=3, alpha=0.75) + 
    scale_color_manual(values = colours) + 
    scale_fill_manual(values = colours) +
    theme_cowplot()
  pca_plot$layers[[1]] <- NULL
  
  if (add_ellipse == TRUE){
    pca_plot <- pca_plot +
      geom_polygon(stat = "ellipse", aes(fill = .data [[ var ]] ), alpha = 0.3)
  }
  
  return(pca_plot)
}

plot_volcano <- function(results_df, pthresh=0.05, FCthresh=0.58, ...) {
  res <- results_df$res
  # replace rownames with highest classified level
  rownames(res) <- as.character(lapply(strsplit(as.character(rownames(res)), split="\\|"),
                                       tail, n=1))
  groups <- results_df$groups
  EnhancedVolcano(res, x="log2FC", y="p.adj", 
                  lab=rownames(res),
                  labSize = 4,
                  pCutoff = pthresh, 
                  FCcutoff = FCthresh, 
                  title = paste(groups[1], 
                                "versus", 
                                groups[2], 
                                sep=" "), 
                  subtitle = NULL,
                  max.overlaps = 30,
                  ...)
}

maaslin2_func <- function(func_profiles, feattype, out, 
                          fixed, abun_thresh=0, prev_thresh=0, ...){
  # write check for filtering and AST transform 
  if (feattype %in% c("UNIREF", "EC", "KO", "pathways")){
    mat_in <- func_profiles[[feattype]]
    metadata_in <- func_profiles$Metadata
  } else {
    stop("Feature type not supported, choose from: 'UNIREF', 'EC', 'KO', 'pathways'")
  }
  
  Maaslin2(input_data = mat_in,
           input_metadata = metadata_in, 
           output = out, 
           normalization = "NONE",
           transform = "AST",
           min_abundance = abun_thresh,
           min_prevalence = prev_thresh,
           fixed_effects = fixed,
           cores=6,
           ...
           )
}


