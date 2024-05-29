if(!require("pacman")) {
  install.packages("pacman", repos = "http://cran.us.r-project.org")
}
pacman::p_load(
  BiocManager,
  tidyverse,
  ggpubr,
  ggsci,
  rstatix,
  phyloseq,
  vegan,
  zoo,
  mixOmics,
  lemon,
  EnhancedVolcano,
  microbiome,
  microViz,
  cowplot,
  Maaslin2
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
  "#52414C",
  "#D81E5B",
  "#331832",
  "#27474E",
  "#573D1C",
  "#404E4D",
  "#DAD4EF",
  "#E86A92",
  "#044389",
  "#6C4B5E",
  "#4E6E58",
  "#826AED",
  "#FF0054",
  "#9E0059",
  "#387D7A",
  "#395E66",
  "#1BE7FF"
)

######################### DATA WRANGLING #########################

taxonomy <- function (ps) {
  return(as.data.frame(tax_table(ps)))
}

meta_to_df <- function(ps, rownames=T) {
  if (rownames == T){
    return(as(sample_data(ps), "data.frame")) 
  } else {
    return(as(sample_data(ps), "data.frame") %>% 
             rownames_to_column("sample_id"))
  }
  
}

ps_to_feattab <- function(ps) {
  return(as.data.frame(ps@otu_table))
}

get_top_n <- function(ps, n, level = "species") {
  if (level != "species") {
    ps <- ps %>% tax_fix %>%
      tax_agg(rank = level)
    ps <- ps_get(ps)
  }
  topn <- ps %>%
    transform(transform = "relative") %>%
    psmelt() %>%
    group_by(OTU) %>%
    summarise(Mean_abund = mean(Abundance)) %>%
    slice_max(Mean_abund, n = n) %>%
    pull(OTU)
  return(topn)
}

get_top_n_group <- function(ps, n, level = "species",var,  
                            group=NULL) {
  if (level != "species") {
    ps <- ps %>% tax_fix %>%
      tax_agg(rank = level)
    ps <- ps_get(ps)
  }
  topn <- ps %>%
    transform(transform = "relative") %>%
    psmelt() %>%
    filter(.data[[ var ]] == group) %>% 
    group_by(OTU) %>%
    summarise(Mean_abund = mean(Abundance)) %>%
    slice_max(Mean_abund, n = n) %>%
    pull(OTU)
  return(topn)
}

filter_rownames <- function(df, filt_vector) {
  # wrapper script around filter to handle rownames
  df_filt <- df %>% 
    rownames_to_column(var="id") %>%
    filter(id %in% filt_vector) %>%
    column_to_rownames(var="id")
  return(df_filt)
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

load_phylo <- function(asvtab, taxa, mapping, tree = NULL) {
  # convert to phyloseq and return list
  phylo_asv <- otu_table(asvtab, taxa_are_rows = T)

  phylo_tax <- tax_table(as.matrix(taxa))

  phylo_map <- sample_data(mapping)

  if (exists("tree")) {
    phylo_tree <- read_tree(tree)
    return(merge_phyloseq(phylo_asv, phylo_tax, phylo_tree, phylo_map))
  }
  else {
    return(merge_phyloseq(phylo_asv, phylo_tax, phylo_map))
  }
}

return_level <- function(level) {
  if (level %in% c("Subspecies", "subspecies", "t")){
    return("t")
  } else if (level %in% c("Species", "species", "s")){
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



######################### NORMALISATION & FILTERING ############################

transform <- function(ps, transform = "relative", offset=1) {
  if (length(is(ps)) == 1 && class(ps) == "phyloseq") {
    x <- ps_to_feattab(ps)
    mean_sum <- mean(colSums(x))
  }
  else {
    print("not a phyloseq object, exiting")
    stop()
  }

  if (transform %in% c("relative", "arcsin", "log10")) {
    if (transform == "relative") {
      ps_t <- t(t(x) / colSums(x))

    } else if (transform == "arcsin") {
        if (mean_sum == 1) {
          ps_t <- asin(sqrt(x))
        } else if (mean_sum == 100) {
          x <- x / 100
          ps_t <- asin(sqrt(x))
        } else {
          stop("not proportion data, exiting")
        } 
    } else if (transform == "log10"){
          if (max_x == 100 | max_x == 1) {
            warning("data are relative abundances", call. = F)
          ps_t <- log10(x + offset)
        }   else {
        ps_t <- log10(x + offset)
      }
    }
    otu_table(ps)@.Data <- as.matrix(ps_t)

    return(ps)

  } else {
      print("Not a valid transform, exiting")
    stop()
  }

}

filter_ps <- function(ps, abun, prev=NULL) {
  
  x <- ps_to_feattab(ps) %>% 
    as.matrix()
  
  
  if (is.null(prev)) {
    ps_filt <- x[rowSums(x)  >= abun, ]
  } else if (!is.null(prev)) {
   ps_filt <- x[rowSums(x >= abun) >= ncol(x) * prev, ]
  }
  
  otu_table(ps)@.Data <- ps_filt
  
  return(ps)
}

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
######################### ALPHA DIVERSITY #########################

# calculate shannon effective
Shannon.E <- function(x) {
  summed <- sum(x)
  shannon.e <-
    round(exp(-sum(x[x > 0] / summed * log(x[x > 0] / summed))), digits = 2)
  return(shannon.e)
}

Shannon <- function(x) {
  summed <- sum(x)
  shannon.e <-
    round(exp(-sum(x[x > 0] / summed * log(x[x > 0] / summed))), digits = 2)
  return(shannon.e)
}

InvSimpson <- function(x) {
  summed <- sum(x)
  shannon.e <-
    round(exp(-sum(x[x > 0] / summed * log(x[x > 0] / summed))), digits = 2)
  return(shannon.e)
}

# calculate richness
Richness <- function(x, detection = 1e-5) {
  observed <- sum(x > detection)
  return(observed)
}



# calculate all alpha diversity matrices and return dataframe
calc_alpha <- function(ps, ...) {
  mat_in <- ps_to_feattab(ps) %>%
    t()

  diversity <-
    setNames(data.frame(matrix(ncol = 2, nrow = nsamples(ps))),
             c("Richness", "Shannon.Effective"))
  rownames(diversity) <- rownames(meta_to_df(ps))

  diversity$Richness <- apply(mat_in, 1, Richness, ...)
  diversity$Shannon.Effective <- apply(mat_in, 1, Shannon.E)
  diversity$Shannon <- apply(mat_in, 1, function(x) diversity(x))
  diversity$InvSimpson <- apply(mat_in, 1, function(x) diversity(x, index = "invsimpson"))

  return(diversity)
}



######################### BETA DIVERSITY #########################

calc_betadiv <- function(ps, dist, ord_method = "NMDS") {
  if (ord_method %in% c("NMDS", "MDS", "PCoA")) {

    if (dist %in% c("bray", "jaccard")) {
      dist_mat <- distance(ps, dist)
      if (ord_method == "NMDS"){
        ord <- ordinate(ps, ord_method, dist_mat, trace=FALSE)
      } else {
        ord <- ordinate(ps, ord_method, dist_mat)
      }
      
      return_list <- list("Distance_Matrix" = dist_mat,
                          "Ordination" = ord)
      return(return_list)
    }
    else {
      print(
        "Distance metric not supported, supported metrics are; bray, jaccard"
      )
    }
  }
  else {
    print("Ordination method not supported, supported methods are: NMDS, MDS, PCoA")
    stop()
  }

}

################# DIFFERENTIAL ABUNDANCE #######################################

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


tax_fold_change <- function(ps, log.n0=1e-05, mult.corr="fdr",
                var, case, ctrl, feat_type="taxa", fc_method="fc") {
  
  seqtab <- ps_to_feattab(ps)
  meta <- meta_to_df(ps, rownames = F)
  
  if (feat_type == "functional"){
    log.n0 = 1e-08
  }
  # infer name of sample column 
  sampleid_col <- colnames(meta)[1]
  
  stopifnot(all(meta[sampleid_col] %>% pull() 
                %in% colnames(seqtab)))
  
  # create pval matrix with row for each feature
  p.val <- matrix(NA, nrow=nrow(seqtab), ncol=1, 
                  dimnames=list(row.names(seqtab), c("pval")))
  # duplicate pval matrix to store gfc 
  fc <- p.val
  colnames(fc) <- c("FC") 
  
  # calculate wilcoxon test and effect size for each feature and study
  for (f in row.names(seqtab)) {
    
    x <- unlist(seqtab[f, meta %>% 
                         filter(.data[[var]]==case) 
                       %>% pull(sampleid_col)])
    
    y <- unlist(seqtab[f, meta %>% 
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
   
   }
  # multiple hypothesis correction
  p.adj <- data.frame(apply(p.val, MARGIN=2, FUN=p.adjust, method=mult.corr),
                      check.names = FALSE)
  ### Return datafram containing pval and log2FC
  
  df_res <- data.frame(p.adj,fc)
  colnames(df_res) <- c("p.adj", "log2FC")
  
  
  comparison <- sort(c(case, ctrl))
  
  return_list <- list(res = df_res, groups = comparison)
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

# calculate adonis R2 and p-value
# better NA handling here
phyloseq_adonis <- function(ps, dist_matrix, group_variable, ...) {
  meta_df <- meta_to_df(ps)
  if (sum(is.na(meta_df[[group_variable]])) > 0) {
    print("metadata contains NAs, remove these samples with subset_samples
          before continuing")
    return (NULL)
  } else {
    # convert distance matix object to string
    dist_str <- deparse(substitute(dist_matrix))
    # define formula
    form <- as.formula(paste(dist_str, group_variable, sep="~"))
    ps_ad <- adonis2(form, data = meta_df, ...)
    return(ps_ad)
  }
}

# calculate betadisper
phyloseq_betadisper <- function(ps, dist_matrix, group_variable, ...) {
  meta_df <- meta_to_df(ps)
  if (sum(is.na(meta_df[[group_variable]])) > 0) {
    print("metadata contains NAs, remove these samples with subset_samples
          before continuing")
    return (NULL)
  } else {

    bd <- betadisper(dist_matrix, meta_df[[group_variable]])
    anova_res <- anova(bd)
    return(anova_res)
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
    } else if (datamax < 1){
      ybreak <- 0.1
      round_any(statmax, accuracy = 0.1)
      ybreakmax <- round_any(datamax, accuracy = 0.1)
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

plot_taxonomic_comp  <- function(ps, tax_level, var, ord=NULL, n_taxa=10, 
                                 per_group=F, groups=NULL) {
  ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  if (!tax_level %in% ranks){
    stop("Provide a valid taxonomic rank")
  }
  
  if (per_group == T){
    top_taxa <- map(groups, function(x) ps %>% 
                     get_top_n_group(n = n_taxa, level = tax_level,var = var,
                                     group =x))
    
    top_taxa <- c(unique(unlist(top_taxa)), "other")
    # redefine n_taxa as in this case it reflects n across groups
    n_taxa <- length(top_taxa)
  } else {
    top_taxa <- c(get_top_n(ps, n=n_taxa, level = tax_level), "other") 
  }
  
  taxa_pal <- c(stacked_bar.palette[1:length(top_taxa)-1], "#DCDCDC", "white")
  names(taxa_pal) <- top_taxa
  
  
  if (!is.null(ord)) {
    ps <- ps %>%
      ps_mutate(plot_var = factor(.data[[var]], levels = ord))
  }
   
  melt_ps <- ps %>% 
    aggregate_taxa(level = tax_level) %>% 
    psmelt() 
  
  # desired taxa - this does not work at present because taxa names only contain species designation and not entire string 
  taxa_plot <- melt_ps %>% filter(.data[[ tax_level ]] %in% top_taxa)
  # group rest into other
  repl_cols <- c("OTU", ranks)
  
  other_plot <- melt_ps %>% filter(!.data[[ tax_level ]] %in% top_taxa) %>% 
    group_by(Sample) %>% 
    mutate(Abundance = sum(Abundance)) %>% 
    ungroup %>% 
    distinct(Sample, .keep_all = T) %>% 
   mutate_at(vars(repl_cols), ~ paste("other"))
  # bind dfs
  plot_df <- bind_rows(taxa_plot, other_plot)
  
  comp_fig <- plot_df %>%
    ggplot(aes(fill=.data[[ tax_level ]], y=Abundance, x=Sample)) + 
    geom_bar(position="fill", stat="identity") + 
    scale_fill_manual(values = taxa_pal) + 
    facet_grid(. ~ factor(plot_var),
               scales = "free", space = "free") +
    theme_cowplot() +
    scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
    theme(
      panel.background = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 22),
      axis.text.y = element_text(size = 18),
      axis.ticks = element_blank()
    )
  return(comp_fig)
}

plot_beta_div <- function(ps,
                          dist_matrix,
                          ordination,
                          group_variable,
                          add_ellipse = FALSE,
                          cols = NULL,
                          shape = NULL) {
  if (is.null(cols)) {
    cols <- calc_pal(ps, group_variable)
  }
  # significance
  ad <- phyloseq_adonis(ps, dist_matrix, group_variable)
  betadisp <- phyloseq_betadisper(ps, dist_matrix, group_variable)
  
  # check sample data dimensions >1 otherwise phyloseq
  # fails to colour samples due to bug:https://github.com/joey711/phyloseq/issues/541
  if (dim(sample_data(ps))[2] < 2) {
    # add repeat of first column as dummy
    sample_data(ps)[, 2] <- sample_data(ps)[, 1]
  }
  plot <-
    plot_ordination(ps, ordination , color = group_variable, shape = shape)
  plot$layers[[1]] <- NULL
  
  
  
  if (is.numeric(meta_to_df(ps)[[group_variable]])) {
    plot_out <- plot + geom_point(size = 3, alpha = 0.75) +
      theme_bw() +
      #scale_color_viridis_c() + 
      #scale_fill_viridis_c() +
    labs(caption = bquote(Adonis ~ R ^ 2 ~ .(round(ad$R2[1], 2)) ~
                            ~ p - value ~ .(ad$`Pr(>F)`[1])))
  } else {
    plot_out <- plot + geom_point(size = 3, alpha = 0.75) +
      theme_bw() +
      scale_fill_manual(values = cols) +
      scale_color_manual(values = cols) +
      labs(caption = bquote(Adonis ~ R ^ 2 ~ .(round(ad$R2[1], 2)) ~
                              ~ p - value ~ .(ad$`Pr(>F)`[1])))
  }
  
  if (add_ellipse == TRUE) {
    plot_out <- plot_out +
      geom_polygon(stat = "ellipse", aes(fill = .data [[group_variable]]), alpha = 0.3)
  }
  
  # ideally this needs to be added to the main plot eventually underneath adonis
  # temporary fix for when betadispersion calculation fails
  # if (!is.null(betadisp)){
  #   if (betadisp$`Pr(>F)`[[1]] < 0.05) {
  #     warning("Group dispersion is not homogenous, interpret results carefully",
  #             call. = F)
  #   }
  #   
  # }

  return(plot_out)
}

plot_volcano <- function(results_df, pthresh=0.05, FCthresh=0.58, ...) {
  res <- results_df$res
  # replace rownames with highest classified level
  rownames(res) <- as.character(lapply(strsplit(as.character(rownames(res)), split="\\|"),
                                       tail, n=1))
  groups <- results_df$groups
  EnhancedVolcano(res, x="log2FC", y="p.adj", 
                  lab=rownames(res),
                  labSize = 3,
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

maaslin2_tax <- function(tax_profile, out, 
                          fixed, abun_thresh=0, prev_thresh=0, ...){
  mat_in <- ps_to_feattab(tax_profile) %>% 
    t_df()
  # this only works for species level
  colnames(mat_in) <- gsub(".*s__","",colnames(mat_in))
  
  metadata_in <- meta_to_df(tax_profile)
  
  Maaslin2(input_data = mat_in,
           input_metadata = metadata_in, 
           output = out, 
           normalization = "NONE",
           transform = "NONE",
           min_abundance = abun_thresh,
           min_prevalence = prev_thresh,
           fixed_effects = fixed,
           cores=6,
           ...
  )
}

