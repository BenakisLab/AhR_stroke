---
title: "AHR_scRNAseq_analysis"
authors: "Alba Simats & Corinne Benakis"
date: "25/04/2026"
---

library(devtools)
library(Seurat)
library(dplyr)
library(Matrix)
library(openxlsx)
library(RColorBrewer)
library(patchwork)
library(readxl)
library(RColorBrewer)
library(ggplot2)

setwd("...")
getwd()

matrix_dir = "/Volumes/.../filtered_feature_bc_matrix/"
list.files(matrix_dir)

barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
mat <- Matrix::readMM(file = matrix.path)

barcode.names = read.delim(barcode.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
feature.names = read.delim(features.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)

GEX <- mat[0:32285,] #check the feature.name file and adjust the rows here
GEX_feature.names <- feature.names[0:32285,]

HTO <- mat[32286:32289,]
HTO_feature.names <- feature.names[32286:32289,]
HTO_feature.names

colnames(GEX) = barcode.names$V1
colnames(HTO) = barcode.names$V1
rownames(GEX) = GEX_feature.names$V2
rownames(HTO) = HTO_feature.names$V1

rownames(GEX) <- make.unique(rownames(GEX))
rownames(HTO) <- paste0(rownames(HTO), "-HTO")


# Create seurat object 
s1 <- CreateSeuratObject (counts = GEX, project = "Ahr")
s1[["HTO"]] <- CreateAssayObject(counts = HTO)
s1 <- NormalizeData(s1, assay = "HTO", normalization.method = "CLR")
s1 <- HTODemux(s1, assay = "HTO", positive.quantile = 0.99)
table(s1$HTO_classification.global)

Idents(s1) <- "HTO_maxID"
RidgePlot(s1, assay = "HTO", features = rownames(s1[["HTO"]])[1:4], ncol = 2)

Idents(s1) <- "HTO_classification.global"
VlnPlot(s1, features = "nCount_RNA", pt.size = 0.1, log = TRUE)
head(s1)

# Extract singlets
Ahr2 <- subset(s1, idents = c("Doublet", "Negative"), invert = TRUE)
table(Ahr2$orig.ident)

# Get batches based on cell names
sample <- sapply(colnames(GetAssayData(object = Ahr2, slot = "counts")),
                 FUN=function(x){substr(x,18,18)})
sample <- as.numeric(as.character(sample))
names(sample) <- colnames(GetAssayData(object = Ahr2, slot = "counts"))
Ahr2 <- AddMetaData(Ahr2, sample, "sample")

new.grouping <- c("organ")
Ahr2[[new.grouping]] <- new.grouping
colnames(Ahr2@meta.data)
Ahr2$organ[Ahr2$HTO_classification == "B0301-HTO" ] <- "Sp" #spleen
Ahr2$organ[Ahr2$HTO_classification == "B0302-HTO" ] <- "Br" #brain
Ahr2$organ[Ahr2$HTO_classification == "B0303-HTO" ] <- "Lp" #ileal lamina propria
Ahr2$organ[Ahr2$HTO_classification == "B0304-HTO" ] <- "Bl" #blood
table(Ahr2$organ)
saveRDS(Ahr2, "Ahr.rds")

Ahr2[["percent.mt"]] <- PercentageFeatureSet(Ahr2, pattern = "^mt-")
VlnPlot(object = Ahr2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(Ahr2, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(Ahr2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
Ahr2 <- subset(Ahr2, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10)
VlnPlot(object = Ahr2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
saveRDS(Ahr2, "Ahr.rds")

library(sctransform)
Ahr2 <- SCTransform(Ahr2, vars.to.regress = "percent.mt", verbose = FALSE)
Ahr2 <- RunPCA(object = Ahr2, features = VariableFeatures(object = Ahr2), verbose = FALSE)
Ahr2 <- ProjectDim(object = Ahr2)
ElbowPlot(object = Ahr2)
DimHeatmap(object = Ahr2, dims = 8:16, cells = 500, balanced = TRUE)

Ahr2 <- FindNeighbors(object = Ahr2, dims = 1:10) 
Ahr2 <- FindClusters(object = Ahr2, resolution = 0.1) 
Ahr2 <- RunUMAP(object = Ahr2, dims = 1:10)
DimPlot(Ahr2, reduction = 'umap', label = TRUE, group.by = "seurat_clusters")
saveRDS(Ahr2, file = "Ahr.rds")

Ahr <-Ahr2
DimPlot(Ahr)
new.grouping <- c("experiment")
Ahr[[new.grouping]] <- new.grouping
colnames(Ahr@meta.data)
Ahr$experiment[Ahr$sample == "1" ] <- "WT"
Ahr$experiment[Ahr$sample == "2" ] <- "WT"
Ahr$experiment[Ahr$sample == "3" ] <- "WT"
Ahr$experiment[Ahr$sample == "4" ] <- "KO"
Ahr$experiment[Ahr$sample == "5" ] <- "KO"
Ahr$experiment[Ahr$sample == "6" ] <- "KO"
DimPlot(object = Ahr, group.by="experiment")


Ahr
Ahr@active.assay = "SCT"
Ahr <- FindVariableFeatures(object = Ahr, selection.method = 'mean.var.plot', mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf))
Ahr[["SCT"]] <- split(Ahr[["SCT"]], f = Ahr$experiment)
Ahr <- IntegrateLayers(object = Ahr, method = HarmonyIntegration, assay = "SCT", orig.reduction = "pca", 
                       new.reduction = 'harmony', verbose = FALSE)
Ahr[["RNA"]] <- JoinLayers(Ahr[["RNA"]])
Ahr <- FindNeighbors(object = Ahr, reduction = "harmony", dims = 1:14) #15
Ahr <- FindClusters(object = Ahr, resolution = 0.1) #0.2
Ahr <- RunUMAP(object = Ahr, dims = 1:14, reduction = "harmony")
saveRDS(Ahr, "Ahr_integrated.rds")

Ahr@active.assay = "RNA"
Ahr <- NormalizeData(object = Ahr, normalization.method = "LogNormalize", scale.factor = 1e4)
Ahr <- FindVariableFeatures(object = Ahr, selection.method = 'mean.var.plot', mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf))
length(x = VariableFeatures(object = Ahr))
Ahr <- ScaleData(object = Ahr, features = rownames(x = Ahr), vars.to.regress = c("nCount_RNA", "percent.mito"))
saveRDS(Ahr, "Ahr_integrated.rds")

## further analyses were performed on organs of interest separately ##

#UMAPs, DEGs were generated following the removal of mitochondrial, ribosomal, hemoglobin, genemodel (Gm) annotated genes, and the sex-specific genes (Tsix, Xist).
DefaultAssay(Ahr_organ) <- "RNA"
unwanted_pattern <- paste0("^mt-|","^Rps|","^Rpl|", "^Hbb|", "^Gm|","Rik$|","^Tsix$|", "^Xist$")
genes_to_remove <- grep(pattern = unwanted_pattern, x = rownames(Ahr_organ),value = TRUE)
genes_to_keep <- setdiff(rownames(Ahr_organ),genes_to_remove)
Ahr_organ_noRibo <- subset(Ahr_organ,features = genes_to_keep)
stopifnot( !any(grepl(unwanted_pattern, rownames(Ahr_organ_noRibo))))# Verify that no unwanted genes remain



# GSEA using the Gene Ontology (GO) database. 
Ahr_lp -> Ahr 

Layers(Ahr[["RNA"]])
options(spe = c("mouse"))
Ahr <- GeneSetAnalysisGO(Ahr, parent = "GO:0002376")
matr <- Ahr@misc$AUCell$GO$"GO:0002376"
matr <- RenameGO(matr)
head(matr, 4:3)

GeneSetAnalysisGO()
SeuratExtend::Heatmap(CalcStats(matr, f = Ahr_lp_noNKnoRibo$seurat_clusters, order = "p", n = 3), lab_fill = "zscore")

stats_cluster <- CalcStats(
  matr,
  f = Ahr$seurat_clusters,
  order = "p",
  n = 3
)

features_keep <- unique(rownames(stats_cluster))
cluster_cond <- interaction(
        Ahr$experiment,
        Ahr$seurat_clusters,
  sep = "_"
)

stats_final <- CalcStats(
  matr[features_keep, ],
  f = cluster_cond
)
SeuratExtend::Heatmap(
  stats_final,
  lab_fill = "zscore"
)

WaterfallPlot(matr, f = Ahr$experiment, ident.1 = "KO", ident.2 = "WT", top.n = 5)
p <- WaterfallPlot(matr, f = Ahr$experiment,ident.1 = "KO",ident.2 = "WT",style = "segment", color_theme = "D", top.n = 5, len.threshold = 2)

p + 
  ggtitle(" Enriched Pathways (immune_system_process) in Ahr lp ") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 9),
    axis.text.y = element_text(size = 8)  # Smaller y-axis font
  )


# module score on a specific gene set #
AddModuleScore(
        Ahr,
  features=gene_DCTolerance_merge_GOandLiterature_indata,
  pool = NULL,
  nbin = 24,
  ctrl = 100,
  k = FALSE,
  assay = "RNA",
  name = "Cluster",
  seed = 1,
  search = FALSE,
  layer = "data"
)

Ahr <- AddModuleScore(
        Ahr,
  features = list(gene_DCTolerance_merge_GOandLiterature_indata),  
  name = "gene_DCTolerance_merge_GOandLiterature_indata_Module"
)

VlnPlot(Ahr, features = "gene_DCTolerance_merge_GOandLiterature_indata_Module1", 
        group.by = "seurat_clusters", cols=color_palette_lp_noNK, split.by = "experiment")


# DELTA AHR module (KO-WT)
df <- FetchData(
        Ahr,
  vars = c(
    "gene_DCTolerance_merge_GOandLiterature_indata_Module1", 
    "seurat_clusters",
    "experiment", 
    "sample"     
  )
)

df_mouse <- df %>%
  group_by(sample, experiment, seurat_clusters) %>%
  summarise(
    module_mean = mean(gene_DCTolerance_merge_GOandLiterature_indata_Module1),
    .groups = "drop"
  )

sum_stats <- df_mouse %>%
  group_by(seurat_clusters, experiment) %>%
  summarise(
    mean = mean(module_mean, na.rm = TRUE),
    sd   = sd(module_mean, na.rm = TRUE),
    n    = n(),
    se   = sd / sqrt(n),
    .groups = "drop"
  ) %>%
  select(seurat_clusters, experiment, mean, se) %>%
  pivot_wider(names_from = experiment, values_from = c(mean, se))

df_effect2 <- sum_stats %>%
  mutate(
    delta_KO_WT = mean_KO - mean_WT,
    se_delta    = sqrt(se_KO^2 + se_WT^2)
  ) %>%
  select(seurat_clusters, delta_KO_WT, se_delta)

df_stats <- df_mouse %>%
  group_by(seurat_clusters) %>%
  wilcox_test(module_mean ~ experiment, exact = FALSE) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p, method = "BH")) %>%
  select(seurat_clusters, p, p_adj)

df_plot <- df_effect2 %>%
  left_join(df_stats, by = "seurat_clusters") %>%
  mutate(
    seurat_clusters = factor(seurat_clusters),
    sig = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE ~ ""
    )
  )
df_plot <- df_plot %>%
  arrange(delta_KO_WT) %>%
  mutate(seurat_clusters = factor(seurat_clusters, levels = seurat_clusters))

ggplot(df_plot, aes(x = seurat_clusters, y = delta_KO_WT)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbar(aes(ymin = delta_KO_WT - se_delta,
                    ymax = delta_KO_WT + se_delta),
                width = 0.25) +
  geom_point(size = 3) +
  geom_text(aes(label = sig),
            vjust = -1.1, size = 4) +
  labs(
    x = "Cluster",
    y = expression(Delta*" AHR module (KO - WT)")
  ) +
  theme_classic()+ labs(
    subtitle = "Mean ± SE of KO−WT difference per cluster; n = 3 mice per genotype\nWilcoxon tests per cluster; BH-FDR across clusters"
  )

# Sample-level stats for module score per cluster (KO vs WT)
obj <- Ahr
score_col   <- "gene_DCTolerance_merge_GOandLiterature_indata_Module1" 
cluster_col <- "seurat_clusters"
cond_col    <- "experiment"   
sample_col  <- "sample"       

md <- obj@meta.data

df_sample <- md %>%
  dplyr::select(all_of(c(sample_col, cond_col, cluster_col, score_col))) %>%
  dplyr::filter(!is.na(.data[[score_col]])) %>%
  dplyr::group_by(
    sample    = .data[[sample_col]],
    condition = .data[[cond_col]],
    cluster   = .data[[cluster_col]]
  ) %>%
  dplyr::summarise(
    score   = mean(.data[[score_col]], na.rm = TRUE),
    n_cells = dplyr::n(),
    .groups = "drop"
  )

sample_counts <- df_sample %>%
  group_by(cluster, condition) %>%
  summarise(n_samples = n_distinct(sample), .groups = "drop") %>%
  pivot_wider(
    names_from = condition,
    values_from = n_samples,
    values_fill = 0
  ) %>%
  arrange(cluster)
print(sample_counts)

eligible_clusters <- sample_counts %>%
  filter(KO >= 3, WT >= 3) %>%
  pull(cluster)

message("Clusters retained for analysis: ",
        paste(eligible_clusters, collapse = ", "))

df_sample_filt <- df_sample %>%
  filter(cluster %in% eligible_clusters)

stats_tbl <- map_dfr(sort(unique(df_sample_filt$cluster)), function(cl) {
  df_cl <- df_sample_filt %>% filter(cluster == cl)
  tst <- wilcox_test(df_cl, score ~ condition)
  tibble(
    cluster = cl,
    p = tst$p
  )
}) %>%
  mutate(p.adj = p.adjust(p, method = "BH")) %>%
  add_significance("p.adj")

print(stats_tbl)

p <- ggplot(
  df_sample_filt,
  aes(x = factor(cluster), y = score, fill = condition)
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.35,
    position = position_dodge(width = 0.8)
  ) +
  geom_point(
    aes(size = n_cells, color = condition),
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.8
    ),
    alpha = 0.85
  ) +
  scale_fill_manual(values = c("KO" = "darkgrey", "WT" = "firebrick")) +
  scale_color_manual(values = c("KO" = "darkgrey", "WT" = "firebrick")) +
  theme_classic() +
  labs(
    x = "Seurat cluster",
    y = "Mean module score per sample",
    size = "Cells in cluster"
  )

print(p)
stats_tbl

#######################  heatmap to calculate scale.data by cluster
# ============================================================
# DIFFERENTIALLY EXPRESSED GENE-SET HEATMAPS BY CLUSTER
#
# For each cluster, this script:
#   1. compares KO versus WT using the normalized RNA data layer;
#   2. tests only genes from a predefined gene set;
#   3. applies BH correction within that predefined gene set;
#   4. retains significant upregulated and downregulated genes;
#   5. runs ScaleData() jointly on WT and KO cells from that cluster;
#   6. plots cell-level scaled expression and the associated log2FC;
#   7. saves one PDF and one significant-gene table per cluster;
#   8. saves a multipage PDF, all DE results, and a cluster summary.
#
# Interpretation:
#   Positive avg_log2FC = higher in KO
#   Negative avg_log2FC = higher in WT
#
# IMPORTANT:
#   - Differential expression and log2FC are calculated from RNA "data".
#   - RNA "scale.data" is used only for heatmap visualization.
#   - WT and KO are scaled together within each cluster, not separately.
# ============================================================

# 1. LOAD PACKAGES
library(Seurat)
library(dplyr)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(grid)

# 2. USER PARAMETERS
# Seurat object
obj
# Metadata columns
cluster_column <- "seurat_clusters"
condition_column <- "experiment"

# Condition order and fold-change direction
reference_condition <- "WT"
test_condition <- "KO"
condition_levels <- c(reference_condition, test_condition)

# Differential-expression thresholds
padj_threshold <- 0.05
minimum_pct <- 0.10
minimum_cells_per_group <- 3
minimum_cells_expressing_gene <- 3
minimum_abs_log2fc <- 0.25


# Heatmap display limits
expression_limit <- 2
fc_limit <- 1.5
# Fixed dimensions of the expression heatmap body
heatmap_body_width <- unit(85, "mm")
heatmap_body_height <- unit(85, "mm")
# PDF dimensions
pdf_width <- 8.5
pdf_height <- 8.5

# Output directory
output_directory <- "ComplexHeatmaps_byCluster"
dir.create(output_directory,showWarnings = FALSE,recursive = TRUE)

# 3. DEFINE THE PREDEFINED GENE SET
dc_gene_set_focused <- unique(c("Zbtb46","Batf3","Irf8","Xcr1","Clec9a","Itgae","Cadm1",
"Irf4","Rorc", "Sirpa","Clec10a","Cd209a","Clec4a2","Cd40","Cd80","Cd83","Cd86","Relb","H2-DMb1","Ifi30","Icosl","Cd274","Pdcd1lg2",
"Ccr7","Fscn1", "Lamp3","Ccl17","Ccl22","Nrp2","Sema7a", "Aldh1a2", "Itgb8", "Tgfb1", "Tgfb2","Tgfbr1", "Tgfbr2",
"Btla","Cd200", "Cd200r1", "Fas", "Havcr2", "Tnfrsf14", "Ido1", "Il4i1", "Il10", "Il10ra","Il27", "Inhba", "Cblb","Foxo3", "Irak3",
 "Itch","Phlpp1","Pten","Socs1", "Socs2","Socs3", "Tnfaip3", "Maf", "Pparg", "Tnfsf4","Tnfsf9"))
module_genes <- unique(dc_gene_set_focused)


# 4. PREPARE AND VALIDATE THE OBJECT
DefaultAssay(obj) <- "RNA"
required_metadata <- c(cluster_column, condition_column)
missing_metadata <- setdiff(required_metadata, colnames(obj@meta.data))
if (length(missing_metadata) > 0) {stop( "Missing metadata columns: ",paste(missing_metadata, collapse = ", "))}

conditions_present <- unique(as.character(obj[[condition_column, drop = TRUE]]))
if (!all(condition_levels %in% conditions_present)) {stop("The condition column '", condition_column, "' must contain both '", reference_condition,"' and '", test_condition, "'.")}

# Join split RNA layers, if present. This does not rerun normalization.
rna_layers <- Layers(obj[["RNA"]])
if (any(grepl("^(counts|data)\\.", rna_layers))) {obj <- JoinLayers(object = obj, assay = "RNA")}

if (!"data" %in% Layers(obj[["RNA"]])) {stop( "The RNA assay does not contain a normalized 'data' layer. ","Run NormalizeData() before this script if appropriate for your workflow." )}

# Determine which predefined genes are available for DE testing.
normalized_matrix <- LayerData(object = obj,assay = "RNA",layer = "data")

module_genes_present <- intersect(module_genes, rownames(normalized_matrix))
module_genes_missing <- setdiff(module_genes, rownames(normalized_matrix))

cat("\nPredefined genes requested: ", length(module_genes), "\n", sep = "")
cat("Predefined genes present in RNA data: ",length(module_genes_present), "\n",sep = "")

if (length(module_genes_missing) > 0) {cat("\nGenes absent from the normalized RNA data layer:\n")print(module_genes_missing)}
if (length(module_genes_present) == 0) {stop("None of the predefined genes are present in the RNA data layer.")}


# 5. DEFINE CLUSTER ORDER AND OPTIONAL LABELS
cluster_labelsLP <- data.frame(cluster = as.character(0:9),
label = c("Trem1+ cDC2","Rorc+ cDC2", "Cx3cr1+ macrophages","Xcr1+ cDC1","Apoe+ macrophages","pre-cDCs","pDCs","pDCs/cDCs", "Ccr7+ DCs","Ccr2+ cells"),
stringsAsFactors = FALSE)

cluster_levels <- unique(as.character(obj[[cluster_column, drop = TRUE]]))
if (all(grepl("^[0-9]+$", cluster_levels))) {cluster_levels <- as.character(sort(as.integer(cluster_levels)))} else {cluster_levels <- sort(cluster_levels)}

# Optional table with columns named "cluster" and "label".
if (exists("cluster_labelsLP") && all(c("cluster", "label") %in% colnames(cluster_labelsLP))) {
cluster_label_map <- setNames(
    as.character(cluster_labelsLP$label),
    as.character(cluster_labelsLP$cluster) )} 
else {cluster_label_map <- setNames( paste0("Cluster ", cluster_levels),cluster_levels)}
get_cluster_label <- function(cluster_id) {label <- unname(cluster_label_map[as.character(cluster_id)])
  if (length(label) == 0 || is.na(label) || label == "") { label <- paste0("Cluster ", cluster_id)  label}
safe_filename <- function(x) { x <- gsub("[^A-Za-z0-9_-]+", "_", x)gsub("_+", "_", x)}

# 6. DEFINE COLOURS
condition_palette <- c( "WT" = "#FF7A00", "KO" = "#005B96")

# Fall back to generated colours if condition names are changed above.
if (!all(condition_levels %in% names(condition_palette))) {condition_palette <- setNames( c("#FF7A00", "#005B96"), condition_levels)}
up_label <- paste0("Higher in ", test_condition)
down_label <- paste0("Higher in ", reference_condition)

direction_palette <- setNames(c("#B2182B", "#2166AC"),c(up_label, down_label))

# Gene-wise scaled expression.
expression_colors <- circlize::colorRamp2(c(-expression_limit, 0, expression_limit),c("#F000D0", "#080808", "#FFFF00"))

# log2FC: negative = higher in reference; positive = higher in test.
fc_colors <- circlize::colorRamp2(c(-fc_limit, 0, fc_limit),c("royalblue", "white", "firebrick"))

# 7. INITIALIZE OUTPUT OBJECTS
de_results_by_cluster <- list()
summary_by_cluster <- list()
complex_heatmaps_by_cluster <- list()

# 8. LOOP OVER CLUSTERS
for (cl in cluster_levels) {
  
  message("\n========================================")
  message("Analyzing cluster ", cl)
  message("========================================")
  
  cluster_cells <- colnames(obj)[
    as.character(obj[[cluster_column, drop = TRUE]]) == as.character(cl)
  ]
   if (length(cluster_cells) == 0) {
    message("Skipping cluster ", cl, ": no cells found.")
    next
  }
  
  obj_cluster <- subset( x = obj,cells = cluster_cells )
  
  DefaultAssay(obj_cluster) <- "RNA"
  Idents(obj_cluster) <- condition_column
  
  condition_counts <- table(Idents(obj_cluster))
  
  n_reference <- if (reference_condition %in% names(condition_counts)) {
    unname(condition_counts[reference_condition])
  } else {
    0
  }
  
  n_test <- if (test_condition %in% names(condition_counts)) {
    unname(condition_counts[test_condition])
  } else {
    0
  }
  
  message( "Cluster ", cl, ": ",reference_condition, " = ", n_reference, "; ", test_condition, " = ", n_test )
  
  if (
    n_reference < minimum_cells_per_group ||
    n_test < minimum_cells_per_group
  ) {
    message(
      "Skipping cluster ", cl,
      ": insufficient cells in one or both conditions."
    )
    
    summary_by_cluster[[as.character(cl)]] <- data.frame(
      cluster = as.character(cl),
      label = get_cluster_label(cl),
      n_reference = n_reference,
      n_test = n_test,
      n_genes_tested = 0,
      n_significant = 0,
      n_higher_in_test = 0,
      n_higher_in_reference = 0,
      status = "Insufficient cells",
      stringsAsFactors = FALSE
    )
    
    next
  }
  

  # 8A. WITHIN-CLUSTER DIFFERENTIAL EXPRESSION
  # The Wilcoxon test and avg_log2FC use normalized RNA data.
  # All predefined genes are supplied for testing.
 de_cluster <- FindMarkers(
    object = obj_cluster,
    assay = "RNA",
    ident.1 = test_condition,
    ident.2 = reference_condition,
    features = module_genes_present,
    test.use = "wilcox",
    slot = "data",
    logfc.threshold = 0,
    min.pct = 0,
    min.cells.feature = minimum_cells_expressing_gene,
    min.cells.group = minimum_cells_per_group,
    only.pos = FALSE,
    verbose = FALSE
  )
  
  if (nrow(de_cluster) == 0) {
    message("Skipping cluster ", cl, ": FindMarkers returned no genes.")
    
    summary_by_cluster[[as.character(cl)]] <- data.frame(
      cluster = as.character(cl),
      label = get_cluster_label(cl),
      n_reference = n_reference,
      n_test = n_test,
      n_genes_tested = 0,
      n_significant = 0,
      n_higher_in_test = 0,
      n_higher_in_reference = 0,
      status = "FindMarkers returned no genes",
      stringsAsFactors = FALSE
    )
    
    next
  }
  
  de_cluster <- de_cluster %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::mutate(
      cluster = as.character(cl),
      
      # BH correction across the predefined genes tested in this cluster.
      p_val_BH = p.adjust(p_val, method = "BH"),
      
      maximum_pct = pmax(pct.1, pct.2),
      
      direction = dplyr::case_when(
        avg_log2FC > 0 ~ up_label,
        avg_log2FC < 0 ~ down_label,
        TRUE ~ "No change"
      )
    )
  
  de_results_by_cluster[[as.character(cl)]] <- de_cluster
  
  de_cluster %>%
    dplyr::filter(
      !is.na(p_val_BH),
      p_val_BH < padj_threshold,
      maximum_pct >= minimum_pct,
      avg_log2FC != 0
    ) %>%
    dplyr::summarise(
      n_significant_before_FC_filter = dplyr::n(),
      n_below_0.25 = sum(abs(avg_log2FC) < 0.25),
      n_at_or_above_0.25 = sum(abs(avg_log2FC) >= 0.25),
      smallest_abs_log2FC = min(abs(avg_log2FC))
    ) %>%
    print()

  # 8B. SELECT SIGNIFICANT GENES FROM THE PREDEFINED GENE SET
  # minimum_pct is applied after testing, so it does not change
  # which predefined genes enter the BH correction.
   significant_genes_tbl <- de_cluster %>%
    dplyr::filter(
      !is.na(p_val_BH),
      p_val_BH < padj_threshold,
      maximum_pct >= minimum_pct,
      !is.na(avg_log2FC),
      abs(avg_log2FC) >= minimum_abs_log2fc,
      avg_log2FC != 0
    ) %>%
    dplyr::mutate(
      direction = factor(
        direction,
        levels = c(up_label, down_label)
      )
    ) %>%
    dplyr::arrange(
      direction,
      dplyr::desc(abs(avg_log2FC))
    )
  
  n_significant <- nrow(significant_genes_tbl)
  n_higher_in_test <- sum(significant_genes_tbl$direction == up_label)
  n_higher_in_reference <- sum(significant_genes_tbl$direction == down_label)
  
  message(
    "Significant genes: ", n_significant,
    " (", up_label, ": ", n_higher_in_test,
    "; ", down_label, ": ", n_higher_in_reference, ")"
  )
  
  if (n_significant == 0) {
    summary_by_cluster[[as.character(cl)]] <- data.frame(
      cluster = as.character(cl),
      label = get_cluster_label(cl),
      n_reference = n_reference,
      n_test = n_test,
      n_genes_tested = nrow(de_cluster),
      n_significant = 0,
      n_higher_in_test = 0,
      n_higher_in_reference = 0,
      status = "No significant predefined genes",
      stringsAsFactors = FALSE
    )
    
    next
  }
  
  cluster_file_id <- safe_filename(as.character(cl))
  write.csv(significant_genes_tbl,file = file.path( output_directory,paste0("cluster_", cluster_file_id, "_significant_gene_set_DE.csv" ) ), row.names = FALSE)
  
  # 8C. SCALE THE SIGNIFICANT GENES WITHIN THIS CLUSTER
  # ScaleData() is run jointly on all WT and KO cells in the
  # cluster. It centres and scales each gene across those cells.
  # It is used only to create the heatmap values.
  genes_for_heatmap <- significant_genes_tbl$gene
  
  obj_cluster <- ScaleData(
    object = obj_cluster,
    assay = "RNA",
    features = genes_for_heatmap,
    do.center = TRUE,
    do.scale = TRUE,
    scale.max = 10,
    verbose = FALSE
  )
  
  expression_scaled <- LayerData(
    object = obj_cluster,
    assay = "RNA",
    layer = "scale.data"
  )
  
  expression_scaled <- as.matrix(
    expression_scaled[
      genes_for_heatmap,
      cluster_cells,
      drop = FALSE
    ]
  )
  
  # Guard against non-finite values from zero-variance rows.
  expression_scaled[!is.finite(expression_scaled)] <- 0
  
  # Clip only for colour display.
  expression_scaled[expression_scaled > expression_limit] <- expression_limit
  expression_scaled[expression_scaled < -expression_limit] <- -expression_limit
  
  # 8D. ORDER CELLS BY CONDITION
cell_metadata <- obj_cluster@meta.data[ cluster_cells, , drop = FALSE ]
cell_metadata$Condition <- factor( as.character(cell_metadata[[condition_column]]),levels = condition_levels )
cell_order <- rownames(cell_metadata)[ order(cell_metadata$Condition) ]
cell_metadata_ordered <- cell_metadata[cell_order, ,drop = FALSE ]
  
expression_scaled <- expression_scaled[ , cell_order,  drop = FALSE ]
  
  # 8E. PREPARE ROW ANNOTATIONS AND log2FC MATRIX
  
  row_direction <- factor( as.character(significant_genes_tbl$direction),levels = c(up_label, down_label) )
  names(row_direction) <- significant_genes_tbl$gene
  
  fc_matrix <- matrix( significant_genes_tbl$avg_log2FC,
    ncol = 1,
    dimnames = list(
    significant_genes_tbl$gene,
      paste0(test_condition, " vs ", reference_condition)
    )
  )
  
  fc_matrix_plot <- fc_matrix
  fc_matrix_plot[fc_matrix_plot > fc_limit] <- fc_limit
  fc_matrix_plot[fc_matrix_plot < -fc_limit] <- -fc_limit
  
  # 8F. CREATE ANNOTATIONS

  top_annotation <- ComplexHeatmap::HeatmapAnnotation(
    Condition = cell_metadata_ordered$Condition,
    col = list(Condition = condition_palette),
    show_annotation_name = TRUE,
    annotation_name_gp = gpar(
      fontsize = 9,
      fontface = "bold"
    ),
    simple_anno_size = unit(4, "mm"),
    annotation_legend_param = list(
      Condition = list(
        title = "Condition",
        at = condition_levels
      )
    )
  )
  
  direction_annotation <- ComplexHeatmap::rowAnnotation(
    Direction = row_direction,
    col = list(Direction = direction_palette),
    simple_anno_size = unit(4, "mm"),
    annotation_legend_param = list(
      Direction = list(
        title = "Differential\nexpression",
        at = c(up_label, down_label)
      )
    )
  )
  

  # 8G. CREATE THE EXPRESSION HEATMAP
 cluster_label <- get_cluster_label(cl)
  
  heatmap_title <- paste0(
    "Cluster ",
    cl,
    ": ",
    cluster_label,
    "\nSignificant genes from predefined gene set"
  )
  
  expression_heatmap <- ComplexHeatmap::Heatmap(
    expression_scaled,
    name = "Scaled\nexpression",
    col = expression_colors,
    
    # Cluster genes within the upregulated/downregulated sections.
    cluster_rows = TRUE,
    clustering_distance_rows = "euclidean",
    clustering_method_rows = "complete",
    row_split = row_direction,
    cluster_row_slices = FALSE,
    row_gap = unit(2, "mm"),
    
    # Cluster cells separately within WT and KO.
    cluster_columns = TRUE,
    column_split = cell_metadata_ordered$Condition,
    cluster_column_slices = FALSE,
    column_gap = unit(2, "mm"),
    
    top_annotation = top_annotation,
    
    show_column_names = FALSE,
    show_row_names = FALSE,
    
    column_title = heatmap_title,
    column_title_gp = gpar(
      fontsize = 11,
      fontface = "bold"
    ),
    
    row_dend_width = unit(15, "mm"),
    width = heatmap_body_width,
    height = heatmap_body_height,
    border = FALSE,
    
    # Avoid Cairo temporary-raster issues.
    use_raster = FALSE,
    
    heatmap_legend_param = list(
      title = "Scaled\nexpression",
      at = c(-expression_limit, -1, 0, 1, expression_limit),
      labels = c(
        paste0("<=-", expression_limit),
        "-1",
        "0",
        "1",
        paste0(">=", expression_limit)
      )
    )
  )
  

  # 8H. CREATE THE log2FC PANEL
fc_heatmap <- ComplexHeatmap::Heatmap(
    
    fc_matrix_plot,
    name = paste0("log2FC\n", test_condition, " vs ", reference_condition),
    col = fc_colors,
    
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    
    show_row_names = TRUE,
    row_names_side = "right",
    row_names_gp = gpar(fontsize = 8),
    
    show_column_names = TRUE,
    column_names_gp = gpar(
      fontsize = 8,
      fontface = "bold"
    ),
    
    width = unit(8, "mm"),
    
    rect_gp = gpar(
      col = "white",
      lwd = 0.5
    ),
    
    heatmap_legend_param = list(
      title = paste0(
        "log2FC\n",
        test_condition,
        " vs ",
        reference_condition
      ),
      at = c(-fc_limit, 0, fc_limit),
      labels = c(
        paste0("<=-", fc_limit),
        "0",
        paste0(">=+", fc_limit)
      )
    )
  )
  

  # 8I. COMBINE AND STORE THE HEATMAP
 combined_heatmap <-
    direction_annotation +
    expression_heatmap +
    fc_heatmap
  
  complex_heatmaps_by_cluster[[as.character(cl)]] <- combined_heatmap
  
  summary_by_cluster[[as.character(cl)]] <- data.frame(
    cluster = as.character(cl),
    label = cluster_label,
    n_reference = n_reference,
    n_test = n_test,
    n_genes_tested = nrow(de_cluster),
    n_significant = n_significant,
    n_higher_in_test = n_higher_in_test,
    n_higher_in_reference = n_higher_in_reference,
    status = "Heatmap created",
    stringsAsFactors = FALSE
  )
}

# 9. SAVE COMBINED RESULT TABLES
if (length(de_results_by_cluster) > 0) {
  all_de_results <- dplyr::bind_rows(de_results_by_cluster)
  
  write.csv(
    all_de_results,
    file = file.path(
      output_directory,
      "all_clusters_gene_set_DE_results_July27.csv"
    ),
    row.names = FALSE
  )
}

if (length(summary_by_cluster) > 0) {
  cluster_summary <- dplyr::bind_rows(summary_by_cluster)
  
  write.csv(
    cluster_summary,
    file = file.path(
      output_directory,
      "cluster_heatmap_summary.csv"
    ),
    row.names = FALSE
  )
}


# 10. SAVE ONE PDF PER CLUSTER
graphics.off()

if (length(complex_heatmaps_by_cluster) > 0) {
  
  for (cl in names(complex_heatmaps_by_cluster)) {
    
    cluster_label <- get_cluster_label(cl)
    
    individual_pdf <- file.path(
      output_directory,
      paste0(
        "cluster_",
        safe_filename(cl),
        "_",
        safe_filename(cluster_label),
        "_significant_gene_set_heatmap_July25.pdf"
      )
    )
    
    grDevices::pdf(
      file = individual_pdf,
      width = pdf_width,
      height = pdf_height,
      useDingbats = FALSE
    )
    
    tryCatch(
      {
        ComplexHeatmap::draw(
          object = complex_heatmaps_by_cluster[[cl]],
          merge_legends = TRUE,
          heatmap_legend_side = "right",
          annotation_legend_side = "right",
          newpage = TRUE
        )
      },
      finally = {
        grDevices::dev.off()
      }
    )
    
    message(
      "Saved: ",
      normalizePath(individual_pdf)
    )
  }
}

# 11. SAVE A MULTIPAGE PDF CONTAINING ALL CLUSTERS
if (length(complex_heatmaps_by_cluster) > 0) {
  
  combined_pdf <- file.path(
    output_directory,
    "all_clusters_significant_gene_set_heatmaps.pdf"
  )
  
  grDevices::pdf(
    file = combined_pdf,
    width = pdf_width,
    height = pdf_height,
    onefile = TRUE,
    useDingbats = FALSE
  )
  
  tryCatch(
    {
      for (cl in names(complex_heatmaps_by_cluster)) {
        
        ComplexHeatmap::draw(
          object = complex_heatmaps_by_cluster[[cl]],
          merge_legends = TRUE,
          heatmap_legend_side = "right",
          annotation_legend_side = "right",
          newpage = TRUE
        )
      }
    },
    finally = {
      grDevices::dev.off()
    }
  )
  
  message("Combined PDF saved to: ", normalizePath(combined_pdf))
}

##################
