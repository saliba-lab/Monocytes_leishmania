################################################
# Script Description ====
# This script describes the analysis of the single cell dataset with Seurat
# for further analysis
# MAIN STEPS:
# 1- Load normalized count matrix and setup a Seurat object
# 2- Perform a clustering analysis including batch correction with Harmony
# 3- 

# Set working directory ====
dir <- "."
setwd(dir)

# Load data ====
norm_mtx <- read.delim("../outputs/scNorm/Normalized_mtx.tsv", check.names = F)
#colnames(norm_mtx) <- sub("\\_S\\d+$", "", colnames(norm_mtx))

metadata <- read.csv("../outputs/metadata_table.csv")
metadata <- metadata[metadata$cell %in% colnames(norm_mtx), ] #keep info only about cells passing QC
rownames(metadata) <- metadata$cell

#sanity check; output 0 if cell are correctly ordered in both mtx and metadata
print(
  paste0(
    "Sanity check. Number of cell names different between normalized count matrix and metadata table : ", 
    sum(colnames(norm_mtx) != metadata$cell)
  )
)

# Setup SeuratObject ====
ds <- Seurat::CreateSeuratObject(
  counts = norm_mtx, 
  project = "Mono_Leishmania", 
  meta.data = metadata
)

# Fill data slot with log-transformed normalized counts
ds <- Seurat::SetAssayData(
  ds, 
  assay = "RNA", 
  slot = "data", 
  new.data = as.matrix(log1p(ds@assays$RNA@counts))
) 

#ds$mouse_per_ear <- paste0(ds$mouse, "_", ds$mouse_ear) # !!!!!!!!!!!!!!!!!!
#ds <- subset(ds, subset = mouse_per_ear != "Mouse_3_R" & mouse_per_ear != "Mouse_3_L")


# Highly variable features selection ====
#HVF <- Seurat::VariableFeatures(ds, assay = "RNA")
ds <- Seurat::FindVariableFeatures(ds, assay = "RNA", selection.method = "vst", nfeatures = 1000) #determine highly variable features
Seurat::VariableFeaturePlot(ds, assay = "RNA")
Seurat::LabelPoints(
  plot = Seurat::VariableFeaturePlot(ds, assay = "RNA"), 
  points = head(Seurat::VariableFeatures(ds, assay = "RNA"), 20), 
  repel = TRUE
)

# Skip scaling step ====
ds <- Seurat::SetAssayData(
  ds, 
  assay = "RNA", 
  slot = "scale.data", 
  new.data = as.matrix(ds@assays$RNA@data)
) 

# PCA ====
ds <- Seurat::RunPCA(ds)
Seurat::ElbowPlot(ds)

# Fix for Harmony function ====
# See Harmony Github -> https://github.com/immunogenomics/harmony/issues/187

RunHarmony.Seurat <- function(
  object,
  group.by.vars,
  reduction = 'pca',
  dims.use = NULL,
  theta = NULL,
  lambda = NULL,
  sigma = 0.1,
  nclust = NULL,
  tau = 0,
  block.size = 0.05,
  max.iter.harmony = 10,
  max.iter.cluster = 20,
  epsilon.cluster = 1e-5,
  epsilon.harmony = 1e-4,
  plot_convergence = FALSE,
  verbose = TRUE,
  reference_values = NULL,
  reduction.save = "harmony",
  assay.use = NULL,
  project.dim = TRUE,
  ...
) {
  if (!requireNamespace('Seurat', quietly = TRUE)) {
    stop("Running Harmony on a Seurat object requires Seurat")
  }
  assay.use <- assay.use %||% Seurat::DefaultAssay(object)
  if (reduction == "pca" && !reduction %in% Seurat::Reductions(object = object)) {
    if (isTRUE(x = verbose)) {
      message("Harmony needs PCA. Trying to run PCA now.")
    }
    object <- tryCatch(
      expr = Seurat::RunPCA(
        object = object,
        assay = assay.use,
        verbose = verbose,
        reduction.name = reduction
      ),
      error = function(...) {
        stop("Harmony needs PCA. Tried to run PCA and failed.")
      }
    )
  }
  if (!reduction %in% Seurat::Reductions(object = object)) {
    stop("Requested dimension reduction is not present in the Seurat object")
  }
  embedding <- Seurat::Embeddings(object, reduction = reduction)
  if (is.null(dims.use)) {
    dims.use <- seq_len(ncol(embedding))
  }
  dims_avail <- seq_len(ncol(embedding))
  if (!all(dims.use %in% dims_avail)) {
    stop("trying to use more dimensions than computed. Rereun dimension reduction
         with more dimensions or run Harmony with fewer dimensions")
  }
  if (length(dims.use) == 1) {
    stop("only specified one dimension in dims.use")
  }
  metavars_df <- Seurat::FetchData(
    object,
    group.by.vars,
    cells = Seurat::Cells(x = object[[reduction]])
  )
  
  harmonyEmbed <- HarmonyMatrix(
    embedding,
    metavars_df,
    group.by.vars,
    FALSE,
    0,
    theta,
    lambda,
    sigma,
    nclust,
    tau,
    block.size,
    max.iter.harmony,
    max.iter.cluster,
    epsilon.cluster,
    epsilon.harmony,
    plot_convergence,
    FALSE,
    verbose,
    reference_values
  )
  
  #reduction.key <- Seurat::Key(reduction.save, quiet = TRUE)
  reduction.key <- reduction.save # HERE IS THE FIX 
  rownames(harmonyEmbed) <- rownames(embedding)
  colnames(harmonyEmbed) <- paste0(reduction.key, seq_len(ncol(harmonyEmbed)))
  
  object[[reduction.save]] <- Seurat::CreateDimReducObject(
    embeddings = harmonyEmbed,
    stdev = as.numeric(apply(harmonyEmbed, 2, stats::sd)),
    assay = Seurat::DefaultAssay(object = object[[reduction]]),
    key = reduction.key
  )
  if (project.dim) {
    object <- Seurat::ProjectDim(
      object,
      reduction = reduction.save,
      overwrite = TRUE,
      verbose = FALSE
    )
  }
  return(object)
}

# Harmony correction ====
#ds$mouse <- factor(ds$mouse, levels = unique(ds$mouse)[order(unique(ds$mouse))])
library(Seurat)
library(harmony)

ds <- harmony::RunHarmony(
  object = ds,
  group.by.vars = "mouse",
  reduction = "pca",
  reduction.save = "harmony"
)

# UMAP ====
ds <- Seurat::RunUMAP(
  ds, reduction = "harmony",
  dims = 1:8,
  n.neighbors = 10,
  min.dist = 0.5
)

# Leiden clustering ====
ds <- Seurat::FindNeighbors(ds, reduction = "harmony", dims = 1:8, k.param = 10)

adj.matrix <- ds@graphs$RNA_snn
adj.matrix <- as(object = adj.matrix, Class = "dgCMatrix")
igraph <- igraph::graph_from_adjacency_matrix(adjmatrix = adj.matrix, weighted = TRUE)

ds@meta.data$Leiden_clusters <- as.factor(
  leidenbase::leiden_find_partition(
    igraph               = igraph,
    resolution_parameter = 0.08,
    seed                 = 29
  )$membership
)

# Some plots to explore dataset ====
Seurat::DimPlot(ds, reduction = "umap", group.by = "Leiden_clusters", pt.size = 2)
Seurat::DimPlot(ds, reduction = "umap", group.by = "infection_status", pt.size = 2)
Seurat::DimPlot(ds, reduction = "umap", group.by = "mouse", pt.size = 2)


Seurat::FeaturePlot(ds, reduction = "umap", features = "Il7r", pt.size = 2) + 
  ggplot2::scale_colour_viridis_c(option = "magma", direction = -1)

Seurat::FeaturePlot(ds, reduction = "umap", features = "H2-Aa", pt.size = 2) + 
  ggplot2::scale_colour_viridis_c(option = "magma", direction = -1)

Seurat::FeaturePlot(ds, reduction = "umap", features = "Ccr2", pt.size = 2) + 
  ggplot2::scale_colour_viridis_c(option = "magma", direction = -1)

Seurat::FeaturePlot(ds, reduction = "umap", features = "Tmem176a", pt.size = 2) + 
  ggplot2::scale_colour_viridis_c(option = "magma", direction = -1)

# Annotation of Leiden clusters ====
lookup_annotation <- tibble::tibble(
  Leiden = c(1, 2, 3),
  Cluster_paper = c("Cluster1", "Cluster3", "Cluster2"),
  Annotation = c("Hi1/U1", "Lo3", "Hi2/U2") 
)

ds$Cluster_number_paper <- sapply(ds$Leiden_clusters, function(x) lookup_annotation$Cluster_paper[match(x, lookup_annotation$Leiden)])
ds$Annotation <- sapply(ds$Leiden_clusters, function(x) lookup_annotation$Annotation[match(x, lookup_annotation$Leiden)])
ds$Annotation <- factor(ds$Annotation, levels = c("Hi1/U1", "Hi2/U2", "Lo3"))

Seurat::DimPlot(ds, reduction = "umap", group.by = "Annotation", pt.size = 2)

# Save annotated Seurat Object ====
saveRDS(ds, "../outputs/Annotated_SeuObj.rds")

# Testing for differentially expressed genes ====
deg_output <- "../outputs/DEG_testing/"
dir.create(deg_output)
Seurat::Idents(ds) <- "Annotation"

# DEG testing with FindAllMarkers
deg <- Seurat::FindAllMarkers(ds)

write.csv( 
  deg[deg$p_val_adj < 0.05, ],
  paste0(deg_output, "DEG_FindAllMarkers.csv")
)

# DEG testing for Lo3 (Cluster3) vs Hi1/U1 and Hi2/U2 (Cluster1 + Cluster2)
deg <- Seurat::FindMarkers(
  ds,
  ident.1 = "Lo3",
  ident.2 = c("Hi2/U2","Hi1/U1"),
  logfc.threshold = 0,
  min.pct = 0,
  min.diff.pct = -Inf,
  test.use = "wilcox",
  only.pos = F,
  return.thresh = Inf
)
deg$gene <- rownames(deg)

write.csv( #Save all DEGs for Lo3
  deg,
  paste0(deg_output, "DEG_all_genes_Lo3_vs_all.csv")
)

write.csv( #Save only the upregulated genes in Lo3
  deg[deg$avg_log2FC > 0 & deg$p_val_adj < 0.05, ],
  paste0(deg_output, "DEG_upregulated_genes_Lo3_vs_all.csv")
)

# DEG testing for Hi1/U1 (Cluster1) vs Hi2/U2 and Lo3 (Cluster2 + Cluster3)
deg <- Seurat::FindMarkers(
  ds,
  ident.1 = "Hi1/U1",
  ident.2 = c("Hi2/U2","Lo3"),
  logfc.threshold = 0,
  min.pct = 0,
  min.diff.pct = -Inf,
  test.use = "wilcox",
  only.pos = F,
  return.thresh = Inf
)
deg$gene <- rownames(deg)

write.csv( #Save all DEGs for Hi1/U1
  deg,
  paste0(deg_output, "DEG_all_genes_Hi1U1_vs_all.csv")
)

write.csv( #Save only the upregulated genes in Hi1/U1
  deg[deg$avg_log2FC > 0 & deg$p_val_adj < 0.05, ],
  paste0(deg_output, "DEG_upregulated_genes_Hi1U1_vs_all.csv")
)

# DEG testing for Hi2/U2 (Cluster2) vs Hi1/U1 and Lo3 (Cluster1 + Cluster3)
deg <- Seurat::FindMarkers(
  ds,
  ident.1 = "Hi2/U2",
  ident.2 = c("Hi1/U1","Lo3"),
  logfc.threshold = 0,
  min.pct = 0,
  min.diff.pct = -Inf,
  test.use = "wilcox",
  only.pos = F,
  return.thresh = Inf
)
deg$gene <- rownames(deg)

write.csv( #Save all DEGs for Lo3
  deg,
  paste0(deg_output, "DEG_all_genes_Hi2U2_vs_all.csv")
)

write.csv( #Save only the upregulated genes in Lo3
  deg[deg$avg_log2FC > 0 & deg$p_val_adj < 0.05, ],
  paste0(deg_output, "DEG_upregulated_genes_Hi2U2_vs_all.csv")
)

# Prepare DEG list compatible with Gene Set Enrichment Analysis (GSEA) ====
# Create function to produce a ranked gene list
rank_DEG_for_GSEA <- function(DEG_table, name, col_pval, col_sign) {
  GSEA_score <- -log10(DEG_table[, col_pval])*sign(DEG_table[, col_sign])
  DEG_table <- cbind(DEG_table, GSEA_score)
  DEG_table <- DEG_table %>% dplyr::select(gene, GSEA_score)
  DEG_table <- DEG_table[order(DEG_table$GSEA_score, decreasing = T), ]
}

# Generate ranked DEG list for Lo3
write.table( 
  rank_DEG_for_GSEA(
    DEG_table = read.csv(
      paste0(deg_output, "DEG_all_genes_Lo3_vs_all.csv")
    ),
    col_pval = "p_val",
    col_sign = "avg_log2FC"
  ),
  file = paste0(deg_output, "GSEA_list_Lo3_vs_all.rnk"),
  sep = "\t",
  col.names = F,
  row.names = F
)

# Generate ranked DEG list for Hi1/U1
write.table( 
  rank_DEG_for_GSEA(
    DEG_table = read.csv(
      paste0(deg_output, "DEG_all_genes_Hi1U1_vs_all.csv")
    ),
    col_pval = "p_val",
    col_sign = "avg_log2FC"
  ),
  file = paste0(deg_output, "GSEA_list_Hi1U1_vs_all.rnk"),
  sep = "\t",
  col.names = F,
  row.names = F
)

# Generate ranked DEG list for Hi2/U2
write.table( 
  rank_DEG_for_GSEA(
    DEG_table = read.csv(
      paste0(deg_output, "DEG_all_genes_Hi2U2_vs_all.csv")
    ),
    col_pval = "p_val",
    col_sign = "avg_log2FC"
  ),
  file = paste0(deg_output, "GSEA_list_Hi2U2_vs_all.rnk"),
  sep = "\t",
  col.names = F,
  row.names = F
)

# Dotplot with marker genes for each cluster ====
# Select marker genes 
genes <- c(
  "Il1b", "Ccl5", "Cxcl9", "Ccr2",
  "Cd84", "Cxcr4", "Cd63", "Pdpn",
  "Il7r", "Cd274", "Gpnmb", 'Ctsb', 
  "Ctsd", "Ctss", "Ctsk","H2-DMb1",
  "H2-Ab1", "H2-Aa", "H2-Eb1"
)

# Prepare data for Dotplot
ids <- rownames(ds@assays$RNA@data) %in% genes
data <- t(as.matrix(ds@assays$RNA@data[ids, ]))
data <- dplyr::as_tibble(data)
data <- dplyr::select(data, genes) #reorder columns
data$Leiden_clusters <- ds@meta.data$Annotation

data$Leiden_clusters <- factor(
  x = data$Leiden_clusters, 
  levels = c("Lo3", "Hi2/U2", "Hi1/U1")
)

data <- tidyr::gather(data, "Gene", "Expression", -Leiden_clusters)

data$Gene <- factor(
  x      = data$Gene,
  levels = unique(data$Gene)
)

data$Pct <- data$Expression > 0
data$Freq <- rep(1, length(data$Pct))

data <- dplyr::group_by(data, Gene)
data <- dplyr::mutate(data, Scaled = scale(as.numeric(Expression))[, 1])
data <- dplyr::group_by(data, Leiden_clusters, Gene)

data <- dplyr::summarise(
  data,
  Mean   = mean(as.numeric(Expression)),
  Scaled = mean(Scaled),
  Pct    = sum(Pct)/sum(Freq)*100
)

# Dotplot
colorscale <- RColorBrewer::brewer.pal(n = 6, name = "RdBu")

dotplot <- ggplot2::ggplot(
  data = data,
  mapping = ggplot2::aes(
    y = Leiden_clusters,
    x = Gene,
    size = Pct,
    col = Scaled
  )
) +
  ggplot2::geom_point() +
  ggplot2::scale_color_gradient2(midpoint=0, 
                                 high=colorscale[length(colorscale)], 
                                 mid="white",
                                 low=colorscale[1]
  ) +
  ggplot2::theme_classic(base_size = 15) +
  ggplot2::theme(
    axis.title       = ggplot2::element_blank(),
    axis.text.x      = ggplot2::element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y      = ggplot2::element_text(
      face  = "italic"
    ),
    legend.position  = "right",
    legend.title     = ggplot2::element_blank()
  ) +
  ggplot2::guides(
    color = ggplot2::guide_colorbar(
      barheight = 11, frame.colour = "black", ticks = FALSE, order = 1
    ),
    size  = ggplot2::guide_legend(label.position = "right")
  ) +
  ggplot2::scale_size_area(max_size = 8)

dotplot
