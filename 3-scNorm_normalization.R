################################################
# Script Description ====
# This script describes the normalization of the count matrix with scNorm
# MAIN STEPS:
# 1- Load count matrix and metadata table
# 2- Exclude Kikume, Leishmania and lowly detected genes
# 3- Normalize count matrix with SCnorm (gene length and plate number are considered during the procedure)

# Set working directory ====
dir <- "."
setwd(dir)

# Load count matrix and metadata from QC ====
mtx <- read.delim("../data/Raw_count_mtx.tsv", check.names = F)
mouse_gene_conversion <- read.csv("../outputs/gene_conversion_table.csv")
metadata <- read.csv("../outputs/metadata_table.csv")

# Extract gene info from count matrix ====
gene_info <- tibble::tibble(
  gene_id = mtx[, "Gene_id"],
  gene_length = mtx[, "Gene_length"]
) 

mtx <- mtx[, -c(1,2)] #remove gene_id and gene length information from the count matrix
colnames(mtx) <- sub("_S\\d+$", "", colnames(mtx)) 

# Sanity check for identical cell names in the count matrix and metadata table
print(
  paste0(
    "Sanity check. Number of cell names different between count matrix and metadata table : ", 
    sum(colnames(mtx) != metadata$cell)
  )
)

# Filter out low quality cells ===
low_quality_cells <- which(metadata$pass_qc == "N")
mtx <- mtx[, -low_quality_cells]
metadata <- metadata[-low_quality_cells, ]

# Filter out mkikume and Leishmania genes ====
mkikume_gene <- grep("^mKiku", gene_info$gene_id)
leishmania_gene <- c(
  grep("^LMJF", gene_info$gene_id), 
  grep("^ENSRNA", gene_info$gene_id)
)

gene_to_filter <- c(mkikume_gene, leishmania_gene)

mtx <- mtx[-gene_to_filter, ]
gene_info <- gene_info[-gene_to_filter, ]

# Filter out undetected genes based on chosen thresholds; but save all ERCC genes ====
min.reads <- 5
min.cells <- 10

ercc_gene <- grep("^ERCC", gene_info$gene_id)
detected_genes <- apply(mtx[-ercc_gene, ] >= min.reads, MARGIN = 1, sum) >= min.cells
detected_genes <- c(detected_genes, rep(TRUE, length(ercc_gene)))

mtx <- mtx[detected_genes, ]
gene_info <- gene_info[detected_genes, ]

# Update gene_info with the mouse gene symbol names and ERCC names ====
gene_info$gene_name <- sapply(
  gene_info$gene_id, 
  function(x) mouse_gene_conversion$symbol[match(x, mouse_gene_conversion$gene_id)]
)

# Normalization with scNorm ====
# Investigate biases based on the count - depth relationship
SCnorm::plotCountDepth(as.matrix(mtx))

# Investigate gene length biases
rownames(mtx) <- gene_info$gene_name

gene_length <- gene_info$gene_length
names(gene_length) <- gene_info$gene_name

SCnorm::plotWithinFactor(Data = mtx, withinSample = gene_length) #plot gene length biases

# Set up a SingleCell Experiment Object with spike-ins information ====
sce <- SingleCellExperiment::SingleCellExperiment(list(counts = mtx))
myspikes <- grepl("^ERCC-", rownames(sce))
sce <- SingleCellExperiment::splitAltExps(sce, ifelse(myspikes, "ERCC", "gene"))

# Update gene_length to remove ERCC genes
gene_length <- gene_info$gene_length[-grep("^ERCC-", gene_info$gene_name)]
names(gene_length) <- gene_info$gene_name[-grep("^ERCC-", gene_info$gene_name)]

# Perform Normalization ====
Conditions <- metadata$plate #set condition based on plate

DataNorm <- SCnorm::SCnorm(
  Data = sce, 
  Conditions = Conditions, 
  withinSample = gene_length,
  useSpikes = T,
  FilterCellNum = 10, 
  PrintProgressPlots = TRUE, 
  NCores=3, 
  reportSF = TRUE, 
  useZerosToScale = FALSE
)
# CAUTION: the useZerosToScale across condition argument should be set to T or F according to the intended DE testing. 
# See SCnorm github / guide.

# Check normalization results
SCnorm::plotCountDepth(sce, NormalizedData = SingleCellExperiment::normcounts(DataNorm))
SCnorm::plotWithinFactor(Data = SingleCellExperiment::normcounts(DataNorm), withinSample = gene_length)

# Save scNorm outputs ====
# Save normalized matrix
dir.create("../outputs/scNorm")
write.table(SingleCellExperiment::normcounts(DataNorm), sep = "\t", "../outputs/scNorm/Normalized_mtx.tsv")

# Get Genes filtered out 
genes_not_normalized <- SCnorm::results(DataNorm, type ="GenesFilteredOut")
genes_not_normalized <- list2DF(genes_not_normalized)
write.csv(genes_not_normalized, "../outputs/scNorm/Genes_not_normalized_during_scNorm.csv", row.names = F)

# End of document 
################################################
