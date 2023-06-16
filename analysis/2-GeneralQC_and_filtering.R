################################################
# Script Description ====
# This script describes the processing of the raw count matrix and general quality control of the dataset
# MAIN STEPS:
# 1- Loading count matrix and extracting gene information (name + length)
# 2- Converting gene names (ensembl IDs) to Symbols
# 3- Assembling the metadata table by computing various metrics and from Info_experiment.csv
# 4- Explore QC metrics to determine which cells can be kept for further analysis

# Set working directory ====
dir <- "."
setwd(dir)

# Load count matrix ====
mtx <- read.delim("../data/Raw_count_mtx.tsv", check.names = F)

gene_info <- dplyr::tibble(
  gene_id = mtx[, "Gene_id"],
  gene_length = mtx[, "Gene_length"]
)

mtx <- mtx[, -c(1,2)] #remove gene_id and gene length information from count matrix

# Get Leishmania / Mouse / ERCC genes ====
unique(stringr::str_sub(gene_info$gene_id, start = 1, end = 5)) #show the different name patterns of genes

mkikume_gene <- grep("^mKiku", gene_info$gene_id)

leishmania_gene <- c(
  grep("^LMJF", gene_info$gene_id), 
  grep("^ENSRNA", gene_info$gene_id)
)

mouse_gene <- grep("^ENSMUS", gene_info$gene_id)

ercc_gene <- grep("^ERCC", gene_info$gene_id)

# Create a lookup table with BiomaRt: Ensembl Ids to gene symbols ====
ensembl <- biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
gene_biomart <- biomaRt::getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"), 
  filters = "ensembl_gene_id", 
  values = sub("\\.\\d+$", "", gene_info$gene_id[mouse_gene]), #need to remove the ".XX" pattern at the end of each ensembl ID
  mart = ensembl
)

gene_info$trimmed_id <- sub("\\.\\d+$", "", gene_info$gene_id)

# Deal with ensembl IDs without symbol names
empty_symbol <- which(gene_biomart$external_gene_name == "")
gene_biomart$external_gene_name[empty_symbol] <- gene_biomart$ensembl_gene_id[empty_symbol]

# Deal with duplicated symbol names by pasting together Symbol and EnsemblID
dup_names <- which(duplicated(gene_biomart$external_gene_name) | duplicated(gene_biomart$external_gene_name, fromLast = TRUE))

gene_biomart$external_gene_name[dup_names] <- paste0(
  gene_biomart$external_gene_name[dup_names],
  "_",
  gene_biomart$ensembl_gene_id[dup_names]
)

# Add gene name symbol based on biomaRt conversion ===
gene_info$symbol <- sapply(gene_info$trimmed_id, function(x) gene_biomart$external_gene_name[match(x, gene_biomart$ensembl_gene_id)])

# Deal with NAs coming from 1) mouse genes not found in biomaRt or 2) non-mouse genes 
gene_info$symbol[is.na(gene_info$symbol)] <- gene_info$gene_id[is.na(gene_info$symbol)]

# Load experimental info ====
exp_info <- read.csv("../data/Info_experiment.csv")

# Create metadata: compute various metrics and add experimental information ====
metadata <- tibble::tibble(
  cell = sub("\\_S\\d+$", "", colnames(mtx)),
  orig.ident = "Mono_Leishmania",
  plate = stringr::str_sub(colnames(mtx), start = 1, end = 1),
  total_counts = colSums(mtx),
  mouse_percent = colSums(mtx[mouse_gene, ])/colSums(mtx) * 100,
  ercc_percent = colSums(mtx[ercc_gene, ])/colSums(mtx) * 100,
  mito_percent = colSums(mtx[grep("^mt-", gene_info$symbol), ])/colSums(mtx) * 100,
  leishmania_percent = colSums(mtx[leishmania_gene, ])/colSums(mtx) * 100,
  leishmania_counts = colSums(mtx[leishmania_gene, ]),
  mkikume_counts = colSums(mtx[mkikume_gene, ]),
  detected_genes = apply(mtx[mouse_gene, ], 2, function(x) sum(x >= 1)) # number of genes with at least 1 count
)

#Add metadata from exp_info
metadata$infection_status <- sapply(metadata$cell, function(x) exp_info$Infection_status[match(x, exp_info$Cell)])
metadata$mouse <- sapply(metadata$cell, function(x) exp_info$Mouse[match(x, exp_info$Cell)])
metadata$mouse_ear <- sapply(metadata$cell, function(x) exp_info$Ear[match(x, exp_info$Cell)])
metadata$Cd11c_fluo <- sapply(metadata$cell, function(x) exp_info$Cd11c_fluorescence_intensity[match(x, exp_info$Cell)])
metadata$Ly6c_fluo <- sapply(metadata$cell, function(x) exp_info$Ly6c_fluorescence_intensity[match(x, exp_info$Cell)])
metadata$KikumeGreen_fluo <- sapply(metadata$cell, function(x) exp_info$KikumeGreen_fluorescence_intensity[match(x, exp_info$Cell)])
metadata$KikumeRed_fluo <- sapply(metadata$cell, function(x) exp_info$KikumeRed_fluorescence_intensity[match(x, exp_info$Cell)])

# Explore QC metrics ====
ggplot2::theme_set(cowplot::theme_cowplot())
data <- metadata

ggplot2::ggplot(data, ggplot2::aes(x = mouse_percent, y = ercc_percent, col = infection_status)) + ggplot2::geom_point()
ggplot2::ggplot(data, ggplot2::aes(x = detected_genes, y = ercc_percent)) + ggplot2::geom_point()
ggplot2::ggplot(data, ggplot2::aes(x = ercc_percent, y = mito_percent, col = infection_status)) + ggplot2::geom_point()
ggplot2::ggplot(data, ggplot2::aes(x = detected_genes, y = mito_percent, col = infection_status)) + ggplot2::geom_point()

ggplot2::ggplot(data, ggplot2::aes(x = orig.ident, y = mouse_percent)) + ggplot2::geom_violin()
ggplot2::ggplot(data, ggplot2::aes(x = orig.ident, y = mito_percent)) + ggplot2::geom_violin()
ggplot2::ggplot(data, ggplot2::aes(x = orig.ident, y = ercc_percent)) + ggplot2::geom_violin()
ggplot2::ggplot(data, ggplot2::aes(x = orig.ident, y = leishmania_percent, col = infection_status)) + ggplot2::geom_boxplot()

# Decision plots
color_pal <- scales::hue_pal()(length(unique(data$infection_status)))

p1 <- ggplot2::ggplot(data, ggplot2::aes(x = infection_status, y = detected_genes, fill = infection_status)) + ggplot2::geom_violin() +
  ggplot2::geom_hline(yintercept = 1000, linetype = "longdash") +
  ggplot2::scale_fill_manual(values = color_pal) +
  ggplot2::xlab("")

p2 <- ggplot2::ggplot(data, ggplot2::aes(x = infection_status, y = mito_percent, fill = infection_status)) + 
  ggplot2::geom_violin() +
  ggplot2::geom_hline(yintercept = 2, linetype = "longdash") +
  ggplot2::scale_fill_manual(values = color_pal) +
  ggplot2::xlab("")

p3 <- ggplot2::ggplot(data, ggplot2::aes(x = infection_status, y = ercc_percent, fill = infection_status)) +
  ggplot2::geom_violin() + 
  ggplot2::geom_hline(yintercept = 20, linetype = "longdash") +
  ggplot2::scale_fill_manual(values = color_pal) +
  ggplot2::xlab("")

dir.create("../outputs")
ggplot2::ggsave(
  filename = "../outputs/QC_metrics.pdf", 
  plot = gridExtra::marrangeGrob(list(p1, p2, p3), nrow=1, ncol=1), 
  width = 7, height = 7
)

# Final decision ====
cells_to_keep <- subset(metadata, subset = mito_percent < 2 & ercc_percent < 20 & detected_genes > 1000)

metadata$pass_qc <- "N"
metadata$pass_qc[data$cell %in% cells_to_keep$cell] <- "Y"

# Save QC/metadata and gene conversion results ====
dir.create("../outputs")
write.csv(metadata, "../outputs/metadata_table.csv", row.names = F)
write.csv(gene_info, "../outputs/gene_conversion_table.csv", row.names = F)

# End of document 
################################################
