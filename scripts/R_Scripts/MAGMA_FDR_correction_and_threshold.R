args <- commandArgs(TRUE)
gene_output <- args[1]  # include path with file stem
path_output <- args[2]

# Define directories for gene and pathway results
gene_results_dir <- "/data/h_vmac/contra2/01_APOE4_Strat/05_MAGMA/gene_results/"
pathway_results_dir <- "/data/h_vmac/contra2/01_APOE4_Strat/05_MAGMA/pathway_results/"

# Load necessary library
library(data.table)

# Read in results
genes <- fread(gene_output, header=TRUE, stringsAsFactors=FALSE)
paths <- fread(path_output, header=TRUE, stringsAsFactors=FALSE, skip=4) # Skip first 4 lines for metadata

# Do FDR correction
print("FDR Correcting...")
genes$P.fdr <- p.adjust(genes$P, method="fdr")
paths$P.fdr <- p.adjust(paths$P, method="fdr")

# Construct output file paths
gene_fdr_output <- paste0(gene_results_dir, basename(gene_output), "_with_FDR")
path_fdr_output <- paste0(pathway_results_dir, basename(path_output), "_with_FDR")

# Write FDR-corrected results
write.table(genes, file=gene_fdr_output, row.names=FALSE, quote=FALSE)
write.table(paths, file=path_fdr_output, row.names=FALSE, quote=FALSE)

# Subset by significance
genes_sig <- subset(genes, P.fdr < 0.05)
paths_sig <- subset(paths, P.fdr < 0.05)

# Construct output file paths for significant results
gene_sig_output <- paste0(gene_results_dir, basename(gene_output), "_sig")
path_sig_output <- paste0(pathway_results_dir, basename(path_output), "_sig")

# Write significant results
write.table(genes_sig, file=gene_sig_output, row.names=FALSE, quote=FALSE)
write.table(paths_sig, file=path_sig_output, row.names=FALSE, quote=FALSE)
print("Done")