args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
root <- normalizePath(args[1], mustWork = TRUE)
mode <- args[2]
suppressPackageStartupMessages(library(Matrix))
dir <- file.path(root, "raw_feature_bc_matrix")
read_tab <- function(name, header = FALSE) {
  con <- gzfile(file.path(dir, paste0(name, ".gz")), "rt")
  on.exit(close(con))
  read.delim(con, header = header, stringsAsFactors = FALSE, check.names = FALSE,
             quote = "", comment.char = "")
}
con <- gzfile(file.path(dir, "matrix.mtx.gz"), "rt")
counts <- as(readMM(con), "CsparseMatrix")
close(con)
features <- read_tab("features.tsv")
stopifnot(ncol(features) == 3L)
names(features) <- c("gene_id", "gene_name", "feature_type")
barcodes <- read_tab("barcodes.tsv")[[1]]
cells <- read.delim(file.path(dir, "cell_metadata.tsv"), check.names = FALSE,
                    stringsAsFactors = FALSE, quote = "", comment.char = "")
stopifnot(nrow(counts) == nrow(features), ncol(counts) == length(barcodes),
          !anyDuplicated(features$gene_id), !anyDuplicated(barcodes),
          identical(barcodes, cells$cell_id),
          all(is.finite(counts@x)), all(counts@x >= 0),
          all(counts@x == round(counts@x)))
if (sum(counts) == 0) stop("No gene-assigned UMIs; inspect alignment/annotation/strand before proceeding")
rownames(counts) <- features$gene_id
colnames(counts) <- barcodes
rownames(cells) <- cells$cell_id
stopifnot(all(Matrix::colSums(counts) == cells$nCount_RNA),
          all(Matrix::colSums(counts > 0) == cells$nFeature_RNA))
mito <- grepl("^MT-", features$gene_name)
cells$percent_mito <- ifelse(cells$nCount_RNA > 0,
  100 * Matrix::colSums(counts[mito, , drop = FALSE]) / cells$nCount_RNA, NA_real_)
provenance <- list(sample = "GW18", run = "SRR17878454", mode = mode,
  counts_type = "Unnormalized integer UMI requantification using retained corrected CB/UB",
  features = "Gene (exonic)", normalized = FALSE, cell_filtered = FALSE,
  method = readLines(file.path(root, "method.txt")),
  input = readLines(file.path(root, "input.txt")),
  reference_SHA256SUMS = readLines(file.path(root, "reference_SHA256SUMS")),
  code_SHA256SUMS = readLines(file.path(root, "code_SHA256SUMS")))
output <- list(counts = counts, features = features, cell_metadata = cells, provenance = provenance)
name <- if (mode == "pilot") "GW18_pilot_NOT_FULL_DATA.rds" else "GW18_for_qc.rds"
saveRDS(output, file.path(root, name), compress = "gzip")
capture.output({
  cat("Mode:", mode, "\nRDS:", name, "\nGenes:", nrow(counts), "\nBarcodes:", ncol(counts),
      "\nTotal UMIs:", sum(counts), "\nNo normalization or cell filtering was applied.\n")
  cat("\nBarcode QC by retained library (diagnostic, not a pass/fail cutoff):\n")
  for (lib in unique(cells$library)) {
    cat("\n", lib, "\n", sep = "")
    print(summary(cells[cells$library == lib, c("nCount_RNA", "nFeature_RNA", "percent_mito")]))
  }
  cat("\nRDS is a list with counts, features, cell_metadata and provenance; it is not a Seurat object.\n")
}, file = file.path(root, "QC_summary.txt"))
writeLines(capture.output(sessionInfo()), file.path(root, "R_sessionInfo.txt"))
cat(readLines(file.path(root, "QC_summary.txt")), sep = "\n")
