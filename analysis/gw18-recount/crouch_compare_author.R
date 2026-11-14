# Run on a compute node: Rscript --vanilla crouch_compare_author.R COUNTS_RDS AUTHOR_GW18_DIR OUTPUT_DIR
# Mapping file must be reviewed first; no sample/library identities are guessed.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3L)
suppressPackageStartupMessages(library(Matrix))
recount <- readRDS(args[1])
if (!identical(recount$provenance$mode, "full")) stop("Compare full-depth counts, not the pilot")
author <- normalizePath(args[2], mustWork = TRUE)
out <- args[3]
if (dir.exists(out)) stop("Use a new comparison output directory")
dir.create(out, recursive = TRUE)
tab <- function(file, header = FALSE) {
  con <- gzfile(file.path(author, file), "rt"); on.exit(close(con))
  read.delim(con, header = header, quote = "", comment.char = "", check.names = FALSE,
             stringsAsFactors = FALSE)
}
af <- tab("features.tsv.gz")
ab <- tab("barcodes.tsv.gz")[[1]]
meta <- read.delim(file.path(author, "meta.tsv"), quote = "", comment.char = "",
                   check.names = FALSE, stringsAsFactors = FALSE)
mapping <- read.delim(file.path(author, "library_mapping_to_review.tsv"), quote = "",
                      comment.char = "", check.names = FALSE, stringsAsFactors = FALSE,
                      colClasses = "character", na.strings = "")
if (!all(c("author_prefix", "recount_library", "recount_gem") %in% names(mapping)) ||
    anyNA(mapping[, c("author_prefix", "recount_library", "recount_gem")]) ||
    anyDuplicated(mapping$author_prefix)) stop("Fill and review all library/GEM mappings first")
parts <- regmatches(ab, regexec("^(.*?)([ACGT]{16})(?:-([1-9][0-9]*))?$", ab, perl = TRUE))
if (any(lengths(parts) < 3L)) stop("Unexpected author cell naming; explicit parsing needed")
prefix <- vapply(parts, `[[`, "", 2L)
cb <- vapply(parts, `[[`, "", 3L)
mi <- match(prefix, mapping$author_prefix)
if (anyNA(mi)) stop("Some author prefixes have no reviewed mapping")
lib <- mapping$recount_library[mi]
gem <- mapping$recount_gem[mi]
if (any(!grepl("^[1-9][0-9]*$", gem))) stop("Invalid GEM mapping")
author_gem <- vapply(parts, function(x) if (length(x) >= 4L) x[4] else "", "")
if (any(nzchar(author_gem) & author_gem != gem)) stop("Author GEM suffix disagrees with mapping")
if (anyDuplicated(meta$cellId) || anyNA(match(ab, meta$cellId))) stop("Author metadata/cell IDs disagree")
meta <- meta[match(ab, meta$cellId), , drop = FALSE]
if (any(toupper(trimws(meta$age)) != "GW18")) stop("Author subset contains non-GW18 cells")
actual <- recount$cell_metadata
# Match explicit fields rather than guessing how encoded display identifiers look.
key <- paste(lib, paste0(cb, "-", gem), sep = "|")
if (anyDuplicated(key)) stop("Author mapping merges different cells; fix mapping")
ri <- match(key, paste(actual$library, actual$CB, sep = "|"))
match_table <- data.frame(author_cell = ab, author_prefix = prefix, library = lib,
  CB = paste0(cb, "-", gem), matched = !is.na(ri),
  recount_cell = ifelse(is.na(ri), NA, actual$cell_id[ri]), stringsAsFactors = FALSE)
write.table(match_table, file.path(out, "cell_matches.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
if (!any(!is.na(ri))) stop("No cells match; review mapping before expression comparison")
unique_pairs <- function(a, b) {
  a <- as.character(a); b <- as.character(b)
  # Exact, one-to-one identifiers only; no silent gene-symbol aggregation.
  aa <- !duplicated(a) & !duplicated(a, fromLast = TRUE) & !is.na(a) & nzchar(a)
  bb <- !duplicated(b) & !duplicated(b, fromLast = TRUE) & !is.na(b) & nzchar(b)
  j <- match(a, b)
  good <- aa & !is.na(j)
  good[good] <- bb[j[good]]
  list(a = which(good), b = j[good])
}
gp <- unique_pairs(af[[1]], recount$features$gene_id)
gene_method <- "Exact unique gene ID"
if (!length(gp$a)) {
  if (ncol(af) < 2L) stop("No matching gene IDs; provide an explicit feature mapping")
  gp <- unique_pairs(af[[2]], recount$features$gene_name)
  gene_method <- "Exact unique gene symbol (no gene IDs matched)"
}
if (!length(gp$a)) stop("No genes match")
write.table(data.frame(author_row = gp$a, recount_row = gp$b,
  gene_id = recount$features$gene_id[gp$b], gene_name = recount$features$gene_name[gp$b]),
  file.path(out, "gene_matches.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
con <- gzfile(file.path(author, "matrix.mtx.gz"), "rt")
am <- as(readMM(con), "CsparseMatrix"); close(con)
stopifnot(nrow(am) == nrow(af), ncol(am) == length(ab), all(is.finite(am@x)))
# Deterministic representatives from each library; maximum 200 total cells.
eligible <- which(!is.na(ri))
by_library <- split(eligible, lib[eligible])
budget <- max(1L, floor(200L / length(by_library)))
chosen <- unlist(lapply(by_library, function(ii) ii[unique(round(seq(1, length(ii), length.out = min(length(ii), budget))))]), use.names = FALSE)
rows <- lapply(chosen, function(j) {
  values <- as.numeric(recount$counts[gp$b, ri[j]])
  # Normalization denominator uses every counted gene, not only the common genes.
  total <- sum(recount$counts[, ri[j]])
  norm <- log1p(values / max(1, total) * 10000)
  av <- as.numeric(am[gp$a, j])
  use <- values != 0 | av != 0
  rho <- if (sum(use) >= 3L && sd(norm[use]) > 0 && sd(av[use]) > 0)
    cor(norm[use], av[use], method = "spearman") else NA_real_
  data.frame(author_cell = ab[j], recount_cell = actual$cell_id[ri[j]], library = lib[j],
    nCount_RNA = total, common_genes = length(values), genes_in_union = sum(use),
    spearman_detected_gene_union = rho)
})
correlations <- do.call(rbind, rows)
write.table(correlations, file.path(out, "cell_expression_correlations.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
if ("cell_type" %in% names(meta)) {
  markers <- c("PECAM1", "VWF", "CLDN5", "PDGFRB", "RGS5", "ACTA2", "TAGLN", "COL1A1", "DCN", "PTPRC", "AIF1", "RBFOX3")
  groups <- split(eligible, as.character(meta$cell_type[eligible]))
  marker_rows <- list()
  for (marker in markers) {
    k <- which(recount$features$gene_name[gp$b] == marker)
    if (length(k) != 1L) next
    for (type in names(groups)) {
      jj <- groups[[type]]
      own <- as.numeric(recount$counts[gp$b[k], ri[jj]])
      totals <- actual$nCount_RNA[ri[jj]]
      marker_rows[[length(marker_rows) + 1L]] <- data.frame(gene = marker, author_cell_type = type,
        matched_cells = length(jj), recount_detection_fraction = mean(own > 0),
        recount_mean_log1p_CP10k = mean(log1p(own / pmax(1, totals) * 10000)),
        author_mean_provided_expression = mean(as.numeric(am[gp$a[k], jj])))
    }
  }
  if (length(marker_rows)) write.table(do.call(rbind, marker_rows),
    file.path(out, "marker_by_author_cell_type.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
}
capture.output({
  cat("Author expression values were kept unchanged; transformation remains to be confirmed.\n")
  cat("Spearman correlation is exploratory concordance, not a proof of raw-count recovery.\n")
  cat("Gene matching:", gene_method, "\nCommon one-to-one genes:", length(gp$a), "\n")
  cat("Author GW18 cells:", length(ab), "\nMatched cells:", sum(!is.na(ri)),
      "\nMatch fraction:", mean(!is.na(ri)), "\n")
  print(aggregate(matched ~ library, match_table, mean))
  cat("\nPer-cell expression rank correlation, excluding shared-zero genes:\n")
  print(summary(correlations$spearman_detected_gene_union))
  cat("\nDo not impose an arbitrary correlation threshold. Inspect mapping, features, gene-count rules, cell-type markers, and the author's normalization method.\n")
}, file = file.path(out, "comparison_summary.txt"))
cat(readLines(file.path(out, "comparison_summary.txt")), sep = "\n")
