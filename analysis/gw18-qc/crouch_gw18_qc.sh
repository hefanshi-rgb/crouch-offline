#!/bin/bash
# GW18 cell calling and QC; launch on login, compute only inside Slurm.
set -Eeuo pipefail
umask 077
BASE=/public4/usersdir/guojunping/guojunping_2026-11-08_001/upload/Crouch
SOURCE="$BASE/recount_slurm/GW18/full_1713860"
SELF=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)/$(basename "${BASH_SOURCE[0]}")
if [[ "${1:-}" != --worker ]]; then
    [[ -s "$SOURCE/GW18_for_qc.rds" && -f "$SOURCE/COMPLETE" ]] || { echo 'Full result missing'; exit 1; }
    if ! type module >/dev/null 2>&1; then
        for init in /etc/profile.d/modules.sh /usr/share/Modules/init/bash /usr/share/lmod/lmod/init/bash; do
            if [[ -r "$init" ]]; then set +u; source "$init"; set -u; break; fi
        done
    fi
    set +u
    module use /public/software/modules
    module load compiler/gcc/7.4.0 system/singularity/3.6.3 system/R/4.3.1
    set -u
    lib=$(gfortran -print-file-name=libgfortran.so.4)
    [[ "$lib" == /* && -r "$lib" ]] || { echo 'GCC runtime missing'; exit 1; }
    export LD_LIBRARY_PATH="${lib%/*}${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
    export GW18_QC_R=$(command -v Rscript)
    export GW18_QC_SINGULARITY=$(command -v singularity)
    mkdir -p "$BASE/logs"
    reply=$(sbatch --parsable --export=ALL -p bigmem --nodelist=comput159 \
        -N 1 -n 1 -c 4 --mem=64G -t 06:00:00 --job-name=gw18-cell-qc \
        --chdir="$BASE" --output="$BASE/logs/gw18-cell-qc-%j.log" "$SELF" --worker)
    job=${reply%%;*}
    [[ "$job" =~ ^[0-9]+$ ]] || { echo "$reply"; exit 1; }
    printf 'Job=%s\nLog=%s/logs/gw18-cell-qc-%s.log\nResult=%s/recount_slurm/GW18/qc_%s\n' "$job" "$BASE" "$job" "$BASE" "$job"
    exit 0
fi
: "${SLURM_JOB_ID:?This worker must run through Slurm}"
OUT="$BASE/recount_slurm/GW18/qc_${SLURM_JOB_ID}"
mkdir "$OUT"
cp "$SELF" "$OUT/submitted_script.sh"
export TZ=UTC OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export R_LIBS="$BASE/software/r-libraries/R-4.3.1-Matrix-1.6.0${R_LIBS:+:$R_LIBS}"
export R_LIBS_USER="$R_LIBS"
export GW18_QC_BASE="$BASE" GW18_QC_OUT="$OUT"
cat > "$OUT/run_STAR.sh" <<'STAR_WRAPPER'
#!/bin/bash
set -Eeuo pipefail
exec "$GW18_QC_SINGULARITY" exec --cleanenv \
    --home "$GW18_QC_BASE/software/container_home:/home/crouch" \
    --bind "$GW18_QC_BASE:$GW18_QC_BASE" --pwd "$GW18_QC_OUT" \
    "$GW18_QC_BASE/software/ubuntu-22.04.sif" \
    "$GW18_QC_BASE/software/crouch-recount/STAR-2.7.10b" "$@"
STAR_WRAPPER
sed -n '/^# BEGIN_R$/,/^# END_R$/p' "$SELF" > "$OUT/qc.R"
"$GW18_QC_R" --vanilla "$OUT/qc.R" "$SOURCE" "$OUT" "$OUT/run_STAR.sh"
touch "$OUT/COMPLETE"
printf '\nSUCCESS: %s\n' "$OUT"
cat "$OUT/QC_summary.txt"
exit 0
: <<'R_EMBEDDED'
# BEGIN_R
args <- commandArgs(TRUE)
stopifnot(length(args) == 3L)
suppressPackageStartupMessages(library(Matrix))
src <- normalizePath(args[1]); out <- normalizePath(args[2]); runner <- normalizePath(args[3])
version <- system2('/bin/bash', shQuote(c(runner, '--version')), stdout=TRUE)
if (!identical(trimws(version), '2.7.10b')) stop('Expected STAR 2.7.10b')
# Exploratory thresholds: flags and a separate candidate object, never final QC.
min.genes <- 200L; min.umi <- 500L; max.mito <- 20
# Official STARsolo EmptyDrops_CR defaults; 3000 is a prior, not a forced cell count.
ed <- c(3000, 0.99, 10, 45000, 90000, 500, 0.01, 20000, 0.01, 10000)
input <- file.path(src, 'GW18_for_qc.rds')
x <- readRDS(input); m <- x$counts; meta <- x$cell_metadata; feat <- x$features
stopifnot(inherits(m, 'dgCMatrix'), ncol(m) == nrow(meta), nrow(m) == nrow(feat),
          !anyDuplicated(colnames(m)), !anyDuplicated(rownames(m)),
          identical(colnames(m), as.character(meta$cell_id)),
          identical(rownames(m), as.character(feat$gene_id)),
          all(is.finite(m@x)), all(m@x >= 0), all(m@x == round(m@x)),
          identical(x$provenance$mode, 'full'), !isTRUE(x$provenance$normalized))
stopifnot(all(grepl('^[ACGT]{16}-[0-9]+$', meta$CB)))
meta$gem <- sub('.*-', '', meta$CB)
meta$nCount_RNA <- as.numeric(Matrix::colSums(m))
meta$nFeature_RNA <- diff(m@p)
stopifnot(all(m@x > 0))
mito <- grepl('^MT-', feat$gene_name)
if (!any(mito)) stop('No MT- genes; mitochondrial QC cannot be computed')
meta$percent_mito <- ifelse(meta$nCount_RNA > 0,
    100 * Matrix::colSums(m[mito, , drop=FALSE]) / meta$nCount_RNA, NA_real_)
meta$cell_called <- FALSE
meta$flag_low_genes <- meta$nFeature_RNA < min.genes
meta$flag_low_UMI <- meta$nCount_RNA < min.umi
meta$flag_high_mito <- !is.finite(meta$percent_mito) | meta$percent_mito > max.mito
meta$qc_candidate <- FALSE
key <- paste(meta$library, meta$gem, sep='|GEM=')
write.tsv <- function(d, p) write.table(d, p, sep='\t', quote=FALSE, row.names=FALSE, na='NA')
med <- function(v) if (length(v)) median(v) else NA_real_
summaries <- list(); sensitivity <- list()
pdf(file.path(out, 'GW18_QC_plots.pdf'), width=11, height=8.5, useDingbats=FALSE)
for (k in sort(unique(key))) {
    ii <- which(key == k); id <- sprintf('library_%02d', length(summaries)+1L)
    folder <- file.path(out, id); raw <- file.path(folder, 'raw'); filtered <- file.path(folder, 'called')
    dir.create(raw, recursive=TRUE); dir.create(filtered)
    cat('\nCELL CALLING:', k, 'barcodes=', length(ii), '\n')
    a <- m[, ii, drop=FALSE]
    stopifnot(!anyDuplicated(meta$CB[ii]), sum(a) > 0)
    # STAR 2.7.10b standalone import omits its last nonzero-column group.
    # An explicit ZERO entry in a sentinel column protects the last real barcode.
    # This adds no molecule; sentinel is never in metadata or exported counts.
    sentinel <- '__ZERO_COUNT_BOUNDARY_SENTINEL__'
    padded <- cbind(a, sparseMatrix(i=1L, j=1L, x=0, dims=c(nrow(a), 1L)))
    stopifnot(length(padded@x) == length(a@x)+1L, tail(padded@x,1L) == 0,
              sum(padded) == sum(a))
    invisible(writeMM(padded, file.path(raw, 'matrix.mtx')))
    write.table(feat[, c('gene_id','gene_name','feature_type')], file.path(raw, 'features.tsv'),
                sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)
    writeLines(c(meta$CB[ii], sentinel), file.path(raw, 'barcodes.tsv'))
    rm(padded); gc(FALSE)
    argv <- c(runner, '--runMode', 'soloCellFiltering', paste0(raw, '/'),
              paste0(filtered, '/'), '--outFileNamePrefix', paste0(folder, '/'),
              '--runThreadN', '4', '--soloCellFilter', 'EmptyDrops_CR', as.character(ed))
    log <- file.path(folder, 'STAR.stdout.txt')
    status <- system2('/bin/bash', shQuote(argv), stdout=log, stderr=log)
    if (status != 0L) stop('STAR cell calling failed; inspect ', log)
    bcfile <- file.path(filtered, 'barcodes.tsv')
    if (!file.exists(bcfile)) stop('No called barcode file; inspect ', folder)
    cb <- readLines(bcfile)
    stopifnot(length(cb) > 0L, !anyDuplicated(cb), !sentinel %in% cb, all(cb %in% meta$CB[ii]))
    jj <- match(cb, meta$CB[ii]); selected <- ii[jj]
    cm <- as(readMM(file.path(filtered, 'matrix.mtx')), 'CsparseMatrix')
    stopifnot(identical(dim(cm), c(nrow(m), length(selected))),
              Matrix::nnzero(cm - m[, selected, drop=FALSE]) == 0)
    meta$cell_called[selected] <- TRUE
    meta$qc_candidate[selected] <- !(meta$flag_low_genes[selected] |
        meta$flag_low_UMI[selected] | meta$flag_high_mito[selected])
    good <- selected[meta$qc_candidate[selected]]
    lines <- readLines(file.path(folder, 'Log.out'))
    ambient <- any(grepl('finished ambient profile', lines, fixed=TRUE))
    summaries[[k]] <- data.frame(group=k, directory=id, raw_barcodes=length(ii),
        positive_barcodes=sum(meta$nCount_RNA[ii]>0), called_cells=length(selected),
        qc_candidates=length(good), raw_UMIs=sum(meta$nCount_RNA[ii]),
        called_UMIs=sum(meta$nCount_RNA[selected]),
        called_UMI_fraction=sum(meta$nCount_RNA[selected])/sum(meta$nCount_RNA[ii]),
        called_median_UMI=med(meta$nCount_RNA[selected]),
        called_median_genes=med(meta$nFeature_RNA[selected]),
        called_median_mito=med(meta$percent_mito[selected]),
        candidate_median_UMI=med(meta$nCount_RNA[good]),
        candidate_median_genes=med(meta$nFeature_RNA[good]),
        ambient_profile_used=ambient)
    sensitivity[[k]] <- do.call(rbind, lapply(c(100,200,300,500), function(g)
        do.call(rbind, lapply(c(10,20,30), function(mt) data.frame(group=k,
            min_genes=g, min_UMIs=min.umi, max_mito=mt,
            cells=sum(meta$nFeature_RNA[selected]>=g & meta$nCount_RNA[selected]>=min.umi &
                      is.finite(meta$percent_mito[selected]) & meta$percent_mito[selected]<=mt))))))
    par(mfrow=c(2,2), mar=c(4,4,3,1), oma=c(0,0,2,0))
    ord <- ii[order(meta$nCount_RNA[ii], decreasing=TRUE)]
    ord <- ord[meta$nCount_RNA[ord]>0]
    thin <- unique(pmax(1L, round(exp(seq(0,log(length(ord)),length.out=min(6000L,length(ord)))))))
    plot(thin, meta$nCount_RNA[ord[thin]], log='xy', type='l', col='grey50',
         xlab='Barcode rank', ylab='UMIs', main='Barcode rank (blue = called)')
    cr <- which(meta$cell_called[ord]); points(cr, meta$nCount_RNA[ord[cr]], pch=16, cex=.25, col='#2166AC')
    bg <- ord[thin]
    plot(log10(1+meta$nCount_RNA[bg]), log10(1+meta$nFeature_RNA[bg]),
         pch=16, cex=.35, col='grey75', xlab='log10(1 + UMIs)', ylab='log10(1 + genes)',
         main='Called blue; QC candidates green',
         xlim=range(log10(1+meta$nCount_RNA[ord])), ylim=range(log10(1+meta$nFeature_RNA[ord])))
    points(log10(1+meta$nCount_RNA[selected]), log10(1+meta$nFeature_RNA[selected]), pch=16, cex=.25, col='#2166AC')
    if(length(good)) points(log10(1+meta$nCount_RNA[good]), log10(1+meta$nFeature_RNA[good]), pch=16, cex=.25, col='#1B9E77')
    abline(v=log10(1+min.umi), h=log10(1+min.genes), lty=2)
    hist(meta$nFeature_RNA[selected], breaks=50, main='Genes in called cells', xlab='Detected genes', col='grey80')
    abline(v=min.genes, col='red', lty=2)
    hist(meta$percent_mito[selected], breaks=seq(0,100,by=2), main='Mitochondrial fraction in called cells', xlab='Percent mitochondrial', col='grey80')
    abline(v=max.mito, col='red', lty=2)
    mtext(k, outer=TRUE, cex=1.2)
    rm(a, cm); gc(FALSE)
}
dev.off()
meta$doublet_status <- 'not_assessed'
rownames(meta) <- meta$cell_id
con <- gzfile(file.path(out, 'all_barcode_QC.tsv.gz'), 'wt')
write.tsv(meta, con); close(con)
summary <- do.call(rbind, summaries); rownames(summary) <- NULL
write.tsv(summary, file.path(out, 'QC_summary.tsv'))
capture.output({
    cat('Exploratory QC: genes >= 200; UMIs >= 500; mitochondrial percent <= 20\n')
    for (r in seq_len(nrow(summary))) {
        z <- summary[r, ]
        cat('\n', z$group, '\n', sep='')
        cat('Raw barcodes:', z$raw_barcodes, '\nCalled cells:', z$called_cells,
            '\nQC candidates:', z$qc_candidates,
            '\nUMI fraction in called cells:', round(z$called_UMI_fraction,4),
            '\nCalled median UMIs / genes / mito%:', z$called_median_UMI,
            z$called_median_genes, round(z$called_median_mito,2),
            '\nAmbient profile used:', z$ambient_profile_used, '\n')
    }
}, file=file.path(out,'QC_summary.txt'))
write.tsv(do.call(rbind, sensitivity), file.path(out, 'threshold_sensitivity.tsv'))
write.tsv(meta[meta$cell_called, , drop=FALSE], file.path(out, 'called_cell_QC.tsv'))
prov <- x$provenance
prov$cell_filtered <- TRUE; prov$normalized <- FALSE
prov$cell_calling <- list(tool='STARsolo 2.7.10b EmptyDrops_CR', parameters=ed,
    grouping='physical library and original GEM suffix',
    input_scope='Only SRA records with recoverable corrected CB/UB; original ambient pool completeness unknown',
    boundary_guard='Explicit zero-count sentinel excluded from output; counts verified against original',
    source_RDS=input, source_MD5=unname(tools::md5sum(input)))
prov$qc <- list(min_genes=min.genes, min_UMIs=min.umi, max_percent_mito=max.mito,
    thresholds_status='Exploratory, requires review', doublets_assessed=FALSE)
for (stage in c('cellcalled', 'qc_candidate')) {
    keep <- if(stage=='cellcalled') meta$cell_called else meta$qc_candidate
    obj <- list(counts=m[, keep, drop=FALSE], features=feat,
                cell_metadata=meta[keep, , drop=FALSE], provenance=prov)
    obj$provenance$qc_applied <- stage=='qc_candidate'
    obj$provenance$mode <- stage
    saveRDS(obj, file.path(out, paste0('GW18_', stage, '_for_qc.rds')), compress='gzip')
}
writeLines(c('Cell calling and exploratory QC, not final integration-ready data.',
    'Counts and all reference gene rows preserved exactly; no normalization or imputation.',
    'Calls depend on retained SRA barcode background; original droplet completeness unknown.',
    'Expected cells = 3000 per library is an initial STAR default, not a forced number.',
    'ambient_profile_used=FALSE indicates knee-only calls without the additional ambient test.',
    'QC candidate: called AND genes >= 200 AND UMIs >= 500 AND mitochondrial percent <= 20.',
    'No upper-count cutoff, doublet removal, ambient-RNA subtraction or author comparison applied.',
    'Review QC_summary.tsv, threshold_sensitivity.tsv and GW18_QC_plots.pdf before finalizing.'),
    file.path(out,'README.txt'))
writeLines(capture.output(sessionInfo()), file.path(out, 'R_sessionInfo.txt'))
if (!any(meta$qc_candidate)) stop('No QC candidates: inspect reports; do not use for integration')
# END_R
R_EMBEDDED
