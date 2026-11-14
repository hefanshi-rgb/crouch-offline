#!/bin/bash
set -Eeuo pipefail
umask 077
BASE=/public4/usersdir/guojunping/guojunping_2026-11-08_001/upload/Crouch
SELF=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)/$(basename "${BASH_SOURCE[0]}")
if [[ "${1:-}" != --worker ]]; then
    if ! type module >/dev/null 2>&1; then
        for init in /etc/profile.d/modules.sh /usr/share/Modules/init/bash /usr/share/lmod/lmod/init/bash; do
            if [[ -r "$init" ]]; then set +u; source "$init"; set -u; break; fi
        done
    fi
    set +u
    module use /public/software/modules
    module load compiler/gcc/7.4.0 system/R/4.3.1
    set -u
    lib=$(gfortran -print-file-name=libgfortran.so.4)
    [[ "$lib" == /* && -r "$lib" ]] || { echo 'GCC runtime missing'; exit 1; }
    export LD_LIBRARY_PATH="${lib%/*}${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
    EXPR_R=$(command -v Rscript);export EXPR_R
    mkdir -p "$BASE/logs"
    reply=$(sbatch --parsable --export=ALL -p bigmem --nodelist=comput159 \
      -N 1 -n 1 -c 2 --mem=32G -t 02:00:00 --job-name=gw18-expression-check \
      --chdir="$BASE" --output="$BASE/logs/gw18-expression-check-%j.log" "$SELF" --worker)
    job=${reply%%;*}; [[ "$job" =~ ^[0-9]+$ ]] || { echo "$reply";exit 1; }
    printf 'Job=%s\nLog=%s/logs/gw18-expression-check-%s.log\nResult=%s/recount_slurm/GW18/expression_%s\n' "$job" "$BASE" "$job" "$BASE" "$job"
    exit 0
fi
: "${SLURM_JOB_ID:?Submit using the launcher}"
OUT="$BASE/recount_slurm/GW18/expression_${SLURM_JOB_ID}"
mkdir "$OUT";cp "$SELF" "$OUT/submitted_script.sh"
export TZ=UTC OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export R_LIBS="$BASE/software/r-libraries/R-4.3.1-Matrix-1.6.0${R_LIBS:+:$R_LIBS}"
export R_LIBS_USER="$R_LIBS"
sed -n '/^# BEGIN_R$/,/^# END_R$/p' "$SELF" > "$OUT/expression_check.R"
"$EXPR_R" --vanilla "$OUT/expression_check.R" \
 "$BASE/recount_slurm/GW18/meta_compare_1717231/GW18_author_cells_raw_counts.rds" \
 "$BASE/author_expression_check" "$OUT"
touch "$OUT/COMPLETE"
printf '\nSUCCESS: %s\n' "$OUT"
exit 0
: <<'R_EMBEDDED'
# BEGIN_R
suppressPackageStartupMessages(library(Matrix))
a<-commandArgs(TRUE); stopifnot(length(a)>=3L)
raw.path<-a[1]; src<-a[2]; out<-a[3]; min.genes<-if(length(a)>=4L) as.integer(a[4]) else 1000L
own<-readRDS(raw.path)
stopifnot(own$provenance$mode=='author_cohort',!isTRUE(own$provenance$normalized),
    all(own$counts@x>=0),all(own$counts@x==round(own$counts@x)))
tab<-function(name) read.delim(gzfile(file.path(src,name),'rt'),header=FALSE,
    quote='',comment.char='',stringsAsFactors=FALSE,check.names=FALSE)
af<-tab('features.tsv.gz'); ab<-tab('barcodes.tsv.gz')[[1]]
meta<-read.delim(file.path(src,'meta.tsv'),quote='',comment.char='',stringsAsFactors=FALSE,check.names=FALSE)
stopifnot(!anyDuplicated(ab),!anyDuplicated(meta$cellId),all(toupper(trimws(meta$age))=='GW18'))
meta<-meta[match(ab,meta$cellId),,drop=FALSE];stopifnot(identical(ab,meta$cellId))
stopifnot(!anyDuplicated(own$cell_metadata$author_cellId),
    identical(colnames(own$counts),own$cell_metadata$cell_id))
ci<-match(ab,own$cell_metadata$author_cellId)
if(anyNA(ci)) stop('Some sampled author cells do not match the final author-cohort RDS')
pairs<-function(x,y) {
    x<-as.character(x);y<-as.character(y)
    ux<-!duplicated(x)&!duplicated(x,fromLast=TRUE)&!is.na(x)&nzchar(x)
    uy<-!duplicated(y)&!duplicated(y,fromLast=TRUE)&!is.na(y)&nzchar(y)
    j<-match(x,y);ok<-ux&!is.na(j);ok[ok]<-uy[j[ok]]
    list(author=which(ok),own=j[ok])
}
p<-pairs(af[[1]],own$features$gene_id);method<-'Exact unique gene IDs'
if(length(p$author)<min.genes) {
    sy<-if(ncol(af)>=2L) af[[2]] else af[[1]]
    q<-pairs(sy,own$features$gene_name)
    if(length(q$author)>length(p$author)) {p<-q;method<-'Exact unique gene symbols'}
}
if(length(p$author)<min.genes) stop('Too few common genes; inspect feature identifiers before comparing')
con<-gzfile(file.path(src,'matrix.mtx.gz'),'rt'); am<-as(readMM(con),'CsparseMatrix');close(con)
stopifnot(nrow(am)==nrow(af),ncol(am)==length(ab),all(is.finite(am@x)))
write.tsv<-function(d,name) write.table(d,file.path(out,name),sep='\t',row.names=FALSE,quote=FALSE)
write.tsv(data.frame(author_row=p$author,recount_row=p$own,
    gene_id=own$features$gene_id[p$own],gene_name=own$features$gene_name[p$own]),'gene_matches.tsv')
raw<-own$counts[p$own,ci,drop=FALSE]
# The denominator includes ALL counted genes, before common-gene subsetting.
totals<-as.numeric(Matrix::colSums(own$counts[,ci,drop=FALSE]))
stopifnot(all(totals>0))
cor.safe<-function(x,y,method) if(length(x)>=3L&&sd(x)>0&&sd(y)>0) cor(x,y,method=method) else NA_real_
rows<-lapply(seq_along(ab),function(j){
    rv<-as.numeric(raw[,j]); nv<-log1p(rv/totals[j]*10000); av<-as.numeric(am[p$author,j])
    use<-rv!=0|av!=0; both<-rv>0&av!=0
    data.frame(author_cell=ab[j],recount_cell=own$cell_metadata$cell_id[ci[j]],
      library=own$cell_metadata$library[ci[j]],
      author_cell_type=if('cell_type'%in%names(meta)) as.character(meta$cell_type[j]) else NA_character_,
      common_genes=length(rv),genes_in_union=sum(use),genes_detected_in_both=sum(both),
      detection_jaccard=if(sum(use)>0) sum(both)/sum(use) else NA_real_,
      raw_UMIs=totals[j],UMI_fraction_in_common_genes=sum(rv)/totals[j],
      spearman_union=cor.safe(nv[use],av[use],'spearman'),
      pearson_logNormalize_vs_author=cor.safe(nv[use],av[use],'pearson'),
      spearman_both_detected=cor.safe(nv[both],av[both],'spearman'))
})
res<-do.call(rbind,rows);write.tsv(res,'cell_expression_correlations.tsv')
summ<-do.call(rbind,lapply(split(res,res$library),function(z){
    data.frame(library=z$library[1],cells=nrow(z),finite_spearman=sum(is.finite(z$spearman_union)),
      spearman_q25=unname(quantile(z$spearman_union,.25,na.rm=TRUE)),
      spearman_median=median(z$spearman_union,na.rm=TRUE),
      spearman_q75=unname(quantile(z$spearman_union,.75,na.rm=TRUE)),
      pearson_median=median(z$pearson_logNormalize_vs_author,na.rm=TRUE),
      detection_jaccard_median=median(z$detection_jaccard,na.rm=TRUE))
}));rownames(summ)<-NULL;write.tsv(summ,'library_expression_summary.tsv')
report<-c('Expression concordance in a stratified sample, not validation of every GW18 cell.',
 'Our counts: log1p(10000 * count / total UMI across ALL original genes).',
 'Author values: retained exactly as supplied; their transformation has not been established.',
 'Shared-zero genes excluded from both correlation metrics.',
 'Spearman evaluates ranking; monotonic per-cell normalization does not change that ranking.',
 'Pearson is exploratory until the author normalization method is confirmed.',
 'High correlation is not proof of exact recovery of author raw counts.',
 paste('Gene matching:',method),paste('Common genes:',length(p$author)),
 paste('Author features:',nrow(af)),paste('Recount features:',nrow(own$features)),
 paste('Matched sampled cells:',nrow(res)),
 paste('Author supplied values: min=',min(am@x),'max=',max(am@x),'negative entries=',sum(am@x<0)),
 capture.output(print(summ,row.names=FALSE,digits=4)))
writeLines(report,file.path(out,'expression_report.txt'))
if(!any(is.finite(res$spearman_union))) stop('No finite correlations; inspect expression_report.txt')
pdf(file.path(out,'expression_concordance.pdf'),width=11,height=8.5,useDingbats=FALSE)
par(mfrow=c(2,2),mar=c(6,4,3,1))
short<-sub(':.*','',res$library)
boxplot(res$spearman_union~short,ylim=c(-1,1),ylab='Spearman (expressed-gene union)',xlab='Library')
hist(res$spearman_union,breaks=30,xlim=c(-1,1),main='Per-cell Spearman',xlab='Spearman')
hist(res$pearson_logNormalize_vs_author,breaks=30,xlim=c(-1,1),main='Pearson: exploratory',xlab='Pearson')
plot(res$detection_jaccard,res$spearman_union,pch=16,cex=.6,xlim=c(0,1),ylim=c(-1,1),
    xlab='Detected-gene Jaccard',ylab='Spearman',main='Detection overlap and expression ranking')
par(mfrow=c(2,2),mar=c(4,4,3,1))
for(lib in sort(unique(res$library))) {
    ii<-which(res$library==lib&is.finite(res$spearman_union))
    if(!length(ii)) next
    j<-ii[which.min(abs(res$spearman_union[ii]-median(res$spearman_union[ii])))]
    rv<-as.numeric(raw[,j]);nv<-log1p(rv/totals[j]*10000);av<-as.numeric(am[p$author,j]);use<-which(rv!=0|av!=0)
    if(length(use)>10000L) use<-use[unique(round(seq(1,length(use),length.out=10000)))]
    plot(nv[use],av[use],pch=16,cex=.3,col=rgb(0,0,0,.25),
      xlab='Our log1p CP10k',ylab='Author supplied expression',
      main=paste(sub(':.*','',lib),'median-rho cell; rho=',round(res$spearman_union[j],3)))
}
dev.off()
writeLines(capture.output(sessionInfo()),file.path(out,'R_sessionInfo.txt'))
cat(report,sep='\n')
# END_R
R_EMBEDDED
