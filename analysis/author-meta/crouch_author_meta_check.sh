#!/bin/bash
# Metadata-only validation. Counts and all genes remain unchanged.
set -Eeuo pipefail
umask 077
BASE=/public4/usersdir/guojunping/guojunping_2026-11-08_001/upload/Crouch
SELF=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)/$(basename "${BASH_SOURCE[0]}")
if [[ "${1:-}" != --worker ]]; then
    mode=${1:-inspect}; mapping=${2:-}
    [[ "$mode" == inspect || "$mode" == compare ]] || { echo 'Usage: bash crouch_author_meta_check.sh inspect | compare REVIEWED_MAPPING.tsv'; exit 1; }
    if [[ "$mode" == compare ]]; then
        [[ -s "$mapping" ]] || { echo 'Provide a reviewed mapping file'; exit 1; }
        mapping=$(cd "$(dirname "$mapping")" && pwd -P)/$(basename "$mapping")
    fi
    [[ -s "$BASE/author_meta_all/meta.tsv" ]] || { echo 'Missing author_meta_all/meta.tsv'; exit 1; }
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
    AUTHOR_META_R=$(command -v Rscript); export AUTHOR_META_R
    mkdir -p "$BASE/logs"
    reply=$(sbatch --parsable --export=ALL -p bigmem --nodelist=comput159 \
        -N 1 -n 1 -c 2 --mem=32G -t 02:00:00 --job-name="gw18-meta-$mode" \
        --chdir="$BASE" --output="$BASE/logs/gw18-meta-$mode-%j.log" "$SELF" --worker "$mode" "$mapping")
    job=${reply%%;*}
    [[ "$job" =~ ^[0-9]+$ ]] || { echo "$reply"; exit 1; }
    printf 'Job=%s\nLog=%s/logs/gw18-meta-%s-%s.log\nResult=%s/recount_slurm/GW18/meta_%s_%s\n' "$job" "$BASE" "$mode" "$job" "$BASE" "$mode" "$job"
    exit 0
fi
: "${SLURM_JOB_ID:?Submit through the launcher}"
mode=${2:?}; mapping=${3:-}
OUT="$BASE/recount_slurm/GW18/meta_${mode}_${SLURM_JOB_ID}"
mkdir "$OUT"
cp "$SELF" "$OUT/submitted_script.sh"
export TZ=UTC OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export R_LIBS="$BASE/software/r-libraries/R-4.3.1-Matrix-1.6.0${R_LIBS:+:$R_LIBS}"
export R_LIBS_USER="$R_LIBS"
sed -n '/^# BEGIN_R$/,/^# END_R$/p' "$SELF" > "$OUT/meta_check.R"
"$AUTHOR_META_R" --vanilla "$OUT/meta_check.R" "$BASE" "$BASE/author_meta_all/meta.tsv" "$mode" "$mapping" "$OUT"
touch "$OUT/COMPLETE"
printf '\nSUCCESS: %s\n' "$OUT"
cat "$OUT/report.txt"
exit 0
: <<'R_EMBEDDED'
# BEGIN_R
suppressPackageStartupMessages(library(Matrix))
a <- commandArgs(TRUE); stopifnot(length(a)==5L)
base<-a[1]; path<-a[2]; mode<-a[3]; mapping.path<-a[4]; out<-a[5]
write.tsv<-function(x,p) write.table(x,p,sep='\t',quote=FALSE,row.names=FALSE,na='NA')
allmeta<-read.delim(path,check.names=FALSE,stringsAsFactors=FALSE,quote='',comment.char='')
stopifnot(all(c('cellId','age') %in% names(allmeta)),!anyDuplicated(allmeta$cellId),
          all(!is.na(allmeta$cellId)&nzchar(allmeta$cellId)))
allmeta$age<-trimws(as.character(allmeta$age))
invcols<-intersect(c('age','sample','area'),names(allmeta))
inventory<-as.data.frame(table(allmeta[,invcols,drop=FALSE]),stringsAsFactors=FALSE)
inventory<-inventory[inventory$Freq>0,,drop=FALSE]
write.tsv(inventory,file.path(out,'all_sample_inventory.tsv'))
meta<-allmeta[toupper(allmeta$age)=='GW18' & !is.na(allmeta$age),,drop=FALSE]
if(!nrow(meta)) stop('No age=GW18 rows; inspect age labels')
metrics<-intersect(c('ncount_rna','nCount_RNA','nfeature_rna','nFeature_RNA','percent_mito','percent_ribo'),names(allmeta))
audit<-list()
for(scope in c('all','GW18')) for(field in metrics) {
    dat<-if(scope=='all') allmeta else meta
    z<-suppressWarnings(as.numeric(as.character(dat[[field]]))); finite<-is.finite(z)
    q<-if(any(finite)) quantile(z[finite],c(0,.25,.5,.75,1),names=FALSE) else rep(NA_real_,5)
    audit[[length(audit)+1L]]<-data.frame(scope=scope,field=field,rows=length(z),
       missing_or_nonfinite=sum(!finite),noninteger=sum(abs(z[finite]-round(z[finite]))>1e-6),
       minimum=q[1],q25=q[2],median=q[3],q75=q[4],maximum=q[5])
}
write.tsv(do.call(rbind,audit),file.path(out,'author_QC_column_audit.tsv'))
parts<-regmatches(meta$cellId,regexec('^(.*?)([ACGT]{16})(?:-([1-9][0-9]*))?$',meta$cellId,perl=TRUE))
if(any(lengths(parts)<3L)) stop('Unexpected author cell IDs; explicit parsing required')
prefix<-vapply(parts,function(z) z[2],''); cb<-vapply(parts,function(z) z[3],'')
gem<-vapply(parts,function(z) if(length(z)>=4L) z[4] else '','')
full<-readRDS(file.path(base,'recount_slurm/GW18/full_1713860/GW18_for_qc.rds'))
actual<-full$cell_metadata; full<-NULL; gc(FALSE)
actual$gem<-sub('.*-','',actual$CB); actual$sequence<-sub('-[1-9][0-9]*$','',actual$CB)
stopifnot(!anyDuplicated(paste(actual$library,actual$CB,sep='|')))
qcdir<-file.path(base,'recount_slurm/GW18/qc_1717229')
stopifnot(file.exists(file.path(qcdir,'COMPLETE')))
called<-readRDS(file.path(qcdir,'GW18_cellcalled_for_qc.rds'))$cell_metadata$cell_id
candidate<-readRDS(file.path(qcdir,'GW18_qc_candidate_for_qc.rds'))
qc<-candidate$cell_metadata$cell_id
actual$called<-actual$cell_id %in% called; actual$qc<-actual$cell_id %in% qc
ag<-split(seq_len(nrow(actual)),paste(actual$library,actual$gem,sep='|GEM='))
rows<-list(); prefix.info<-list()
for(p in sort(unique(prefix))) {
    jj<-which(prefix==p)
    z<-data.frame(author_prefix=p,cells=length(jj))
    for(f in intersect(c('sample','area','cell_type'),names(meta)))
        z[[f]]<-paste(sort(unique(as.character(meta[jj,f]))),collapse='; ')
    prefix.info[[p]]<-z
    for(ix in ag) {
        thisgem<-actual$gem[ix[1]]; mi<-match(cb[jj],actual$sequence[ix])
        ok<-(!nzchar(gem[jj])|gem[jj]==thisgem)&!is.na(mi); hit<-ix[mi[ok]]
        rows[[length(rows)+1L]]<-data.frame(author_prefix=p,author_cells=length(jj),
            recount_library=actual$library[ix[1]],recount_gem=thisgem,
            raw_matches=sum(ok),raw_fraction=mean(ok),called_matches=sum(actual$called[hit]),
            qc_matches=sum(actual$qc[hit]),qc_fraction=sum(actual$qc[hit])/length(jj))
    }
}
overlap<-do.call(rbind,rows)
overlap<-overlap[order(overlap$author_prefix,-overlap$raw_fraction),]
write.tsv(overlap,file.path(out,'prefix_library_overlap.tsv'))
write.tsv(do.call(rbind,prefix.info),file.path(out,'author_prefix_inventory.tsv'))
write.tsv(data.frame(author_prefix=sort(unique(prefix)),recount_library='',recount_gem=''),
          file.path(out,'library_mapping_to_review.tsv'))
write.tsv(meta,file.path(out,'author_GW18_meta.tsv'))
report<-c('Metadata-only comparison: no gene-expression matrix was supplied.',
    paste('All author metadata cells:',nrow(allmeta)),paste('Author GW18 cells:',nrow(meta)),
    paste('Our called cells:',length(called)),paste('Our QC candidates:',length(qc)),
    'Mapping candidates below are not adopted automatically; inspect labels and barcode overlaps.',
    capture.output(print(overlap,row.names=FALSE,digits=4)),
    '\nAuthor prefix labels:',capture.output(print(do.call(rbind,prefix.info),row.names=FALSE)),
    '\nAuthor QC-column audit (GW18):',
    capture.output(print(do.call(rbind,audit)[do.call(rbind,audit)$scope=='GW18',],row.names=FALSE,digits=5)))
if(mode=='compare') {
    mapping<-read.delim(mapping.path,check.names=FALSE,stringsAsFactors=FALSE,
        colClasses='character',quote='',comment.char='',na.strings='')
    req<-c('author_prefix','recount_library','recount_gem')
    stopifnot(all(req %in% names(mapping)),!anyNA(mapping[,req]),!anyDuplicated(mapping$author_prefix))
    mm<-match(prefix,mapping$author_prefix)
    if(anyNA(mm)) stop('Some author prefixes have no reviewed mapping')
    lib<-mapping$recount_library[mm]; mg<-mapping$recount_gem[mm]
    stopifnot(all(grepl('^[1-9][0-9]*$',mg)),all(!nzchar(gem)|gem==mg))
    stopifnot(all(paste(lib,mg,sep='|GEM=') %in% names(ag)))
    keys<-paste(lib,paste0(cb,'-',mg),sep='|')
    if(anyDuplicated(keys)) stop('Mapping merges different author cells; review library identities')
    ri<-match(keys,paste(actual$library,actual$CB,sep='|'))
    if(!any(!is.na(ri))) stop('No cell matches; review mapping')
    comparison<-data.frame(author_cell=meta$cellId,author_prefix=prefix,recount_library=lib,
        recount_gem=mg,CB=paste0(cb,'-',mg),raw_match=!is.na(ri),
        called_match=!is.na(ri)&actual$called[ri],qc_match=!is.na(ri)&actual$qc[ri],
        recount_cell=actual$cell_id[ri],recount_UMIs=actual$nCount_RNA[ri],
        recount_genes=actual$nFeature_RNA[ri],recount_percent_mito=actual$percent_mito[ri])
    for(f in metrics) comparison[[paste0('author_',f)]]<-meta[[f]]
    write.tsv(comparison,file.path(out,'cell_metadata_comparison.tsv'))
    mi<-match(candidate$cell_metadata$cell_id,comparison$recount_cell)
    original.names<-names(candidate$cell_metadata)
    if(any(paste0('author_',names(meta)) %in% original.names)) stop('Author annotation columns already exist')
    candidate$cell_metadata$author_match<-!is.na(mi)
    for(f in names(meta)) candidate$cell_metadata[[paste0('author_',f)]]<-meta[[f]][mi]
    candidate$provenance$author_metadata<-list(source=path,source_MD5=unname(tools::md5sum(path)),
        mapping=mapping,comparison='Cell identity, annotation and QC metrics only; no gene-expression validation',
        counts_changed=FALSE,features_removed=FALSE)
    saveRDS(candidate,file.path(out,'GW18_qc_candidate_with_author_metadata.rds'),compress='gzip')
    corr<-list()
    for(l in sort(unique(lib))) for(f in intersect(c('ncount_rna','nCount_RNA','nfeature_rna','nFeature_RNA'),metrics)) {
        ix<-which(lib==l & comparison$qc_match)
        av<-suppressWarnings(as.numeric(as.character(meta[[f]][ix])))
        own<-if(grepl('feature',f,ignore.case=TRUE)) comparison$recount_genes[ix] else comparison$recount_UMIs[ix]
        ok<-is.finite(av)&is.finite(own)
        rho<-if(sum(ok)>=3 && sd(av[ok])>0 && sd(own[ok])>0) cor(av[ok],own[ok],method='spearman') else NA_real_
        corr[[length(corr)+1L]]<-data.frame(library=l,author_field=f,cells=sum(ok),spearman_QC_metric=rho)
    }
    write.tsv(do.call(rbind,corr),file.path(out,'QC_metric_correlations.tsv'))
    coverage<-do.call(rbind,lapply(sort(unique(lib)),function(l){ii<-which(lib==l)
        data.frame(library=l,author_cells=length(ii),raw_matches=sum(comparison$raw_match[ii]),
            called_matches=sum(comparison$called_match[ii]),qc_matches=sum(comparison$qc_match[ii]))}))
    write.tsv(coverage,file.path(out,'author_cell_coverage.tsv'))
    report<-c(report,'\nREVIEWED MAPPING APPLIED:',capture.output(print(coverage,row.names=FALSE)),
       paste('QC candidate cells with author annotation:',sum(!is.na(mi)),'/',length(mi)),
       'All QC candidate cells and all original features are retained. Counts were not normalized or changed.',
       'QC-metric correlation is not gene-expression concordance or proof of exact count recovery.')
}
writeLines(report,file.path(out,'report.txt'))
writeLines(capture.output(sessionInfo()),file.path(out,'R_sessionInfo.txt'))
# END_R
R_EMBEDDED
