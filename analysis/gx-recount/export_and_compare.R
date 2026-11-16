args<-commandArgs(TRUE);stopifnot(length(args)>=2L)
suppressPackageStartupMessages(library(Matrix));suppressPackageStartupMessages(library(jsonlite))
dir<-normalizePath(args[1]);author_path<-normalizePath(args[2]);old_path<-if(length(args)>=3L)args[3]else''
readtab<-function(n)read.delim(gzfile(file.path(dir,n)),header=FALSE,quote='',comment.char='',stringsAsFactors=FALSE,check.names=FALSE)
f<-readtab('features.tsv.gz');names(f)<-c('gene_id','gene_name','feature_type');bc<-readtab('barcodes.tsv.gz')[[1]]
con<-gzfile(file.path(dir,'matrix.mtx.gz'),'rt');m<-as(readMM(con),'CsparseMatrix');close(con)
md<-read.delim(file.path(dir,'cell_metadata.tsv'),quote='',comment.char='',check.names=FALSE,stringsAsFactors=FALSE)
p<-fromJSON(file.path(dir,'provenance.json'),simplifyVector=FALSE)
stopifnot(isTRUE(p$full_BAM_processed),!isTRUE(p$normalized),nrow(m)==33694L,ncol(m)==1759L,identical(bc,md$author_cellId),!anyDuplicated(f$gene_id),!anyDuplicated(bc),all(is.finite(m@x)),all(m@x>=0),all(m@x==round(m@x)),all(colSums(m)==md$nCount_GX),all(colSums(m>0)==md$nFeature_GX),sum(m)>0)
rownames(m)<-f$gene_id;colnames(m)<-bc;rownames(md)<-bc
x<-list(counts=m,features=f,cell_metadata=md,provenance=p)
saveRDS(x,file.path(dir,'GW18_R7CE_BAM_GX_author_cells_raw_counts.rds'),compress='gzip')
safe<-function(a,b,method)if(length(a)>=3L&&sd(a)>0&&sd(b)>0)cor(a,b,method=method)else NA_real_
unique_ids<-function(z)z[!is.na(z)&nzchar(z)&!duplicated(z)&!duplicated(z,fromLast=TRUE)]
a<-readRDS(author_path);ab<-colnames(a$data);aj<-which(ab%in%bc);stopifnot(length(aj)==50L)
asy<-rownames(a$data);shared<-intersect(unique_ids(f$gene_name),unique_ids(asy));stopifnot(length(shared)>10000L)
gi<-match(shared,f$gene_name);ai<-match(shared,asy);tot<-as.numeric(colSums(m))
ac<-do.call(rbind,lapply(aj,function(j){k<-match(ab[j],bc);g<-as.numeric(m[gi,k]);ad<-as.numeric(a$data[ai,j]);ac<-as.numeric(a$counts[ai,j]);u<-g>0|ad>0;gn<-log1p(10000*g/tot[k]);data.frame(author_cellId=ab[j],shared_genes=length(shared),expressed_union_genes=sum(u),GX_UMIs=tot[k],spearman_GX_vs_author_data=safe(gn[u],ad[u],'spearman'),pearson_logCP10k_GX_vs_author_data=safe(gn[u],ad[u],'pearson'),jaccard=if(any(u))sum(g>0&ad>0)/sum(u)else NA_real_,exploratory_pearson_integer_GX_vs_fractional_author_counts=safe(g[u],ac[u],'pearson'))}))
write.table(ac,file.path(dir,'GX_vs_author_50_cells.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
summary<-data.frame(cells=nrow(ac),spearman_median=median(ac$spearman_GX_vs_author_data,na.rm=TRUE),pearson_median=median(ac$pearson_logCP10k_GX_vs_author_data,na.rm=TRUE),jaccard_median=median(ac$jaccard,na.rm=TRUE))
write.table(summary,file.path(dir,'GX_vs_author_summary.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
report<-c('BAM GX reconstruction finished: author cohort R7CE 1759 cells, all 33694 legacy reference genes.',paste('Total integer UMIs:',sum(m)),'Raw matrices are never normalized on disk. Author RNA/counts are processed fractional values, not original integer UMI truth.','Author comparison uses exact unique symbols; excludes shared zeros per cell; GX library totals include all reference genes.','MAPQ255 single-GX CB/UB reconstruction is not guaranteed identical to original Cell Ranger molecule filtering.',capture.output(print(summary)))
if(nzchar(old_path)&&file.exists(old_path)){
 o<-readRDS(old_path);stopifnot(is.list(o),!isTRUE(o$provenance$normalized),identical(o$provenance$reference_label,'GRCh38-1.2.0'),!anyDuplicated(o$features$gene_id),!anyDuplicated(o$cell_metadata$author_cellId))
 if(grepl('GeneFull',paste(unlist(o$provenance$features),collapse=' '),fixed=TRUE))stop('Expected exonic Gene STARsolo output')
 oi<-match(f$gene_id,o$features$gene_id);oj<-match(bc,o$cell_metadata$author_cellId);stopifnot(!anyNA(oi),!anyNA(oj),identical(as.character(o$cell_metadata$CB[oj]),md$CB),all(grepl('^R7CE:',o$cell_metadata$library[oj])))
 om<-o$counts[oi,oj,drop=FALSE];stopifnot(all(om@x>=0),all(om@x==round(om@x)))
 ot<-as.numeric(colSums(o$counts[,oj,drop=FALSE]));d<-drop0(m-om)
 pcs<-do.call(rbind,lapply(seq_along(bc),function(j){g<-as.numeric(m[,j]);s<-as.numeric(om[,j]);u<-g>0|s>0;gn<-log1p(10000*g/tot[j]);sn<-log1p(10000*s/ot[j]);data.frame(author_cellId=bc[j],GX_UMIs=tot[j],STARsolo_UMIs=ot[j],differing_gene_entries=sum(g!=s),absolute_UMI_difference=sum(abs(g-s)),identical_fraction_in_expressed_union=if(any(u))mean(g[u]==s[u])else NA_real_,spearman=safe(g[u],s[u],'spearman'),raw_pearson=safe(g[u],s[u],'pearson'),logCP10k_pearson=safe(gn[u],sn[u],'pearson'),jaccard=if(any(u))sum(g>0&s>0)/sum(u)else NA_real_)}))
 write.table(pcs,file.path(dir,'GX_vs_STARsolo_per_cell.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
 sm<-data.frame(cells=ncol(m),genes=nrow(m),GX_total_UMIs=sum(m),STARsolo_total_UMIs=sum(om),differing_entries=length(d@x),sum_absolute_difference=sum(abs(d@x)),spearman_median=median(pcs$spearman,na.rm=TRUE),raw_pearson_median=median(pcs$raw_pearson,na.rm=TRUE),logCP10k_pearson_median=median(pcs$logCP10k_pearson,na.rm=TRUE),jaccard_median=median(pcs$jaccard,na.rm=TRUE))
 write.table(sm,file.path(dir,'GX_vs_STARsolo_summary.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
 genes<-data.frame(gene_id=f$gene_id,gene_symbol=f$gene_name,GX_UMIs=rowSums(m),STARsolo_UMIs=rowSums(om),absolute_difference=rowSums(abs(d)));genes<-genes[order(-genes$absolute_difference),]
 write.table(genes,file.path(dir,'GX_vs_STARsolo_gene_differences.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
 report<-c(report,'','STARsolo comparison: same exact gene IDs and full author cell IDs, raw equality first, correlations exclude shared-zero genes.',capture.output(print(sm)))
 writeLines('COMPLETE',file.path(dir,'COMPARISON_COMPLETE'))
}else{report<-c(report,'','WAITING_FOR_STARSOLO_INPUT: put reviewed old-reference R7CE subset at inputs/STARsolo_GW18_R7CE_GRCh38_1.2.0.rds, then rerun this script.');writeLines('WAITING_FOR_STARSOLO_INPUT',file.path(dir,'COMPARISON_STATUS'))}
writeLines(report,file.path(dir,'report.txt'));writeLines(capture.output(sessionInfo()),file.path(dir,'R_sessionInfo.txt'));cat(paste(report,collapse='\n'),'\n')
