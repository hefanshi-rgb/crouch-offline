compare_gx_starsolo <- function(gx_path, star_path, author_path, output=NULL, strict_source=TRUE) {
 suppressPackageStartupMessages(library(Matrix))
 g<-readRDS(gx_path); s<-readRDS(star_path); a<-readRDS(author_path)
 stopifnot(is.list(g),is.list(s),!isTRUE(g$provenance$normalized),!isTRUE(s$provenance$normalized),
  isTRUE(g$provenance$full_BAM_processed),identical(s$provenance$reference_label,'GRCh38-1.2.0'),
  nrow(g$counts)==33694L,ncol(g$counts)==1759L,nrow(s$counts)==33694L,
  !anyDuplicated(g$features$gene_id),!anyDuplicated(s$features$gene_id),
  identical(colnames(g$counts),g$cell_metadata$author_cellId),
  !anyDuplicated(g$cell_metadata$author_cellId),!anyDuplicated(s$cell_metadata$author_cellId),
  all(is.finite(g$counts@x)),all(g$counts@x>=0),all(g$counts@x==round(g$counts@x)),
  all(is.finite(s$counts@x)),all(s$counts@x>=0),all(s$counts@x==round(s$counts@x)))
 if(grepl('GeneFull',paste(unlist(s$provenance$features),collapse=' '),fixed=TRUE))stop('Select the exonic Gene STARsolo matrix')
 if(strict_source){
  full<-ncol(s$counts)==13906L && sum(s$counts)==28758442
  exported<-ncol(s$counts)==1759L && identical(as.numeric(s$provenance$GX_comparison_export$source_cells),13906) && identical(as.numeric(s$provenance$GX_comparison_export$source_total_UMIs),28758442)
  if(!full && !exported)stop('STARsolo input must be the reviewed full old-reference GW18 matrix, or its validated R7CE export')
 }
 ids<-as.character(g$features$gene_id); cells<-as.character(g$cell_metadata$author_cellId)
 si<-match(ids,s$features$gene_id);sj<-match(cells,s$cell_metadata$author_cellId)
 stopifnot(!anyNA(si),!anyNA(sj),identical(as.character(s$cell_metadata$CB[sj]),as.character(g$cell_metadata$CB)),all(grepl('^R7CE:',s$cell_metadata$library[sj])))
 st<-as.numeric(colSums(s$counts[,sj,drop=FALSE])); sm<-s$counts[si,sj,drop=FALSE];rm(s);gc(FALSE)
 m<-g$counts;gt<-as.numeric(colSums(m));stopifnot(all(gt>0),all(st>0));d<-drop0(m-sm)
 safe<-function(x,y,method)if(length(x)>=3L&&sd(x)>0&&sd(y)>0)cor(x,y,method=method)else NA_real_
 pcs<-do.call(rbind,lapply(seq_along(cells),function(j){
  x<-as.numeric(m[,j]);y<-as.numeric(sm[,j]);u<-x>0|y>0
  data.frame(author_cellId=cells[j],GX_UMIs=gt[j],STARsolo_UMIs=st[j],GX_detected_genes=sum(x>0),STARsolo_detected_genes=sum(y>0),
   differing_gene_entries=sum(x!=y),absolute_UMI_difference=sum(abs(x-y)),identical_fraction_in_detected_union=mean(x[u]==y[u]),
   spearman=safe(x[u],y[u],'spearman'),raw_pearson=safe(x[u],y[u],'pearson'),
   logCP10k_pearson=safe(log1p(10000*x[u]/gt[j]),log1p(10000*y[u]/st[j]),'pearson'),detection_jaccard=sum(x>0&y>0)/sum(u))
 }))
 med<-function(x)if(any(is.finite(x)))median(x[is.finite(x)])else NA_real_
 summary<-data.frame(cells=ncol(m),genes=nrow(m),GX_total_UMIs=sum(m),STARsolo_total_UMIs=sum(sm),
  differing_entries=length(d@x),sum_absolute_UMI_difference=sum(abs(d@x)),spearman_median=med(pcs$spearman),
  raw_pearson_median=med(pcs$raw_pearson),logCP10k_pearson_median=med(pcs$logCP10k_pearson),
  detection_jaccard_median=med(pcs$detection_jaccard),identical_fraction_in_detected_union_median=med(pcs$identical_fraction_in_detected_union))
 unique.symbols<-function(z)z[!is.na(z)&nzchar(z)&!duplicated(z)&!duplicated(z,fromLast=TRUE)]
 sym<-as.character(g$features$gene_name);asy<-rownames(a$data);common<-intersect(unique.symbols(sym),unique.symbols(asy))
 stopifnot(length(common)>10000L,identical(colnames(a$counts),colnames(a$data)))
 ai<-match(common,asy);gi<-match(common,sym);aj<-which(colnames(a$data)%in%cells);stopifnot(length(aj)==50L)
 au<-do.call(rbind,lapply(aj,function(j){k<-match(colnames(a$data)[j],cells);x<-as.numeric(m[gi,k]);y<-as.numeric(sm[gi,k]);z<-as.numeric(a$data[ai,j]);u<-x>0|y>0|z>0
  nx<-log1p(10000*x/gt[k]);ny<-log1p(10000*y/st[k]); ux<-x>0|z>0;uy<-y>0|z>0
  data.frame(author_cellId=cells[k],shared_unique_gene_symbols=length(common),fixed_three_way_expressed_union=sum(u),
   GX_spearman=safe(nx[u],z[u],'spearman'),STARsolo_spearman=safe(ny[u],z[u],'spearman'),
   GX_pearson=safe(nx[u],z[u],'pearson'),STARsolo_pearson=safe(ny[u],z[u],'pearson'),
   GX_detection_jaccard=sum(x>0&z>0)/sum(ux),STARsolo_detection_jaccard=sum(y>0&z>0)/sum(uy))
 }))
 au_summary<-data.frame(cells=nrow(au),GX_spearman=med(au$GX_spearman),STARsolo_spearman=med(au$STARsolo_spearman),
  GX_pearson=med(au$GX_pearson),STARsolo_pearson=med(au$STARsolo_pearson),
  GX_detection_jaccard=med(au$GX_detection_jaccard),STARsolo_detection_jaccard=med(au$STARsolo_detection_jaccard))
 genes<-data.frame(gene_id=ids,gene_symbol=sym,GX_UMIs=as.numeric(rowSums(m)),STARsolo_UMIs=as.numeric(rowSums(sm)),absolute_UMI_difference=as.numeric(rowSums(abs(d))))
 genes<-genes[order(-genes$absolute_UMI_difference),]
 if(is.null(output))output<-file.path(dirname(normalizePath(gx_path)),paste0('GX_STARsolo_validation_',format(Sys.time(),'%Y%m%d_%H%M%S')))
 if(dir.exists(output))stop('Output directory exists; preserving it');dir.create(output,recursive=TRUE)
 tab<-function(x,name)write.table(x,file.path(output,name),sep='\t',quote=FALSE,row.names=FALSE)
 tab(summary,'GX_vs_STARsolo_summary.tsv');tab(pcs,'GX_vs_STARsolo_per_cell.tsv');tab(genes,'gene_differences.tsv');tab(au,'GX_STARsolo_vs_author_50_cells.tsv');tab(au_summary,'GX_STARsolo_vs_author_summary.tsv')
 grDevices::pdf(file.path(output,'validation_plots.pdf'),width=10,height=8)
 par(mfrow=c(2,2));hist(pcs$spearman,main='GX vs STARsolo: per-cell Spearman',xlab='Spearman',col='grey80');hist(pcs$logCP10k_pearson,main='GX vs STARsolo: logCP10k Pearson',xlab='Pearson',col='grey80')
 plot(st,gt,pch=16,cex=.4,xlab='STARsolo UMIs per cell',ylab='BAM GX UMIs per cell');abline(0,1,col='red')
 boxplot(au[,c('GX_spearman','STARsolo_spearman')],names=c('GX','STARsolo'),ylab='Spearman vs author RNA/data',main='Same 50 cells, genes and expressed-gene mask');grDevices::dev.off()
 report<-c('Full-depth GW18 R7CE: original-BAM GX UMI reconstruction versus reviewed legacy-reference STARsolo Gene counts.',
  'All 1759 author cells and 33694 genes retained. No new cell QC, gene filtering or normalization of saved raw counts.',
  'First assess exact raw-count differences; correlations exclude shared-zero genes per cell.',
  'Author-expression comparison uses the SAME 50 cells, unique common symbols, and fixed THREE-way expressed-gene union.',
  'All original reference genes contribute to each per-cell normalization denominator.',
  'Author RNA/data are processed expression and RNA/counts are fractional; neither is original integer-UMI ground truth.',
  'Simple CB/UB/GX deduplication is not asserted to exactly reproduce original Cell Ranger molecule filtering.',
  paste('GX input:',normalizePath(gx_path)),paste('STARsolo input:',normalizePath(star_path)),paste('Author input:',normalizePath(author_path)),
  capture.output(print(summary,row.names=FALSE,digits=6)),capture.output(print(au_summary,row.names=FALSE,digits=6)))
 writeLines(report,file.path(output,'report.txt'));writeLines(capture.output(sessionInfo()),file.path(output,'R_sessionInfo.txt'));cat(paste(report,collapse='\n'),'\nOutput:',output,'\n')
 invisible(list(summary=summary,author_summary=au_summary,output=output))
}
if(interactive()){
 cat('1. Select GW18_R7CE_BAM_GX_author_cells_raw_counts.rds\n');gx_path<-file.choose()
 cat('2. Select reviewed GW18_author_cells_GRCh38_1.2.0_raw_counts.rds\n');star_path<-file.choose()
 cat('3. Select author_50_R7CE_slots.rds from this bundle\n');author_path<-file.choose()
 gx_validation<-compare_gx_starsolo(gx_path,star_path,author_path)
}
