args<-commandArgs(TRUE);suppressPackageStartupMessages(library(Matrix))
g<-readRDS(args[1]);a<-readRDS(args[2]);out<-args[3]
unique.symbols<-function(z)z[!is.na(z)&nzchar(z)&!duplicated(z)&!duplicated(z,fromLast=TRUE)]
shared<-intersect(unique.symbols(g$features$gene_name),unique.symbols(rownames(a$counts)))
gi<-match(shared,g$features$gene_name);ai<-match(shared,rownames(a$counts));j<-match(colnames(a$counts),g$cell_metadata$author_cellId);stopifnot(length(j)==50L,!anyNA(j))
x<-as.matrix(g$counts[gi,j,drop=FALSE]);y<-as.matrix(a$counts[ai,,drop=FALSE]);joint<-x>0&y>0
r<-matrix(NA_real_,nrow(x),ncol(x));r[joint]<-y[joint]/x[joint]
measure<-function(v){v<-v[is.finite(v)&v>0];if(length(v)<5L)return(c(n=length(v),median=NA,log2_IQR=NA,relative_range=NA,near_median_fraction=NA));m<-median(v);c(n=length(v),median=m,log2_IQR=IQR(log2(v)),relative_range=(max(v)-min(v))/m,near_median_fraction=mean(abs(v/m-1)<1e-6))}
pc<-as.data.frame(t(apply(r,2,measure)));pc$author_cellId<-colnames(a$counts)
pg<-as.data.frame(t(apply(r,1,measure)));pg$gene_symbol<-shared
write.table(pc,file.path(out,'author_counts_over_GX_per_cell_ratios.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(pg,file.path(out,'author_counts_over_GX_per_gene_ratios.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
valid<-pg$n>=5;summary<-data.frame(cells=ncol(r),genes_tested=sum(valid),per_cell_log2_ratio_IQR_median=median(pc$log2_IQR,na.rm=TRUE),per_gene_log2_ratio_IQR_median=median(pg$log2_IQR[valid],na.rm=TRUE),genes_with_nearly_constant_ratio=sum(pg$relative_range[valid]<1e-6,na.rm=TRUE),per_gene_near_median_fraction_median=median(pg$near_median_fraction[valid],na.rm=TRUE))
write.table(summary,file.path(out,'author_counts_ratio_summary.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
writeLines(c('Exploratory ratios: fractional author RNA/counts divided by original-BAM reconstructed integer GX UMIs, jointly positive exact unique shared symbols, the same 50 R7CE cells.',
 'Constant ratios within a cell would support simple cell scaling IF integer molecules were identical; constant ratios within a gene would support gene scaling under the same condition.',
 'BAM-derived counts may differ from original Cell Ranger molecule filtering; these diagnostics do not prove a normalization formula.',capture.output(print(summary,row.names=FALSE))),file.path(out,'author_counts_ratio_notes.txt'))
cat(readLines(file.path(out,'author_counts_ratio_notes.txt')),sep='\n')
