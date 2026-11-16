args<-commandArgs(TRUE);suppressPackageStartupMessages(library(Matrix));suppressPackageStartupMessages(library(jsonlite));w<-args[1];code<-args[2];author<-args[3];d<-file.path(w,'comparison_test');dir.create(d,showWarnings=FALSE)
f<-read.delim(file.path(w,'inputs/features.tsv'),header=FALSE,stringsAsFactors=FALSE);md<-read.delim(file.path(w,'inputs/author_cells.tsv'),stringsAsFactors=FALSE,check.names=FALSE)
set.seed(75);m<-sparseMatrix(i=rep(seq_len(5),nrow(md)),j=rep(seq_len(nrow(md)),each=5),x=sample(1:9,5*nrow(md),replace=TRUE),dims=c(nrow(f),nrow(md)));rownames(m)<-f[[1]];colnames(m)<-md$author_cellId
writeMM(m,file.path(d,'matrix.mtx'));zz<-gzfile(file.path(d,'matrix.mtx.gz'),'wb');writeBin(readBin(file.path(d,'matrix.mtx'),'raw',n=file.info(file.path(d,'matrix.mtx'))$size),zz);close(zz)
zz<-gzfile(file.path(d,'features.tsv.gz'),'wt');write.table(f,zz,sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE);close(zz);zz<-gzfile(file.path(d,'barcodes.tsv.gz'),'wt');writeLines(md$author_cellId,zz);close(zz)
md$nCount_GX<-colSums(m);md$nFeature_GX<-colSums(m>0);write.table(md,file.path(d,'cell_metadata.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write_json(list(full_BAM_processed=TRUE,normalized=FALSE,fixture=TRUE),file.path(d,'provenance.json'),auto_unbox=TRUE)
names(f)<-c('gene_id','gene_name','feature_type');revrow<-rev(seq_len(nrow(f)));revcol<-rev(seq_len(nrow(md)));oldmd<-md[revcol,,drop=FALSE];oldmd$cell_id<-oldmd$author_cellId
old<-list(counts=m[revrow,revcol],features=f[revrow,],cell_metadata=oldmd,provenance=list(normalized=FALSE,reference_label='GRCh38-1.2.0',features='Gene/exonic',fixture=TRUE));saveRDS(old,file.path(d,'STARsolo_fixture.rds'))
run<-function()system2('/usr/local/bin/Rscript',shQuote(c(code,d,author,file.path(d,'STARsolo_fixture.rds'))),stdout=file.path(d,'run.log'),stderr=file.path(d,'run.err'))
stopifnot(run()==0L);s<-read.delim(file.path(d,'GX_vs_STARsolo_summary.tsv'));stopifnot(s$differing_entries==0,s$sum_absolute_difference==0,s$GX_total_UMIs==s$STARsolo_total_UMIs)
old$counts[nrow(old$counts),ncol(old$counts)]<-old$counts[nrow(old$counts),ncol(old$counts)]+3;saveRDS(old,file.path(d,'STARsolo_fixture.rds'));stopifnot(run()==0L);s<-read.delim(file.path(d,'GX_vs_STARsolo_summary.tsv'));stopifnot(s$differing_entries==1,s$sum_absolute_difference==3)
cat('PASS: all 33694 gene IDs and 1759 cells retained; reversed gene/cell orders matched exactly; introduced 3-UMI difference detected as one differing entry; author 50-cell comparison exported.\n')
