# Run in Windows RStudio. Choose the reviewed full old-reference GW18 RDS.
suppressPackageStartupMessages(library(Matrix))
cat('Select GW18_author_cells_GRCh38_1.2.0_raw_counts.rds\n')
p<-file.choose();x<-readRDS(p)
stopifnot(is.list(x),nrow(x$counts)==33694L,ncol(x$counts)==13906L,sum(x$counts)==28758442,!isTRUE(x$provenance$normalized),identical(x$provenance$reference_label,'GRCh38-1.2.0'),identical(colnames(x$counts),x$cell_metadata$cell_id),!anyDuplicated(x$cell_metadata$author_cellId))
i<-which(grepl('^all_R7_ctx_endo_',x$cell_metadata$author_cellId))
stopifnot(length(i)==1759L,all(grepl('^R7CE:',x$cell_metadata$library[i])))
x$counts<-x$counts[,i,drop=FALSE];x$cell_metadata<-x$cell_metadata[i,,drop=FALSE]
x$provenance$GX_comparison_export<-list(source_file=basename(p),source_MD5=unname(tools::md5sum(p)),source_cells=13906L,source_total_UMIs=28758442,selection='All exact author R7 ctx endo cells; no additional QC, filtering or normalization')
out<-file.path(path.expand('~/Desktop'),'STARsolo_GW18_R7CE_GRCh38_1.2.0.rds')
if(!dir.exists(dirname(out)))out<-file.path(dirname(p),basename(out))
saveRDS(x,out,compress='gzip')
cat('\nSaved:',normalizePath(out),'\nCells:',ncol(x$counts),'\nGenes:',nrow(x$counts),'\nUMIs:',sum(x$counts),'\nSize MB:',round(file.info(out)$size/1024^2,2),'\n')
