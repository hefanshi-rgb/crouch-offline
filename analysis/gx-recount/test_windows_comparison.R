args<-commandArgs(TRUE);source(args[1]);w<-args[2];author<-args[3]
a<-readRDS(author);aj<-which(grepl('^all_R7_ctx_endo_',colnames(a$data)));stopifnot(length(aj)==50L)
a$counts<-a$counts[,aj,drop=FALSE];a$data<-a$data[,aj,drop=FALSE];a$meta<-a$meta[aj,,drop=FALSE]
saveRDS(a,file.path(w,'inputs/author_50_R7CE_slots.rds'),compress='gzip')
fx<-file.path(w,'comparison_test');out<-file.path(w,paste0('windows_comparison_test_',format(Sys.time(),'%Y%m%d_%H%M%S')))
r<-compare_gx_starsolo(file.path(fx,'GW18_R7CE_BAM_GX_author_cells_raw_counts.rds'),file.path(fx,'STARsolo_fixture.rds'),file.path(w,'inputs/author_50_R7CE_slots.rds'),out,strict_source=FALSE)
stopifnot(r$summary$differing_entries==1,r$summary$sum_absolute_UMI_difference==3)
failed<-inherits(try(compare_gx_starsolo(file.path(fx,'GW18_R7CE_BAM_GX_author_cells_raw_counts.rds'),file.path(fx,'STARsolo_fixture.rds'),file.path(w,'inputs/author_50_R7CE_slots.rds'),paste0(out,'_must_not_exist'),strict_source=TRUE),silent=TRUE),'try-error')
stopifnot(failed,!dir.exists(paste0(out,'_must_not_exist')))
cat('PASS: column/gene permutation aligned; 3-UMI single-entry difference detected; unified author 50-cell comparisons computed; unreviewed fixture blocked by default source validation.\n')
