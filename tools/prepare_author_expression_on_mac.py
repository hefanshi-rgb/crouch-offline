#!/usr/bin/env python3
"""Stream the author matrix locally; relay a stratified GW18 expression sample only."""
import argparse,collections,csv,gzip,hashlib,importlib.util,io,json,math,pathlib,random,re,shutil,subprocess,sys,tarfile
CLUSTER_SCRIPT="#!/bin/bash\nset -Eeuo pipefail\numask 077\nBASE=/public4/usersdir/guojunping/guojunping_2026-11-08_001/upload/Crouch\nSELF=$(cd \"$(dirname \"${BASH_SOURCE[0]}\")\" && pwd -P)/$(basename \"${BASH_SOURCE[0]}\")\nif [[ \"${1:-}\" != --worker ]]; then\n    if ! type module >/dev/null 2>&1; then\n        for init in /etc/profile.d/modules.sh /usr/share/Modules/init/bash /usr/share/lmod/lmod/init/bash; do\n            if [[ -r \"$init\" ]]; then set +u; source \"$init\"; set -u; break; fi\n        done\n    fi\n    set +u\n    module use /public/software/modules\n    module load compiler/gcc/7.4.0 system/R/4.3.1\n    set -u\n    lib=$(gfortran -print-file-name=libgfortran.so.4)\n    [[ \"$lib\" == /* && -r \"$lib\" ]] || { echo 'GCC runtime missing'; exit 1; }\n    export LD_LIBRARY_PATH=\"${lib%/*}${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}\"\n    EXPR_R=$(command -v Rscript);export EXPR_R\n    mkdir -p \"$BASE/logs\"\n    reply=$(sbatch --parsable --export=ALL -p bigmem --nodelist=comput159 \\\n      -N 1 -n 1 -c 2 --mem=32G -t 02:00:00 --job-name=gw18-expression-check \\\n      --chdir=\"$BASE\" --output=\"$BASE/logs/gw18-expression-check-%j.log\" \"$SELF\" --worker)\n    job=${reply%%;*}; [[ \"$job\" =~ ^[0-9]+$ ]] || { echo \"$reply\";exit 1; }\n    printf 'Job=%s\\nLog=%s/logs/gw18-expression-check-%s.log\\nResult=%s/recount_slurm/GW18/expression_%s\\n' \"$job\" \"$BASE\" \"$job\" \"$BASE\" \"$job\"\n    exit 0\nfi\n: \"${SLURM_JOB_ID:?Submit using the launcher}\"\nOUT=\"$BASE/recount_slurm/GW18/expression_${SLURM_JOB_ID}\"\nmkdir \"$OUT\";cp \"$SELF\" \"$OUT/submitted_script.sh\"\nexport TZ=UTC OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1\nexport R_LIBS=\"$BASE/software/r-libraries/R-4.3.1-Matrix-1.6.0${R_LIBS:+:$R_LIBS}\"\nexport R_LIBS_USER=\"$R_LIBS\"\nsed -n '/^# BEGIN_R$/,/^# END_R$/p' \"$SELF\" > \"$OUT/expression_check.R\"\n\"$EXPR_R\" --vanilla \"$OUT/expression_check.R\" \\\n \"$BASE/recount_slurm/GW18/meta_compare_1717231/GW18_author_cells_raw_counts.rds\" \\\n \"$BASE/author_expression_check\" \"$OUT\"\ntouch \"$OUT/COMPLETE\"\nprintf '\\nSUCCESS: %s\\n' \"$OUT\"\nexit 0\n: <<'R_EMBEDDED'\n# BEGIN_R\nsuppressPackageStartupMessages(library(Matrix))\na<-commandArgs(TRUE); stopifnot(length(a)>=3L)\nraw.path<-a[1]; src<-a[2]; out<-a[3]; min.genes<-if(length(a)>=4L) as.integer(a[4]) else 1000L\nown<-readRDS(raw.path)\nstopifnot(own$provenance$mode=='author_cohort',!isTRUE(own$provenance$normalized),\n    all(own$counts@x>=0),all(own$counts@x==round(own$counts@x)))\ntab<-function(name) read.delim(gzfile(file.path(src,name),'rt'),header=FALSE,\n    quote='',comment.char='',stringsAsFactors=FALSE,check.names=FALSE)\naf<-tab('features.tsv.gz'); ab<-tab('barcodes.tsv.gz')[[1]]\nmeta<-read.delim(file.path(src,'meta.tsv'),quote='',comment.char='',stringsAsFactors=FALSE,check.names=FALSE)\nstopifnot(!anyDuplicated(ab),!anyDuplicated(meta$cellId),all(toupper(trimws(meta$age))=='GW18'))\nmeta<-meta[match(ab,meta$cellId),,drop=FALSE];stopifnot(identical(ab,meta$cellId))\nstopifnot(!anyDuplicated(own$cell_metadata$author_cellId),\n    identical(colnames(own$counts),own$cell_metadata$cell_id))\nci<-match(ab,own$cell_metadata$author_cellId)\nif(anyNA(ci)) stop('Some sampled author cells do not match the final author-cohort RDS')\npairs<-function(x,y) {\n    x<-as.character(x);y<-as.character(y)\n    ux<-!duplicated(x)&!duplicated(x,fromLast=TRUE)&!is.na(x)&nzchar(x)\n    uy<-!duplicated(y)&!duplicated(y,fromLast=TRUE)&!is.na(y)&nzchar(y)\n    j<-match(x,y);ok<-ux&!is.na(j);ok[ok]<-uy[j[ok]]\n    list(author=which(ok),own=j[ok])\n}\np<-pairs(af[[1]],own$features$gene_id);method<-'Exact unique gene IDs'\nif(length(p$author)<min.genes) {\n    sy<-if(ncol(af)>=2L) af[[2]] else af[[1]]\n    q<-pairs(sy,own$features$gene_name)\n    if(length(q$author)>length(p$author)) {p<-q;method<-'Exact unique gene symbols'}\n}\nif(length(p$author)<min.genes) stop('Too few common genes; inspect feature identifiers before comparing')\ncon<-gzfile(file.path(src,'matrix.mtx.gz'),'rt'); am<-as(readMM(con),'CsparseMatrix');close(con)\nstopifnot(nrow(am)==nrow(af),ncol(am)==length(ab),all(is.finite(am@x)))\nwrite.tsv<-function(d,name) write.table(d,file.path(out,name),sep='\\t',row.names=FALSE,quote=FALSE)\nwrite.tsv(data.frame(author_row=p$author,recount_row=p$own,\n    gene_id=own$features$gene_id[p$own],gene_name=own$features$gene_name[p$own]),'gene_matches.tsv')\nraw<-own$counts[p$own,ci,drop=FALSE]\n# The denominator includes ALL counted genes, before common-gene subsetting.\ntotals<-as.numeric(Matrix::colSums(own$counts[,ci,drop=FALSE]))\nstopifnot(all(totals>0))\ncor.safe<-function(x,y,method) if(length(x)>=3L&&sd(x)>0&&sd(y)>0) cor(x,y,method=method) else NA_real_\nrows<-lapply(seq_along(ab),function(j){\n    rv<-as.numeric(raw[,j]); nv<-log1p(rv/totals[j]*10000); av<-as.numeric(am[p$author,j])\n    use<-rv!=0|av!=0; both<-rv>0&av!=0\n    data.frame(author_cell=ab[j],recount_cell=own$cell_metadata$cell_id[ci[j]],\n      library=own$cell_metadata$library[ci[j]],\n      author_cell_type=if('cell_type'%in%names(meta)) as.character(meta$cell_type[j]) else NA_character_,\n      common_genes=length(rv),genes_in_union=sum(use),genes_detected_in_both=sum(both),\n      detection_jaccard=if(sum(use)>0) sum(both)/sum(use) else NA_real_,\n      raw_UMIs=totals[j],UMI_fraction_in_common_genes=sum(rv)/totals[j],\n      spearman_union=cor.safe(nv[use],av[use],'spearman'),\n      pearson_logNormalize_vs_author=cor.safe(nv[use],av[use],'pearson'),\n      spearman_both_detected=cor.safe(nv[both],av[both],'spearman'))\n})\nres<-do.call(rbind,rows);write.tsv(res,'cell_expression_correlations.tsv')\nsumm<-do.call(rbind,lapply(split(res,res$library),function(z){\n    data.frame(library=z$library[1],cells=nrow(z),finite_spearman=sum(is.finite(z$spearman_union)),\n      spearman_q25=unname(quantile(z$spearman_union,.25,na.rm=TRUE)),\n      spearman_median=median(z$spearman_union,na.rm=TRUE),\n      spearman_q75=unname(quantile(z$spearman_union,.75,na.rm=TRUE)),\n      pearson_median=median(z$pearson_logNormalize_vs_author,na.rm=TRUE),\n      detection_jaccard_median=median(z$detection_jaccard,na.rm=TRUE))\n}));rownames(summ)<-NULL;write.tsv(summ,'library_expression_summary.tsv')\nreport<-c('Expression concordance in a stratified sample, not validation of every GW18 cell.',\n 'Our counts: log1p(10000 * count / total UMI across ALL original genes).',\n 'Author values: retained exactly as supplied; their transformation has not been established.',\n 'Shared-zero genes excluded from both correlation metrics.',\n 'Spearman evaluates ranking; monotonic per-cell normalization does not change that ranking.',\n 'Pearson is exploratory until the author normalization method is confirmed.',\n 'High correlation is not proof of exact recovery of author raw counts.',\n paste('Gene matching:',method),paste('Common genes:',length(p$author)),\n paste('Author features:',nrow(af)),paste('Recount features:',nrow(own$features)),\n paste('Matched sampled cells:',nrow(res)),\n paste('Author supplied values: min=',min(am@x),'max=',max(am@x),'negative entries=',sum(am@x<0)),\n capture.output(print(summ,row.names=FALSE,digits=4)))\nwriteLines(report,file.path(out,'expression_report.txt'))\nif(!any(is.finite(res$spearman_union))) stop('No finite correlations; inspect expression_report.txt')\npdf(file.path(out,'expression_concordance.pdf'),width=11,height=8.5,useDingbats=FALSE)\npar(mfrow=c(2,2),mar=c(6,4,3,1))\nshort<-sub(':.*','',res$library)\nboxplot(res$spearman_union~short,ylim=c(-1,1),ylab='Spearman (expressed-gene union)',xlab='Library')\nhist(res$spearman_union,breaks=30,xlim=c(-1,1),main='Per-cell Spearman',xlab='Spearman')\nhist(res$pearson_logNormalize_vs_author,breaks=30,xlim=c(-1,1),main='Pearson: exploratory',xlab='Pearson')\nplot(res$detection_jaccard,res$spearman_union,pch=16,cex=.6,xlim=c(0,1),ylim=c(-1,1),\n    xlab='Detected-gene Jaccard',ylab='Spearman',main='Detection overlap and expression ranking')\npar(mfrow=c(2,2),mar=c(4,4,3,1))\nfor(lib in sort(unique(res$library))) {\n    ii<-which(res$library==lib&is.finite(res$spearman_union))\n    if(!length(ii)) next\n    j<-ii[which.min(abs(res$spearman_union[ii]-median(res$spearman_union[ii])))]\n    rv<-as.numeric(raw[,j]);nv<-log1p(rv/totals[j]*10000);av<-as.numeric(am[p$author,j]);use<-which(rv!=0|av!=0)\n    if(length(use)>10000L) use<-use[unique(round(seq(1,length(use),length.out=10000)))]\n    plot(nv[use],av[use],pch=16,cex=.3,col=rgb(0,0,0,.25),\n      xlab='Our log1p CP10k',ylab='Author supplied expression',\n      main=paste(sub(':.*','',lib),'median-rho cell; rho=',round(res$spearman_union[j],3)))\n}\ndev.off()\nwriteLines(capture.output(sessionInfo()),file.path(out,'R_sessionInfo.txt'))\ncat(report,sep='\\n')\n# END_R\nR_EMBEDDED\n"
PREFIXES=('all_R7_ctx_endo_','all_R7_ctx_peri_','all_R7_GE_endo_','all_R7_GE_peri_')
def sha(b):return hashlib.sha256(b).hexdigest()
def opened(p):return gzip.open(p,'rt',encoding='utf-8-sig',newline='')
def main():
 ap=argparse.ArgumentParser(description=__doc__)
 ap.add_argument('--directory',type=pathlib.Path)
 ap.add_argument('--repo',type=pathlib.Path,default=pathlib.Path.home()/'Documents/GitHub/crouch-offline')
 ap.add_argument('--output',type=pathlib.Path,default=pathlib.Path.home()/'Desktop/crouch-expression-check-transfer')
 ap.add_argument('--per-library',type=int,default=50);ap.add_argument('--seed',type=int,default=2718)
 a=ap.parse_args()
 if a.per_library<1:raise SystemExit('per-library must be positive')
 if a.directory is None:
  result=subprocess.run(['osascript','-e','POSIX path of (choose folder with prompt "选择同时包含作者 meta.tsv、matrix.mtx.gz、features.tsv.gz、barcodes.tsv.gz 的文件夹")'],check=True,capture_output=True,text=True)
  a.directory=pathlib.Path(result.stdout.rstrip('\r\n'))
 for name in ('meta.tsv','matrix.mtx.gz','features.tsv.gz','barcodes.tsv.gz'):
  if not (a.directory/name).is_file():raise SystemExit('Missing '+str(a.directory/name))
 helper=a.repo/'tools/prepare_author_meta_on_mac.py'
 if not helper.is_file():raise SystemExit('Missing previous metadata relay helper: '+str(helper))
 spec=importlib.util.spec_from_file_location('metadata_relay',helper);relay=importlib.util.module_from_spec(spec);spec.loader.exec_module(relay)
 with (a.directory/'meta.tsv').open(encoding='utf-8-sig',newline='') as f:
  r=csv.DictReader(f,delimiter='\t');fields=r.fieldnames or []
  if not {'cellId','age'}.issubset(fields):raise SystemExit('Metadata needs cellId and age')
  rows=[x for x in r if x['age'].strip().upper()=='GW18']
 if len({x['cellId'] for x in rows})!=len(rows):raise SystemExit('Duplicate author GW18 cell IDs')
 groups={p:[] for p in PREFIXES}
 for row in rows:
  m=re.fullmatch(r'(.*?)([ACGT]{16})(?:-([1-9][0-9]*))?',row['cellId'])
  if not m or m[1] not in groups:raise SystemExit('Unexpected GW18 prefix; inspect '+row['cellId'])
  groups[m[1]].append(row)
 rng=random.Random(a.seed); selected={}; selection={}
 for p in PREFIXES:
  types=collections.defaultdict(list)
  for row in sorted(groups[p],key=lambda x:x['cellId']):types[row.get('cell_type','unknown')].append(row)
  if not types:raise SystemExit('No author cells for '+p)
  for pool in types.values():rng.shuffle(pool)
  chosen=[];goal=min(a.per_library,len(groups[p]))
  while len(chosen)<goal:
   for name in sorted(types):
    if types[name] and len(chosen)<goal:chosen.append(types[name].pop())
  for row in chosen:selected[row['cellId']]=row
  selection[p]={'available':len(groups[p]),'selected':len(chosen),'cell_types':dict(collections.Counter(x.get('cell_type','unknown') for x in chosen))}
 with opened(a.directory/'barcodes.tsv.gz') as f:barcodes=[line.rstrip('\r\n').split('\t')[0] for line in f]
 if len(set(barcodes))!=len(barcodes):raise SystemExit('Duplicate author matrix cell IDs; explicit mapping required')
 if not set(selected).issubset(barcodes):raise SystemExit('Metadata and matrix cell IDs do not match EXACTLY; do not strip prefixes or infer positional correspondence')
 columns={};kept=[]
 for j,cell in enumerate(barcodes,1):
  if cell in selected:kept.append(cell);columns[j]=len(kept)
 release=a.repo/'author-expression-check-v1'
 if a.output.exists() or release.exists():raise SystemExit('Existing output/release preserved. Use a fresh output or inspect the previous run; do not replace a published release silently.')
 a.output.mkdir(parents=True);sample=a.output/'author_expression_check';sample.mkdir()
 with opened(a.directory/'features.tsv.gz') as f:features=sum(1 for line in f)
 matrix=a.directory/'matrix.mtx.gz';before=matrix.stat();seen=kept_nz=0
 with opened(matrix) as f,(sample/'.entries.tmp').open('w',encoding='utf-8') as body:
  header=f.readline().strip().lower().split()
  if len(header)!=5 or header[:3]!=['%%matrixmarket','matrix','coordinate'] or header[3] not in ('real','integer') or header[4]!='general':raise SystemExit('Expected a real/integer general coordinate MatrixMarket file')
  for line in f:
   if line.strip() and not line.startswith('%'):nr,nc,nz=map(int,line.split());break
  else:raise SystemExit('Matrix dimensions missing')
  if nc!=len(barcodes) or nr!=features:raise SystemExit('Features/barcodes disagree with matrix dimensions')
  print('Scanning author matrix once; selecting',len(kept),'cells from',nc,'cells;',nz,'stored entries',flush=True)
  for line in f:
   if not line.strip():continue
   i,j,v=line.split();ii=int(i);jj=int(j)
   if not(1<=ii<=nr and 1<=jj<=nc):raise SystemExit('Invalid matrix coordinate')
   seen+=1
   if jj in columns:
    if not math.isfinite(float(v)):raise SystemExit('Non-finite selected author value')
    body.write(i+' '+str(columns[jj])+' '+v+'\n');kept_nz+=1
   if seen%10000000==0:print('Scanned {:,}/{:,} entries'.format(seen,nz),flush=True)
 if seen!=nz:raise SystemExit('Truncated author matrix')
 after=matrix.stat()
 if (before.st_size,before.st_mtime_ns)!=(after.st_size,after.st_mtime_ns):raise SystemExit('Source matrix changed during extraction')
 with gzip.open(sample/'matrix.mtx.gz','wt',compresslevel=1) as dest,(sample/'.entries.tmp').open() as body:
  dest.write('%%MatrixMarket matrix coordinate real general\n% Author values unchanged; selected GW18 cells only\n'+f'{nr} {len(kept)} {kept_nz}\n')
  shutil.copyfileobj(body,dest)
 (sample/'.entries.tmp').unlink()
 shutil.copyfile(a.directory/'features.tsv.gz',sample/'features.tsv.gz')
 with gzip.open(sample/'barcodes.tsv.gz','wt') as f:f.write('\n'.join(kept)+'\n')
 with (sample/'meta.tsv').open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(selected[cell] for cell in kept)
 (sample/'check_expression.sh').write_text(CLUSTER_SCRIPT,encoding='utf-8')
 info={'seed':a.seed,'strategy':'Round-robin across author cell types within each library; not population-weighted',
  'selection':selection,'cells':len(kept),'features':nr,'stored_entries':kept_nz,
  'author_values_changed':False,'source_matrix_bytes':before.st_size,'source_matrix_mtime_ns':before.st_mtime_ns}
 (sample/'SAMPLE_INFO.json').write_text(json.dumps(info,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
 content=['meta.tsv','matrix.mtx.gz','features.tsv.gz','barcodes.tsv.gz','check_expression.sh','SAMPLE_INFO.json']
 (sample/'EXPRESSION_CHECK_SHA256SUMS').write_text(''.join(sha((sample/n).read_bytes())+'  '+n+'\n' for n in content),encoding='utf-8')
 archive='crouch_author_expression_check.tar.gz'
 with tarfile.open(a.output/archive,'w:gz',compresslevel=1,format=tarfile.USTAR_FORMAT) as tar:
  for name in content+['EXPRESSION_CHECK_SHA256SUMS']:
   p=sample/name;ti=tarfile.TarInfo('author_expression_check/'+name);ti.size=p.stat().st_size;ti.mode=0o644
   with p.open('rb') as f:tar.addfile(ti,f)
 release.mkdir();lines=[];whole=hashlib.sha256();total=0
 with (a.output/archive).open('rb') as f:
  while True:
   b=f.read(2*1024*1024)
   if not b:break
   total+=1;whole.update(b);name=archive+'.part.'+str(total).zfill(5)
   (release/name).write_bytes(b);lines.append(sha(b)+'  '+name+'\n')
 manifest=''.join(lines).encode();(release/'PARTS_SHA256SUMS').write_bytes(manifest)
 digest=whole.hexdigest();size=(a.output/archive).stat().st_size
 (release/'ARCHIVE_SHA256SUMS').write_text(digest+'  '+archive+'\n',encoding='utf-8')
 receiver=relay.receiver_for(archive,digest,sha(manifest),size,total)
 receiver=receiver.replace('author-meta-v1','author-expression-check-v1').replace('crouch-author-meta-transfer','crouch-expression-check-transfer').replace('crouch_author_meta','crouch_author_expression')
 cloud=a.repo/'cloud-transfer';cloud.mkdir(exist_ok=True)
 (cloud/'crouch_author_expression_receive.sh').write_text(receiver,encoding='utf-8')
 installer=relay.BASE_INSTALLER.replace('crouch_cloud_receive.sh','crouch_author_expression_receive.sh').replace('crouch-cloud-transfer','crouch-expression-check-transfer')
 installer,n=re.subn(r'^script_hash=[0-9a-f]{64}$','script_hash='+sha(receiver.encode()),installer,flags=re.M);assert n==1
 (a.output/'Windows_Git_Bash.txt').write_text(installer,encoding='utf-8')
 (a.output/'crouch_author_expression_SHA256SUMS').write_text(digest+'  '+archive+'\n',encoding='utf-8')
 attrs=a.repo/'.gitattributes';old=attrs.read_text() if attrs.exists() else '';rule='author-expression-check-v1/*.part.* -text'
 if rule not in old.splitlines():attrs.write_text(old+('' if not old or old.endswith('\n') else '\n')+rule+'\n')
 analysis=a.repo/'analysis/author-expression';analysis.mkdir(parents=True,exist_ok=True)
 (analysis/'check_expression.sh').write_text(CLUSTER_SCRIPT,encoding='utf-8')
 target=a.repo/'tools/prepare_author_expression_on_mac.py'
 if pathlib.Path(__file__).resolve()!=target.resolve():shutil.copyfile(__file__,target)
 print('READY:',len(kept),'cells;',nr,'features;',kept_nz,'stored entries;',size,'archive bytes;',total,'parts')
 print('Selection:',json.dumps(selection,ensure_ascii=False))
 print('Commit author-expression-check-v1/, new cloud-transfer receiver, analysis/author-expression/, tools/prepare_author_expression_on_mac.py and .gitattributes; Push origin.')
 print('Windows Git Bash installer:',a.output/'Windows_Git_Bash.txt')
if __name__=='__main__':main()
