#!/usr/bin/env python3
"""Disk-backed original-BAM UMI reconstruction; no realignment."""
import argparse,collections,csv,gzip,hashlib,json,os,pathlib,re,sqlite3,time
import pysam
METHOD='Primary mapped alignments, MAPQ=255, exactly one GX, complete valid CB/UB, fixed author cohort; unique (physical library, full CB incl GEM, GX, UB). No normalization or extra cell/gene QC.'
def write_json(p,obj):
 q=p.with_suffix(p.suffix+'.tmp');q.write_text(json.dumps(obj,indent=2)+'\n');os.replace(q,p)
def count(bam_path,input_dir,out,max_records=0,fixture=False):
 out.mkdir(parents=True,exist_ok=True)
 if (out/'molecules.sqlite').exists():raise FileExistsError('Refuse to overwrite molecule database: '+str(out))
 cells=list(csv.DictReader((input_dir/'author_cells.tsv').open(),delimiter='\t'))
 features=[x.rstrip('\n').split('\t') for x in (input_dir/'features.tsv').open()]
 cbindex={r['CB']:i for i,r in enumerate(cells)};bare_author={x.split('-')[0] for x in cbindex};geneindex={r[0]:i for i,r in enumerate(features)}
 assert len(cbindex)==len(cells) and len(geneindex)==len(features)
 if not fixture:assert len(cells)==1759 and len(features)==33694
 db=sqlite3.connect(out/'molecules.sqlite');db.execute('PRAGMA journal_mode=WAL');db.execute('PRAGMA synchronous=NORMAL');db.execute('PRAGMA cache_size=-65536');db.execute('CREATE TABLE molecules(cell INTEGER,gene INTEGER,umi TEXT,exonic INTEGER,hits INTEGER,PRIMARY KEY(cell,gene,umi)) WITHOUT ROWID')
 sql='INSERT INTO molecules VALUES(?,?,?,?,1) ON CONFLICT(cell,gene,umi) DO UPDATE SET hits=hits+1,exonic=MAX(exonic,excluded.exonic)'
 stats=collections.Counter();regions=collections.Counter();umis=collections.Counter();rgs=collections.Counter();batch=[];started=time.time();limit_reached=False
 def report(stage):
  elapsed=time.time()-started
  write_json(out/'status.json',dict(stage=stage,updated_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),elapsed_seconds=round(elapsed,2),records_per_second=round(stats['records']/max(elapsed,.001),1),counters=dict(stats),regions_of_eligible_reads=dict(regions),UB_lengths=dict(umis)))
 with pysam.AlignmentFile(str(bam_path),'rb',threads=1 if max_records else 3,ignore_truncation=bool(max_records)) as bam:
  header=bam.header.to_dict();(out/'original_header.sam').write_text(str(bam.header))
  for rg in header.get('RG',[]):
   if not fixture and not rg.get('ID','').startswith('R7CE:MissingLibrary:1:'):raise ValueError('Unexpected physical library/GEM in BAM header '+str(rg))
  for r in bam.fetch(until_eof=True):
   stats['records']+=1
   if stats['records']%500000==0:
    db.executemany(sql,batch);db.commit();batch.clear();report('COUNTING');print('PROGRESS',stats['records'],'eligible',stats['eligible_reads'],'seconds',round(time.time()-started),flush=True)
   if not r.is_secondary and not r.is_supplementary:
    stats['primary_records']+=1
    if r.has_tag('CB'):
     cb=r.get_tag('CB');ci=cbindex.get(cb)
     if ci is None and cb.split('-')[0] in bare_author:raise ValueError('Author barcode with unexpected GEM '+cb)
     if ci is not None:
      stats['author_cell_records']+=1
      if not r.is_unmapped and r.mapping_quality==255:
       stats['author_mapped_MAPQ255']+=1
       gx=r.get_tag('GX') if r.has_tag('GX') else ''
       if gx and len(gx.split(';'))==1:
        stats['author_single_GX_MAPQ255']+=1
        if gx not in geneindex:raise ValueError('GX absent from reference '+gx)
        if r.has_tag('UB'):
         ub=r.get_tag('UB');umis[len(ub)]+=1
         if len(ub)!=10:raise ValueError('Unexpected complete UB length '+str(len(ub)))
         if re.fullmatch('[ACGT]{10}',ub):
          if not re.fullmatch('[ACGT]{16}-1',cb):raise ValueError('Unexpected CB '+cb)
          rg=r.get_tag('RG') if r.has_tag('RG') else ''
          if not fixture and not rg.startswith('R7CE:MissingLibrary:1:'):raise ValueError('Unexpected RG '+rg)
          rt=r.get_tag('RE') if r.has_tag('RE') else 'MISSING';regions[rt]+=1;rgs[rg]+=1
          stats['eligible_reads']+=1;stats['eligible_TX_present']+=int(r.has_tag('TX') and bool(r.get_tag('TX')));stats['eligible_duplicate_flag']+=int(r.is_duplicate);stats['eligible_xf_present']+=int(r.has_tag('xf'));stats['eligible_mm_present']+=int(r.has_tag('mm'))
          batch.append((ci,geneindex[gx],ub,int(rt=='E')))
          if len(batch)>=50000:db.executemany(sql,batch);db.commit();batch.clear()
         else:stats['author_single_GX_UB_contains_N_or_invalid']+=1
        else:stats['author_single_GX_missing_UB']+=1
       elif gx:stats['author_multigene_GX']+=1
       else:stats['author_missing_GX']+=1
    else:stats['primary_missing_CB']+=1
   if max_records and stats['records']>=max_records:limit_reached=True;break
 db.executemany(sql,batch);db.commit();report('EXPORTING')
 ncell=len(cells);ngene=len(features);totals=[0]*ncell;ngenes=[0]*ncell
 rows=db.execute('SELECT gene,cell,COUNT(*) FROM molecules GROUP BY gene,cell').fetchall()
 with gzip.open(out/'matrix.mtx.gz','wt') as f:
  f.write('%%MatrixMarket matrix coordinate integer general\n% BAM-derived CB UB GX\n');f.write('%d %d %d\n'%(ngene,ncell,len(rows)))
  for gi,ci,n in rows:f.write('%d %d %d\n'%(gi+1,ci+1,n));totals[ci]+=n;ngenes[ci]+=1
 with gzip.open(out/'features.tsv.gz','wt') as f:
  for row in features:f.write('\t'.join(row)+'\n')
 with gzip.open(out/'barcodes.tsv.gz','wt') as f:
  for c in cells:f.write(c['author_cellId']+'\n')
 for i,c in enumerate(cells):c['nCount_GX']=totals[i];c['nFeature_GX']=ngenes[i]
 with (out/'cell_metadata.tsv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(cells[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(cells)
 stats['unique_UMIs']=sum(totals);stats['nonzero_gene_cell_entries']=len(rows);stats['cells_with_zero_UMIs']=totals.count(0)
 stats['unique_UMIs_with_exonic_RE_support']=db.execute('SELECT COUNT(*) FROM molecules WHERE exonic=1').fetchone()[0]
 stats['CB_UB_pairs_assigned_to_multiple_GX']=db.execute('SELECT COUNT(*) FROM (SELECT cell,umi FROM molecules GROUP BY cell,umi HAVING COUNT(*)>1)').fetchone()[0]
 stats['reads_per_UMI_max']=db.execute('SELECT MAX(hits) FROM molecules').fetchone()[0]
 db.execute('PRAGMA wal_checkpoint(TRUNCATE)');db.close()
 prov=dict(method=METHOD,sample='GW18',library='R7CE:MissingLibrary:1',reference_label='GRCh38-1.2.0',normalized=False,counts_type='BAM-derived integer UMI reconstruction; exact original Cell Ranger matrix is not established',author_cohort_cells=ncell,all_reference_genes=ngene,source_BAM=str(bam_path),pysam=pysam.__version__,full_BAM_processed=not max_records and not fixture,benchmark_record_limit=max_records,code_SHA256=hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(),eligible_read_RG_counts=dict(rgs),eligible_read_regions=dict(regions),counters=dict(stats),elapsed_seconds=round(time.time()-started,2))
 write_json(out/'provenance.json',prov);report('BENCHMARK_COMPLETE' if max_records else 'COUNT_COMPLETE');print(json.dumps(prov,indent=2),flush=True);return prov
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('bam',type=pathlib.Path);p.add_argument('inputs',type=pathlib.Path);p.add_argument('output',type=pathlib.Path);p.add_argument('--max-records',type=int,default=0);a=p.parse_args();count(a.bam,a.inputs,a.output,a.max_records)
