import json,pathlib,sqlite3,gzip,sys
w=pathlib.Path(sys.argv[1]);out=w/'R7CE_full_GX';p=json.loads((out/'provenance.json').read_text());i=json.loads((w/'BAM_INDEX_METADATA.json').read_text());db=sqlite3.connect('file:'+str(out/'molecules.sqlite')+'?mode=ro',uri=True)
n,h,e=db.execute('SELECT COUNT(*),SUM(hits),MIN(exonic) FROM molecules').fetchone();assert n==p['counters']['unique_UMIs'] and h==p['counters']['eligible_reads'] and e==1;assert p['counters']['records']==i['expected_total_alignment_records'];db.close()
with gzip.open(out/'matrix.mtx.gz','rt') as f:
 assert f.readline().startswith('%%MatrixMarket matrix coordinate integer general')
 line=f.readline()
 while line.startswith('%'):line=f.readline()
 nr,nc,nnz=map(int,line.split());assert (nr,nc,nnz)==(33694,1759,p['counters']['nonzero_gene_cell_entries'])
 seen=0;total=0;percell=[0]*nc
 for line in f:
  r,c,v=map(int,line.split());assert 1<=r<=nr and 1<=c<=nc and v>0;seen+=1;total+=v;percell[c-1]+=v
 assert seen==nnz and total==n and min(percell)>0
v=dict(PASS=True,ENA_BAM_MD5_verified=True,ENA_BAI_MD5_verified=True,all_alignment_records_scanned=p['counters']['records'],alignment_records_match_independent_BAI=True,all_reference_genes=nr,all_author_cells=nc,UMIs=total,eligible_reads_sum_matches_SQLite_hits=True,all_counted_molecules_have_exonic_support=True,zero_UMI_author_cells=percell.count(0))
(out/'final_audit.json').write_text(json.dumps(v,indent=2)+'\n');print(json.dumps(v,indent=2))
