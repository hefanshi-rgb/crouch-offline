#!/usr/bin/env python3
import csv,gzip,hashlib,io,json,pathlib,re,sys,tarfile
root=pathlib.Path(sys.argv[1]);meta_archive=pathlib.Path(sys.argv[2]);partdir=pathlib.Path(sys.argv[3]);out=root/'inputs';out.mkdir(exist_ok=True)
with tarfile.open(meta_archive,'r:gz') as t:
 rows=list(csv.DictReader(io.TextIOWrapper(t.extractfile('author_meta_all/meta.tsv')),delimiter='\t'))
rows=[r for r in rows if r['cellId'].startswith('all_R7_ctx_endo_')]
assert len(rows)==1759
for r in rows:r['author_cellId']=r['cellId'];r['CB']=r['cellId'].removeprefix('all_R7_ctx_endo_')+'-1';r['library']='R7CE:MissingLibrary:1'
assert len({r['CB'] for r in rows})==1759
with (out/'author_cells.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
checks={}
for line in (partdir/'PARTS_SHA256SUMS').read_text().splitlines():
 h,name=line.split(maxsplit=1);checks[pathlib.Path(name.strip().lstrip('*')).name]=h
parts=sorted(partdir.glob('*.part.*'));assert len(parts)==519
class Parts(io.RawIOBase):
 def __init__(self):self.i=0;self.f=None
 def readable(self):return True
 def readinto(self,b):
  while True:
   if self.f is None:
    if self.i==len(parts):return 0
    p=parts[self.i];self.i+=1
    assert hashlib.sha256(p.read_bytes()).hexdigest()==checks[p.name],str(p)
    self.f=p.open('rb')
   n=self.f.readinto(b)
   if n:return n
   self.f.close();self.f=None
features={};found=False
with io.BufferedReader(Parts(),buffer_size=1024*1024) as inp,tarfile.open(fileobj=inp,mode='r|gz') as t:
 for m in t:
  if m.name.endswith('/genes/genes.gtf'):
   for raw in t.extractfile(m):
    line=raw.decode()
    if line.startswith('#'):continue
    fields=line.rstrip('\n').split('\t');attr=fields[8]
    gid=re.search(r'gene_id "([^"]+)"',attr);gn=re.search(r'gene_name "([^"]+)"',attr)
    if gid:features.setdefault(gid[1],gn[1] if gn else gid[1])
   found=True;break
assert found and len(features)==33694,len(features)
with (out/'features.tsv').open('w') as f:
 for gid,gn in features.items():f.write(gid+'\t'+gn+'\tGene Expression\n')
(out/'input_provenance.json').write_text(json.dumps(dict(reference='GRCh38-1.2.0',reference_archive=json.loads((partdir/'TRANSFER_INFO.json').read_text()),author_cells=1759,all_reference_genes=33694,meta_archive_SHA256=hashlib.sha256(meta_archive.read_bytes()).hexdigest(),author_CB_rule='Exact author all_R7_ctx_endo_ barcode + -1, checked against original BAM RG GEM=1',source_partdir=str(partdir)),indent=2)+'\n')
print('READY author cells',len(rows),'all genes',len(features),flush=True)
