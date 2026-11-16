#!/usr/bin/env python3
import hashlib,json,pathlib,sys,time,os
import pysam
work=pathlib.Path(sys.argv[1]);bam=pathlib.Path(sys.argv[2]);src=json.loads((work/'source.json').read_text())
assert bam.stat().st_size==src['bytes'],'BAM size differs from ENA'
md5=hashlib.md5();sha=hashlib.sha256();t=time.time()
with bam.open('rb') as f:
 while data:=f.read(8*1024*1024):md5.update(data);sha.update(data)
assert md5.hexdigest()==src['md5'],'Full ENA MD5 mismatch'
pysam.quickcheck(str(bam))
with pysam.AlignmentFile(str(bam),'rb') as f:
 h=f.header.to_dict();rgs=h.get('RG',[])
 assert rgs and all(x['ID'].startswith('R7CE:MissingLibrary:1:') for x in rgs),'Wrong physical library/GEM'
 assert len(h['SQ'])==194,'Unexpected contig dictionary'
assert (work/'inputs/features.tsv').read_text().count('\n')==33694
assert (work/'inputs/author_cells.tsv').read_text().count('\n')==1760
out=dict(**src,local_path=str(bam),SHA256=sha.hexdigest(),ENA_MD5_verified=True,BAM_quickcheck_pass=True,original_download_reused=True,verified_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),verification_elapsed_seconds=round(time.time()-t,2))
p=work/'BAM_VERIFIED.json';tmp=p.with_suffix('.tmp');tmp.write_text(json.dumps(out,indent=2)+'\n');os.replace(tmp,p)
print(json.dumps(out,indent=2),flush=True)
