import json,pathlib,hashlib,pysam,sys
w=pathlib.Path(sys.argv[1]);b=pathlib.Path(sys.argv[2]);idx=pathlib.Path(str(b)+'.bai')
s=json.loads((w/'ENA_source.json').read_text())[0];urls=s['submitted_ftp'].split(';');name=idx.name;i=next(i for i,x in enumerate(urls) if x.endswith('/'+name));expected=s['submitted_md5'].split(';')[i];actual=hashlib.md5(idx.read_bytes()).hexdigest();assert expected==actual,'BAI MD5 mismatch'
with pysam.AlignmentFile(str(b),'rb',index_filename=str(idx)) as f:
 stats=[dict(contig=x.contig,mapped=x.mapped,unmapped=x.unmapped) for x in f.get_index_statistics()];nocoordinate=f.nocoordinate
v=dict(BAI=str(idx),ENA_MD5_verified=True,MD5=actual,contigs=stats,no_coordinate_records=nocoordinate,expected_total_alignment_records=sum(x['mapped']+x['unmapped'] for x in stats)+nocoordinate)
(w/'BAM_INDEX_METADATA.json').write_text(json.dumps(v,indent=2)+'\n');print('BAI_MD5_PASS expected alignment records:',v['expected_total_alignment_records']);print('no-coordinate records:',nocoordinate)
