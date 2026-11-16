#!/usr/bin/env python3
import importlib.util,json,pathlib,tempfile
import pysam
spec=importlib.util.spec_from_file_location('gx',pathlib.Path(__file__).with_name('count_gx.py'));gx=importlib.util.module_from_spec(spec);spec.loader.exec_module(gx)
root=pathlib.Path(tempfile.mkdtemp(prefix='gx_scientific_test_'));inp=root/'inputs';inp.mkdir()
(inp/'author_cells.tsv').write_text('author_cellId\tCB\tlibrary\ncellA\tAAAAAAAAAAAAAAAA-1\tL1\ncellZero\tCCCCCCCCCCCCCCCC-1\tL1\n');(inp/'features.tsv').write_text('G1\tOne\tGene Expression\nG2\tTwo\tGene Expression\nGzero\tZero\tGene Expression\n')
header={'HD':{'VN':'1.6'},'SQ':[{'SN':'chr1','LN':10000}]}
def record(i,cb='AAAAAAAAAAAAAAAA-1',ub='AAAAAAAAAA',gene='G1',flag=0,mapq=255,region='E'):
 r=pysam.AlignedSegment();r.query_name='r'+str(i);r.query_sequence='A'*30;r.query_qualities=pysam.qualitystring_to_array('I'*30);r.flag=flag;r.reference_id=0;r.reference_start=10;r.mapping_quality=mapq;r.cigarstring='30M'
 for k,v in [('CB',cb),('UB',ub),('GX',gene),('RE',region)]:
  if v is not None:r.set_tag(k,v)
 return r
bam=root/'test.bam'
with pysam.AlignmentFile(str(bam),'wb',header=header) as f:
 for i,kwargs in enumerate([{}, {'flag':1024}, {'ub':'CCCCCCCCCC'}, {'gene':'G2','region':None}, {'flag':256}, {'flag':2048}, {'mapq':3}, {'gene':'G1;G2'}, {'cb':None}, {'gene':None}, {'ub':None}, {'ub':'NNNNNNNNNN'}, {'cb':'GGGGGGGGGGGGGGGG-1'}, {'flag':4}]):f.write(record(i,**kwargs))
p=gx.count(bam,inp,root/'out',fixture=True)
assert p['counters']['unique_UMIs']==3
assert p['counters']['cells_with_zero_UMIs']==1
assert p['counters']['unique_UMIs_with_exonic_RE_support']==2
assert p['counters']['eligible_reads']==4
assert p['counters']['CB_UB_pairs_assigned_to_multiple_GX']==1
assert p['counters']['eligible_duplicate_flag']==1
import gzip
text=gzip.open(root/'out/matrix.mtx.gz','rt').read();assert text.startswith('%%MatrixMarket');assert '3 2 2\n' in text;assert '1 1 2\n' in text and '2 1 1\n' in text
for label,kw in [('wrong_UMI_length',{'ub':'A'*12}),('wrong_GEM',{'cb':'AAAAAAAAAAAAAAAA-2'}),('unknown_gene',{'gene':'Gunknown'})]:
 bad=root/(label+'.bam')
 with pysam.AlignmentFile(str(bad),'wb',header=header) as f:f.write(record(100,**kw))
 try:gx.count(bad,inp,root/label,fixture=True)
 except ValueError:pass
 else:raise AssertionError('Expected failure '+label)
print('PASS: PCR duplicates deduplicated; distinct UMIs retained; gene identity preserved; low-MAPQ/multigene/secondary/supplementary/unmapped/missing tags excluded; zero genes and cells retained; GEM, UMI-length and feature mismatches stop the run.');print('fixture',root)
