#!/usr/bin/env python3
"""Run full verified-BAM GX counting and comparisons after resumable download."""
import fcntl,hashlib,json,os,pathlib,subprocess,sys,time,traceback
work=pathlib.Path(sys.argv[1]);code=pathlib.Path(sys.argv[2]);author=pathlib.Path(sys.argv[3]);lock=(work/'pipeline.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
python=work/'venv/bin/python';out=work/'R7CE_full_GX';star=work/'inputs/STARsolo_GW18_R7CE_GRCh38_1.2.0.rds'
def state(stage,**extra):
 p=work/'status.json';q=p.with_suffix('.tmp');q.write_text(json.dumps(dict(stage=stage,updated_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),**extra),indent=2)+'\n');os.replace(q,p)
def command(args):
 print('RUN',json.dumps([str(x) for x in args]),flush=True);subprocess.run([str(x) for x in args],check=True)
try:
 source=json.loads((work/'source.json').read_text());start=time.time();initial=sum(p.stat().st_size for p in (work/'download_parts').glob('*.part'))
 while not (work/'BAM_VERIFIED.json').exists():
  ds=json.loads((work/'download_status.json').read_text())
  if ds['stage']=='FAILED':raise RuntimeError('BAM download failed: '+ds.get('error','unknown'))
  n=sum(p.stat().st_size for p in (work/'download_parts').glob('*.part'));rate=(n-initial)/max(.001,time.time()-start)
  state('DOWNLOADING',downloaded_bytes=n,total_bytes=source['bytes'],fraction_complete=round(n/source['bytes'],5),recent_average_MB_per_second=round(rate/1e6,3),estimated_download_seconds_remaining=round((source['bytes']-n)/rate) if rate>0 else None)
  if time.time()-start>86400:raise RuntimeError('Download wait exceeded 24 hours; partial files retained')
  time.sleep(10)
 verified=json.loads((work/'BAM_VERIFIED.json').read_text());bam=pathlib.Path(verified['local_path'])
 assert verified['ENA_MD5_verified'] is True and verified['md5']==source['md5'] and bam.stat().st_size==source['bytes']
 assert (work/'inputs/features.tsv').is_file() and (work/'inputs/author_cells.tsv').is_file()
 state('COUNTING',BAM=str(bam),output=str(out));command([python,code/'count_gx.py',bam,work/'inputs',out])
 prov=json.loads((out/'provenance.json').read_text());assert prov['full_BAM_processed'] is True
 state('EXPORTING_AND_COMPARING',output=str(out));command(['/usr/local/bin/Rscript',code/'export_and_compare.R',out,author,star])
 prov['verified_source']=verified;(out/'source_provenance.json').write_text(json.dumps(verified,indent=2)+'\n')
 state('COMPLETE' if (out/'COMPARISON_COMPLETE').exists() else 'GX_COMPLETE_WAITING_FOR_STARSOLO_INPUT',output=str(out),GX_integer_UMIs=prov['counters']['unique_UMIs'],required_STARsolo_input=str(star))
 print('PIPELINE_FINISHED',flush=True)
except BaseException as e:
 state('FAILED',error=repr(e));traceback.print_exc();raise
