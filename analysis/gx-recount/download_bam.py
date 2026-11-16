#!/usr/bin/env python3
"""Resumable bounded-memory ENA BAM download with strict Range and ENA MD5 checks."""
import concurrent.futures,hashlib,json,os,pathlib,re,sys,threading,time,urllib.request
root=pathlib.Path(sys.argv[1]);src=json.loads((root/'source.json').read_text());parts=root/'download_parts';parts.mkdir(exist_ok=True)
lock=threading.Lock();started=time.time();partsize=32*1024*1024;total=src['bytes']
def atomic(name,value):
 p=root/name;t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(value,indent=2)+'\n');os.replace(t,p)
def progress(stage,**extra):
 with lock:
  n=sum(min(p.stat().st_size,max(0,total-int(p.stem))) for p in parts.glob('*.part'))
  v=dict(stage=stage,updated_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),downloaded_bytes=n,total_bytes=total,elapsed_seconds=round(time.time()-started),**extra);atomic('download_status.json',v)
def part(start):
 end=min(total,start+partsize)-1;size=end-start+1;p=parts/('%012d.part'%start)
 for attempt in range(20):
  have=p.stat().st_size if p.exists() else 0
  if have==size:return
  if have>size:raise ValueError('Oversize partial '+str(p))
  try:
   req=urllib.request.Request(src['url'],headers={'Range':'bytes=%d-%d'%(start+have,end),'User-Agent':'Crouch-GX-audit/1.0','Accept-Encoding':'identity'})
   with urllib.request.urlopen(req,timeout=60) as r:
    cr=r.headers.get('Content-Range','');want='bytes %d-%d/%d'%(start+have,end,total)
    if r.status!=206 or cr!=want:raise ValueError('Invalid Range response '+str(r.status)+' '+cr)
    with p.open('ab') as f:
     while True:
      b=r.read(min(1024*1024,size-have)) if have<size else b''
      if not b:break
      f.write(b);have+=len(b)
   if have!=size:raise IOError('Incomplete segment '+str(have))
   print('PART_DONE',start,size,flush=True);progress('DOWNLOADING');return
  except Exception as e:
   print('RETRY',start,attempt+1,type(e).__name__,str(e),flush=True);time.sleep(min(2**attempt,30))
 raise RuntimeError('Segment failed after retries: '+str(start))
try:
 progress('DOWNLOADING')
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:list(ex.map(part,range(0,total,partsize)))
 progress('ASSEMBLING_AND_MD5')
 target=root/'R7ctxendo.original.bam';tmp=root/'R7ctxendo.original.bam.assembling';md5=hashlib.md5();sha=hashlib.sha256()
 with tmp.open('wb') as out:
  for start in range(0,total,partsize):
   p=parts/('%012d.part'%start)
   if p.stat().st_size!=min(partsize,total-start):raise ValueError('Part size')
   with p.open('rb') as inp:
    while b:=inp.read(4*1024*1024):out.write(b);md5.update(b);sha.update(b)
 if tmp.stat().st_size!=total or md5.hexdigest()!=src['md5']:raise ValueError('ENA MD5 mismatch')
 os.replace(tmp,target);atomic('BAM_VERIFIED.json',dict(**src,local_path=str(target),SHA256=sha.hexdigest(),ENA_MD5_verified=True))
 progress('BAM_VERIFIED',md5=md5.hexdigest(),SHA256=sha.hexdigest());print('BAM_VERIFIED',str(target),flush=True)
except Exception as e:
 progress('FAILED',error=repr(e));raise
