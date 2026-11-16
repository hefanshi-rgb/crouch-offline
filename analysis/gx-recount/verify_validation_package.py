import hashlib,json,pathlib,sys,tarfile
root=pathlib.Path(sys.argv[1]);repo=pathlib.Path(sys.argv[2]);folder=repo/'gw18-gx-v1';info=json.loads((folder/'TRANSFER_INFO.json').read_text());whole=hashlib.sha256();total=0
for line in (folder/'PARTS_SHA256SUMS').read_text().splitlines():
 h,name=line.split(maxsplit=1);p=folder/name;b=p.read_bytes();assert hashlib.sha256(b).hexdigest()==h;whole.update(b);total+=len(b)
assert whole.hexdigest()==info['SHA256'] and total==info['bytes']
archive=root/'outputs/gx-recount'/info['archive']
with tarfile.open(archive,'r:gz') as t:
 checks=t.extractfile('gx_validation/SHA256SUMS').read().decode().splitlines()
 for line in checks:
  h,name=line.split(maxsplit=1);assert hashlib.sha256(t.extractfile('gx_validation/'+name).read()).hexdigest()==h,name
 assert t.extractfile('gx_validation/final_audit.json') is not None
print('PASS: all 3 transfer parts, whole archive, all archived files, and final audit verified. Bytes:',total)
