#!/bin/bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
sha256sum -c PARTS_SHA256SUMS
cat crouch_refseq_payload.tar.gz.part.* > crouch_refseq_payload.tar.gz.partial
gzip -t crouch_refseq_payload.tar.gz.partial
mv -f crouch_refseq_payload.tar.gz.partial crouch_refseq_payload.tar.gz
echo 'RESTORE COMPLETE: upload crouch_refseq_payload.tar.gz to Crouch/software on the cluster.'
