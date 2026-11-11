#!/bin/bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
sha256sum -c PARTS_SHA256SUMS
cat crouch_refseq_payload.tar.gz.part.* > crouch_refseq_payload.tar.gz.partial
read -r hash name < ARCHIVE_SHA256SUMS
[[ "$name" == crouch_refseq_payload.tar.gz && "$hash" =~ ^[0-9a-fA-F]{64}$ ]]
printf '%s  %s\n' "$hash" crouch_refseq_payload.tar.gz.partial | sha256sum -c -
gzip -t crouch_refseq_payload.tar.gz.partial
mv -f crouch_refseq_payload.tar.gz.partial crouch_refseq_payload.tar.gz
echo 'RESTORE COMPLETE: upload crouch_refseq_payload.tar.gz to Crouch/software on the cluster.'
