#!/bin/bash
# Paste this complete block into Windows Git Bash once.
(
set -euo pipefail
export GIT_TERMINAL_PROMPT=0 GCM_INTERACTIVE=Never
unset GIT_SSL_NO_VERIFY
receiver="$HOME/Desktop/crouch_cloud_receive.sh"
script_hash=8f3597307326a6aa2a4a85046e02d38e375535d4650afaa19e3e162c2538947e
net=(-c http.https://github.com/.curloptResolve=
     -c http.https://github.com/.curloptResolve=github.com:443:172.30.100.20
     -c http.https://github.com/.proxy=
     -c http.https://github.com/.sslBackend=schannel
     -c http.https://github.com/.sslVerify=true
     -c http.version=HTTP/1.1 -c protocol.version=2
     -c http.lowSpeedLimit=512 -c http.lowSpeedTime=60)
g() { git "${net[@]}" "$@"; }
repo="$HOME/Desktop/crouch-cloud-transfer/repository"
if [[ -s "$receiver" ]] && printf '%s  %s\n' "$script_hash" "$receiver" | sha256sum -c -; then
    bash "$receiver"
    exit
fi
mkdir -p "$HOME/Desktop/crouch-cloud-transfer"
if [[ ! -d "$repo/.git" ]]; then
    [[ ! -e "$repo" ]] || { printf 'Inspect existing path: %s\n' "$repo" >&2; exit 1; }
    git init -q "$repo"
    g -C "$repo" remote add origin https://github.com/hefanshi-rgb/crouch-offline.git
fi
[[ "$(g -C "$repo" remote get-url origin)" == https://github.com/hefanshi-rgb/crouch-offline.git ]] || exit 1
g -C "$repo" config remote.origin.promisor true
g -C "$repo" config remote.origin.partialclonefilter blob:none
for attempt in {1..10}; do
    printf 'Installing cloud receiver: attempt %s/10\n' "$attempt"
    if g -C "$repo" fetch --progress --depth=1 --no-tags --filter=blob:none origin main &&
       g -C "$repo" show FETCH_HEAD:cloud-transfer/crouch_cloud_receive.sh > "$receiver.partial" &&
       printf '%s  %s\n' "$script_hash" "$receiver.partial" | sha256sum -c -; then
        mv -f "$receiver.partial" "$receiver"
        bash "$receiver"
        exit
    fi
    sleep 10
done
printf 'Installation stopped after 10 attempts; existing downloads were retained.\n' >&2
exit 1
)
