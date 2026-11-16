#!/bin/bash
# Run manually in Windows Git Bash. Pushes only a harmless test on a NEW branch.
set -euo pipefail
umask 077
parent=${CROUCH_RETURN_PARENT:-$HOME/Desktop}
root=$(mktemp -d "$parent/crouch-return-test.XXXXXX")
branch="cloud-return-test-$(date -u +%Y%m%dT%H%M%SZ)-$(basename "$root" | cut -d. -f2)"
url=${CROUCH_TEST_REMOTE:-https://github.com/hefanshi-rgb/crouch-offline.git}
git init -q "$root"
git -C "$root" config core.autocrlf false
git -C "$root" checkout -q -b "$branch"
printf 'Crouch Windows-to-Mac return test\nUTC: %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" > "$root/cloud-return-test.txt"
(cd "$root" && sha256sum cloud-return-test.txt > SHA256SUMS)
git -C "$root" add -- cloud-return-test.txt SHA256SUMS
git -C "$root" -c user.name='Cloud return test' -c user.email='cloud-return-test@users.noreply.github.com' commit -q -m 'Test Windows file return on independent branch'
git -C "$root" remote add origin "$url"
opts=(-c credential.helper= -c credential.interactive=false -c http.version=HTTP/1.1 -c http.lowSpeedLimit=1 -c http.lowSpeedTime=30)
if [[ -z "${CROUCH_TEST_REMOTE:-}" ]]; then
 opts+=(-c http.https://github.com/.curloptResolve= -c http.https://github.com/.curloptResolve=github.com:443:172.30.100.20 -c http.https://github.com/.proxy= -c http.https://github.com/.sslVerify=true)
 case "$(uname -s)" in MINGW*|MSYS*|CYGWIN*)opts+=(-c http.https://github.com/.sslBackend=schannel);;esac
fi
export GIT_TERMINAL_PROMPT=0 GCM_INTERACTIVE=Never
child=;guard=
cleanup(){ [[ -z "$guard" ]] || kill "$guard" 2>/dev/null || true; [[ -z "$child" ]] || kill "$child" 2>/dev/null || true; }
trap cleanup EXIT
printf 'TEST_BRANCH=%s\n' "$branch"
git "${opts[@]}" -C "$root" push --porcelain origin "HEAD:refs/heads/$branch" > "$root/push.log" 2>&1 &
child=$!
(sleep 60; kill -TERM "$child" 2>/dev/null || true) &
guard=$!
if wait "$child"; then rc=0; else rc=$?; fi
child=
kill "$guard" 2>/dev/null || true;wait "$guard" 2>/dev/null || true;guard=
cat "$root/push.log"
if [[ "$rc" == 0 ]]; then
 printf '\nRETURN_TEST_PUSH_SUCCEEDED\nBranch: %s\nSHA256: ' "$branch";cat "$root/SHA256SUMS"
else
 printf '\nRETURN_TEST_PUSH_FAILED (exit %s). No credentials requested. This can mean write authentication is required, the mirror rejects writes, or a network failure.\nLog: %s/push.log\n' "$rc" "$root"
fi
exit "$rc"
