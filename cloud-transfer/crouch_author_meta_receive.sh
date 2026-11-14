#!/bin/bash
# Windows Git Bash: apply the approved internal GitHub route on every Git call.
# Does not depend on hosts/DNS cache edits or require administrator privileges.
set -euo pipefail
set -m
umask 077
root=${CROUCH_CLOUD_ROOT:-$HOME/Desktop/crouch-author-meta-transfer}
repo="$root/repository"
files="$root/download"
url=${CROUCH_REPO_URL:-https://github.com/hefanshi-rgb/crouch-offline.git}
mirror=${CROUCH_MIRROR_IP:-172.30.100.20}
remote_dir=author-meta-v1
archive=crouch_author_meta.tar.gz
archive_hash=b9bb55f9539a49806997c4bd27652b6c41a2a729ef12d252f4f4fbb02e3aa922
manifest_hash=57a864f87aaeabd8486246fa2156cf29a5e8f07a48bc16ccfe4cff2f48f7403e
archive_bytes=4909962
max_attempts=${CROUCH_MAX_ATTEMPTS:-20}
attempt_seconds=${CROUCH_ATTEMPT_SECONDS:-600}
no_progress_seconds=${CROUCH_NO_PROGRESS_SECONDS:-1800}
retry_base=${CROUCH_RETRY_BASE_SECONDS:-10}
retry_max=${CROUCH_RETRY_MAX_SECONDS:-60}
poll_seconds=${CROUCH_POLL_SECONDS:-1}
for value in "$max_attempts" "$attempt_seconds" "$no_progress_seconds" "$retry_base" "$retry_max" "$poll_seconds"; do
    [[ "$value" =~ ^[1-9][0-9]*$ && ${#value} -le 6 ]] || { printf 'Invalid retry setting: %s\n' "$value" >&2; exit 1; }
done
[[ "$mirror" =~ ^([0-9]{1,3}\.){3}[0-9]{1,3}$ ]] || { printf 'Invalid mirror IPv4 address\n' >&2; exit 1; }
mkdir -p "$files"
export GIT_TERMINAL_PROMPT=0 GCM_INTERACTIVE=Never
unset GIT_SSL_NO_VERIFY
git_options=(-c core.autocrlf=false -c protocol.version=2 -c http.version=HTTP/1.1
    -c http.https://github.com/.proxy= -c http.https://github.com/.sslVerify=true
    -c http.https://github.com/.curloptResolve=
    -c "http.https://github.com/.curloptResolve=github.com:443:$mirror"
    -c http.lowSpeedLimit=512 -c http.lowSpeedTime=60)
case "$(uname -s)" in
    MINGW*|MSYS*|CYGWIN*) git_options+=(-c http.https://github.com/.sslBackend=schannel);;
esac
g() { git "${git_options[@]}" "$@"; }
log() { printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" | tee -a "$root/transfer.log"; }
die() { log "STOP: $*" >&2; exit 1; }
active_pid=
lock_owned=false
lock="$root/.running"
last_progress=$(date +%s)
network_shown=0
stop_active() {
    local n
    if [[ -n "$active_pid" ]]; then
        kill -TERM -- "-$active_pid" 2>/dev/null || kill -TERM "$active_pid" 2>/dev/null || true
        for n in 1 2 3; do kill -0 -- "-$active_pid" 2>/dev/null || break; sleep 1; done
        kill -KILL -- "-$active_pid" 2>/dev/null || true
        wait "$active_pid" 2>/dev/null || true
        active_pid=
    fi
}
cleanup() {
    stop_active
    if [[ "$lock_owned" == true ]]; then rm -f "$lock/pid"; rmdir "$lock" 2>/dev/null || true; fi
}
trap cleanup EXIT
trap 'log "Stopped by Ctrl+C; verified parts are retained."; exit 130' INT
trap 'log "Stopped by signal; verified parts are retained."; exit 143' TERM
guard() {
    [[ ! -e "$root/STOP" ]] || die "STOP file found: $root/STOP. Remove it before resuming."
    (( $(date +%s) - last_progress < no_progress_seconds )) || die "No verified progress for ${no_progress_seconds}s. Inspect transfer.log and last-network.log."
}
pause_retry() {
    local seconds=$1 n
    for ((n=0; n<seconds; n++)); do guard; sleep 1; done
}
flush_network_output() {
    local bytes delta
    bytes=$(wc -c < "$root/last-network.log")
    delta=$((bytes - network_shown))
    if (( delta > 0 )); then
        dd if="$root/last-network.log" bs=1 skip="$network_shown" count="$delta" 2>/dev/null | cat >&2
        network_shown=$bytes
    fi
}
run_network() {
    local input=$1 started status
    shift
    guard
    started=$(date +%s)
    : > "$root/last-network.log"
    network_shown=0
    git "${git_options[@]}" "$@" < "$input" > "$root/last-network.log" 2>&1 &
    active_pid=$!
    while kill -0 "$active_pid" 2>/dev/null; do
        flush_network_output
        guard
        if (( $(date +%s) - started >= attempt_seconds )); then
            stop_active
            flush_network_output
            log "Attempt timed out after ${attempt_seconds}s."
            return 124
        fi
        sleep "$poll_seconds"
    done
    if wait "$active_pid"; then status=0; else status=$?; fi
    active_pid=
    flush_network_output
    if (( status != 0 )) && grep -Eiq 'Authentication failed|could not read Username|Repository not found|returned error: (401|403|407)|not our ref|does not allow request|filtering not recognized' "$root/last-network.log"; then
        die 'Repository, permission, or server-capability error; see last-network.log. This requires correction rather than endless retries.'
    fi
    return "$status"
}
retry_network() {
    local label=$1 input=$2 n delay=$retry_base
    shift 2
    for ((n=1; n<=max_attempts; n++)); do
        guard
        log "$label, attempt $n/$max_attempts"
        if run_network "$input" "$@"; then return 0; fi
        (( n < max_attempts )) || die "$label failed $max_attempts times. Verified downloads are retained."
        log "Retrying in ${delay}s."
        pause_retry "$delay"
        delay=$((delay * 2))
        if (( delay > retry_max )); then delay=$retry_max; fi
    done
}
verified() { [[ -s "$2" ]] && printf '%s  %s\n' "$1" "$2" | sha256sum -c - >/dev/null 2>&1; }
publish_desktop() {
    local desktop="$HOME/Desktop/$archive"
    if [[ ! -e "$desktop" ]]; then
        local temp
        temp=$(mktemp "$HOME/Desktop/.crouch-author-meta.XXXXXX")
        cp "$files/$archive" "$temp"
        mv "$temp" "$desktop"
    elif ! verified "$archive_hash" "$desktop"; then
        log "An existing different Desktop file was preserved. Use $files/$archive"
        return
    fi
    printf '%s  %s\n' "$archive_hash" "$archive" > "$HOME/Desktop/crouch_author_meta_SHA256SUMS"
    log "READY ON DESKTOP: $desktop"
}
install_shortcut() {
    case "$(uname -s)" in
        MINGW*|MSYS*|CYGWIN*)
            local bash_path script_path
            bash_path=$(cygpath -m "$BASH")
            script_path=$(cygpath -m "$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)/$(basename -- "${BASH_SOURCE[0]}")")
            # Percent signs in a Windows path must be doubled in a batch file.
            bash_path=${bash_path//%/%%}; script_path=${script_path//%/%%}
            printf '@echo off\r\n"%s" --login "%s"\r\nset "CROUCH_EXIT_CODE=%%ERRORLEVEL%%"\r\necho.\r\npause\r\nexit /b %%CROUCH_EXIT_CODE%%\r\n' \
                "$bash_path" "$script_path" > "$HOME/Desktop/crouch_author_meta_receive.cmd"
            ;;
    esac
}

guard
if ! mkdir "$lock" 2>/dev/null; then
    old_pid=$(cat "$lock/pid" 2>/dev/null || true)
    if [[ "$old_pid" =~ ^[0-9]+$ ]] && kill -0 "$old_pid" 2>/dev/null; then die "Another receiver is running: PID $old_pid"; fi
    rm -f "$lock/pid"
    rmdir "$lock" || die "Inspect stale lock: $lock"
    mkdir "$lock"
fi
lock_owned=true
printf '%s\n' "$$" > "$lock/pid"
install_shortcut
log "GitHub route configured for this run: github.com:443 -> $mirror; HTTPS verification enabled."
log "Stop: Ctrl+C, or create $root/STOP. Retry limit $max_attempts; no-progress limit ${no_progress_seconds}s."
if verified "$archive_hash" "$files/$archive"; then
    gzip -t "$files/$archive"
    publish_desktop
    log 'Already downloaded and verified; no network transfer needed.'
    exit 0
fi
if [[ ! -d "$repo/.git" ]]; then
    [[ ! -e "$repo" ]] || die "Inspect existing non-repository path: $repo"
    git init -q "$repo"
    g -C "$repo" remote add origin "$url"
fi
[[ "$(g -C "$repo" remote get-url origin)" == "$url" ]] || die 'Existing repository has a different origin; preserving it.'
g -C "$repo" config remote.origin.promisor true
g -C "$repo" config remote.origin.partialclonefilter blob:none
retry_network 'Refreshing repository metadata only' /dev/null -C "$repo" fetch --progress --depth=1 --no-tags --filter=blob:none origin main
snapshot=$(g -C "$repo" rev-parse --verify FETCH_HEAD)
log "Pinned snapshot: $snapshot"

get_file() {
    local name=$1 destination=$2 expected=$3 oid
    oid=$(g -C "$repo" rev-parse --verify "$snapshot:$remote_dir/$name") || die "Release path missing: $remote_dir/$name. The internal mirror may not have synchronized yet."
    printf '%s\n' "$oid" > "$root/object-want.txt"
    retry_network "Fetching $name" "$root/object-want.txt" -C "$repo" -c fetch.negotiationAlgorithm=noop fetch --progress --no-tags \
        --no-write-fetch-head --recurse-submodules=no --filter=blob:none --stdin origin
    g -C "$repo" cat-file blob "$oid" > "$destination"
    verified "$expected" "$destination" || die "Checksum mismatch for $name; preserving .partial for inspection."
}
get_file PARTS_SHA256SUMS "$files/PARTS_SHA256SUMS.partial" "$manifest_hash"
mv -f "$files/PARTS_SHA256SUMS.partial" "$files/PARTS_SHA256SUMS"
total=0
while read -r hash name extra; do
    [[ "$hash" =~ ^[0-9a-f]{64}$ && "$name" =~ crouch_author_meta\.tar\.gz\.part\.[0-9]{5}$ && -z "${extra:-}" ]] || die 'Invalid part manifest'
    total=$((total + 1))
done < "$files/PARTS_SHA256SUMS"
[[ "$total" == 3 ]] || die 'Unexpected metadata part count'
i=0
while read -r hash name; do
    guard
    i=$((i + 1))
    if verified "$hash" "$files/$name"; then
        log "[$i/$total] SKIP verified $name"
    else
        log "[$i/$total] DOWNLOAD $name"
        get_file "$name" "$files/$name.partial" "$hash"
        (( $(wc -c < "$files/$name.partial") <= 2097152 )) || die 'Part exceeds 2 MiB'
        mv -f "$files/$name.partial" "$files/$name"
        log "[$i/$total] VERIFIED $name"
    fi
    last_progress=$(date +%s)
done < "$files/PARTS_SHA256SUMS"
(
    while read -r hash name; do guard; cat "$files/$name"; done < "$files/PARTS_SHA256SUMS"
) > "$files/$archive.partial"
verified "$archive_hash" "$files/$archive.partial" || die 'Reassembled archive checksum mismatch'
[[ "$(wc -c < "$files/$archive.partial" | tr -d ' ')" == "$archive_bytes" ]] || die 'Archive size mismatch'
gzip -t "$files/$archive.partial"
mv -f "$files/$archive.partial" "$files/$archive"
printf '%s  %s\n' "$archive_hash" "$archive" > "$files/ARCHIVE_SHA256SUMS"
publish_desktop
log "DOWNLOAD COMPLETE: $files/$archive"
printf '\nUpload the archive AND crouch_author_meta_SHA256SUMS from Desktop to the cluster Crouch/ directory.\n'
