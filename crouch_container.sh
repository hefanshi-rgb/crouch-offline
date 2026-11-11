#!/bin/bash
# SSH usage: bash crouch_container.sh [gw18|remaining] [SCRATCH_ROOT]
set -euo pipefail
umask 077
die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }

load_tools() {
    if ! type module >/dev/null 2>&1; then
        local init
        for init in /etc/profile.d/modules.sh /usr/share/Modules/init/bash /usr/share/lmod/lmod/init/bash; do
            if [[ -r "$init" ]]; then
                set +u
                source "$init"
                set -u
                break
            fi
        done
    fi
    type module >/dev/null 2>&1 || die "Environment Modules initialization was not found."
    set +u
    module load system/singularity/3.6.3 system/anaconda/3/2024.06 system/R/4.3.1 ||
        die "Could not load the Singularity, Anaconda and R modules."
    set -u
    command -v singularity >/dev/null || die "Singularity is unavailable after module load."
    python3 -c 'import sys; assert sys.version_info >= (3, 6); print(sys.version)'
    Rscript --vanilla -e 'stopifnot(requireNamespace("Matrix", quietly = TRUE)); print(packageVersion("Matrix"))'
}

work=
cleanup() {
    local status=$? name
    if (( status != 0 )) && [[ -n "$work" ]]; then
        for name in _log _errors; do
            if [[ -f "$work/compat_test/$name" ]]; then
                cp -- "$work/compat_test/$name" "$base/logs/container-${SLURM_JOB_ID}${name}.txt" || true
            fi
        done
    fi
    if [[ -n "$work" ]]; then rm -rf -- "$work" || true; fi
}
trap cleanup EXIT
trap 'exit 130' INT
trap 'exit 143' TERM

run_container() {
    local python_root
    python_root=$(dirname -- "$(dirname -- "$(readlink -f -- "$(command -v python3)")")")
    mkdir -p -- "$base/software/container_home"
    SINGULARITYENV_PATH="$python_root/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin" \
    SINGULARITYENV_SLURM_JOB_ID="$SLURM_JOB_ID" \
    SINGULARITYENV_SLURM_CPUS_PER_TASK="${SLURM_CPUS_PER_TASK:-1}" \
    SINGULARITYENV_SLURM_MEM_PER_NODE="${SLURM_MEM_PER_NODE:-2048}" \
    SINGULARITYENV_SLURM_ARRAY_TASK_ID="${SLURM_ARRAY_TASK_ID:-}" \
    SINGULARITYENV_CREATE_BAM="${CREATE_BAM:-false}" \
    SINGULARITYENV_PYTHONNOUSERSITE=1 \
    singularity exec --cleanenv \
        --home "$base/software/container_home:/home/crouch" --pwd "$base" \
        --bind "$base:$base" --bind "$scratch:$scratch" --bind "$python_root:$python_root" \
        "$image" "$@"
}

if [[ "${1:-}" == --worker ]]; then
    : "${SLURM_JOB_ID:?The worker must be submitted with sbatch.}"
    stage=${2:?Missing worker stage}
    base=${3:?Missing project directory}
    scratch=${4:?Missing scratch directory}
    scope=${5:?Missing sample scope}
    pipeline="$base/crouch_slurm.sh"
    image="$base/software/ubuntu-22.04.sif"
    [[ "$scratch" == /* && -d "$scratch" && -w "$scratch" ]] || die "Scratch directory is unavailable: $scratch"
    printf 'Stage=%s Scope=%s Job=%s Node=%s\n' "$stage" "$scope" "$SLURM_JOB_ID" "$(hostname)"
    load_tools
    case "$stage" in
        setup)
            work=$(mktemp -d "$scratch/crouch-container.${SLURM_JOB_ID}.XXXXXX")
            export SINGULARITY_CACHEDIR="$work/cache"
            export SINGULARITY_TMPDIR="$work/tmp"
            mkdir -p -- "$SINGULARITY_CACHEDIR" "$SINGULARITY_TMPDIR"
            [[ -s "$image" ]] ||
                die "Offline cluster: upload a Linux amd64 Ubuntu image to $image before submitting."
            sha256sum -- "$image" > "$base/software/container_image.sha256"
            run_container /bin/bash -c 'getconf GNU_LIBC_VERSION; python3 --version'
            run_container /bin/bash "$pipeline" prepare "$base"
            cr=$(readlink -f -- "$base/software/crouch_runtime/cellranger")
            run_container "$cr" sitecheck > "$base/logs/container-${SLURM_JOB_ID}-sitecheck.txt"
            run_container /bin/bash -c 'cd -- "$1" && exec "$2" testrun --id=compat_test --jobmode=local --localcores=8 --localmem=48 --disable-ui' bash "$work" "$cr"
            [[ -s "$work/compat_test/outs/filtered_feature_bc_matrix.h5" ]] ||
                die "The bundled-data test did not produce its expression matrix."
            printf 'Container and Cell Ranger test completed successfully.\n'
            ;;
        count)
            [[ -s "$image" ]] || die "The container image is missing."
            run_container /bin/bash "$pipeline" count "$base" "$scratch"
            ;;
        qc)
            if [[ "$scope" == gw18 ]]; then bash "$pipeline" qc "$base" GW18
            else bash "$pipeline" qc "$base"; fi
            ;;
        *) die "Unknown worker stage: $stage" ;;
    esac
    exit 0
fi

# Reuse a complete installation or recognize standard names of the same versions.
# Alias creation lets the existing main script keep its original input names.
archive_ready() {
    local expected=$1
    shift
    if [[ -s "$software/$expected" ]]; then return 0; fi
    if [[ -e "$software/$expected" || -L "$software/$expected" ]]; then
        printf 'Empty file or broken link: %s\n' "$software/$expected" >&2
        return 1
    fi
    local candidate
    local matches=()
    for candidate in "$@"; do
        if [[ -s "$software/$candidate" ]]; then matches+=("$candidate"); fi
    done
    if (( ${#matches[@]} == 1 )); then
        ln -s -- "${matches[0]}" "$software/$expected"
        printf 'Recognized archive: %s -> %s\n' "$expected" "${matches[0]}"
        return 0
    fi
    printf 'Missing required archive: %s\n' "$software/$expected" >&2
    if (( ${#matches[@]} > 1 )); then
        printf 'Multiple candidate files; choose the intended package:\n' >&2
        printf '  %s\n' "${matches[@]}" >&2
    fi
    return 1
}

check_project_inputs() {
    local runtime="$software/crouch_runtime"
    if [[ -e "$runtime" || -L "$runtime" ]]; then
        local missing=0 required
        for required in sra-bin/fasterq-dump sra-bin/vdb-validate cellranger; do
            if [[ ! -x "$runtime/$required" ]]; then
                printf 'Missing executable in installed runtime: %s\n' "$runtime/$required" >&2
                missing=1
            fi
        done
        if [[ ! -s "$runtime/reference/reference.json" || ! -d "$runtime/reference/star" ]]; then
            printf 'Incomplete installed reference: %s/reference\n' "$runtime" >&2
            missing=1
        fi
        (( missing == 0 )) || die "The installed runtime is incomplete. Preserve it and inspect it before retrying."
        printf 'Reusing installed software and reference: %s\n' "$runtime"
        return
    fi
    [[ -d "$software" ]] || die "Software directory is missing: $software"
    local missing=0
    archive_ready sratoolkit-3.4.1.tar.gz \
        sratoolkit.3.4.1-ubuntu64.tar.gz sratoolkit.3.4.1-almalinux64.tar.gz || missing=1
    archive_ready cellranger-10.0.0.tar cellranger-10.0.0.tar.gz cellranger-10.0.0.tar.xz || missing=1
    archive_ready refdata-gex-GRCh38-2024-A.tar.gz \
        refdata-gex-GRCh38-2024-A.tar refdata-gex-GRCh38-2024-A.tar.xz || missing=1
    if (( missing != 0 )); then
        printf '\nActual contents of %s:\n' "$software" >&2
        ls -lah -- "$software" >&2
        die "Restore the files listed above, or use a complete software/crouch_runtime installation."
    fi
}

(( $# <= 2 )) || die "Usage: bash crouch_container.sh [gw18|remaining] [SCRATCH_ROOT]"
scope=${1:-gw18}
[[ "$scope" == gw18 || "$scope" == remaining ]] || die "Choose gw18 or remaining."
scratch=${2:-/public4/usersdir/guojunping/scratch}
[[ "$scratch" == /* ]] || die "Provide an absolute scratch directory."
base=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
pipeline="$base/crouch_slurm.sh"
software="$base/software"
[[ -r "$pipeline" ]] || die "Place the current crouch_slurm.sh beside this launcher."
check_project_inputs
[[ -s "$software/ubuntu-22.04.sif" ]] ||
    die "Offline cluster: upload a Linux amd64 Ubuntu image to $software/ubuntu-22.04.sif before submitting."
mkdir -p -- "$base/logs"

submit() {
    local reply job
    reply=$(sbatch --parsable --chdir="$base" --output="$base/logs/container-%j.log" "$@") || return
    job=${reply%%;*}
    [[ "$job" =~ ^[0-9]+$ ]] || die "Unexpected sbatch output: $reply"
    printf '%s\n' "$job"
}

setup_job=$(submit --job-name=crouch-container-test -p bigmem -c 8 --mem=64G -t 02:00:00 \
    "$base/crouch_container.sh" --worker setup "$base" "$scratch" "$scope")
printf 'Environment preparation and bundled-data test: %s\n' "$setup_job"
array=2
if [[ "$scope" == remaining ]]; then array='0,1,3-6%1'; fi
count_job=$(submit --job-name=crouch-container-count -p bigmem -c 24 --mem=256G -t 24:00:00 \
    --array="$array" --dependency="afterok:$setup_job" --kill-on-invalid-dep=yes \
    "$base/crouch_container.sh" --worker count "$base" "$scratch" "$scope")
printf 'Cell Ranger count: %s\n' "$count_job"
qc_job=$(submit --job-name=crouch-container-rds -p general -c 4 --mem=64G -t 04:00:00 \
    --dependency="afterok:$count_job" --kill-on-invalid-dep=yes \
    "$base/crouch_container.sh" --worker qc "$base" "$scratch" "$scope")
printf 'Separate sample RDS files: %s\n' "$qc_job"
printf 'Monitor: sacct -j %s,%s,%s --format=JobID,JobName,NodeList,State,ExitCode\n' "$setup_job" "$count_job" "$qc_job"
printf 'Start by reading: %s/logs/container-%s.log\n' "$base" "$setup_job"
printf 'RDS output: %s/qc_slurm/%s/\n' "$base" "$qc_job"
