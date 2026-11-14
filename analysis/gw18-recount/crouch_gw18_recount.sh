#!/bin/bash
# Submit from the login node; all data processing runs in a Slurm allocation.
set -Eeuo pipefail
umask 077
die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }

if [[ "${1:-}" == --inside ]]; then
    shift
    exec "$@"
fi

if [[ "${1:-}" == --worker ]]; then
    : "${SLURM_JOB_ID:?Submit this script through its launcher.}"
    mode=${2:?}; base=${3:?}; scratch=${4:?}
    cd "$base"
    tool="$base/software/crouch_runtime/sra-bin"
    reference="$base/software/crouch_runtime/reference"
    helper="$base/crouch_recount.py"
    star="$base/software/crouch-recount/STAR-2.7.10b"
    image="$base/software/ubuntu-22.04.sif"
    sra="$base/sra_files/SRR17878454.sra"
    for exe in "$CROUCH_SINGULARITY_BIN" "$CROUCH_PYTHON_BIN" "$CROUCH_RSCRIPT_BIN" "$star"; do
        [[ -x "$exe" ]] || die "Executable unavailable on $(hostname): $exe"
    done
    [[ -d "$scratch" && -w "$scratch" ]] || die "Scratch is unavailable: $scratch"
    mkdir -p "$base/recount_slurm/GW18" "$base/software/container_home"
    result="$base/recount_slurm/GW18/${mode}_${SLURM_JOB_ID}"
    mkdir "$result"
    work=$(mktemp -d "$scratch/crouch-recount-${SLURM_JOB_ID}.XXXXXX")
    index_stage=
    cleanup() {
        status=$?
        if [[ -n "${index_stage:-}" && -d "$index_stage" ]]; then rm -rf -- "$index_stage"; fi
        if [[ -d "$work" ]]; then rm -rf -- "$work"; fi
        if (( status != 0 )); then
            printf 'FAILED status=%s; diagnostic output retained in %s\n' "$status" "$result" >&2
        fi
    }
    trap cleanup EXIT
    trap 'exit 130' INT
    trap 'exit 143' TERM
    cp "$base/crouch_gw18_recount.sh" "$helper" "$base/crouch_recount_export.R" "$result/"
    printf 'Mode=%s Job=%s Node=%s\nResult=%s\n' "$mode" "$SLURM_JOB_ID" "$(hostname)" "$result"
    if [[ "$mode" == compare ]]; then
        full_result=${5:?Missing completed full-count directory}
        author_result=${6:?Missing author GW18 directory}
        [[ -f "$full_result/COMPLETE" && -s "$full_result/GW18_for_qc.rds" ]] || die 'Full-count result is incomplete'
        export R_LIBS="$base/software/r-libraries/R-4.3.1-Matrix-1.6.0${R_LIBS:+:$R_LIBS}"
        export R_LIBS_USER="$R_LIBS"
        cp "$base/crouch_compare_author.R" "$result/"
        "$CROUCH_RSCRIPT_BIN" --vanilla "$base/crouch_compare_author.R" \
            "$full_result/GW18_for_qc.rds" "$author_result" "$result/comparison"
        touch "$result/COMPLETE"
        printf 'Comparison complete: %s/comparison/comparison_summary.txt\n' "$result"
        exit 0
    fi
    "$CROUCH_PYTHON_BIN" - "$work" "$mode" <<'SPACE'
import shutil, sys
free = shutil.disk_usage(sys.argv[1]).free
needed = (600 if sys.argv[2] == "full" else 80) * 1024**3
print("Shared scratch free space: %.1f GiB (does not verify user quota)" % (free / 1024**3))
if free < needed:
    raise SystemExit("Insufficient shared scratch space")
SPACE
    export NCBI_SETTINGS="$work/ncbi/user-settings.mkfg"
    mkdir -p "$work/ncbi"
    container() {
        SINGULARITYENV_NCBI_SETTINGS="$NCBI_SETTINGS" \
            "$CROUCH_SINGULARITY_BIN" exec --cleanenv \
            --home "$base/software/container_home:/home/crouch" \
            --bind "$base:$base" --bind "$scratch:$scratch" --pwd "$base" "$image" "$@"
    }
    container "$tool/vdb-config" --set "/repository/user/main/public/root=$base/software/sra-cache"
    container "$tool/vdb-config" --set '/repository/user/main/public/apps/refseq/volumes/refseq=refseq'
    container "$tool/vdb-config" --set '/repository/remote/disabled=true'
    container "$star" --version > "$result/STAR.version.txt"
    [[ "$(cat "$result/STAR.version.txt")" == 2.7.10b ]] || die 'Unexpected STAR version'
    sha256sum "$star" "$helper" "$base/crouch_gw18_recount.sh" "$base/crouch_recount_export.R" > "$result/code_SHA256SUMS"
    # Input fingerprint: large SRA is identified by path, size and mtime; it is
    # not rehashed here. The reference files are hashed on the compute node.
    stat -Lc 'SRA path=%n bytes=%s mtime=%Y' "$sra" > "$result/input.txt"
    fasta="$reference/fasta/genome.fa"
    gtf="$reference/genes/genes.gtf"
    if [[ ! -s "$gtf" ]]; then gtf="$reference/genes/genes.gtf.gz"; fi
    [[ -s "$fasta" && -s "$gtf" ]] || die 'Existing genome.fa or genes.gtf[.gz] not found'
    sha256sum "$fasta" "$gtf" > "$result/reference_SHA256SUMS"
    cat > "$result/method.txt" <<'METHOD'
GW18 / SRR17878454; SRA corrected-CB/UB requantification, not raw barcode re-correction.
SEQUENCE rows are exported once, with original biological-read quality scores.
Separate SAM and STARsolo run for each physical library and original GEM suffix.
Final cell identifier = percent-encoded first three SPOT_GROUP fields + | + full CB including GEM suffix.
Flowcell/lane are retained in read_groups.tsv; they are not separate GEM libraries.
STARsolo 2.7.10b: original corrected 16-base CB, no barcode correction, Exact UMI deduplication.
MultiGeneUMI filtering; unique-gene integer counts only. Gene and GeneFull are both retained.
Original CR/CY/UR/UY are not reconstructed. STAR's missing barcode-quality placeholder is not a measured quality score; no whitelist matching or barcode correction is performed.
Source UB was corrected under the author's original annotation; remapping with GRCh38-2024-A can change gene assignments. Exact reproduction of author counts is not claimed.
No cell calling, cell QC filtering, expression normalization, log transformation or Harmony is applied.
Raw matrix covers recoverable barcodes with retained CB/UB; completeness of original droplets is not assumed.
GW18_for_qc.rds uses exonic Gene counts; GeneFull counts are separately available as MEX.
Pilot counts are a technical test only and must not be used as full-depth biological data.
METHOD
    cp "$result/code_SHA256SUMS" "$work/current_fingerprint"
    cat "$result/reference_SHA256SUMS" "$result/input.txt" >> "$work/current_fingerprint"
    if [[ "$mode" == full ]]; then
        [[ -s "$base/recount_slurm/GW18/pilot_fingerprint.txt" ]] || die 'Run and complete pilot first'
        cmp -s "$work/current_fingerprint" "$base/recount_slurm/GW18/pilot_fingerprint.txt" ||
            die 'Code/reference/SRA changed since pilot; rerun pilot'
    fi
    if [[ "$mode" == pilot ]]; then
        # Spread the technical pilot across the full archive instead of four
        # long contiguous blocks, which can overrepresent repeated/local loci.
        ranges=
        for ((block_idx=0; block_idx<1000; block_idx++)); do
            block_start=$((1 + (1063261690 - 100) * block_idx / 999))
            ranges="${ranges}${ranges:+,}${block_start}-$((block_start + 99))"
        done
        row_args=(-R "$ranges")
        expected=100000
    else
        row_args=()
        expected=1063261690
    fi
    printf 'Exporting %s SEQUENCE rows, preserving read/tag correspondence.\n' "$mode"
    # CSV quoting preserves all original Phred33 characters. The row ID comes
    # from vdb-dump itself, not an assumed order or a second independently read file.
    container "$tool/vdb-dump" "$sra" -T SEQUENCE -I -f csv \
        -C 'READ,SPOT_GROUP,LINKAGE_GROUP,READ_LEN,READ_START,(INSDC:quality:text:phred_33)QUALITY' \
        ${row_args[@]+"${row_args[@]}"} 2> "$result/vdb-dump.stderr" | \
        "$CROUCH_PYTHON_BIN" "$helper" convert "$work/converted" --expected "$expected" \
        > "$result/conversion.log"
    cp "$work/converted/conversion.json" "$work/converted/read_groups.tsv" "$work/converted/input_preview.csv" "$work/converted/library_inputs.tsv" "$result/"
    cat "$result/conversion.log"
    # Reuse a STARsolo-specific index only when reference fingerprints match.
    index="$base/software/crouch-recount/index-GRCh38-2024-A-2.7.10b"
    if [[ -d "$index" ]]; then
        [[ -s "$index/REFERENCE_SHA256SUMS" ]] && cmp -s "$result/reference_SHA256SUMS" "$index/REFERENCE_SHA256SUMS" ||
            die 'Existing STARsolo index has a different/incomplete reference fingerprint'
    else
        index_stage=$(mktemp -d "$base/software/crouch-recount/.index.XXXXXX")
        gtf_input="$gtf"
        if [[ "$gtf" == *.gz ]]; then
            gzip -cd -- "$gtf" > "$work/genes.gtf"
            gtf_input="$work/genes.gtf"
        fi
        printf 'Building the reusable STARsolo genome index on the compute node.\n'
        container "$star" --runMode genomeGenerate --runThreadN "${SLURM_CPUS_PER_TASK:-24}" \
            --genomeDir "$index_stage" --genomeFastaFiles "$fasta" --sjdbGTFfile "$gtf_input" \
            --sjdbOverhang 90 --limitGenomeGenerateRAM 180000000000 \
            --outFileNamePrefix "$result/index." > "$result/index.stdout" 2>&1
        cp "$result/reference_SHA256SUMS" "$index_stage/REFERENCE_SHA256SUMS"
        [[ -s "$index_stage/Genome" && -s "$index_stage/SA" ]] || die 'Index generation incomplete'
        [[ ! -e "$index" ]] || die 'Another job published an index; rerun this job'
        mv -T -- "$index_stage" "$index"
        index_stage=
    fi
    printf 'Mapping and counting with STARsolo; no barcode re-correction or cell filtering.\n'
    mkdir "$work/star"
    while IFS=$'\t' read -r file_id library gem reads; do
        [[ "$file_id" != file_id ]] || continue
        mkdir -p "$work/star/$file_id" "$result/STARsolo/$file_id"
        container "$star" --runThreadN "${SLURM_CPUS_PER_TASK:-24}" --genomeDir "$index" \
            --readFilesIn "$work/converted/$file_id.sam" --readFilesType SAM SE \
            --readFilesSAMattrKeep None --outSAMtype None --outFileNamePrefix "$work/star/$file_id/" \
            --soloType CB_UMI_Simple --soloCBwhitelist None \
            --soloInputSAMattrBarcodeSeq CB UB --soloInputSAMattrBarcodeQual - \
            --soloCBlen 16 --soloUMIstart 17 --soloBarcodeReadLength 26 --soloUMIlen 10 --soloCBmatchWLtype Exact \
            --soloUMIdedup Exact --soloUMIfiltering MultiGeneUMI \
            --soloStrand Forward --soloFeatures Gene GeneFull --soloMultiMappers Unique \
            --soloCellFilter None --clipAdapterType CellRanger4 \
            > "$result/STARsolo/$file_id/stdout.txt" 2>&1 || {
                cp "$work/star/$file_id"/Log* "$result/STARsolo/$file_id/" 2>/dev/null || true
                die "STARsolo failed for $library GEM $gem; inspect STARsolo/$file_id/"
            }
        cp "$work/star/$file_id"/Log* "$result/STARsolo/$file_id/"
        for feature in Gene GeneFull; do
            mkdir -p "$result/STARsolo/$file_id/$feature"
            for name in Summary.csv Barcodes.stats Features.stats UMIperCellSorted.txt; do
                file="$work/star/$file_id/Solo.out/$feature/$name"
                if [[ -f "$file" ]]; then cp "$file" "$result/STARsolo/$file_id/$feature/"; fi
            done
        done
    done < "$work/converted/library_inputs.tsv"
    for feature in Gene GeneFull; do
        "$CROUCH_PYTHON_BIN" "$helper" merge "$work/converted/library_inputs.tsv" "$work/star" "$work/merged-$feature" --feature "$feature"
    done
    "$CROUCH_PYTHON_BIN" "$helper" package "$work/merged-Gene" "$result/raw_feature_bc_matrix"
    "$CROUCH_PYTHON_BIN" "$helper" package "$work/merged-GeneFull" "$result/raw_gene_full_matrix"
    export R_LIBS="$base/software/r-libraries/R-4.3.1-Matrix-1.6.0${R_LIBS:+:$R_LIBS}"
    export R_LIBS_USER="$R_LIBS"
    "$CROUCH_RSCRIPT_BIN" --vanilla "$base/crouch_recount_export.R" "$result" "$mode"
    touch "$result/COMPLETE"
    if [[ "$mode" == pilot ]]; then
        cp "$work/current_fingerprint" "$base/recount_slurm/GW18/pilot_fingerprint.txt.tmp"
        mv -f "$base/recount_slurm/GW18/pilot_fingerprint.txt.tmp" "$base/recount_slurm/GW18/pilot_fingerprint.txt"
    fi
    printf '\nSUCCESS: %s\n' "$result"
    cat "$result/raw_feature_bc_matrix/matrix_summary.json"
    [[ "$mode" != pilot ]] || printf 'PILOT ONLY: inspect read_groups.tsv, conversion.json, STARsolo/lib_*/Log.final.out and QC_summary.txt before full counting.\n'
    exit 0
fi

mode=${1:-pilot}
[[ "$mode" == pilot || "$mode" == full || "$mode" == compare ]] || die 'Usage: bash crouch_gw18_recount.sh [pilot|full] [SCRATCH_ROOT] | compare FULL_RESULT AUTHOR_GW18_DIR'
base=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
extra_args=()
if [[ "$mode" == compare ]]; then
    full_result=$(cd -- "${2:?Provide FULL_RESULT}" && pwd -P)
    author_result=$(cd -- "${3:?Provide AUTHOR_GW18_DIR}" && pwd -P)
    extra_args=("$full_result" "$author_result")
    scratch=/public4/usersdir/guojunping/scratch
else
    scratch=${2:-/public4/usersdir/guojunping/scratch}
fi
for file in crouch_recount.py crouch_recount_export.R crouch_compare_author.R software/crouch-recount/STAR-2.7.10b software/ubuntu-22.04.sif sra_files/SRR17878454.sra; do
    [[ -s "$base/$file" ]] || die "Missing: $base/$file"
done
[[ -d "$scratch" ]] || die "Scratch directory does not exist: $scratch"
if ! type module >/dev/null 2>&1; then
    for init in /etc/profile.d/modules.sh /usr/share/Modules/init/bash /usr/share/lmod/lmod/init/bash; do
        if [[ -r "$init" ]]; then set +u; source "$init"; set -u; break; fi
    done
fi
type module >/dev/null 2>&1 || die 'Environment Modules unavailable'
set +u
module use /public/software/modules
module load compiler/gcc/7.4.0 system/singularity/3.6.3 system/anaconda/3/2024.06 system/R/4.3.1
set -u
fortran_lib=$(gfortran -print-file-name=libgfortran.so.4)
[[ "$fortran_lib" == /* && -r "$fortran_lib" ]] || die 'GCC 7.4 libgfortran.so.4 unavailable'
export LD_LIBRARY_PATH="${fortran_lib%/*}${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
CROUCH_SINGULARITY_BIN=$(command -v singularity)
CROUCH_PYTHON_BIN=$(command -v python3)
CROUCH_RSCRIPT_BIN=$(command -v Rscript)
export CROUCH_SINGULARITY_BIN CROUCH_PYTHON_BIN CROUCH_RSCRIPT_BIN
chmod u+x "$base/software/crouch-recount/STAR-2.7.10b"
mkdir -p "$base/logs"
node=${CROUCH_NODE:-comput159}
[[ "$node" =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ ]] || die 'Invalid node name'
cpus=24; memory=256G
if [[ "$mode" == pilot ]]; then
    duration=${CROUCH_TIME:-06:00:00}
elif [[ "$mode" == full ]]; then
    duration=${CROUCH_TIME:-24:00:00}
else
    cpus=4; memory=64G; duration=${CROUCH_TIME:-02:00:00}
fi
reply=$(sbatch --parsable --export=ALL --job-name="gw18-recount-$mode" -p bigmem --nodelist="$node" \
    -N 1 -n 1 -c "$cpus" --mem="$memory" -t "$duration" --chdir="$base" \
    --output="$base/logs/gw18-recount-$mode-%j.log" \
    "$base/crouch_gw18_recount.sh" --worker "$mode" "$base" "$scratch" ${extra_args[@]+"${extra_args[@]}"})
job=${reply%%;*}
[[ "$job" =~ ^[0-9]+$ ]] || die "Unexpected sbatch reply: $reply"
printf 'Job=%s\nLog=%s/logs/gw18-recount-%s-%s.log\nResults=%s/recount_slurm/GW18/%s_%s/\n' "$job" "$base" "$mode" "$job" "$base" "$mode" "$job"
printf 'Monitor: sacct -j %s --format=JobID,State,ExitCode,Elapsed\n' "$job"
