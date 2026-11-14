GW18: existing SRA -> unnormalized integer UMI counts, features and RDS
====================================================================

This bundle adds a new workflow. It does not replace crouch_container.sh,
crouch_slurm.sh, the source SRA, the existing container, or the reference cache.
The Linux STAR 2.7.10b static executable and its MIT license are included.
Upstream source: https://github.com/alexdobin/STAR/releases/tag/2.7.10b
No further genomic-data download is performed on the offline cluster.

1. Upload and unpack
--------------------
Upload crouch_gw18_recount_bundle.tar.gz into:
/public4/usersdir/guojunping/guojunping_2026-11-08_001/upload/Crouch/

In the SSH shell (directory management/submission only on login):

cd /public4/usersdir/guojunping/guojunping_2026-11-08_001/upload/Crouch
tar -xzf crouch_gw18_recount_bundle.tar.gz
sha256sum -c RECOUNT_SHA256SUMS
bash crouch_gw18_recount.sh pilot

This submits 24 CPUs / 256 GiB / 6 hours on bigmem comput159.
Actual processing runs under Slurm, using the existing Ubuntu 22.04 container
and existing host Python/R module installations. The GCC 7.4 Fortran library
and the existing shared R Matrix snapshot are reused.

Pilot: 100,000 SEQUENCE rows from four distant blocks, not a representative
biological downsampling. It verifies extraction, grouping, alignment and counts.
Index generation uses the existing GRCh38-2024-A FASTA/GTF and can take a while;
the persistent STAR-specific index needs roughly tens of GB additional space.
Index construction, SRA extraction, mapping and R work all run on the compute node.
Shared scratch free-space checks do not verify the account's individual quota.

The command prints Job, Log and Results. Use the actual printed job ID:
sacct -j JOB_ID --format=JobID,State,ExitCode,Elapsed
tail -n 100 logs/gw18-recount-pilot-JOB_ID.log

Pilot reports:
recount_slurm/GW18/pilot_JOB_ID/conversion.json
recount_slurm/GW18/pilot_JOB_ID/read_groups.tsv
recount_slurm/GW18/pilot_JOB_ID/library_inputs.tsv
recount_slurm/GW18/pilot_JOB_ID/STARsolo/lib_*/Log.final.out
recount_slurm/GW18/pilot_JOB_ID/QC_summary.txt

Review these reports before full counting. A successful pilot is necessary but
not sufficient to establish full-depth data quality or completeness. Do not use
GW18_pilot_NOT_FULL_DATA.rds for downstream integration.

2. Full-depth counting, after the pilot is checked
--------------------------------------------------
bash crouch_gw18_recount.sh full

Default allocation: 24 CPUs / 256 GiB / 24 hours, comput159.
Only if needed and permitted by the partition, override walltime at submission:
CROUCH_TIME=48:00:00 bash crouch_gw18_recount.sh full

The full run checks that code, source SRA size/mtime and reference fingerprints
match a completed pilot. It expects 1,063,261,690 SEQUENCE rows, as measured in
the user's previous full dump. Unexpected layouts or truncated exports fail.
Full conversion/mapping can take many hours. At least 600 GiB free shared
scratch is required conservatively; final outputs and the index also need quota.
Do not start multiple pilot/full jobs concurrently. No automatic infinite retry.
Stop a job with: scancel JOB_ID
Temporary SAMs and STAR work files are removed when the job exits normally,
fails, or receives a catchable cancellation signal. A hard node/Slurm kill may
leave scratch files. The source SRA and completed published results are retained.

Completed full-depth output:
recount_slurm/GW18/full_JOB_ID/
  COMPLETE
  GW18_for_qc.rds
  raw_feature_bc_matrix/{matrix.mtx.gz,features.tsv.gz,barcodes.tsv.gz}
  raw_feature_bc_matrix/cell_metadata.tsv
  raw_gene_full_matrix/{matrix.mtx.gz,features.tsv.gz,barcodes.tsv.gz}
  QC_summary.txt
  conversion.json, read_groups.tsv, library_inputs.tsv
  method.txt, input.txt, reference_SHA256SUMS, code_SHA256SUMS
  STARsolo/ (per-library alignment and counting reports)

GW18_for_qc.rds is an R LIST (not an already integrated Seurat object):
  counts        sparse unnormalized integer UMI matrix, Gene/exonic counts
  features      gene_id, gene_name, feature_type
  cell_metadata sample/run/library/CB and count-based QC metrics
  provenance    method, inputs and fingerprints

No normalization, cell calling, doublet removal, mitochondrial filtering or
Harmony is applied. All recoverable corrected CB/UB labels are considered.
The matrix includes ambient/empty droplets if they are retained in the SRA;
cell selection is still needed before integration. An author-cell subset can
be made later after an explicit library+barcode match.

3. Prepare only the author's GW18 reference on the Mac
------------------------------------------------------
Unpack this bundle locally too. From a LOCAL Mac terminal run:

python3 /PATH/TO/UNPACKED/crouch_recount.py author-subset \
  /PATH/TO/scHE10 /PATH/TO/NEW/GW18_author_reference

Quote paths containing spaces. Replace /PATH/... with actual local paths.
The input directory must contain meta.tsv, barcodes.tsv.gz, features.tsv.gz,
matrix.mtx.gz. The script selects age=GW18, scans the large matrix as a stream,
and copies matching columns with expression values UNCHANGED. It does not load
the full matrix into RAM and never modifies the author's input files.
Progress prints every 10 million matrix entries. Run on the online Mac, not the
cluster login node. A new output directory is required to avoid mixing attempts.

The smaller GW18_author_reference directory contains the three MEX files,
meta.tsv, subset_summary.json and library_mapping_to_review.tsv.
The library mapping is deliberately left blank. Fill recount_library and
recount_gem using the real SRA read-group report and author sample definitions.
Do not guess matches by removing every prefix or matching barcode alone.
Send the metadata/prefix summary for review if mapping is unclear. Only after
review, upload the GW18 subset to Crouch/author_reference/GW18/ on the cluster.
You do not need to transfer the 10 GB RDS or the full 2 GB author matrix.

4. Compare full-depth counts with the author reference
-------------------------------------------------------
After library_mapping_to_review.tsv has been filled and reviewed:

bash crouch_gw18_recount.sh compare \
  "$PWD/recount_slurm/GW18/full_JOB_ID" \
  "$PWD/author_reference/GW18"

This submits a separate 4 CPU / 64 GiB / 2 hour Slurm job. Reports are in:
recount_slurm/GW18/compare_JOB_ID/comparison/
  comparison_summary.txt
  cell_matches.tsv
  gene_matches.tsv
  cell_expression_correlations.tsv
  marker_by_author_cell_type.tsv (when matching markers/labels are available)

The comparison checks exact library+GEM+barcode matching, exact unique gene IDs
(or unique symbols only if no gene IDs match), and exploratory per-cell Spearman
expression agreement on at most 200 matched cells. Genes zero in both cells are
excluded from the per-cell correlation to avoid shared-zero inflation.
Recount values are log1p(CP10k) only for comparison; stored raw counts are intact.
Author values are not normalized again. Their original transformation must be
confirmed before a claim about equivalent numerical expression scales.
Marker profiles are summarized under the author's matched cell-type labels.

High correlation supports expression concordance, not exact reconstruction of
the author's counts. Inspect mapping, library sizes, detected genes, markers,
reference/annotation versions, exonic/intronic policies and source tag retention.
No arbitrary correlation threshold is treated as a scientific pass/fail rule.

Method choices and limits
--------------------------
* The read and its CB/UB are exported from the SAME SEQUENCE row, with actual
  row ID and original read quality. No primary/secondary-alignment duplication.
* This audited SRA stores one biological read/spot. A different layout stops.
* First three SPOT_GROUP components identify the library; flowcell/lane do not
  make separate libraries. Different GEM suffixes are kept separate too.
* Mapping/counting is independent per library/GEM, with the original corrected
  16-base CB. Final matrix identifiers contain library and full CB/GEM suffix.
* Exact UMI collapsing avoids re-correcting an already corrected UB. Standard
  STARsolo MultiGeneUMI handling removes lower-supported cross-gene UMIs.
* No fractional multi-mapper allocation; integer unique-gene counts only.
* Missing CB/UB and invalid UMI exclusions are reported, not filled or guessed.
* Barcode quality tags are absent. STAR emits a warning about its default 'H'
  placeholder. No whitelist matching/barcode correction uses that placeholder;
  it is not claimed to be measured sequencer barcode quality.
* Reusing UB retains the author's gene-dependent correction decisions. Mapping
  with GRCh38-2024-A can change gene assignments relative to the old annotation.
  These are documented UMI RECOUNTS, not a claim to reproduce original counts.
* Recovery of all original low-count droplets cannot be inferred from a 20-row
  audit; inspect full conversion statistics and barcode distributions.

Validation performed before delivery
--------------------------------------
Python conversion tests: original quality preservation including CSV quotes,
missing tags, unexpected layout, repeated rows, library/lane/GEM identities,
non-integer/invalid matrix rejection, exact-value author subsetting.
Actual STAR 2.7.10b synthetic tests on Mac: cross-lane duplicate UMI collapse,
cross-library identical CB/UB separation, exonic/intronic counts and strand rules.
R sparse-count export and explicit author-library mapping tests.
These local tests do not replace the cluster pilot on the actual SRA.
