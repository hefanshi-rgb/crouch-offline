# GW18 R7CE original-BAM GX reconstruction

This is an independent no-realignment baseline. It does not replace the reviewed STARsolo output and is not asserted to exactly reproduce Cell Ranger's original molecule filtering.

Inputs:
- ENA SRR17878454 / R7ctxendo_possorted_genome_bam.bam.1, 11724258461 bytes, ENA MD5 9b978e60107ee907e55c12659d4f5a52.
- All 1759 exact author all_R7_ctx_endo_ cells. Full corrected CB suffix -1 is retained. Original RG is checked as R7CE:MissingLibrary:1; sequencing lanes are pooled within that physical library.
- All 33694 gene IDs from the checked GRCh38-1.2.0 GTF, including zero-count genes.

Count rules: primary, mapped, MAPQ255, exactly one nonempty GX, complete valid 16bp CB plus GEM and 10bp corrected UB. Deduplicate (library,CB,GX,UB). Read duplicate flag alone is NOT an exclusion: reads supporting one UMI are counted once. Secondary/supplementary alignments, ambiguous GX, absent tags and N-containing UBs are recorded/excluded. Unexpected UMI lengths, GEMs, library identity or gene IDs stop the run. RE, TX and available xf/mm tag presence are reported; modern xf is not assumed to exist in this old BAM. Multiple genes sharing CB/UB are reported, not silently reassigned or removed.

The SQLite molecule table keeps bounded Python memory usage. No cell QC, gene filtering, scaling or normalization is applied to saved integer counts. Exonic RE support is audited separately. Cross-gene UMI and other original Cell Ranger filters may make this reconstruction differ from the original matrix.

Files:
- download_bam.py: eight-way 32MiB resumable ranges, validated HTTP206 Content-Range, full ENA MD5 and SHA256 before BAM acceptance.
- prepare_inputs.py: exact author cohort and reference features; validates source reference-part SHA256s.
- count_gx.py: streaming SQLite UMI reconstruction.
- export_and_compare.R: raw count RDS integrity checks; author processed-expression comparison in the same existing 50-cell sample; optional full 1759-cell STARsolo comparison by exact gene IDs/cell IDs.
- export_STARsolo_R7CE_Windows.R: extracts the real reviewed full-depth old-reference R7CE subset in Windows RStudio without changing its counts.
- test_count_gx.py and test_comparison.R: scientific edge cases and deliberate row/column permutation / introduced count difference.
- run_pipeline.py: waits for fully verified BAM then counts and compares automatically.

Current runtime and outputs: work/gx-recount under the project root. status.json is current pipeline state; download.log and pipeline.log contain details. R7CE_full_GX is the formal full-depth output; benchmark and comparison_test are explicitly non-biological test directories and must not be interpreted as full results. A caffeinate process prevents idle system sleep while the pipeline is active. Closing the app does not intentionally terminate this detached process, but shutdown/manual sleep/network loss may interrupt it. Completed and partially downloaded segments are retained for restart.

Current agreed workflow: complete GX counts and the fixed author 50-cell expression comparison on Mac, then transfer the small GX validation bundle to Windows using the established GitHub mirror download route. Run compare_GX_STARsolo_RStudio.R in Windows RStudio and select the GX RDS, the reviewed full-depth GW18_author_cells_GRCh38_1.2.0_raw_counts.rds (33694 genes,13906 cells,28758442 UMIs), and bundled author_50_R7CE_slots.rds. No Windows-to-Mac matrix upload is required; return the resulting summary tables or screenshots. export_STARsolo_R7CE_Windows.R remains an optional alternative if an approved return channel becomes available.

Raw comparison reports exact differing entries and absolute UMI differences before correlations. Correlations exclude shared-zero genes per cell, logCP10k uses all original reference genes as denominators. Author RNA/counts are fractional processed expression and are not integer raw UMI ground truth. No arbitrary expected correlation threshold is assumed.

No GitHub push, credential access, cluster submission or existing-matrix replacement is performed.

Return-channel probe: cloud_return_test_Windows.sh creates only a timestamped harmless text and SHA256SUMS, manually pushes to a new cloud-return-test branch through the existing Windows mirror route. No main update, force push, existing-file upload or credential-helper access. A local bare-remote + separate Mac receiver roundtrip verified SHA256. Actual cloud mirror upload capability is unverified until the user runs it in Windows; an authentication error does not prove the mirror is read-only.

Supplied original BAM reused on 2026-11-16: full ENA MD5 and SHA256 verified, pysam quickcheck passed, original library/GEM and 194-reference dictionary checked. Adjacent original BAI independently matched ENA MD5; expected complete scan is 273774264 alignment records. No whole-BAM download is needed.
