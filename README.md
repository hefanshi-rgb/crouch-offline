# Crouch analysis scripts

Completed transfer payloads have been removed from the current branch. Analysis source code is retained under `analysis/`; existing cluster software, input SRA files, counts and RDS outputs are not changed by this repository cleanup.

- `analysis/gw18-recount/`: SRA corrected-CB/UB UMI requantification (STARsolo 2.7.10b).
- `analysis/gw18-qc/`: per-library cell calling and exploratory QC; all reference genes retained.
- `analysis/author-meta/`: metadata-only cell overlap, QC-column audit and author annotation transfer after explicit library/GEM mapping review.
- `tools/prepare_author_meta_on_mac.py`: select and relay the complete author meta.tsv through 2 MiB Git blobs. No matrix.mtx.gz or RDS is read.

Metadata-only checks cannot validate gene-expression concordance or reproduce a gene-by-cell counts matrix. Both lowercase and uppercase count/feature metadata columns are audited independently.

## Completed transfer cleanup

Removed legacy SRA Toolkit, Ubuntu container, RefSeq (20 MiB and 2 MiB variants), and GW18 recount archive parts and obsolete receiver entry points. Previous commit `42fca0ea8e3aa25d23e96a1a702ad1db3931d3e5` retains the historical transfer manifests. This is a normal deletion commit; Git history was not rewritten.

The retained source assumes the established cluster paths and existing offline software. Do not overwrite running cluster scripts from a repository update.
