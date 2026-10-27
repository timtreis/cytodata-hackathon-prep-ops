# AI Agent Instructions

Data prep for the CytoData hackathon at SBI2 (Nov 1-2 2026): Biohub's "Altair - A549 1k" optical pooled screen (Liu, Leonetti et al. 2026, doi:10.64898/2026.06.01.728087), QC'd and packaged for participants.

## Dataset

- A549, inducible Cas9, 1,000-gene KO library (4 sgRNAs/gene + 211 NTCs). Guide identity via in situ sequencing.
- Public bucket `s3://ops-explorer-public/leonetti_ops/ops_data_portal_submission/` (anonymous HTTPS listing). Use `v2.0.20260724`: 88 per-screen datasets `Biohub_OPS0001..0088`.
- Each screen is its own set of cells: Phase2D + 0-2 live markers. Marker names live in the plate `zarr.json` under `attributes.channels_metadata`, not in the generic channel labels (GFP/mCherry).
- Multiplexed screens: OPS0043 and OPS0088 (Cell Painting, 7 stains, 2 cycles), OPS0077 (4i, 6 antibodies).
- CROP-seq is a separate experiment; only perturbation-level h5ad is public. Imaging and RNA join on gene, never on cell.
- Analysis code: github.com/czbiohub-sf/cyclops-monorepo, github.com/czbiohub-sf/ops-paper-analysis.

## Staged data (`scripts/fetch_altair.py`, features only, no images)

- `v2.0.20260724/atlas/cell_features/all_cells_fluor_<marker>.h5ad`: 59 markers, per-cell hand-crafted features (`op_*` OrganelleProfiler, `cp_single_object_*` CellProfiler-style), ~1k features each.
- `v2.0.20260724/atlas/cell_features/all_cells_phase.h5ad`: 64.6M cells x 3,317, 870 GB, all screens pooled. Subsample before handing out.
- `v1.0.20260521/atlas/cell_features/all_cells_phase.h5ad`: 8.4M-cell phase subset, 62 GB.
- `obs` carries `gene`, `sgRNA`, `barcode`, `experiment`, `well`, `segmentation_id`.
- `processed/phase_v2_n1000_per_sgRNA.h5ad` (`scripts/subsample_phase.py`, seed 0): 4,208,364 x 3,317, 56.7 GB. 4,211 sgRNAs (4,208 at the 1,000 cap, min 52), 1,001 genes incl. NTC, all 88 screens.
- Also staged: gene-level `aggregated_data.h5ad`, CROP-seq, and per-screen `cell_data.parquet` (cell index only: well, x/y, barcode, perturbation, bbox) + metadata + plate `zarr.json`.

## Rules

- Terse output. Lead with results. Say plainly when something failed.
- Data never goes in git. It lives on the cluster.
- Cluster hostnames, paths, accounts and login details live only in `CLAUDE.local.md` (gitignored). Never write them into tracked files, commits, PRs or issues.
- Use scverse tooling (spatialdata, anndata, scanpy). No new dependencies without asking.
- Never commit or push without being asked. Stage only files you touched.
- `tasks/` is local scratch and gitignored.
