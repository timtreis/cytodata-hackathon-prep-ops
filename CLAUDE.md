# AI Agent Instructions

Data prep for the CytoData hackathon at SBI2 (Nov 1-2 2026): Biohub's "Altair - A549 1k" optical pooled screen (Liu, Leonetti et al. 2026, doi:10.64898/2026.06.01.728087), QC'd and packaged for participants.

## Dataset

- A549, inducible Cas9, 1,000-gene KO library (4 sgRNAs/gene + 211 NTCs). Guide identity via in situ sequencing.
- Public bucket `s3://ops-explorer-public/leonetti_ops/ops_data_portal_submission/` (anonymous HTTPS listing). Use `v2.0.20260724`: 88 per-screen datasets `Biohub_OPS0001..0088`.
- Each screen is its own set of cells: Phase2D + 0-2 live markers. Marker names live in the plate `zarr.json` under `attributes.channels_metadata`, not in the generic channel labels (GFP/mCherry).
- Multiplexed screens: OPS0043 and OPS0088 (Cell Painting, 7 stains, 2 cycles), OPS0077 (4i, 10 antibodies). `screen_map.csv` maps `experiment` (ops0031...) to Biohub_OPS00xx.
- CROP-seq is a separate experiment; only perturbation-level h5ad is public. Imaging and RNA join on guide sequence, never on cell.
- A cell is (experiment, well, tile_pheno, segmentation_id). segmentation_id alone repeats across wells and adjacent tiles.
- Analysis code: github.com/czbiohub-sf/cyclops-monorepo, github.com/czbiohub-sf/ops-paper-analysis.

## Pipeline

1. `scripts/fetch_altair.py`: mirror the per-cell feature h5ads (59 markers + phase, 1.38 TB), atlas tables, CROP-seq and per-screen metadata. No images.
2. `scripts/build_altair.py phase|marker`: per screen, |x| > 1e20 to NaN, shared obs columns + mantispy schema, flag-only QC (`qc_pass`), mantispy robust z vs that screen's NTC cells passing QC. Writes cells_n1000/cells_n100 (hash-ranked per guide), guide_x_screen and guide (median of per-screen medians). Selection = mantispy defaults + drop_outliers(500) on cells_n100 passing QC. Resumable: per-screen checkpoints in `tiers/_ckpt_<name>/`.
3. `scripts/unexpressed_genes.py`: the 96 target genes A549 doesn't express (HPA 25.1), the de facto cutting control.
4. `scripts/final_check.py`: shapes, sizes, max |z| of selected features, mantispy validation, `markers.csv`, SHA-256.
5. `notebooks/make_starter.py` writes `notebooks/01_starter.ipynb`.

mantispy is pinned at commit 3552713 (with the `map` extra) in one shared cluster venv; don't let parallel jobs install it.

## Rules

- Terse output. Lead with results. Say plainly when something failed.
- Data never goes in git. It lives on the cluster.
- Cluster hostnames, paths, accounts and login details live only in `CLAUDE.local.md` (gitignored). Never write them into tracked files, commits, PRs or issues.
- Use scverse tooling (spatialdata, anndata, scanpy). No new dependencies without asking.
- Never commit or push without being asked. Stage only files you touched.
- `tasks/` is local scratch and gitignored.
