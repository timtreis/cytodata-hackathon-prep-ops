# SBI2 CytoData hackathon: optical pooled screen dataset

The dataset is Biohub's Altair screen, an optical pooled CRISPR screen of 1,000 gene knockouts in A549 cells (Liu, Leonetti et al. 2026, [doi:10.64898/2026.06.01.728087](https://doi.org/10.64898/2026.06.01.728087)). Each cell was imaged label-free (quantitative phase) and with one or two fluorescent markers, and its guide was read out by in situ sequencing. We ship per-cell feature tables and guide-level profiles, and we don't ship images.

Biohub published the original data on the [OPS Explorer](https://biohub.ai/ops-explorer) and in the public bucket `s3://ops-explorer-public/leonetti_ops/ops_data_portal_submission/v2.0.20260724/`.

## Experimental design

- The cells are A549 lung adenocarcinoma cells with doxycycline-inducible Cas9. The library targets 1,000 genes with 4 guides each, plus 211 non-targeting guides (`NTC`), so there are 4,211 guides.
- The screen is split into 88 separate screens (`batch`). Every screen images phase, and most screens add one or two live markers, e.g. an endogenously tagged protein or a dye. Two screens are Cell Painting (7 stains) and one screen is 4i (10 antibodies).
- Each screen is its own set of cells. Phase and the markers of one screen come from the same cells, but no cell appears in two screens, so you can only compare markers from different screens at the guide level.

## Download

| File | Contents | Size |
|---|---|---|
| `phase_guide.h5ad` | One phase profile per guide (median of the per-screen medians) | 0.05 GB |
| `phase_guide_x_screen.h5ad` | One phase profile per guide and screen (median over the cells) | 3.1 GB |
| `phase_cells_n100.h5ad` | Up to 100 cells per guide, 3,317 phase features | 6.3 GB |
| `phase_cells_n1000_partNN.h5ad` | Up to 1,000 cells per guide, 3,317 phase features, split by screen into parts below 25 GB (see `n1000_shards.csv`) | 4 parts, 63 GB in total |
| `markers_guide.zip` | `<marker>_guide.h5ad` and `<marker>_guide_x_screen.h5ad` for each of the 59 markers | 2.0 GB |
| `<marker>_cells_n100.h5ad` | Up to 100 cells per guide for one marker, with the selected phase features of the same cells, one file per marker (see `markers.csv`) | 0.9 to 1.7 GB each, 76 GB in total |
| `unexpressed_genes.csv` | The 96 target genes that A549 doesn't express, with their A549 `ntpm` and a `multicopy` flag (by HGNC gene family), derived from the Human Protein Atlas version 25.1 (proteinatlas.org, CC BY-SA 3.0, downloaded 2026-10-27) | 2 KB |
| `n1000_shards.csv` | The screens, cell count, size and SHA-256 of each `phase_cells_n1000` part | 1 KB |
| `markers.csv` | The screens, cell count and size of each `<marker>_cells_n100.h5ad` | 5 KB |
| `screen_map.csv` | Maps each `batch` to Biohub's per-screen dataset id (`Biohub_OPS0001` to `Biohub_OPS0088`) | 9 KB |

Most analyses only need the guide-level files. Start with `phase_guide.h5ad` and `phase_guide_x_screen.h5ad`, add `phase_cells_n100.h5ad` for single-cell work, and use the `phase_cells_n1000` parts and the marker files when you need more cells or the fluorescent markers. For the full depth of the fluorescence data, use Biohub's bucket.

```bash
URL=https://exampledata.scverse.org/sbi2-cytodata/ops
for f in phase_guide.h5ad phase_guide_x_screen.h5ad phase_cells_n100.h5ad screen_map.csv unexpressed_genes.csv; do curl -O $URL/$f; done
shasum -a 256 -c <<EOF
9cf51ab18ac710e88583a2e4ccdae5e28aa367834eb28a2d983bbd81cbf9839f  phase_guide.h5ad
56b295abe4c57a5f2668ed369965a2ad694272232b0126c4da809bfa53386b48  phase_guide_x_screen.h5ad
32f0f289d95acdc7ba44fc90cd0cda1621a37b8b610d1b7cf8eb5ee1861375dd  phase_cells_n100.h5ad
9dbd0486fbb685403ec26eff87dc12bb06e7ac773a2a0428ad3ff3482ad3d509  screen_map.csv
7ca185b0bc4b6f4c0f50f517846e9bd783ca0f61ad19fde04ac6ea4367b16acc  unexpressed_genes.csv
EOF
```

The other files are at the same URL, e.g. `$URL/er_golgi_cop-ii_sec23a_cells_n100.h5ad`, and `manifest.csv` lists the SHA-256 of every file. All h5ad files are gzip-compressed, and anndata reads them as they are. Download with curl, wget, requests or pooch, because the CDN blocks Python's default urllib user agent. To put the `phase_cells_n1000` parts back together:

```python
import anndata as ad, pandas as pd

parts = pd.read_csv("n1000_shards.csv")["shard"]
cells = ad.concat([ad.read_h5ad(p) for p in parts], merge="same", uns_merge="first")  # about 115 GB in memory
```

```python
import anndata as ad, h5py, numpy as np
from anndata.io import read_elem

guides = ad.read_h5ad("phase_guide.h5ad")
with h5py.File("phase_cells_n100.h5ad") as f:  # backed mode would load layers["raw"], 5.5 GB
    obs, var = read_elem(f["obs"]), read_elem(f["var"])
    rows, cols = ~obs["qc_fail"].values, np.flatnonzero(var["selected"].values)
    X = np.vstack([f["X"][i : i + 10_000][rows[i : i + 10_000]][:, cols] for i in range(0, len(obs), 10_000)])
cells = ad.AnnData(X, obs=obs[rows], var=var.iloc[cols])  # about 0.5 GB
```

The notebook `notebooks/01_starter.ipynb` downloads the guide-level files and `phase_cells_n100.h5ad`, scores each knockout against the cutting control with mantispy, builds a first map of the knockouts, joins Biohub's CROP-seq embedding by guide, and shows which knockouts only the SEC23A marker picks up. It runs in about 2 minutes once the files are downloaded and needs about 4 GB of memory. It needs more on a large server than on a laptop, because copairs, which computes the mAP scores, uses one thread per core.

## Contents of the h5ad files

The cell files keep every sampled cell except the ones listed under Known issues. The main `obs` columns are:

- `target_gene` (`NTC` for non-targeting guides), `guide` (the 20 nt guide sequence) and `batch` (the screen)
- the cell key (`experiment`, `well`, `tile_pheno`, `segmentation_id`), which identifies a cell uniquely and joins a marker file to the phase file
- the QC flags, which are true for cells with a problem: `qc_high_nan` (more than half of the features missing), `qc_tile_duplicate` (the same cell on two neighboring tiles, often with two different guide calls) and mantispy's `qc_area_outlier`
- `qc_fail`, which is true when any of the QC flags is set, so `mantispy.pp.filter_cells` (mantispy 0.1.0a2 or newer) works on the files directly

We processed the files with [mantispy](https://github.com/scverse/mantispy) at commit `3552713` and renamed its QC column to `qc_fail` for mantispy 0.1.0a2, which you need to run the starter notebook (`pip install "mantispy[map]>=0.1.0a2"`). Values above 1e20 in absolute size are float overflow in Biohub's Hu moment features, and we set them to missing. `X` holds robust z-scores against the non-targeting cells of the same screen, and the raw values are in `layers["raw"]`. `var["selected"]` marks the features that mantispy's default feature selection keeps, with its outlier filter added (no cell above 500 in absolute z-score among the up to 100 cells per guide). The guide-level files use the same features and the same `selected` column. Clip the z-scores, e.g. to 50 in absolute value, before you use distance-based methods, because a few rare cells have features whose spread among the non-targeting cells is close to zero. The marker files hold normalized values only, and you can get the raw values from Biohub's bucket with the cell key. `qc_fail` reflects the features of each file, so the same cell can pass in the phase file and fail in a marker file, and you should join files on the cell key and not on `qc_fail`.

## Embeddings and the CROP-seq data

Biohub's gene embeddings and a CROP-seq (single-cell RNA) screen of the same library are on Zenodo ([doi:10.5281/zenodo.20495191](https://doi.org/10.5281/zenodo.20495191), CC-BY-4.0). The CROP-seq file `figures/figure_5/cropseq_svaeplus.h5ad` has one embedding per guide, and its `sgRNA` column is the same guide sequence as our `guide` column, so you can join it to `phase_guide.h5ad`. The starter notebook shows the join. The CROP-seq cells are not the imaged cells, so the two modalities only line up at the guide level.

## Known issues

- Biohub's phase features have about 29% missing values. The CellProfiler features are missing for the 20% of cells that weren't matched to a Cell Painting round, and some organelle features are missing in most cells because the organelle is rarely found.
- For about 1,100 features the non-targeting cells of at least one screen have no spread, so the z-score is undefined in that screen. mantispy flags these features in `var["degenerate_scale"]` and divides by a tiny constant (1e-18, as pycytominer does), so their values in that screen are either 0 or very large. Feature selection removes them, so use `var["selected"]` or drop the flagged features.
- The library has no cutting control, so every targeting guide cuts the DNA and the non-targeting guides don't, and the DNA damage from the cut alone can make a knockout separate from `NTC`. In phase, 99.9% of the knockouts separate from `NTC`, but 75% separate from the cutting control described next. 96 of the 1,000 target genes aren't detected in A549 by the Human Protein Atlas (nTPM below 1, listed in `unexpressed_genes.csv`, e.g. CA7, CCR2, CD5 and CDH15), and you can use their knockouts as a cutting control instead of `NTC`, leaving out the 22 genes with `multicopy` true, because their guides match many near-identical copies of the gene and cut the DNA more often.
- Screen `ops0174_20260702` (`Biohub_OPS0088`, a Cell Painting screen) has no phase files, because 2,549 of the 3,317 phase features are missing for all its cells in Biohub's table. Its Cell Painting marker files are shipped, and the phase columns appended to them are missing (NaN) for its cells.
- We left out 165,668 cells (0.26%) that Biohub couldn't place on a well or an imaging tile, 111,282 in `ops0094_20251217`, 54,382 in `ops0174_20260702` (the two Cell Painting screens) and 4 in other screens. They are kept in Biohub's source files, which you can join on `obs_names`.
- Biohub's phase table has 97,371 rows (0.15%) that share a screen, well and `segmentation_id` with another row. In the pairs we checked, the two rows sit on neighboring tiles of the same well, so they look like one cell that was segmented twice where the tiles overlap. In 49.6% of these groups the rows carry different guides, e.g. `segmentation_id` 2049460 in well A/1 of `ops0031_20250424` has guide `GGGTTCGAGTGGTCTCAAGT` on tile `A/1/024005` and guide `GCCGAAGGAGTACCTACCGT` on tile `A/1/023005`. We set `qc_tile_duplicate` on all of these rows, which makes `qc_fail` true, so they stay in the cell files but are left out of the normalization and the guide-level files.

## Reproduce the files

```bash
python scripts/fetch_altair.py raw
python scripts/build_altair.py phase raw/v2.0.20260724/atlas/cell_features/all_cells_phase.h5ad tiers --abs-max 1e20
python scripts/build_altair.py marker raw/v2.0.20260724/atlas/cell_features/all_cells_fluor_<marker>.h5ad tiers --abs-max 1e20 \
    --phase raw/v2.0.20260724/atlas/cell_features/all_cells_phase.h5ad --phase-selected tiers/phase_cells_n100.h5ad
```
