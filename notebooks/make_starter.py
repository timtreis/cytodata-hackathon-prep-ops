"""Write notebooks/01_starter.ipynb from the cells below (run once, then execute the notebook)."""

import json
from pathlib import Path

CELLS = [
    ("md", """# Altair optical pooled screen: starter notebook

This notebook downloads the guide-level phase profiles and up to 100 cells per guide, scores how well each knockout separates from the non-targeting controls with mantispy, draws a map of the knockouts, and joins Biohub's CROP-seq embedding by guide sequence. It runs on a laptop in a few minutes and needs about 4 GB of memory."""),
    ("code", """import os
import time
from pathlib import Path

import anndata as ad
import mantispy as mt
import numpy as np
import pandas as pd
import pooch
import scanpy as sc

T0 = time.time()
DATA = Path(os.environ.get("ALTAIR_DATA", "data"))
URL = "https://exampledata.scverse.org/sbi2-cytodata/ops"
# SHA-256 from manifest.csv; pooch downloads each file once and checks it.
FILES = {
    "phase_guide.h5ad": "8cffb86078d6ed9839821b9176c4ef744f1b9b35386b7710794e607188d7b985",
    "phase_guide_x_screen.h5ad": "3f1ae51ff5958b62c16c539160ae0cc5f51c91b5445659cfdf1dc6142f36e976",
    "phase_cells_n100.h5ad": "0db0b620ef927c56e8b1293c02b443abd19974eb4518d0f0c89da8005339dec1",
    "screen_map.csv": "9dbd0486fbb685403ec26eff87dc12bb06e7ac773a2a0428ad3ff3482ad3d509",
    "unexpressed_genes.csv": "7ca185b0bc4b6f4c0f50f517846e9bd783ca0f61ad19fde04ac6ea4367b16acc",
    "markers_guide.zip": "3dd193c7c24890cdddc56036f24c55ade987e7069fd7ff3890ed84f02486fcd9",
}
for f, sha in FILES.items():
    pooch.retrieve(f"{URL}/{f}", known_hash=f"sha256:{sha}", fname=f, path=DATA, progressbar=False)
{f: f"{(DATA / f).stat().st_size / 1e9:.2f} GB" for f in FILES}"""),
    ("md", """## Guide profiles

`phase_guide.h5ad` holds one robust z-scored phase profile per guide, the median of the guide's per-screen medians. We keep the features that mantispy's feature selection kept, and we clip the z-scores at 50 because a few features have a very small spread among the control cells."""),
    ("code", """guides = ad.read_h5ad(DATA / "phase_guide.h5ad")
guides = guides[:, guides.var["selected"].values].copy()
guides.X = np.clip(np.nan_to_num(guides.X), -50, 50)
print(guides)
guides.obs[["guide", "target_gene", "n_cells", "n_screens"]].head()"""),
    ("md", """## Which knockouts have a phenotype

`mt.tl.map` with `mode="activity"` asks whether the four guides of a gene retrieve each other better than they retrieve the control guides, and it reports a mean average precision (`score`) with a permutation q-value per gene. All guide profiles are pooled across screens, so we put them in one batch.

The choice of control matters. Every targeting guide cuts the DNA and the non-targeting guides don't, so against `NTC` almost every knockout separates, including knockouts of genes that A549 cells don't express. A fairer control is the knockouts of the 74 single-copy genes that A549 doesn't express (`unexpressed_genes.csv`), because those guides cut the DNA without removing a gene the cells use. We run both."""),
    ("code", """unexpressed = pd.read_csv(DATA / "unexpressed_genes.csv")
cutting_control = set(unexpressed.query("not multicopy")["gene"])
guides.obs["Metadata_Plate"] = "all"
guides.obs["cutting_control"] = guides.obs["target_gene"].isin(cutting_control).values
t = time.time()
mt.tl.map(guides, mode="activity", label="Metadata_Gene", batch="Metadata_Plate", reference="negcon", null_size=20_000, key_added="map_ntc")
mt.tl.map(guides, mode="activity", label="Metadata_Gene", batch="Metadata_Plate", reference="cutting_control", null_size=20_000, key_added="map_cut")
vs_ntc, activity = guides.uns["mantispy"]["map_ntc"], guides.uns["mantispy"]["map_cut"]
activity = activity[~activity["Metadata_Gene"].isin(cutting_control | {"NTC"})]
print(f"{time.time() - t:.0f} s")
print(f"against NTC: {vs_ntc['is_significant'].mean():.1%} of {len(vs_ntc)} knockouts separate (q < 0.05)")
print(f"against the cutting control: {activity['is_significant'].mean():.1%} of {len(activity)} knockouts separate (q < 0.05)")
multicopy = set(unexpressed.query("multicopy")["gene"])
rank = vs_ntc.sort_values("score", ascending=False).reset_index(drop=True)
print("rank against NTC of unexpressed multicopy genes (their guides cut many copies):", {g: int(i) + 1 for i, g in rank["Metadata_Gene"].items() if g in multicopy})
activity.sort_values("score", ascending=False).head(10)"""),
    ("md", """## A map of the knockouts

A UMAP of the guide profiles, with the non-targeting guides in grey. Genes in the same complex tend to sit together."""),
    ("code", """sc.pp.pca(guides, n_comps=50)
sc.pp.neighbors(guides, n_neighbors=15, transformer="sklearn")
sc.tl.umap(guides, random_state=0)
guides.obs["is_ntc"] = (guides.obs["target_gene"] == "NTC").map({True: "NTC", False: "knockout"})
sc.pl.umap(guides, color="is_ntc", palette={"NTC": "lightgrey", "knockout": "tab:blue"})"""),
    ("md", """## Single cells

`phase_cells_n100.h5ad` holds up to 100 cells per guide. Reading it in backed mode and loading only the selected features and the cells that pass QC keeps the memory use small."""),
    ("code", """import gc

import h5py
from anndata.io import read_elem

# Read obs, var and the selected columns of X straight from the file. A backed AnnData would load
# layers["raw"] (about 5.5 GB) into memory even though we never use it.
with h5py.File(DATA / "phase_cells_n100.h5ad", "r") as f:
    obs, var = read_elem(f["obs"]), read_elem(f["var"])
    keep_cells = obs["qc_pass"].values
    keep_features = np.flatnonzero(var["selected"].values)
    blocks = [f["X"][lo : lo + 10_000][keep_cells[lo : lo + 10_000]][:, keep_features] for lo in range(0, len(obs), 10_000)]
cells = ad.AnnData(np.clip(np.nan_to_num(np.vstack(blocks)), -50, 50), obs=obs[keep_cells], var=var.iloc[keep_features])
del blocks
gc.collect()
print(cells.shape, f"{cells.X.nbytes / 1e9:.2f} GB")
cells.obs.groupby("target_gene", observed=True).size().describe()"""),
    ("md", """## CROP-seq of the same library

Biohub's Zenodo record has a CROP-seq embedding for each guide. Its `sgRNA` column is the guide sequence, so it joins to our `guide` column. The CROP-seq cells are not the imaged cells, so the two modalities only line up at the guide level."""),
    ("code", """import zipfile

zpath = pooch.retrieve(
    "https://zenodo.org/api/records/21970895/files/altair_v2.zip/content",
    known_hash="md5:6b55793b6532138bf52d9a39035c8d62",  # from the Zenodo record
    fname="altair_v2.zip", path=DATA, progressbar=False,
)
with zipfile.ZipFile(zpath) as z:
    z.extract("figures/figure_5/cropseq_svaeplus.h5ad", DATA)
rna = ad.read_h5ad(DATA / "figures/figure_5/cropseq_svaeplus.h5ad")
rna_emb = pd.DataFrame(rna.obsm["emb_guide"], index=rna.obs["sgRNA"].values)
shared = guides.obs["guide"].astype(str).isin(rna_emb.index)
print(f"{shared.sum()} of {guides.n_obs} guides have a CROP-seq embedding")
joint = guides[shared.values].copy()
joint.obsm["X_rna"] = rna_emb.loc[joint.obs["guide"].astype(str)].values
joint"""),
    ("md", """## How much of the morphology neighbourhood is also an RNA neighbourhood

For each knockout guide we take its 10 nearest guides by phase profile and its 10 nearest guides by CROP-seq embedding (cosine distance), and we count how many of the 10 are shared. We compare the mean overlap with the overlap after shuffling which RNA embedding belongs to which guide."""),
    ("code", """from sklearn.neighbors import NearestNeighbors

ko = joint[joint.obs["target_gene"] != "NTC"]
k = 10


def knn(x):
    return NearestNeighbors(n_neighbors=k + 1, metric="cosine").fit(x).kneighbors(x, return_distance=False)[:, 1:]


def overlap(a, b):
    return np.mean([len(set(x) & set(y)) / k for x, y in zip(a, b)])


phase_nn = knn(ko.X)
rna = ko.obsm["X_rna"]
observed = overlap(phase_nn, knn(rna))
rng = np.random.default_rng(0)
shuffled = [overlap(phase_nn, knn(rna[rng.permutation(len(rna))])) for _ in range(5)]
print(f"{observed:.1%} of a guide's 10 phase neighbours are also among its 10 RNA neighbours, against {np.mean(shuffled):.1%} after shuffling")"""),
    ("md", """## What a fluorescent marker adds

The marker files have the same tiers for each of the 59 markers. Phase alone already separates most knockouts, so a marker adds a few, specific ones. Here we score the COP-II coat marker SEC23A, which labels where cargo leaves the endoplasmic reticulum, against the same cutting control, and we list the knockouts that separate with SEC23A but not with phase. They should be genes of the secretory pathway."""),
    ("code", """with zipfile.ZipFile(DATA / "markers_guide.zip") as z:
    z.extract("er_golgi_cop-ii_sec23a_guide.h5ad", DATA)
marker = ad.read_h5ad(DATA / "er_golgi_cop-ii_sec23a_guide.h5ad")
marker = marker[:, marker.var["selected"].values].copy()
marker.X = np.clip(np.nan_to_num(marker.X), -50, 50)
marker.obs["Metadata_Plate"] = "all"
marker.obs["cutting_control"] = marker.obs["target_gene"].isin(cutting_control).values
mt.tl.map(marker, mode="activity", label="Metadata_Gene", batch="Metadata_Plate", reference="cutting_control", null_size=20_000)
marker_map = marker.uns["mantispy"]["map"]
marker_hits = set(marker_map[~marker_map["Metadata_Gene"].isin(cutting_control | {"NTC"})].query("is_significant")["Metadata_Gene"])
phase_hits = set(activity.query("is_significant")["Metadata_Gene"])
only_marker = sorted(marker_hits - phase_hits)
print(f"{len(marker_hits)} knockouts separate with SEC23A, {len(phase_hits)} with phase, {len(only_marker)} only with SEC23A")
only_marker[:30]"""),
    ("code", """print(f"total runtime {time.time() - T0:.0f} s")"""),
]


def cell(i, kind, src):
    lines = src.splitlines(keepends=True)
    if kind == "md":
        return {"cell_type": "markdown", "id": f"c{i:02d}", "metadata": {}, "source": lines}
    return {"cell_type": "code", "id": f"c{i:02d}", "execution_count": None, "metadata": {}, "outputs": [], "source": lines}


nb = {
    "cells": [cell(i, k, s) for i, (k, s) in enumerate(CELLS)],
    "metadata": {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"}},
    "nbformat": 4,
    "nbformat_minor": 5,
}
Path(__file__).with_name("01_starter.ipynb").write_text(json.dumps(nb, indent=1))
