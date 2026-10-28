"""Package the built tiers for the CDN: gzip every h5ad, shard cells_n1000 by screen, zip the marker guide tiers.

Usage: python package.py TIERS_DIR SHIP_DIR
CloudFront won't serve objects above ~30 GB, so no shipped file may exceed MAX_GB. Writes SHA256 lines to
SHIP_DIR/manifest.csv and the shard table to SHIP_DIR/n1000_shards.csv. Finished outputs are skipped on rerun.
"""

import hashlib
import shutil
import sys
import zipfile
from pathlib import Path

import anndata as ad
import pandas as pd

MAX_GB = 25
SHARD_TARGET_GB = 22
RATIO = 0.7  # gzip level 4 gave 0.64 on these files; plan with a margin


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        while chunk := fh.read(1 << 24):
            h.update(chunk)
    return h.hexdigest()


def gz(src: Path, dst: Path) -> None:
    if not dst.exists():
        a = ad.read_h5ad(src)
        if "guide_x_screen" in src.name:  # rows from per-screen aggregates restart their numbering
            a.obs_names = (a.obs["Metadata_Plate"].astype(str) + "_" + a.obs["Metadata_sgRNA"].astype(str)).to_numpy()
        assert a.obs_names.is_unique, src
        a.write_h5ad(dst.with_suffix(".tmp"), compression="gzip", compression_opts=4)
        dst.with_suffix(".tmp").rename(dst)
    print(f"{dst.name}: {dst.stat().st_size / 1e9:.2f} GB", flush=True)


tiers, ship = Path(sys.argv[1]), Path(sys.argv[2])
ship.mkdir(parents=True, exist_ok=True)

for name in ["phase_guide", "phase_guide_x_screen", "phase_cells_n100"]:
    gz(tiers / f"{name}.h5ad", ship / f"{name}.h5ad")

# cells_n1000: whole screens per shard, packed in screen order up to SHARD_TARGET_GB of compressed size.
src = tiers / "phase_cells_n1000.h5ad"
big = ad.read_h5ad(src, backed="r")
per_row = src.stat().st_size * RATIO / big.n_obs
batch = big.obs["batch"].astype(str).to_numpy()
order = list(dict.fromkeys(batch))  # screens in file order; each screen's rows are contiguous
counts = pd.Series(batch).value_counts()
shards, cur = [], []
for s in order:
    if cur and (sum(counts[c] for c in cur) + counts[s]) * per_row > SHARD_TARGET_GB * 1e9:
        shards.append(cur)
        cur = []
    cur.append(s)
shards.append(cur)
rows = []
for i, screens in enumerate(shards, 1):
    dst = ship / f"phase_cells_n1000_part{i:02d}.h5ad"
    if not dst.exists():
        big[pd.Series(batch).isin(screens).to_numpy()].to_memory().write_h5ad(dst.with_suffix(".tmp"), compression="gzip", compression_opts=4)
        dst.with_suffix(".tmp").rename(dst)
    rows.append({"shard": dst.name, "screens": ",".join(screens), "cells": int(sum(counts[s] for s in screens)), "gb": round(dst.stat().st_size / 1e9, 2), "sha256": sha256(dst)})
    print(rows[-1]["shard"], rows[-1]["cells"], rows[-1]["gb"], "GB", flush=True)
pd.DataFrame(rows).to_csv(ship / "n1000_shards.csv", index=False)

# Markers: one gzipped n100 file each; all guide-level files in one zip (the h5ads inside are gzipped already).
mdir = ship / "markers_tmp"
mdir.mkdir(exist_ok=True)
for p in sorted((tiers / "markers").glob("*_cells_n100.h5ad")):
    gz(p, ship / p.name)
guide_files = sorted((tiers / "markers").glob("*_guide*.h5ad"))
for p in guide_files:
    gz(p, mdir / p.name)
if guide_files:
    with zipfile.ZipFile(ship / "markers_guide.zip", "w", compression=zipfile.ZIP_STORED) as z:
        for p in sorted(mdir.glob("*.h5ad")):
            z.write(p, p.name)
    shutil.rmtree(mdir)
for csv in ["screen_map.csv", "unexpressed_genes.csv", "markers.csv"]:
    if (tiers / csv).exists():
        shutil.copy(tiers / csv, ship / csv)

manifest = pd.DataFrame([{"file": p.name, "gb": p.stat().st_size / 1e9, "sha256": sha256(p)} for p in sorted(ship.iterdir()) if p.is_file() and p.suffix != ".tmp" and p.name != "manifest.csv"])
too_big = manifest[manifest["gb"] > MAX_GB]
assert too_big.empty, f"files above {MAX_GB} GB: {too_big['file'].tolist()}"
manifest.to_csv(ship / "manifest.csv", index=False)
print(manifest.to_string(float_format=lambda v: f"{v:.2f}"))
