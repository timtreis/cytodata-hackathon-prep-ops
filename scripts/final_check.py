"""Check every shipped tier and write markers.csv and the SHA-256 list.

Usage: python final_check.py TIERS_DIR
Reports, per file, the shape, the size and max |z| over selected features (cells passing QC and all cells),
and runs mantispy's validator. Reads cell files in backed mode, a block of rows at a time.
"""

import hashlib
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

import mantispy as mt

BLOCK = 50_000


def max_abs(a: ad.AnnData, cols: np.ndarray, rows: np.ndarray | None) -> float:
    m = 0.0
    for lo in range(0, a.n_obs, BLOCK):
        X = np.asarray(a.X[lo : lo + BLOCK])[:, cols]
        if rows is not None:
            X = X[rows[lo : lo + BLOCK]]
        if X.size:
            m = max(m, float(np.nanmax(np.abs(X), initial=0.0)))
    return m


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        while chunk := fh.read(1 << 24):
            h.update(chunk)
    return h.hexdigest()


d = Path(sys.argv[1])
which = sys.argv[2] if len(sys.argv) > 2 else "all"  # "phase", "markers" or "all"
files = (sorted(d.glob("phase_*.h5ad")) if which != "markers" else []) + (sorted((d / "markers").glob("*.h5ad")) if which != "phase" else [])
rows, markers = [], {}
for p in files:
    a = ad.read_h5ad(p, backed="r")
    sel = a.var["selected"].to_numpy(dtype=bool)
    ok = a.obs["qc_pass"].to_numpy(dtype=bool) if "qc_pass" in a.obs else None
    mt.io.validate(a.to_memory() if a.n_obs < 50_000 else a[:5000].to_memory())
    rows.append(
        {
            "file": p.name,
            "n_obs": a.n_obs,
            "n_vars": a.n_vars,
            "selected": int(sel.sum()),
            "gb": p.stat().st_size / 1e9,
            "max_abs_z_selected_qc_pass": max_abs(a, sel, ok) if ok is not None else np.nan,
            "max_abs_z_selected_all": max_abs(a, sel, None),
            "sha256": sha256(p),
        }
    )
    if p.name.endswith("_cells_n100.h5ad") and p.parent.name == "markers":
        m = p.name.removesuffix("_cells_n100.h5ad")
        markers[m] = {"marker": m, "screens": ",".join(sorted(a.obs["batch"].astype(str).unique())), "cells": a.n_obs, "file": p.name, "gb": round(p.stat().st_size / 1e9, 2)}
    print(rows[-1]["file"], rows[-1]["n_obs"], rows[-1]["selected"], f"{rows[-1]['max_abs_z_selected_qc_pass']:.3g}", f"{rows[-1]['max_abs_z_selected_all']:.3g}", flush=True)
report = pd.DataFrame(rows)
report.to_csv(d / f"final_check_{which}.csv", index=False)
if markers:
    pd.DataFrame(markers.values()).to_csv(d / "markers.csv", index=False)
print(report.drop(columns="sha256").to_string(float_format=lambda v: f"{v:.3g}"))
