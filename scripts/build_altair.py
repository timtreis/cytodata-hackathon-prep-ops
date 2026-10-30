"""Build the hackathon tiers of one Altair per-cell feature table with the standard mantispy workflow.

Usage:
    python build_altair.py phase  PHASE.h5ad OUTDIR --abs-max T
    python build_altair.py marker MARKER.h5ad OUTDIR --abs-max T --phase PHASE.h5ad --phase-selected OUTDIR/phase_cells_n100.h5ad

Per screen (one Biohub experiment = one batch): values with |x| > abs-max become NaN, obs gets the shared
cytodata columns and mantispy's schema, cells are flagged (never dropped), and features are robust-z-scored
against that screen's NTC cells that pass QC. Tiers:
    cells_n1000 / cells_n100  at most N cells per guide, chosen by a hash of the cell key (experiment, well, tile, segmentation_id) so the
                              choice is deterministic, nested, and shared by markers imaged in the same cells.
                              X = normalized, layers["raw"] = Biohub's values untouched (phase only).
    guide_x_screen            median per (sgRNA, experiment) over all cells passing QC, n_cells in obs.
    guide                     median of the guide's guide_x_screen rows (median of per-screen medians), n_screens in obs.
Marker mode keeps only the marker's own features and appends the selected phase features of the same cells.
"""

import argparse
import re
import time
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd

import mantispy as mt
from mantispy.pp import DEFAULT_OPERATIONS

try:
    from anndata.io import read_elem
except ImportError:  # anndata < 0.11
    from anndata.experimental import read_elem

TIERS = {"cells_n1000": 1000, "cells_n100": 100}
# segmentation_id repeats across adjacent tiles of a well (the same cell at a tile overlap), so a cell is this 4-tuple.
CELL_KEY = ["experiment", "well", "tile_pheno", "segmentation_id"]
NAN_FLAG = 0.5  # cells / features missing more than this share of values are flagged


def read_obs(f: h5py.File) -> pd.DataFrame:
    g = f["obs"]
    return pd.DataFrame({c: read_elem(g[c]) for c in g if c != g.attrs["_index"]}, index=read_elem(g[g.attrs["_index"]]).astype(str))


def tier_rank(obs: pd.DataFrame) -> np.ndarray:
    """Rank of each cell within its guide by a hash of its CELL_KEY; -1 for unassigned cells."""
    h = pd.util.hash_pandas_object(obs[CELL_KEY].astype({c: str for c in CELL_KEY[:3]}), index=False).to_numpy()
    rank = pd.Series(h).groupby(obs["sgRNA"].to_numpy(), observed=True).rank(method="first").to_numpy() - 1
    return np.where(obs["sgRNA"].isna().to_numpy(), -1, np.nan_to_num(rank, nan=-1)).astype(np.int64)


def screen_rows(obs: pd.DataFrame) -> dict[str, np.ndarray]:
    """Sorted row positions of each screen (Biohub's files hold each screen in a few contiguous blocks)."""
    codes = obs["experiment"].cat.codes.to_numpy()
    return {e: np.flatnonzero(codes == i) for i, e in enumerate(obs["experiment"].cat.categories) if (codes == i).any()}


def read_rows(ds: h5py.Dataset, rows: np.ndarray) -> np.ndarray:
    """Read sorted rows as a few contiguous slices instead of one slow fancy index."""
    cuts = np.flatnonzero(np.diff(rows) != 1) + 1
    return np.concatenate([ds[r[0] : r[-1] + 1] for r in np.split(rows, cuts)])


# The only screen expected to have no usable phase data (2,549 of 3,317 phase features all-NaN in Biohub's table).
EXPECTED_PHASE_EXCLUSIONS = {"ops0174_20260702"}


class NoReference(Exception):
    """A screen has no non-targeting cell passing QC for this modality, so it can't be normalized."""


def to_mantispy(X: np.ndarray, obs: pd.DataFrame, var: pd.DataFrame, abs_max: float) -> ad.AnnData:
    """Shared cytodata columns + mantispy schema; extreme values to NaN; QC flags; normalized X."""
    clean = np.where(np.abs(X) > abs_max, np.nan, X).astype(np.float32)
    obs = obs.copy()
    obs["target_gene"] = obs["gene"]
    obs["guide"] = obs["sgRNA"]
    obs["batch"] = obs["experiment"]
    rc = obs["well"].astype(str).str.extract(r"^([A-Z])/(\d+)/")
    obs["Metadata_Gene"], obs["Metadata_sgRNA"], obs["Metadata_Barcode"] = obs["gene"], obs["sgRNA"], obs["barcode"]
    obs["Metadata_Plate"] = obs["experiment"]
    obs["Metadata_Well"] = (rc[0] + rc[1].str.zfill(2)).astype("category")
    obs["Metadata_ImageID"] = (obs["experiment"].astype(str) + "_" + obs["tile_pheno"].astype(str)).astype("category")
    obs["Metadata_ObjectNumber"] = obs["segmentation_id"].astype("int64")
    obs["Metadata_ObjectType"] = "cell"
    obs["Metadata_Perturbation"] = obs["gene"]
    obs["Metadata_Perturbation_Type"] = "crispr"
    obs["Metadata_CellLine"] = "A549"
    a = ad.AnnData(clean, obs=obs.drop(columns=["model_confidence", "pma_attention", "predicted_class"], errors="ignore"), var=var.copy())
    mt.io.stamp(a, resolution="object", feature_kind="measurement")
    for col, value in {"object": "Cells", "feature": "Area"}.items():  # stamp names no area for non-CellProfiler names
        a.var[col] = a.var[col].astype(object)
        a.var.loc["op_cell_area", col] = value
        a.var[col] = a.var[col].astype("category")  # an all-NaN object column (marker subsets) can't be written to h5ad
    mt.pp.annotate_controls(a, negcon=["NTC"], perturbation_key="Metadata_Gene")
    mt.pp.calculate_qc_metrics(a, max_nan_fraction=NAN_FLAG)
    # Every qc_* flag is True for a problem; mantispy's qc_fail is their union and what mt.pp.filter_cells reads.
    # mantispy's own qc_fail already covers the NaN share and area outliers.
    a.obs["qc_high_nan"] = (a.obs["qc_nan_fraction"] > NAN_FLAG).to_numpy()
    # The same segmentation_id on adjacent tiles: one cell counted twice, often with conflicting guide calls.
    a.obs["qc_tile_duplicate"] = obs.duplicated(["well", "segmentation_id"], keep=False).to_numpy()
    a.obs["qc_fail"] = a.obs["qc_fail"] | a.obs["qc_high_nan"] | a.obs["qc_area_outlier"] | a.obs["qc_tile_duplicate"]
    a.obs["norm_reference"] = (a.obs["Metadata_Control_Type"] == "negcon") & ~a.obs["qc_fail"]
    if not a.obs["norm_reference"].any():
        raise NoReference(f"{obs['experiment'].iloc[0]}: no NTC cell passes QC (median NaN share {np.isnan(clean).mean(1).mean():.2f})")
    mt.pp.normalize(a, method="mad_robustize", by="Metadata_Plate", reference="norm_reference")
    return a


def flag_var(a: ad.AnnData, n_extreme: np.ndarray) -> None:
    raw = a.layers["raw"]
    nan_share = np.isnan(a.X).mean(0)
    a.var["n_extreme"] = n_extreme
    a.var["qc_all_nan"] = np.isnan(raw).all(0)
    a.var["qc_high_nan"] = nan_share > NAN_FLAG
    a.var["qc_extreme"] = n_extreme > 0
    a.var["qc_degenerate"] = a.var["degenerate_scale"].astype(bool) if "degenerate_scale" in a.var else False


def select_on_cells(cells: ad.AnnData) -> tuple[pd.Series, dict]:
    """mantispy's default selection plus drop_outliers (|z| > 500 in any cell), judged on cells_n100 passing QC.

    On guide profiles drop_outliers sees nothing: medians hide the single-cell blow-ups it exists to catch.
    """
    mt.pp.feature_select(cells, operations=(*DEFAULT_OPERATIONS, "drop_outliers"), outlier_cutoff=500.0)
    return cells.var["selected"], cells.uns["mantispy"].get("feature_select")


def rss_gb() -> float:
    with open("/proc/self/status") as fh:
        return next(int(line.split()[1]) for line in fh if line.startswith("VmRSS")) / 1e6


def build(kind: str, src: str, outdir: Path, abs_max: float, phase: str | None, phase_selected: str | None, max_screens: int | None, max_rss_gb: float | None) -> None:
    t0 = time.time()
    outdir.mkdir(parents=True, exist_ok=True)
    # Only [A-Za-z0-9._-] in file names: the CDN reads "+" in a URL as a space (fe2+_ferhonox -> fe2_ferhonox).
    name = "phase" if kind == "phase" else re.sub(r"[^A-Za-z0-9._-]", "", Path(src).stem.removeprefix("all_cells_fluor_"))
    f = h5py.File(src, "r")
    obs, var = read_obs(f), read_elem(f["var"])
    rank = tier_rank(obs)
    cap = TIERS["cells_n100"] if kind == "marker" else max(TIERS.values())  # marker pilot: n100 only
    keep = np.flatnonzero((rank >= 0) & (rank < cap))
    if kind == "marker":
        with h5py.File(phase, "r") as pf:
            phase_var = read_elem(pf["var"]).index
        own = ~var.index.isin(phase_var)  # shape features duplicated from phase are dropped
        sel_phase = ad.read_h5ad(phase_selected, backed="r").var.query("selected").index
        pf = h5py.File(phase, "r")
        pobs = read_obs(pf)
        pscreens = screen_rows(pobs)
    if kind == "marker":
        var_out = pd.concat([var[own].assign(modality="marker"), read_elem(pf["var"]).loc[sel_phase].assign(modality="phase")])
        var_out.index = list(var.index[own]) + [f"phase__{v}" for v in sel_phase]
    else:
        var_out = var.assign(modality="phase")
    n_own = int(own.sum()) if kind == "marker" else var.shape[0]
    # One checkpoint per finished screen, so a crash or a time limit costs one screen; a rerun skips finished screens.
    ckpt = outdir / f"_ckpt_{name}"
    ckpt.mkdir(exist_ok=True)
    screens = list(screen_rows(obs).items())[:max_screens]  # max_screens only for test runs
    for screen, idx in screens:
        done = ckpt / f"{screen}.h5ad"
        if (done.exists() and (ckpt / f"{screen}.gx.h5ad").exists()) or (ckpt / f"{screen}.excluded").exists():
            continue
        # Cells without a well or tile (0.26%, mostly the Cell Painting screens) can't be placed or held by mantispy's
        # schema, so they are left out. Ranks stay as computed over all cells, so finished checkpoints stay valid.
        idx = idx[~obs["well"].iloc[idx].astype(str).str.startswith("NAN").to_numpy()]
        X = read_rows(f["X"], idx)
        try:
            a = to_mantispy(X, obs.iloc[idx], var, abs_max)
        except NoReference as e:
            if kind == "phase" and screen not in EXPECTED_PHASE_EXCLUSIONS:
                raise  # only ops0174 may be skipped in phase; anything else needs a look
            (ckpt / f"{screen}.excluded").write_text(str(e))
            print(f"{name} {screen}: excluded, {e}", flush=True)
            continue
        own_cols = slice(None) if kind == "phase" else np.flatnonzero(own)
        gx = mt.tl.aggregate(a[~a.obs["qc_fail"].to_numpy(), own_cols], by=("Metadata_Plate", "Metadata_sgRNA"), min_cells=1)
        rows = np.flatnonzero(np.isin(idx, keep))
        # copy: a view would keep this screen's whole object alive through part.obs (OOM after ~30 screens)
        part, raw = a[rows].copy(), X[rows].astype(np.float32)
        if kind == "marker":
            # Phase features of the same cells: normalize the screen's phase rows the same way, then look the cells up.
            pidx = pscreens[screen]
            pidx = pidx[~pobs["well"].iloc[pidx].astype(str).str.startswith("NAN").to_numpy()]  # as for the marker rows
            pX_raw = read_rows(pf["X"], pidx)
            pvar = read_elem(pf["var"])
            try:
                pa = to_mantispy(pX_raw, pobs.iloc[pidx], pvar, abs_max)
            except NoReference:
                pa = None  # phase can't be normalized in this screen (ops0174): its phase columns stay NaN
            key = CELL_KEY[1:]
            pos = pd.Series(np.arange(len(pidx)), index=pd.MultiIndex.from_frame(pobs.iloc[pidx][key].astype(str)))
            prow = pos.reindex(pd.MultiIndex.from_frame(part.obs[key].astype(str))).to_numpy()
            assert not np.isnan(prow).any(), f"{screen}: marker cells missing from phase"
            prow = prow.astype(np.int64)
            pcols = pvar.index.get_indexer(sel_phase)
            pX = np.full((len(prow), len(sel_phase)), np.nan, np.float32) if pa is None else np.asarray(pa.X[prow][:, pcols], dtype=np.float32)
            pR = pX_raw[prow][:, pcols]
            part = ad.AnnData(np.hstack([part.X[:, own], pX]), obs=part.obs)
            raw = np.hstack([raw[:, own], pR]).astype(np.float32)
            del pa, pX_raw
        out = ad.AnnData(np.asarray(part.X, dtype=np.float32), obs=part.obs.assign(_row=idx[rows]), var=pd.DataFrame(index=var_out.index))
        out.layers["raw"] = raw
        out.uns["n_extreme"] = (np.abs(X) > abs_max).sum(0)
        out.uns["degenerate"] = a.var["degenerate_scale"].to_numpy(dtype=bool)
        out.uns["mantispy"] = dict(a.uns["mantispy"])
        del a, X, part, raw
        gx.write_h5ad(ckpt / f"{screen}.gx.tmp")
        (ckpt / f"{screen}.gx.tmp").rename(ckpt / f"{screen}.gx.h5ad")
        out.write_h5ad(ckpt / f"{screen}.tmp")
        (ckpt / f"{screen}.tmp").rename(done)
        del out, gx
        mem = rss_gb()
        print(f"{name} {screen}: {len(idx):,} cells, kept {len(rows):,}, rss {mem:.1f} GB ({time.time() - t0:.0f}s)", flush=True)
        if max_rss_gb and mem > max_rss_gb:
            raise SystemExit(f"RSS {mem:.1f} GB above --max-rss-gb {max_rss_gb} after {screen}; finished screens are checkpointed in {ckpt}, rerun to resume")
    excluded = {s: (ckpt / f"{s}.excluded").read_text() for s, _ in screens if (ckpt / f"{s}.excluded").exists()}
    screens = [(s, i) for s, i in screens if s not in excluded]
    parts = [ad.read_h5ad(ckpt / f"{s}.h5ad") for s, _ in screens]
    gxs = [ad.read_h5ad(ckpt / f"{s}.gx.h5ad") for s, _ in screens]
    n_extreme = np.sum([p.uns["n_extreme"] for p in parts], axis=0)
    degenerate = np.any([p.uns["degenerate"] for p in parts], axis=0)
    uns = {"mantispy": dict(parts[0].uns["mantispy"])}
    big = ad.AnnData(np.concatenate([p.X for p in parts]), obs=pd.concat([p.obs for p in parts]), var=var_out, uns=uns)
    big.layers["raw"] = np.concatenate([p.layers["raw"] for p in parts])
    del parts
    big.obs["tier_rank"] = rank[big.obs.pop("_row").to_numpy()]
    big.obs = big.obs.drop(columns=["qc_no_well"], errors="ignore")  # early checkpoints carry it, always False
    if kind == "marker":
        n_extreme = np.r_[n_extreme[own], np.zeros(len(sel_phase), dtype=np.int64)]
        degenerate = np.r_[degenerate[own], np.zeros(len(sel_phase), dtype=bool)]  # phase flags live in the phase tiers
    big.var["degenerate_scale"] = degenerate  # True if the scale had no spread in any screen
    flag_var(big, n_extreme)
    if kind == "marker":
        del big.layers["raw"]  # raw marker values stay in Biohub's bucket, reachable through the cell key
    # guide_x_screen: median over all of a guide's cells passing QC in that screen.
    # guide: median of the guide's guide_x_screen rows, so every screen weighs the same.
    own_var = big.var.iloc[:n_own]
    gx = ad.concat(gxs, merge="same", uns_merge="first")
    gx.obs_names = (gx.obs["Metadata_Plate"].astype(str) + "_" + gx.obs["Metadata_sgRNA"].astype(str)).to_numpy()
    guide = mt.tl.aggregate(gx, by=("Metadata_sgRNA",), min_cells=1)
    guide.obs["n_screens"] = gx.obs.groupby("Metadata_sgRNA", observed=True).size().reindex(guide.obs["Metadata_sgRNA"]).to_numpy()
    gx.obs["batch"] = gx.obs["Metadata_Plate"]
    guide.var = own_var.loc[guide.var_names]
    selected, counts = select_on_cells(big[(big.obs["tier_rank"] < TIERS["cells_n100"]).to_numpy() & ~big.obs["qc_fail"].to_numpy(), :n_own].copy())
    # Marker files' appended phase columns are the phase selection already.
    big.var["selected"] = selected.reindex(big.var_names).fillna(big.var["modality"].eq("phase") & (kind == "marker")).astype(bool).to_numpy()
    big.uns["mantispy"]["feature_select"] = counts
    for tier, n in TIERS.items():
        if kind == "marker" and tier != "cells_n100":
            continue  # option C pilot builds n100 only
        sub = big if n == max(TIERS.values()) else big[(big.obs["tier_rank"] < n).to_numpy()].copy()
        sub.uns["cytodata"] = {"dataset": "altair", "modality": name, "tier": tier, "abs_max": abs_max, "excluded_screens": excluded}
        sub.write_h5ad(outdir / f"{name}_{tier}.h5ad")
        print(f"wrote {name}_{tier}: {sub.shape} ({time.time() - t0:.0f}s)", flush=True)
    own_var = big.var.iloc[:n_own]
    del big, sub
    for tier, t in {"guide_x_screen": gx, "guide": guide}.items():
        t.var = own_var.loc[t.var_names]
        t.uns["mantispy"]["feature_select"] = counts
        t.obs["guide"] = t.obs["Metadata_sgRNA"]
        t.obs["n_cells"] = t.obs["Metadata_CellCount"]
        t.uns["cytodata"] = {"dataset": "altair", "modality": name, "tier": tier, "abs_max": abs_max, "excluded_screens": excluded}
        t.write_h5ad(outdir / f"{name}_{tier}.h5ad")
        print(f"wrote {name}_{tier}: {t.shape} ({time.time() - t0:.0f}s)", flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("kind", choices=["phase", "marker"])
    p.add_argument("src")
    p.add_argument("outdir", type=Path)
    p.add_argument("--abs-max", type=float, required=True)
    p.add_argument("--phase")
    p.add_argument("--phase-selected")
    p.add_argument("--max-screens", type=int)
    p.add_argument("--max-rss-gb", type=float, help="abort after a screen if resident memory exceeds this")
    a = p.parse_args()
    build(a.kind, a.src, a.outdir, a.abs_max, a.phase, a.phase_selected, a.max_screens, a.max_rss_gb)
