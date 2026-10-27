"""Cap a huge per-cell h5ad at N cells per sgRNA, streaming X so it never loads whole.

Usage: python subsample_phase.py IN.h5ad OUT.h5ad [--n 1000] [--seed 0]
       python subsample_phase.py --selftest
"""

import argparse
import tempfile

import anndata as ad
import h5py
import numpy as np
import pandas as pd

try:
    from anndata.io import read_elem
except ImportError:  # anndata < 0.11
    from anndata.experimental import read_elem

BLOCK = 100_000  # rows per X read; 3,317 float32 features -> ~1.3 GB


def pick(guides: pd.Categorical, n: int, seed: int) -> np.ndarray:
    """Sorted row positions keeping at most n random cells per guide; unassigned (NaN) dropped."""
    codes = np.asarray(guides.codes)
    order = np.random.default_rng(seed).permutation(len(codes))
    order = order[codes[order] >= 0]
    rank = pd.Series(codes[order]).groupby(codes[order]).cumcount().to_numpy()
    return np.sort(order[rank < n])


def subsample(src: str, dst: str, n: int, seed: int) -> None:
    with h5py.File(src, "r") as f:
        obs_g = f["obs"]
        idx = pick(read_elem(obs_g["sgRNA"]), n, seed)
        print(f"{len(idx):,} of {f['X'].shape[0]:,} cells kept", flush=True)
        obs = pd.DataFrame(
            {c: read_elem(obs_g[c])[idx] for c in obs_g if c != obs_g.attrs["_index"]},
            index=pd.Index(read_elem(obs_g[obs_g.attrs["_index"]])[idx].astype(str), name=None),
        )
        var = read_elem(f["var"])
        uns = read_elem(f["uns"]) if "uns" in f else {}
        X = f["X"]
        out = np.empty((len(idx), X.shape[1]), dtype=X.dtype)
        w = 0
        for lo in range(0, X.shape[0], BLOCK):
            sel = idx[(idx >= lo) & (idx < lo + BLOCK)]
            if len(sel):
                out[w : w + len(sel)] = X[lo : lo + BLOCK][sel - lo]
                w += len(sel)
            if lo // BLOCK % 50 == 0:
                print(f"read {lo:,} rows", flush=True)
        assert w == len(idx)
    ad.AnnData(out, obs=obs, var=var, uns=uns).write_h5ad(dst)
    print("done", flush=True)


def selftest():
    tmp = tempfile.mkdtemp() + "/t"
    guides = pd.Categorical(["a"] * 5 + ["b"] * 2 + [None] * 3)
    idx = pick(guides, 3, 0)
    kept = pd.Series(np.asarray(guides)[idx]).value_counts()
    assert kept.to_dict() == {"a": 3, "b": 2}, kept
    assert (np.diff(idx) > 0).all()
    a = ad.AnnData(np.arange(40, dtype="f4").reshape(10, 4), obs=pd.DataFrame({"sgRNA": guides}, index=[f"c{i}" for i in range(10)]))
    a.write_h5ad(f"{tmp}_in.h5ad")
    global BLOCK
    BLOCK = 3  # force multiple blocks
    subsample(f"{tmp}_in.h5ad", f"{tmp}_out.h5ad", 3, 0)
    b = ad.read_h5ad(f"{tmp}_out.h5ad")
    np.testing.assert_array_equal(b.X, a[b.obs_names].X)
    print("selftest ok")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("src", nargs="?")
    p.add_argument("dst", nargs="?")
    p.add_argument("--n", type=int, default=1000)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--selftest", action="store_true")
    a = p.parse_args()
    selftest() if a.selftest else subsample(a.src, a.dst, a.n, a.seed)
