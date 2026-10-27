"""Mirror the per-cell feature tables of Biohub's Altair OPS release from the public bucket.

Usage: python3 fetch_altair.py OUTDIR [--dry-run]

Fetches v2 atlas tables + all per-cell feature h5ads (fluorescence + 870 GB phase), the v1
phase h5ad and small atlas tables, CROP-seq, and per-screen metadata/cell index. No images.
Files are fetched as 1 GB byte ranges in parallel; rerun the same command to resume.
"""

import os
import re
import sys
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

BUCKET = "https://ops-explorer-public.s3.amazonaws.com/"
ROOT = "leonetti_ops/ops_data_portal_submission/"
V1, V2 = ROOT + "v1.0.20260521/", ROOT + "v2.0.20260724/"
CHUNK = 1 << 30
# Listing all of ROOT walks millions of zarr chunks, so list only what we take.
PREFIXES = [ROOT + "CROPseq/", V2 + "experimental_metadata.yaml", V2 + "datasets/"] + [
    v + d for v in (V1, V2) for d in ("atlas/", "atlas/cell_features/")
]


def ls(prefix, delimiter=None):
    tok = None
    while True:
        q = {"list-type": "2", "prefix": prefix} | ({"delimiter": delimiter} if delimiter else {})
        if tok:
            q["continuation-token"] = tok
        xml = urllib.request.urlopen(BUCKET + "?" + urllib.parse.urlencode(q)).read().decode()
        yield from zip(re.findall(r"<Key>([^<]+)</Key>", xml), map(int, re.findall(r"<Size>(\d+)</Size>", xml)))
        m = re.search(r"<NextContinuationToken>([^<]+)", xml)
        if not m:
            return
        tok = m.group(1)


def listing():
    for p in PREFIXES:
        if p == V2 + "datasets/":
            # ponytail: three small listings per screen instead of walking each zarr
            for screen in sorted(_screens(p)):
                yield from ls(f"{p}{screen}/cell_data.parquet")
                yield from ls(f"{p}{screen}/metadata/")
                yield from ls(f"{p}{screen}/{screen}.zarr/zarr.json")
        else:
            yield from ls(p, "/")


def _screens(prefix):
    xml = urllib.request.urlopen(BUCKET + "?" + urllib.parse.urlencode({"list-type": "2", "prefix": prefix, "delimiter": "/"})).read().decode()
    return {s.removeprefix(prefix).strip("/") for s in re.findall(r"<Prefix>([^<]+)</Prefix>", xml) if s != prefix}


def wanted(key):
    if key.endswith("/") or ".zarr/" in key.replace(".zarr/zarr.json", ""):
        return False
    if key.startswith(V1 + "atlas/"):
        rel = key.removeprefix(V1 + "atlas/")
        return "/" not in rel or rel == "cell_features/all_cells_phase.h5ad"  # v1 fluor is superseded by v2
    return True


def fetch_range(key, size, dst, i):
    done = Path(f"{dst}.part.done/{i}")
    if done.exists():
        return 0
    lo, hi = i * CHUNK, min(size, (i + 1) * CHUNK) - 1
    req = urllib.request.Request(BUCKET + urllib.parse.quote(key), headers={"Range": f"bytes={lo}-{hi}"})
    for attempt in range(5):
        try:
            with urllib.request.urlopen(req, timeout=300) as r:
                data = r.read()
            if len(data) != hi - lo + 1:
                raise OSError(f"short read {len(data)}")
            fd = os.open(f"{dst}.part", os.O_WRONLY)
            try:
                os.pwrite(fd, data, lo)
            finally:
                os.close(fd)
            done.touch()
            return len(data)
        except OSError as e:
            print(f"retry {attempt} {key}[{i}]: {e}", flush=True)
    raise RuntimeError(f"failed: {key}[{i}]")


def main():
    outdir, dry = Path(sys.argv[1]), "--dry-run" in sys.argv
    items = [kv for kv in listing() if wanted(kv[0])]
    print(f"{len(items)} files, {sum(s for _, s in items) / 1e9:.1f} GB", flush=True)
    if dry:
        for k, s in items:
            print(f"{s / 1e9:9.3f} GB  {k.removeprefix(ROOT)}")
        return
    jobs, finals = [], []
    for key, size in items:
        dst = outdir / key.removeprefix(ROOT)
        if dst.exists() and dst.stat().st_size == size:
            continue
        dst.parent.mkdir(parents=True, exist_ok=True)
        Path(f"{dst}.part.done").mkdir(exist_ok=True)
        with open(f"{dst}.part", "ab") as f:
            f.truncate(size)
        n = max(1, -(-size // CHUNK))
        jobs += [(key, size, dst, i) for i in range(n)]
        finals.append((dst, n))
    print(f"{len(jobs)} ranges to check", flush=True)
    got = 0
    with ThreadPoolExecutor(16) as pool:
        for j, n in enumerate(pool.map(lambda a: fetch_range(*a), jobs), 1):
            got += n
            if j % 50 == 0:
                print(f"{j}/{len(jobs)} ranges, {got / 1e9:.1f} GB new", flush=True)
    for dst, n in finals:
        marks = Path(f"{dst}.part.done")
        assert len(list(marks.iterdir())) == n, dst
        os.rename(f"{dst}.part", dst)
        for m in marks.iterdir():
            m.unlink()
        marks.rmdir()
    print("done", flush=True)


if __name__ == "__main__":
    main()
