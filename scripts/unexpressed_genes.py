"""List the Altair target genes A549 doesn't express, from the Human Protein Atlas cell line RNA table.

Usage: python unexpressed_genes.py RNA_CELLLINE_ZIP PERTURBATION_LIBRARY_CSV OUT_CSV
(rna_celline.tsv.zip from https://www.proteinatlas.org/download/tsv/rna_celline.tsv.zip, version 25.1)
"""

import sys
import zipfile

import pandas as pd

# Families of near-identical paralogs, where one guide matches many loci. Named by HGNC family, not by alignment.
MULTICOPY = r"^(?:GOLGA6|GOLGA8|USP17L|SPATA31A|PRAMEF|TBC1D3|FOXD4L|KRTAP|PPIAL4|TRIM43|H[1-4]C|PLGLB|GGTLC|SPRR2|LCE3)"

hpa_zip, library_csv, out = sys.argv[1:4]
with zipfile.ZipFile(hpa_zip) as z:
    chunks = pd.read_csv(z.open(z.namelist()[0]), sep="\t", chunksize=2_000_000)
    ntpm = pd.concat(c[c["Cell line"] == "A-549"] for c in chunks).groupby("Gene name")["nTPM"].max()
lib = pd.read_csv(library_csv)
genes = sorted(set(lib.loc[lib["role"] == "targeting", "gene_symbol"]) & set(ntpm.index))
u = pd.DataFrame({"gene": genes, "ntpm": ntpm[genes].to_numpy()}).query("ntpm < 1")  # HPA's "not detected"
u["multicopy"] = u["gene"].str.contains(MULTICOPY)
u.to_csv(out, index=False)
print(f"{len(u)} of {len(genes)} target genes not detected in A549, {int(u['multicopy'].sum())} multicopy")
