import nbformat as nbf

nb = nbf.v4.new_notebook()
cells = []

cells.append(nbf.v4.new_markdown_cell(
"""# Phase 2 — Quality Control (QC)

Goal: understand the incoming AnnData object and apply a defensible,
minimal cell/gene filtering step before any normalization or transformation
(handled separately in `03_preprocessing.ipynb`)."""
))

cells.append(nbf.v4.new_code_cell(
"""import pertpy as pt
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt

adata = pt.dt.kang_2018()
adata"""
))

cells.append(nbf.v4.new_markdown_cell(
"""## Understanding the input object

This dataset arrives already processed via `pertpy` (it is not raw FASTQ
data). Before applying any filtering or transformation, we inspect what is
actually stored in the AnnData object, rather than assuming a standard
raw-counts workflow applies unmodified."""
))

cells.append(nbf.v4.new_code_cell(
"""print("adata.X dtype:", adata.X.dtype)
print("adata.X min/max:", adata.X.min(), adata.X.max())
print()
print("Layers:", list(adata.layers.keys()))
print()
print("obs columns:", adata.obs.columns.tolist())
print()
print("obsm keys:", list(adata.obsm.keys()))
print()
print("adata.raw:", adata.raw)"""
))

cells.append(nbf.v4.new_markdown_cell(
"""## Inspecting mitochondrial content

Mitochondrial percentage is a standard scRNA-seq QC metric used to flag
dying or stressed cells. Before relying on it, we check whether
mitochondrial genes are present in this dataset at all."""
))

cells.append(nbf.v4.new_code_cell(
"""adata.var["mt"] = adata.var["name"].str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
adata.obs[["n_genes_by_counts", "total_counts", "pct_counts_mt"]].describe()"""
))

cells.append(nbf.v4.new_code_cell(
"""adata.var["name"][adata.var["name"].str.contains("MT", case=False, na=False)].head(20)"""
))

cells.append(nbf.v4.new_markdown_cell(
"""No genes matching the `MT-` prefix are present in this dataset, and no
gene names containing "MT" resemble mitochondrial genes either. This
processed release of the Kang dataset appears to have had mitochondrial
genes removed upstream, before distribution via pertpy.

**Mitochondrial QC is therefore not applicable here** and is documented as
a limitation/property of this processed dataset, rather than a step we can
perform. We rely on `total_counts` and `n_genes_by_counts` instead."""
))

cells.append(nbf.v4.new_code_cell(
"""fig, axes = plt.subplots(1, 2, figsize=(12, 4))
axes[0].hist(adata.obs["n_genes_by_counts"], bins=100)
axes[0].set_title("Genes per cell")
axes[1].hist(adata.obs["total_counts"], bins=100)
axes[1].set_title("Total counts per cell")
plt.tight_layout()
plt.show()"""
))

cells.append(nbf.v4.new_code_cell(
"""fig, axes = plt.subplots(1, 2, figsize=(12, 4))
axes[0].hist(adata.obs["n_genes_by_counts"], bins=100, range=(0, 300))
axes[0].set_title("Genes per cell (zoomed, low end)")
axes[1].hist(adata.obs["total_counts"], bins=100, range=(0, 1000))
axes[1].set_title("Total counts per cell (zoomed, low end)")
plt.tight_layout()
plt.show()

for t in [100, 150, 200, 250, 300]:
    print(f"genes < {t}: {(adata.obs['n_genes_by_counts'] < t).sum()} cells")"""
))

cells.append(nbf.v4.new_markdown_cell(
"""## Choosing a filtering threshold

The distribution of genes detected per cell shows a small population of
likely debris/empty droplets (fewer than ~150 genes), a quiet gap, then the
main cell population starting around 250+ genes. We use `min_genes=200`, a
standard scRNA-seq floor that falls within this gap.

`total_counts` requires no additional low-end filtering: the dataset
already has no cells below ~560 counts, indicating an upstream minimum-count
filter was already applied before this processed release."""
))

cells.append(nbf.v4.new_code_cell(
"""print(f"Before filtering: {adata.n_obs} cells")
sc.pp.filter_cells(adata, min_genes=200)
print(f"After filtering: {adata.n_obs} cells")"""
))

cells.append(nbf.v4.new_code_cell(
"""print(f"Before filtering: {adata.n_vars} genes")
sc.pp.filter_genes(adata, min_cells=3)
print(f"After filtering: {adata.n_vars} genes")"""
))

cells.append(nbf.v4.new_markdown_cell(
"""## Save checkpoint

Save the QC-filtered object so `03_preprocessing.ipynb` can load it
directly, without re-running ingestion and filtering."""
))

cells.append(nbf.v4.new_code_cell(
"""adata.write("../results/adata_qc_filtered.h5ad")
print("Saved: results/adata_qc_filtered.h5ad")"""
))

nb['cells'] = cells

with open("notebooks/02_qc.ipynb", "w") as f:
    nbf.write(nb, f)

print("02_qc.ipynb written successfully.")
