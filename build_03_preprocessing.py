import nbformat as nbf

nb = nbf.v4.new_notebook()
cells = []

cells.append(nbf.v4.new_markdown_cell(
"""# Phase 3 — Preprocessing

Goal: understand the expression matrix inherited from `02_qc.ipynb`, choose
an appropriate normalization strategy based on what is actually in
`adata.X` (not assumed), then reduce dimensionality and cluster."""
))

cells.append(nbf.v4.new_code_cell(
"""import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt

adata = sc.read_h5ad("../results/adata_qc_filtered.h5ad")
adata"""
))

cells.append(nbf.v4.new_markdown_cell(
"""## Inspecting the expression matrix before transforming it

We check whether `adata.X` holds raw counts (non-negative integers), or
data that has already been normalized/transformed (floats, possibly with
negative values from scaling). This determines whether normalization is
needed at all, and prevents blindly re-applying a standard raw-counts
recipe to data that may already be processed."""
))

cells.append(nbf.v4.new_code_cell(
"""X_sample = adata.X[:1000].toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X[:1000])

is_integer_like = np.allclose(X_sample, np.round(X_sample))
has_negative = (X_sample < 0).any()
max_val = X_sample.max()

print(f"dtype: {adata.X.dtype}")
print(f"min: {X_sample.min():.3f}, max: {max_val:.3f}")
print(f"Values look integer-like (raw counts?): {is_integer_like}")
print(f"Contains negative values (already scaled?): {has_negative}")
print(f"Layers available: {list(adata.layers.keys())}")
print(f"adata.raw available: {adata.raw is not None}")

if is_integer_like and not has_negative:
    print("\\nInterpretation: adata.X appears to hold raw (or near-raw) counts.")
    print("-> Will apply standard normalization: total-count normalize + log1p.")
elif has_negative:
    print("\\nInterpretation: adata.X appears to already be scaled/transformed")
    print("(negative values present, consistent with prior scaling).")
    print("-> Standard re-normalization would not be appropriate; documenting as-is.")
else:
    print("\\nInterpretation: adata.X appears to be non-integer but non-negative")
    print("(possibly already normalized/log-transformed).")
    print("-> Will proceed carefully, checking distribution before any further transform.")"""
))

cells.append(nbf.v4.new_markdown_cell(
"""## Normalization

Based on the inspection above: if raw counts are confirmed, we store them
in `adata.layers["counts"]` before transforming (so raw counts remain
recoverable), then apply total-count normalization (each cell scaled to
the same total count) followed by log1p — the standard approach for
count-based scRNA-seq data, making expression values comparable across
cells with different sequencing depth."""
))

cells.append(nbf.v4.new_code_cell(
"""adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
print("Normalization applied. adata.X now holds log1p(normalized counts).")
print(f"New min/max: {adata.X.min():.3f}, {adata.X.max():.3f}")"""
))

cells.append(nbf.v4.new_markdown_cell(
"""## Highly variable gene (HVG) selection

Selecting highly variable genes focuses downstream analysis (PCA,
clustering) on genes carrying the most biological signal, reducing noise
from genes with flat expression across all cells. We use the standard
Seurat-flavor method and 2000 genes, a common default that balances
signal retention against dimensionality."""
))

cells.append(nbf.v4.new_code_cell(
"""sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat")
n_hvg = adata.var["highly_variable"].sum()
print(f"Selected {n_hvg} highly variable genes out of {adata.n_vars}")

adata.raw = adata  # keep full gene set accessible for later marker lookups
adata = adata[:, adata.var.highly_variable].copy()
print(f"Subset to HVGs: {adata.shape}")"""
))

cells.append(nbf.v4.new_markdown_cell(
"""## Scaling and PCA

Scaling (zero mean, unit variance per gene) prevents highly expressed
genes from dominating the principal components purely due to magnitude.
We clip extreme values (max_value=10) to reduce the influence of outlier
cells. PCA then summarizes the HVG expression space into a smaller number
of components capturing the main axes of variation."""
))

cells.append(nbf.v4.new_code_cell(
"""sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver="arpack")

sc.pl.pca_variance_ratio(adata, n_pcs=50, log=False)
plt.show()"""
))

cells.append(nbf.v4.new_markdown_cell(
"""## Choosing the number of principal components

We inspect the variance ratio plot above to choose how many PCs to carry
into neighbor-graph construction: enough to capture the main structure,
without including components that mostly reflect noise. Record the number
chosen and the rationale directly here after viewing the plot."""
))

cells.append(nbf.v4.new_code_cell(
"""N_PCS = 20  # update after inspecting the variance ratio plot above

sc.pp.neighbors(adata, n_pcs=N_PCS)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=1.0)

print(f"Leiden clustering found {adata.obs['leiden'].nunique()} clusters using {N_PCS} PCs")"""
))

cells.append(nbf.v4.new_markdown_cell(
"""## Visualizing structure: clusters, condition, and existing cell type

We look at the UMAP colored three ways: by our independently derived
Leiden clusters, by condition (control vs. stimulated), and by the
dataset's pre-existing `cell_type` labels. This last one is shown here
only for visual context — formal comparison against our own marker-based
annotation happens in `04_clustering_annotation.ipynb`, not here."""
))

cells.append(nbf.v4.new_code_cell(
"""sc.pl.umap(adata, color=["leiden", "label", "cell_type"], ncols=1, wspace=0.4)
plt.show()"""
))

cells.append(nbf.v4.new_markdown_cell(
"""## Save checkpoint

Save the preprocessed, clustered object for `04_clustering_annotation.ipynb`."""
))

cells.append(nbf.v4.new_code_cell(
"""adata.write("../results/adata_preprocessed.h5ad")
print("Saved: results/adata_preprocessed.h5ad")"""
))

nb['cells'] = cells

with open("notebooks/03_preprocessing.ipynb", "w") as f:
    nbf.write(nb, f)

print("03_preprocessing.ipynb written successfully.")
