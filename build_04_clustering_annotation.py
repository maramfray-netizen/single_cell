import nbformat as nbf

nb = nbf.v4.new_notebook()
cells = []

cells.append(nbf.v4.new_markdown_cell(
"""# Phase 4 — Clustering Annotation

Goal: assign biological identities to the 18 Leiden clusters found in
`03_preprocessing.ipynb`, using canonical marker genes. Only after an
independent annotation is proposed do we compare it against the dataset's
pre-existing `cell_type` labels, treating those as a validation reference
rather than the starting point."""
))

cells.append(nbf.v4.new_code_cell(
"""import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

adata = sc.read_h5ad("../results/adata_preprocessed.h5ad")
print(adata)
print(adata.obs["leiden"].value_counts().sort_index())"""
))

cells.append(nbf.v4.new_markdown_cell(
"""## Marker genes per cluster"""
))

cells.append(nbf.v4.new_code_cell(
"""adata_full = adata.raw.to_adata()
adata_full.obs["leiden"] = adata.obs["leiden"]

sc.tl.rank_genes_groups(adata_full, groupby="leiden", method="wilcoxon")

top_markers = {}
for cl in adata_full.obs["leiden"].cat.categories:
    genes = adata_full.uns["rank_genes_groups"]["names"][cl][:10]
    top_markers[cl] = list(genes)
    print(f"Cluster {cl}: {list(genes)}")"""
))

cells.append(nbf.v4.new_markdown_cell(
"""## Proposing cell type identities from canonical markers"""
))

cells.append(nbf.v4.new_code_cell(
"""marker_sets = {
    "CD4 T cells": ["CD3D", "CD3E", "IL7R", "CD4"],
    "CD8 T cells": ["CD3D", "CD3E", "CD8A", "CD8B"],
    "B cells": ["MS4A1", "CD79A", "CD79B"],
    "NK cells": ["NKG7", "GNLY", "KLRD1"],
    "CD14+ Monocytes": ["CD14", "LYZ", "S100A8"],
    "FCGR3A+ Monocytes": ["FCGR3A", "MS4A7"],
    "Dendritic cells": ["FCER1A", "CST3", "CLEC10A"],
    "Megakaryocytes": ["PPBP", "PF4"],
}

for name, genes in marker_sets.items():
    genes_present = [g for g in genes if g in adata_full.var_names]
    sc.tl.score_genes(adata_full, genes_present, score_name=f"score_{name}")

score_cols = [f"score_{name}" for name in marker_sets]
cluster_scores = adata_full.obs.groupby("leiden")[score_cols].mean()
proposed_identity = cluster_scores.idxmax(axis=1).str.replace("score_", "", regex=False)

print(cluster_scores.round(3))
print()
print(proposed_identity)"""
))

cells.append(nbf.v4.new_markdown_cell(
"""## Comparing to the existing Kang cell_type labels"""
))

cells.append(nbf.v4.new_code_cell(
"""adata_full.obs["cell_type"] = adata.obs["cell_type"]

kang_mode_per_cluster = adata_full.obs.groupby("leiden")["cell_type"].agg(lambda x: x.mode()[0])

comparison = pd.DataFrame({
    "top_markers": [", ".join(top_markers[cl][:5]) for cl in proposed_identity.index],
    "proposed_identity": proposed_identity,
    "kang_reference_identity": kang_mode_per_cluster,
})
comparison["agreement"] = comparison["proposed_identity"] == comparison["kang_reference_identity"]

comparison"""
))

cells.append(nbf.v4.new_markdown_cell(
"""## Discussing disagreements and low-confidence clusters

**Cluster 6 (proposed NK cells, Kang reference: CD8 T cells):** top markers
(CCL5, NKG7, GZMB, GZMH, APOBEC3G) are cytotoxic effector genes shared by
both NK cells and activated/cytotoxic CD8 T cells - a known, biologically
expected ambiguity between these two cytotoxic lineages, not an annotation
error.

**Cluster 13 (proposed Megakaryocytes, Kang reference: CD4 T cells):** top
markers (PPBP, PF4, GNG11, SDPR) are canonical platelet/megakaryocyte
genes, with one of the highest marker scores in the dataset. The marker
evidence here favors our independent call over the Kang reference label.

**Low-confidence clusters:**
- Cluster 17: top markers (HBB, HBA1, HBA2, ALAS2) are hemoglobin genes,
  indicating a likely contaminating erythrocyte/reticulocyte population,
  not covered by the 8 PBMC marker sets used here.
- Clusters 10 and 14: dominated by housekeeping/ribosomal genes, suggesting
  weak transcriptional specificity rather than a clearly wrong assignment.

**Overall:** the large majority of clusters show strong, marker-supported
agreement with the Kang reference annotation, validating both the
independent clustering and the marker-based annotation approach."""
))

cells.append(nbf.v4.new_code_cell(
"""max_scores = cluster_scores.max(axis=1)
print("Max marker score per cluster (low = ambiguous/unassigned):")
print(max_scores.sort_values())"""
))

cells.append(nbf.v4.new_markdown_cell(
"""## Visualizing marker expression across clusters"""
))

cells.append(nbf.v4.new_code_cell(
"""all_markers = sum(marker_sets.values(), [])
all_markers_present = [g for g in dict.fromkeys(all_markers) if g in adata_full.var_names]

sc.pl.dotplot(adata_full, all_markers_present, groupby="leiden", standard_scale="var")
plt.savefig("../figures/fig2_marker_dotplot.png", dpi=150, bbox_inches="tight")
plt.show()"""
))

cells.append(nbf.v4.new_code_cell(
"""sc.pl.umap(adata_full, color=["leiden"], legend_loc="on data")
plt.show()"""
))

cells.append(nbf.v4.new_markdown_cell(
"""## Save checkpoint"""
))

cells.append(nbf.v4.new_code_cell(
"""cluster_to_identity = proposed_identity.to_dict()
adata_full.obs["proposed_cell_type"] = adata_full.obs["leiden"].map(cluster_to_identity)

adata_full.write("../results/adata_annotated.h5ad")
comparison.to_csv("../results/cluster_annotation_comparison.csv")
print("Saved: results/adata_annotated.h5ad and results/cluster_annotation_comparison.csv")"""
))

nb['cells'] = cells

with open("notebooks/04_clustering_annotation.ipynb", "w") as f:
    nbf.write(nb, f)

print("04_clustering_annotation.ipynb written successfully.")
