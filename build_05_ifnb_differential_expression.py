import nbformat as nbf

nb = nbf.v4.new_notebook()
cells = []

cells.append(nbf.v4.new_markdown_cell(
"""# Phase 5 - IFN-b Differential Expression

Central biological question: how does IFN-b stimulation alter transcriptional
programs across different PBMC cell populations? For each cell type, we
compare stimulated vs. control cells and check whether known
interferon-stimulated genes (ISGs) are recovered.

Note: UMAP visualization was skipped in this environment due to a
package-level conflict (torch/sympy) affecting UMAP's numba-based backend;
clustering (Leiden) and differential expression do not depend on UMAP and
are unaffected."""
))

cells.append(nbf.v4.new_code_cell(
"""import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

adata = sc.read_h5ad("../results/adata_annotated.h5ad")
print(adata.obs.groupby(["proposed_cell_type", "label"]).size().unstack())"""
))

cells.append(nbf.v4.new_code_cell(
"""canonical_isgs = ["ISG15", "IFI6", "MX1", "OAS1", "IFIT1", "IFIT3", "STAT1", "IRF7", "IFITM3", "RSAD2"]
canonical_isgs_present = [g for g in canonical_isgs if g in adata.var_names]
print(f"{len(canonical_isgs_present)} of {len(canonical_isgs)} canonical ISGs found")
print(canonical_isgs_present)"""
))

cells.append(nbf.v4.new_code_cell(
"""cell_types = adata.obs["proposed_cell_type"].unique()
de_results = {}
MIN_CELLS = 30

for ct in cell_types:
    sub = adata[adata.obs["proposed_cell_type"] == ct].copy()
    n_ctrl = (sub.obs["label"] == "ctrl").sum()
    n_stim = (sub.obs["label"] == "stim").sum()
    if n_ctrl < MIN_CELLS or n_stim < MIN_CELLS:
        print(f"Skipping {ct}: insufficient cells")
        continue
    sc.tl.rank_genes_groups(sub, groupby="label", groups=["stim"], reference="ctrl", method="wilcoxon")
    result = sc.get.rank_genes_groups_df(sub, group="stim").set_index("names")
    de_results[ct] = result
    n_sig = (result["pvals_adj"] < 0.05).sum()
    print(f"{ct}: ctrl={n_ctrl}, stim={n_stim}, {n_sig} genes significant")"""
))

cells.append(nbf.v4.new_code_cell(
"""isg_summary = pd.DataFrame(index=canonical_isgs_present, columns=list(de_results.keys()))
for ct, result in de_results.items():
    for gene in canonical_isgs_present:
        if gene in result.index:
            lfc = result.loc[gene, "logfoldchanges"]
            padj = result.loc[gene, "pvals_adj"]
            isg_summary.loc[gene, ct] = f"{lfc:.2f} (p={padj:.1e})"
        else:
            isg_summary.loc[gene, ct] = "n/a"
isg_summary"""
))

cells.append(nbf.v4.new_code_cell(
"""top_up_genes = {}
for ct, result in de_results.items():
    sig_up = result[(result["pvals_adj"] < 0.05) & (result["logfoldchanges"] > 1)].sort_values("pvals_adj").head(50)
    top_up_genes[ct] = set(sig_up.index)
    print(f"{ct}: {len(sig_up)} genes (adj p<0.05, logFC>1)")

from functools import reduce
shared_across_all = reduce(lambda a, b: a & b, top_up_genes.values())
print(f"\\nShared across ALL cell types ({len(shared_across_all)}): {sorted(shared_across_all)}")"""
))

cells.append(nbf.v4.new_code_cell(
"""response_magnitude = pd.Series({ct: (r["pvals_adj"] < 0.05).sum() for ct, r in de_results.items()}).sort_values(ascending=False)
print(response_magnitude)

fig, ax = plt.subplots(figsize=(8, 5))
response_magnitude.plot(kind="barh", ax=ax)
ax.set_xlabel("Number of significant DE genes (adj p < 0.05)")
ax.set_title("IFN-b response magnitude by cell type")
plt.tight_layout()
plt.savefig("../figures/fig3a_response_magnitude.png", dpi=150, bbox_inches="tight")
plt.show()"""
))

cells.append(nbf.v4.new_code_cell(
"""adata.obs["cell_type_condition"] = adata.obs["proposed_cell_type"].astype(str) + "_" + adata.obs["label"].astype(str)
sc.pl.dotplot(adata, canonical_isgs_present, groupby="cell_type_condition", standard_scale="var")
plt.savefig("../figures/fig3b_isg_dotplot.png", dpi=150, bbox_inches="tight")
plt.show()"""
))

cells.append(nbf.v4.new_code_cell(
"""strongest_ct = response_magnitude.index[0]
result = de_results[strongest_ct].copy()
result["neg_log10_padj"] = -np.log10(result["pvals_adj"].clip(lower=1e-300))

fig, ax = plt.subplots(figsize=(8, 6))
ax.scatter(result["logfoldchanges"], result["neg_log10_padj"], s=5, alpha=0.4, color="gray")
highlight = result.loc[result.index.isin(canonical_isgs_present)]
ax.scatter(highlight["logfoldchanges"], highlight["neg_log10_padj"], s=30, color="red")
for gene, row in highlight.iterrows():
    ax.annotate(gene, (row["logfoldchanges"], row["neg_log10_padj"]), fontsize=8)
ax.set_xlabel("log2 fold change (stim vs ctrl)")
ax.set_ylabel("-log10(adjusted p-value)")
ax.set_title(f"IFN-b response volcano plot: {strongest_ct}")
plt.tight_layout()
plt.savefig("../figures/fig3c_volcano.png", dpi=150, bbox_inches="tight")
plt.show()"""
))

cells.append(nbf.v4.new_code_cell(
"""import os
os.makedirs("../results", exist_ok=True)
for ct, result in de_results.items():
    safe_name = ct.replace(" ", "_").replace("+", "plus")
    result.to_csv(f"../results/de_{safe_name}.csv")
isg_summary.to_csv("../results/isg_summary_by_celltype.csv")
response_magnitude.to_csv("../results/response_magnitude_by_celltype.csv")
print("Saved all Phase 5 results")"""
))

nb['cells'] = cells
with open("notebooks/05_ifnb_differential_expression.ipynb", "w") as f:
    nbf.write(nb, f)
print("05_ifnb_differential_expression.ipynb written successfully.")
