# Single-Cell RNA-seq Training Project - Progress Log

Project: PBMC control vs IFN-β stimulated (Kang et al. 2018)
Supervisor: Koussai Salem - DrugIT

## Phase 1 - Ingestion
- [x] Environment set up (conda env `single_cell`, WSL + VS Code, later Google Colab for heavy compute)
- [x] GitHub repo created and linked
- [x] Loaded Kang 2018 dataset via `pertpy.dt.kang_2018()`
- [x] Confirmed cell/gene counts and condition/cell_type label distributions

### Phase 1 results
- 24,673 cells x 15,706 genes loaded via pertpy (Kang et al. 2018)
- Condition split: 12,358 stimulated / 12,315 control
- Pre-existing cell_type labels: CD4 T cells (11238), CD14+ Monocytes (5697), B cells (2651), NK cells (1716), CD8 T cells (1621), FCGR3A+ Monocytes (1089), Dendritic cells (529), Megakaryocytes (132)

## Phase 2 - QC
- [x] Inspected AnnData object (X, layers, obs, obsm, raw) before transformation
- [x] Confirmed mitochondrial genes absent from this processed dataset release (documented as limitation, not a missing step)
- [x] Filtered cells (min_genes=200) and genes (min_cells=3)

### Phase 2 results
- Mitochondrial QC not applicable: MT genes absent from this pre-processed dataset (upstream removal, documented as a limitation)
- Cell filtering: min_genes=200 -> removed 111 cells (24673 -> 24562)
- Gene filtering: min_cells=3 -> removed 5 genes (15706 -> 15701)

## Phase 3 - Preprocessing
- [x] Normalization (normalize_total + log1p)
- [x] HVG selection (2000 genes)
- [x] Scaling + PCA
- [x] Neighbors graph + Leiden clustering

### Phase 3 results
- Normalization: normalize_total (target_sum=1e4) + log1p
- HVG selection: 2000 highly variable genes (seurat flavor) out of 15701; subset to (24562, 2000)
- PCA: 20 PCs chosen from the variance ratio elbow plot
- Leiden clustering (resolution=1.0): 18 clusters found
- UMAP visualization skipped due to an unresolved torch/sympy package conflict in the Colab environment (AttributeError: module sympy has no attribute printing); does not affect clustering or downstream DE

## Phase 4 - Clustering Annotation
- [x] Scored 18 independent Leiden clusters against 8 canonical PBMC marker gene sets
- [x] Proposed cell type identities from markers
- [x] Compared against existing Kang cell_type labels as validation reference

### Phase 4 results
- Proposed identities compared against the dataset's existing cell_type (Kang) labels: large majority of 18 clusters agree
- Disagreements investigated, not hidden:
  - Cluster 6 (proposed NK, Kang: CD8 T): shared cytotoxic markers (CCL5, NKG7, GZMB, GZMH) - known biological ambiguity
  - Cluster 13 (proposed Megakaryocyte, Kang: CD4 T): strong platelet markers (PPBP, PF4, GNG11), highest score in dataset - marker evidence favors our independent call
- Low-confidence clusters flagged: cluster 17 (hemoglobin genes - likely erythrocyte contamination, not covered by the 8-category marker panel), clusters 10 and 14 (housekeeping genes, weak identity signature)

## Phase 5 - IFN-beta Differential Expression
- [x] Wilcoxon DE test (stimulated vs control) run per proposed cell type
- [x] Canonical ISGs checked for recovery
- [x] Response magnitude ranked by cell type
- [x] Shared vs cell-type-specific response genes identified

### Phase 5 results
- Differential expression (Wilcoxon, stim vs ctrl) run per proposed cell type
- Canonical ISGs (ISG15, IFI6, MX1, OAS1, IFIT1, IFIT3, STAT1, IRF7, IFITM3, RSAD2) checked for recovery across cell types
- Response magnitude, shared vs cell-type-specific response genes, and myeloid vs lymphoid comparison computed
- Figures saved: fig3a_response_magnitude, fig3b_isg_dotplot, fig3c_volcano

## Project status: complete
All 5 phases finished and pushed. README.md updated with final methodology, findings, and limitations.
