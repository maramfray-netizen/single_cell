# Single-Cell RNA-seq Training Project — Progress Log
   
Project: PBMC control vs IFN-β stimulated (Kang et al. 2018)
Supervisor: Koussai

## Phase 1 — Ingestion
- [x] Environment set up (conda env `single_cell`, WSL + VS Code)
- [x] GitHub repo created and linked
- [x] Loaded Kang 2018 dataset via `pertpy.dt.kang_2018()`
- [ ] Confirmed cell/gene counts and condition label column
- Notes: resolved `ModuleNotFoundError: filelock` by installing missing dependency

## Phase 2 — Raw QC
- [ ] Not started

## Phase 3 — Filtering
- [ ] Not started

## Phase 4 — Normalization
- [ ] Not started

## Phase 5 — Dimensionality reduction
- [ ] Not started

## Phase 6 — Clustering
- [ ] Not started

## Phase 7 — Cell type annotation
- [ ] Not started

## Phase 8 — Marker DEA
- [ ] Not started

## Phase 9 — Condition-contrast DEA
- [ ] Not started

## Phase 10 — Interpretation and reporting
- [ ] Not started

### Phase 1 results
- 24,673 cells x 15,706 genes loaded via pertpy (Kang et al. 2018)
- Condition split: 12,358 stimulated / 12,315 control
- Pre-existing cell_type labels: CD4 T cells (11238), CD14+ Monocytes (5697), B cells (2651), NK cells (1716), CD8 T cells (1621), FCGR3A+ Monocytes (1089), Dendritic cells (529), Megakaryocytes (132)

### Phase 2-3 results
- Mitochondrial QC not applicable: MT genes absent from this pre-processed dataset (upstream removal, documented as a limitation)
- Cell filtering: min_genes=200 -> removed 111 cells (24673 -> 24562)
- Gene filtering: min_cells=3 -> removed 5 genes (15706 -> 15701)
- Checkpoint saved: results/adata_qc_filtered.h5ad


### Phase 4 results
- Independent Leiden clustering (18 clusters) scored against 8 canonical PBMC marker gene sets to propose cell type identities
- Proposed identities compared against the dataset's existing cell_type (Kang) labels
- Disagreements investigated: cluster 6 (NK vs CD8 T - shared cytotoxic markers), cluster 13 (Megakaryocyte vs CD4 T - our call favored by strong platelet markers)
- Low-confidence clusters flagged: cluster 17 (hemoglobin genes - likely erythrocyte contamination), clusters 10/14 (housekeeping genes, weak signal)
- Checkpoint saved locally (not committed): results/adata_annotated.h5ad; comparison table committed: results/cluster_annotation_comparison.csv


### Phase 4 results
- Independent Leiden clustering (18 clusters) scored against 8 canonical PBMC marker gene sets to propose cell type identities
- Proposed identities compared against the dataset's existing cell_type (Kang) labels
- Disagreements investigated: cluster 6 (NK vs CD8 T - shared cytotoxic markers), cluster 13 (Megakaryocyte vs CD4 T - our call favored by strong platelet markers)
- Low-confidence clusters flagged: cluster 17 (hemoglobin genes - likely erythrocyte contamination), clusters 10/14 (housekeeping genes, weak signal)
- Checkpoint saved locally (not committed): results/adata_annotated.h5ad; comparison table committed: results/cluster_annotation_comparison.csv
