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
