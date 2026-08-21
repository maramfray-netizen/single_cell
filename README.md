# Single-Cell RNA-seq Reanalysis: IFN-β Response in PBMCs

## Project question

How does IFN-β (interferon-beta) stimulation alter transcriptional programs
across different PBMC (peripheral blood mononuclear cell) populations?

This is a reproducible reanalysis of an already-processed single-cell
RNA-seq perturbation dataset — not a from-scratch pipeline starting from
raw sequencing reads.

## Dataset

Kang et al. 2018 (Nature Biotechnology), PBMC control vs. IFN-β
stimulated, loaded via pertpy.dt.kang_2018(). Source: GEO GSE96583.

- 24,673 cells x 15,706 genes as distributed
- Two conditions: control / stimulated (roughly balanced)
- Includes pre-existing cell_type labels, used later in this project
  strictly as a validation reference, not as the primary annotation
  method - clustering and marker-based annotation are performed
  independently first (see 04_clustering_annotation.ipynb).

## Environment

conda env create -f environment.yml
conda activate single_cell
python -m ipykernel install --user --name single_cell

## Repository structure

single_cell/
- README.md
- environment.yml
- PROGRESS.md
- notebooks/
  - 01_ingestion.ipynb
  - 02_qc.ipynb
  - 03_preprocessing.ipynb (in progress)
  - 04_clustering_annotation.ipynb (planned)
  - 05_ifnb_differential_expression.ipynb (planned)
- figures/
- results/ (generated .h5ad checkpoints, not committed - see .gitignore)

Large generated AnnData objects are not committed to this repository.
Each notebook contains the code needed to regenerate them from the
previous checkpoint or from the original dataset.

## Methodology so far

01 - Ingestion: load the Kang 2018 dataset via pertpy; confirm cell/gene
counts and condition/cell_type label distributions.

02 - QC: inspect the AnnData object structure (X, layers, obs, obsm, raw)
before any transformation. Mitochondrial genes are absent from this
processed release (documented as a dataset limitation, not treated as a
missing step). Cells are filtered on min_genes=200 (justified from the
gene-count distribution); genes filtered on min_cells=3. Result: 24,673
to 24,562 cells, 15,706 to 15,701 genes.

(Sections below will be filled in as each phase is completed.)

## Main findings

(To be completed once the IFN-beta differential expression analysis,
phase 5, is done.)

## Limitations

- Mitochondrial-based QC was not possible: mitochondrial genes are absent
  from this processed dataset release.
- This is an educational reanalysis, not a validated clinical or
  diagnostic result.

## Reproducing this analysis

1. Clone this repository.
2. Create the environment as described above.
3. Run notebooks in numerical order (01 through 05); each depends on the
   checkpoint saved by the previous one in results/.
