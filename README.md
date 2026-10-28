# SAVI

SAVI links ATAC peaks to the genes they regulate in **spatial multi-omic data**, and lets the peak–gene
relationship **change across space and pseudotime**. Cells are grouped into spatiotemporal bins, and a
**fused lasso** estimates one set of peak coefficients per group while pulling the coefficients of
neighbouring groups towards each other. A group-wise marginal correlation with a permutation p-value is
computed alongside it.

This repository contains the full workflow: data preparation and imputation, SAVI itself, automated
selection of its grouping parameters, baselines and evaluation.

```
Seurat objects ──► SAVI_data_prep.R ──► GLEAM_impute.py ──► SAVI_data_prep.R ──► SAVI.R ──► SAVI_evaluation.R
 (sc multiome +     (stage 1: write     (python, own        (stage 3: load       (fused     (baselines +
  spatial)           GLEAM inputs)       environment)        imputation)          lasso)     evaluation)
                                                                  SAVI_tuning.R chooses sp_num / time_num
```

## Files

| File | What it does |
|---|---|
| `SAVI_data_prep.R` | Stage 1 `SAVI_prepare_gleam_inputs()`: selects genes and candidate peaks and writes the inputs of GLEAM. Stage 3 `SAVI_load_imputation()`: reads GLEAM's output and returns the SAVI input (raw and imputed data, coordinates, genes, peaks). Self-contained. |
| `GLEAM_impute.py` | Trains GLEAM and imputes RNA and ATAC for the spatial cells. Run on its own, in its own conda environment (`run()` or the command line). |
| `SAVI.R` | The three SAVI steps, in `fusing_bench()`: (1) spatiotemporal grouping, (2) group-specific fused lasso with a within-group elastic-net baseline, (3) group-wise Spearman correlation with permutation p-values. |
| `SAVI_tuning.R` | Module 1: recommends `sp_num` and `time_num` by simulation on your own data. |
| `SAVI_evaluation.R` | Module 2: extra baselines, cross-modality and split-half reproducibility, validation against known gene–CRE pairs. |
| `SAVI_data_prep_example.ipynb` | Example: from the two Seurat objects to `SAVI_input.rds`. |
| `GLEAM_example.ipynb` | Example: running `GLEAM_impute.py`. |
| `SAVI_example.ipynb` | Example: tuning, SAVI, baselines and evaluation (R kernel). |

## Requirements

**R** (4.x)
- Data preparation: `Seurat` (>= 5), `Signac`, `GenomicRanges`, `Matrix`, `data.table`
- SAVI: `dplyr`, `glmnet`, `future`, `furrr`, and **`myFuser`** (provides `fusedLassoProximal`, the fused-lasso solver). Install with `devtools::install_local("code/myFuser/",force = T)`
- Evaluation: `ggplot2`; `pROC` and `PRROC` (AUROC / AUPR); `IRanges` and `GenomicRanges` (gene–CRE overlaps)
- Only for the Signac / Cicero baselines: `cicero` and `monocle3`

**Python** (separate conda environment, GPU recommended)
- `torch`, `torch_geometric`, `torch_sparse`, `torchvision`, `numpy`, `pandas`, `scipy`, `scikit-learn`

## Quick start

All steps are shown in the notebooks; this is the short version.

**1. Prepare the inputs for GLEAM (R).**
```r
source("SAVI_data_prep.R")
SAVI_prepare_gleam_inputs(
  sc_data = e18_brain,          # single-cell multiome Seurat object (RNA + ATAC)
  sp_data = P22mouse,           # spatial Seurat object (RNA + peaks)
  fragment_path = "spatial_atac_fragments.tsv.gz",
  work_dir = "work_dir"
)
```

**2. Impute with GLEAM (shell, GLEAM environment).**
```bash
python GLEAM_impute.py --work_dir work_dir --gpu 0
```
Options: `--epochs` (600), `--hidden_channels` (140), `--lr` (0.001), `--seed` (0), `--skip_atac`.

**3. Load the imputation (R).**
```r
prep <- SAVI_load_imputation("work_dir")
saveRDS(prep, "work_dir/SAVI_input.rds")
```
`prep` holds `imputed` and `raw` (`rna`, `atac`), `spatial_coord`, `cell_types`, `use_genes`, `use_peaks`
(candidate peaks per gene), `sc_peaks`, `sp_peaks` and `sp_cells`.

**4. Choose the grouping parameters (optional).**
```r
source("SAVI_tuning.R")
tuning <- SAVI_tune_grid(prep$imputed, prep$spatial_coord, pseudoT, prep$use_genes, prep$use_peaks,
                         sp_nums = c(3, 5, 7, 10), time_nums = c(1, 3, 5, 7), n_rep = 2)
tuning$recommended        # sp_num, time_num, sp_interval, time_interval
```

**5. Run SAVI.**
```r
library(future); library(furrr); library(glmnet); library(dplyr); library(myFuser)
source("SAVI.R")
plan(multisession, workers = 8)
res <- fusing_bench(
  sp_interval   = diff(range(prep$spatial_coord[, 1])) / sp_num,
  time_interval = diff(range(pseudoT)) / time_num,
  spatial_coord = prep$spatial_coord, pseudoT = pseudoT,
  use_peaks = prep$use_peaks, use_genes = prep$use_genes,
  imputation = prep$imputed, cell_types = prep$cell_types
)
```
`pseudoT` is **one** pseudotime value per spatial cell (same order as `spatial_coord`), shared by all genes.

**6. Baselines and evaluation (optional).**
```r
source("SAVI_evaluation.R")
baselines  <- SAVI_baselines(prep$imputed, names(res$coef), prep$use_peaks, prep$spatial_coord)
evaluation <- SAVI_evaluate(res, names(res$coef), prep$use_peaks, baselines = baselines,
                            gene_cre = gene_cre)   # gene_cre optional
```

## SAVI output

`res$coef[[gene]]` contains, per gene:

| Element | Content |
|---|---|
| `fusion` | peak × group coefficients of the final fused lasso (lambda and gamma chosen by grouped 5-fold cross-validation) |
| `fusion_old` | fused lasso with fixed lambda = gamma = 0.001 |
| `lasso` | within-group elastic net (baseline, no sharing across groups) |
| `marginal`, `marginal_p` | within-group Spearman correlation of each peak with the gene, and its permutation p-value |
| `fusion_cv` | selected lambda and gamma, and the full cross-validation grid |

The grouping is shared by all genes and stored once: `res$groups`, `res$group_all`, `res$group_count`.
Genes that fail return `NULL` and are logged in `fuser_region_error.log`.

## Module 1: choosing `sp_num` and `time_num` (`SAVI_tuning.R`)

`SAVI_tune_grid()` takes the first gene and its imputed ATAC as real data, declares the first 5 peaks the true
enhancers and simulates the gene's expression from them, with enhancer effects that vary smoothly over space and
pseudotime. It then fits SAVI's fused lasso for every `sp_num` × `time_num` in the grid and scores the
estimates against the known truth (AUROC, AUPR, precision, recall, F1, averaged over groups). The recommended
pair has the best mean of AUROC, AUPR and F1 among settings that keep at least 80% of the cells; near-ties go to
the coarser grouping.

- `SAVI_grouping_summary()` shows how many groups and cells each setting keeps, without fitting anything.
- `SAVI_tune_grid(..., dry_run = TRUE)` returns only that grouping summary.
- `plot_tuning_grid()` draws the metrics against `time_num`, one panel per `sp_num`.

## Module 2: baselines and evaluation (`SAVI_evaluation.R`)

**Baselines** (`SAVI_baselines()`): global lasso `Lasso(g)` and Moran's R (computed from the imputed data);
optionally eSpatial (from the eSpatial tables), Signac links and Cicero (from the spatial Seurat object, its
fragment file and a gene annotation, via `build_link_object()`). The within-group lasso `Lasso(w)` and the
`Marginal` correlation come from SAVI itself.

**Evaluations** (`SAVI_evaluate()` runs what its inputs allow and reports what it skipped):

| Evaluation | Needs | Function |
|---|---|---|
| Validation against known gene–CRE pairs: per-gene AUROC and AUPR, precision and recall of the top-k peaks | a gene–CRE table (`Gene`, `Coordinate1`, `Coordinate2`) | `eval_validation_pairs()`, `eval_validation_metrics()`, `eval_topk_precision()` |
| Cross-modality reproducibility (Spearman correlation and Jaccard index), e.g. ATAC vs H3K27ac | a SAVI result for a second modality (optional) | `eval_reproducibility()`, `eval_link_reproducibility()` |
| Split-half reproducibility: the same analysis on two random halves of the cells | nothing extra | `SAVI_split_reproducibility()` |

## Notes

- Moran's R uses a bandwidth of range(x) / 20 by default. An earlier analysis script used
  `max(x) - min(y)` / 20, which mixes the two coordinates; pass `moran_bandwidth` to reproduce those numbers.
- `additional_gene` in `SAVI_prepare_gleam_inputs()` adds genes to both GLEAM and the SAVI analysis. A gene is
  still kept only if it is measured in the spatial data and has more than `min_peaks_per_gene` (default 50)
  variable peaks in its window; lower `min_peaks_per_gene` or raise `n_sc_var_peaks` to keep more genes.
- SAVI is slow: each gene runs a grouped cross-validation over a lambda × gamma grid and 1000 permutations.
  Use `future::plan(multisession, workers = n)` and, for many genes, run them in batches and save each batch.

## Citation and license

Add the citation of the SAVI paper and the license of the repository here.
