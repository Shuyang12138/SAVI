# SAVI
SAVI is a graph-guided fused lasso framework that integrates spatial context with latent cellular state to infer gene-enhancer interactions from spatial multi-modal data.

## Preparation
### Install/load required packages
For the first step, we can download the code by
```bash
git clone https://github.com/Shuyang12138/SAVI.git
```
### Input data
The input data should be an imputed data list `imputation` with one element `imputation$rna` containing the cell-by-gene expression imputation matrix and another element `imputation$atac` containing the cell-by-peak accessibility matrix, a spatial coordinate matrix containing spatial coordinates of cells and a vector (Monocle3) or list (MultiVelo) of pseudo times for each cell.
```R
load('AdultMB_multiplex/real_application.Rdata')
```
## GE identification
```R
devtools::install_local("code/myFuser/",force = T)
source('code/savi.R')
library(furrr)
plan(multicore, workers = 10)
example_gene = names(P22_multipseudoT)
MB5Matac_result <- fusing_bench_multiveloT(
        sp_num = 6,
        time_num = 5,
        spatial_coord = atac_spatial_coord,
        pseudoT = P22_multipseudoT,
        use_peaks = atac_use_peaks[example_gene],
        use_genes = example_gene ,
        imputation = P22_atac_imputation,
        cell_types = P22_celltype,
          intercept=T
      )
```
## Output
The output includes SAVI inferred gene-enhancer regulatory score for each gene as an element named `fusion`

## Requirement
R packages: Seurat, Signac, data.table, GenomicRanges
