################################################################################
# SAVI.R  --  annotated version
#
# OVERVIEW OF THE THREE MAIN STEPS OF SAVI (all implemented inside
# `fusing_bench()` below; Step 1 is run once, Steps 2-3 once per gene):
#
#   STEP 1. Spatiotemporal grouping of cells
#           Cells are binned on a (x, y, pseudotime) grid into groups; groups
#           with too few cells are dropped, and a similarity graph (G_simu)
#           between the remaining groups is built from their grid distances.
#           The grouping uses ONE pseudotime shared by all genes (e.g. the
#           ATAC-based pseudotime), so the groups are the same for every gene.
#
#   STEP 2. Group-specific peak-gene regression with a fused lasso
#           Gene expression is regressed on candidate peak accessibility, with
#           one coefficient vector per group. The fusion penalty (via G_simu)
#           shrinks coefficients of neighbouring groups towards each other.
#           lambda (sparsity) and gamma (fusion) are tuned by grouped 5-fold CV.
#           A within-group elastic-net fit (glmnet) is also computed as a
#           baseline that does not share information across groups.
#
#   STEP 3. Group-wise marginal correlation and permutation-based significance
#           Spearman correlation between each peak and the gene is computed
#           within each group; a null distribution is built by permuting gene
#           expression within groups (1000 times) and converted to a one-sided
#           p-value with `signac_p()`.
#
# Sections before the main function are helper functions.
################################################################################

# ------------------------------------------------------------------------------
# 0. Libraries
# ------------------------------------------------------------------------------
library(Seurat)
library(Signac)
library(ggplot2)
library(dplyr)
library(data.table)
library(GenomicRanges)
library(Matrix)
library(glmnet)   # within-group elastic-net baseline in Step 2

# ------------------------------------------------------------------------------
# 0a. Helper: find_peak()
#     Builds a binary peak x gene sparse matrix indicating which peaks overlap
#     which gene ranges (used to define the candidate peaks for each gene).
#       peaks    : GRanges of ATAC peaks
#       genes    : GRanges of gene bodies (or gene bodies + extended windows)
#       distance : (unused in the body below)
#       sep      : separators used to write ranges as strings, e.g. chr-start-end
#     Returns: sparse matrix, rows = peaks, cols = genes, 1 = overlap.
# ------------------------------------------------------------------------------
find_peak <- function(peaks, genes, distance = 0, sep = c("-", "-")) {
  # All (peak, gene) pairs whose ranges overlap by any amount
  overlaps <- findOverlaps(
    query = peaks,
    subject = genes,
    type = 'any',
    select = 'all'
  )
  # Convert overlap hits to a binary sparse peak x gene matrix
  hit_matrix <- sparseMatrix(
    i = queryHits(x = overlaps),
    j = subjectHits(x = overlaps),
    x = 1,
    dims = c(length(x = peaks), length(x = genes))
  )
  # Name rows/columns by genomic coordinates
  rownames(x = hit_matrix) <- GRangesToString(grange = peaks, sep = sep)
  colnames(x = hit_matrix) <- GRangesToString(grange = genes, sep = sep)
  return(hit_matrix)
}

# ------------------------------------------------------------------------------
# 0b. Helper: CollapseToLongestTranscript()
#     Collapses transcript-level annotations to one range per gene, spanning
#     from the minimum start to the maximum end across all its transcripts.
#       ranges : GRanges with gene_id, gene_biotype and gene_name metadata
#     Returns: GRanges with one entry per gene_id.
# ------------------------------------------------------------------------------
CollapseToLongestTranscript <- function(ranges) {
  seqnames=start=end=strand=gene_biotype=gene_name=NULL
  range.df <- as.data.table(x = ranges)
  range.df$strand <- as.character(x = range.df$strand)
  # Unstranded ("*") ranges are treated as "+"
  range.df$strand <- ifelse(
    test = range.df$strand == "*",
    yes = "+",
    no = range.df$strand
  )
  # One row per gene: outermost start/end, first strand/biotype/name
  collapsed <- range.df[
    , list(unique(seqnames),
        min(start),
        max(end),
        strand[[1]],
        gene_biotype[[1]],
        gene_name[[1]]),
    "gene_id"
  ]
  colnames(x = collapsed) <- c(
    "gene_id", "seqnames", "start", "end", "strand", "gene_biotype", "gene_name"
  )
  # Ensure gene names are unique
  collapsed$gene_name <- make.unique(names = collapsed$gene_name)
  gene.ranges <- makeGRangesFromDataFrame(
    df = collapsed,
    keep.extra.columns = TRUE
  )
  return(gene.ranges)
}

# ------------------------------------------------------------------------------
# 0c. Helper: gaussian_kernel()
#     Gaussian distance-decay weights: exp(-d^2 / (2 * bandwidth^2)).
#     (Only used by the disabled geographically weighted correlation below.)
# ------------------------------------------------------------------------------
gaussian_kernel <- function(d, bandwidth) {
  exp(- (d^2) / (2 * bandwidth^2))
}

# ------------------------------------------------------------------------------
# 0d. Helper: gw_correlation()
#     Geographically weighted (local) Pearson correlation between X and Y at
#     location i, using row i of the weight matrix W.
#     NOTE: the calls to this function inside `fusing_bench()` are commented
#     out; the Spearman correlation per group is used instead.
# ------------------------------------------------------------------------------
gw_correlation <- function(i, X, Y, W) {
  w <- W[i, ]       # Extract weights for location i
  w <- w / sum(w)   # Normalize weights

  # Compute weighted means
  X_mean <- sum(w * X)
  Y_mean <- sum(w * Y)

  # Compute weighted covariance and variances
  cov_xy <- sum(w * (X - X_mean) * (Y - Y_mean))
  var_x <- sum(w * (X - X_mean)^2)
  var_y <- sum(w * (Y - Y_mean)^2)

  # Compute local correlation
  r_local <- cov_xy / sqrt(var_x * var_y)
  return(r_local)
}

# ------------------------------------------------------------------------------
# 0e. Helper: signac_p()
#     Converts observed correlations into p-values using a permutation null.
#       background_corr : list of null correlation matrices (peaks x groups),
#                         one per permutation
#       corre           : observed correlation matrix (peaks x groups)
#     For every peak/group entry, a Z-score is computed against the mean and SD
#     of the null, and a one-sided (upper-tail) normal p-value is returned.
#     Returns: p-value matrix with the same dimensions as `corre`.
# ------------------------------------------------------------------------------
signac_p = function(background_corr, corre){
  # Stack the list of null matrices into a 3D array (peaks x groups x permutations)
  background_array <- simplify2array(background_corr)

  # Mean and standard deviation across the third dimension (i.e., over the list)
  mean_mat <- apply(background_array, c(1, 2), mean)
  sd_mat <- apply(background_array, c(1, 2), sd)

  # Compute Z-score
  z_score <- (corre - mean_mat) / sd_mat
  # One-sided p-value: probability of a correlation at least this large under the null
  pval_mat <- pnorm(z_score, lower.tail = FALSE)
  return(pval_mat)
}

# ------------------------------------------------------------------------------
# 1. MAIN FUNCTION: fusing_bench()
#    Runs the three SAVI steps for every gene in `use_genes`.
#
#    Arguments
#      sp_interval   : width of a spatial bin (same width for x and y), in the
#                      units of `spatial_coord`. E.g. with `num` bins along x:
#                      (max(x) - min(x)) / num
#      time_interval : width of a pseudotime bin, e.g.
#                      (max(pseudoT) - min(pseudoT)) / t_num
#      spatial_coord : cells x 2 matrix/data.frame of spatial coordinates
#      pseudoT       : numeric vector, one pseudotime value per cell (same cell
#                      order as `spatial_coord`). The SAME pseudotime is used
#                      for every gene, so the grouping is shared by all genes.
#      use_peaks     : list (one element per gene) of candidate peak names
#      use_genes     : vector of gene names to analyse
#      imputation    : list with imputed matrices `$rna` (cells x genes) and
#                      `$atac` (cells x peaks)
#      cell_types    : vector of cell-type labels (one per cell)
#      intercept     : whether to fit an intercept (default FALSE)
#
#    Returns: list(coef, groups, group_all, group_count)
#      coef        : named list (one entry per gene) with the fused lasso fits,
#                    within-group lasso, marginal correlation + p-values
#      groups      : integer group ID of every retained cell (shared by genes)
#      group_all   : group label of EVERY cell (incl. cells in dropped groups)
#      group_count : retained groups, their (x, y, t) bin indices and cell counts
#
#    Requires furrr/future (future_map, furrr_options), glmnet, the custom
#    `myFuser` package (fusedLassoProximal), and dplyr to be available.
#    Call `future::plan(multisession, workers = n)` beforehand to run genes
#    in parallel (the default plan runs them sequentially).
# ------------------------------------------------------------------------------
fusing_bench = function(sp_interval,time_interval,spatial_coord,pseudoT,use_peaks,use_genes,imputation,cell_types,intercept=F){
set.seed(2025, kind = "L'Ecuyer-CMRG")  # reproducible random numbers across parallel workers

################################################################################
# STEP 1: Spatiotemporal grouping of cells + group similarity graph
#         (gene-independent: computed once, outside the gene loop)
################################################################################
# Combine spatial coordinates with the (global) pseudotime
grid = data.frame(spatial_coord)
grid$t = pseudoT
# Integer bin index of every cell along x, y and pseudotime
cell_x_simu = floor((grid[,1]-min(grid[,1]))/sp_interval)
cell_y_simu = floor((grid[,2]-min(grid[,2]))/sp_interval)
cell_t_simu = floor((grid$t-min(grid$t))/time_interval)
library("dplyr")
# Count cells per (x, y, t) bin and keep only bins with more than 5 cells
group_count = data.frame(cell_x_simu,cell_y_simu,cell_t_simu)%>%count(cell_x_simu,cell_y_simu,cell_t_simu)
group_count = group_count[group_count$n>5,]
# Label each retained group ("x-y-t") and its spatial region ("x-y")
group_count$group = paste0(group_count$cell_x_simu,'-',group_count$cell_y_simu,'-',group_count$cell_t_simu)
group_count$region = paste0(group_count$cell_x_simu,'-',group_count$cell_y_simu)
rownames(group_count) = 1:nrow(group_count)
# Pairwise Euclidean distance between groups in (x, y, t) bin space
group_dist = as.matrix(dist(scale(group_count[,1:3],center = F)))
# Map from group label to integer group ID (used to index the coefficient matrix)
level = 1:nrow(group_count)
names(level) = group_count$group
# Group label for every cell (including cells in groups that will be dropped)
group_x = data.frame(cell_x_simu,cell_y_simu,cell_t_simu)
group_x$group = paste0(group_x$cell_x_simu,'-',group_x$cell_y_simu,'-',group_x$cell_t_simu)
group_x$region = paste0(group_x$cell_x_simu,'-',group_x$cell_y_simu)
# Group similarity graph: Gaussian kernel of group distance; weak links
# (similarity < exp(-1)) are removed and self-links are set to 0.
# This matrix defines which groups are "fused" in Step 2.
G_simu = exp(-group_dist^2)
G_simu[G_simu<exp(-1)] = 0
diag(G_simu) = 0

# Loop over genes (in parallel if a future plan is set); each iteration handles one gene
result = future_map(1:length(use_genes),function(i){
tryCatch({
        # Response (gene expression) and predictors (peak accessibility)
        y_simu = imputation$rna[, use_genes[i]]
        X_simu = imputation$atac[, use_peaks[[i]]]
        rownames(X_simu) = 1:nrow(X_simu)

        gene_exp = y_simu
        # Keep only cells that belong to retained groups (n > 5)
        X_simu = X_simu[is.element(group_x$group, group_count$group), ]
        X_simu = as.matrix(X_simu)
        y_simu = y_simu[is.element(group_x$group, group_count$group)]
        # Add tiny noise to break ties / avoid constant expression within groups
        y_simu = y_simu+rnorm(n = length(y_simu),mean=0,sd = 1e-3)
        cell_types =  cell_types[is.element(group_x$group, group_count$group)]
        group_x = group_x[is.element(group_x$group, group_count$group), ]

        # Integer group ID for every retained cell
        groups = level[group_x$group]
        names(groups) = NULL

    #######################################################
    # STEP 2: Group-specific peak-gene regression (fused lasso)
    #######################################################
        # 2a. Initial fused lasso fit with fixed lambda / gamma (no tuning).
        #     Returned as `fusion_old`; its column names (group IDs) are also
        #     reused later to loop over groups.
        beta.estimate_simu = myFuser::fusedLassoProximal(
            X_simu, matrix(y_simu, ncol=1), groups, lambda=0.001, 
            tol=1e-6, gamma=0.001, G_simu, intercept=intercept, scaling=T, num.it=5000
        )

        # 2b. Hyperparameter grid for grouped cross-validation
        #     lambda = sparsity penalty, gamma = fusion penalty across groups
        lambda_seq <- c(seq(0, 0.001, 0.0002),0.002,0.005,0.008,0.01)
        gamma_seq <- c(seq(0.0001, 0.001, 0.0002),0.002,0.005,0.008,0.01)
        param_grid <- expand.grid(lam = lambda_seq, gam = gamma_seq)
        set.seed(123)
        # Stratified 5-fold split: cells are assigned to folds *within each
        # group*, so every group is represented in every training/test fold
        group_ids <- unique(groups)
        group_folds <- lapply(group_ids, function(gid) {
              idx <- which(groups == gid)
              fold_assignments <- sample(rep(1:5, length.out = length(idx)))
              split(idx, fold_assignments)
        })
        names(group_folds) <- group_ids

        # CV routine for one (lambda, gamma) pair; returns the CV RMSE.
        # Errors are written to fuser_error_log.txt and NULL is returned.
        process_combination_cv_grouped <- function(params) {
            tryCatch({
              lam <- params$lam
              gam <- params$gam
              all_fold_errors <- numeric(5)

        for (fold in 1:5) {
            train_idx <- c()
            test_idx <- c()

        # Collect the test cells (this fold) and training cells (other folds) from every group
        for (gid in group_ids) {
              folds_for_group <- group_folds[[as.character(gid)]]
              test_idx <- c(test_idx, folds_for_group[[fold]])
              train_idx <- c(train_idx, unlist(folds_for_group[-fold]))
        }

    # Subset
    X_train <- X_simu[train_idx, , drop = FALSE]
    y_train <- y_simu[train_idx]
      y_train = matrix(y_train,ncol=1)
    groups_train <- groups[train_idx]

    X_test <- X_simu[test_idx, , drop = FALSE]
    y_test <- y_simu[test_idx]
      y_test = matrix(y_test,ncol=1)
    groups_test <- groups[test_idx]

    # Fit model
    beta_estimate <- myFuser::fusedLassoProximal(
      X_train, y_train, groups_train,
      lambda = lam, gamma = gam,
      G_simu, intercept = intercept, scaling = TRUE, tol = 1e-6, num.it = 1000
    )

    # Prediction: extract matching columns
    # (each test cell is predicted with the coefficients of its own group)
    col_names <- colnames(X_test)
    y_pred <- diag(X_test[, col_names, drop = FALSE] %*% beta_estimate[col_names, groups_test])
    if(intercept) y_pred = y_pred+beta_estimate[nrow(beta_estimate),groups_test]

    # Compute squared error
    fold_error <- mean((y_test - y_pred)^2)
    all_fold_errors[fold] <- fold_error
  }

  # Average the fold errors and report the CV RMSE for this (lambda, gamma)
  cv_rmse <- sqrt(mean(all_fold_errors))
  list(lambda = lam, gamma = gam, cv_rmse = cv_rmse)
            }, error = function(e) {
    msg <- paste0(
      "[", Sys.time(), "] ",
      "Error for lambda=", lam, ", gamma=", gam, ": ", conditionMessage(e), "\n"
    )
    cat(msg, file = "fuser_error_log.txt", append = TRUE)
    return(NULL)
  })
}
# 2c. Evaluate every (lambda, gamma) combination on the grid
param_grid <- expand.grid(lam = lambda_seq, gam = gamma_seq)

param_list <- split(param_grid, seq(nrow(param_grid)))

# (Run sequentially with lapply: the gene loop is already parallel, so
#  parallelising here as well would oversubscribe the workers.)
cv_results <- lapply(
  param_list,
  function(params) suppressWarnings(process_combination_cv_grouped(params))
)

# 2d. Select the combination with the lowest CV RMSE and refit on all cells
cv_df <- do.call(rbind, lapply(cv_results, as.data.frame))
best_params <- cv_df[which.min(cv_df$cv_rmse), ]
best_lambda <- best_params$lambda
best_gamma <- best_params$gamma
# Final fused lasso estimate (peaks x groups coefficient matrix) -> `fusion`
beta_estimate_final <- myFuser::fusedLassoProximal(
  X_simu, y_simu, groups,
  lambda = best_lambda,
  gamma = best_gamma,
  G_simu,
  intercept = intercept,
  scaling = TRUE,
  tol = 1e-6,
  num.it = 1000
)
    ########################################################
    # STEP 2 (baseline): within-group elastic net (no sharing across groups)
    #   For each group, fit glmnet (alpha = 0.5) on that group's cells only.
    #   lambda is chosen by 5-fold CV when the group has > 3 cells,
    #   otherwise a fixed value (0.002) is used.
    ########################################################
        all_group_result = c()
        for(g in colnames(beta.estimate_simu)){
                x_group = X_simu[groups==g,]
                y_group = y_simu[groups==g]
                if(length(y_group)>3){
                cv_model = cv.glmnet(x_group, y_group, alpha = 0.5,nfolds = 5)
                # NOTE: named `glmnet_lambda` (not `best_lambda`) so it does not
                # overwrite the fused-lasso lambda reported in `fusion_cv`.
                glmnet_lambda <- cv_model$lambda.min
                    }
            else{
                glmnet_lambda = 0.002
            }
                lasso_withingroup = glmnet(x_group, y_group, alpha = 0.5,lambda = glmnet_lambda,intercept = intercept)
                # Column-bind coefficients: one column per group
                all_group_result = cbind(all_group_result,coef(lasso_withingroup))
        }

        ########################################################
        # STEP 3: Group-wise marginal correlation + permutation p-values
        ########################################################
        Y = y_simu[1:length(y_simu)]
        # (Disabled alternative) geographically weighted correlation per peak,
        # averaged within each group:
        #all_peak_cor = mclapply(1:ncol(X_simu),function(j){
        #local_correlations <- sapply(1:length(Y), FUN=function(k) gw_correlation(k, X_simu[,j], Y, W))
        
        #return(local_correlations)
#},mc.cores=2)
        #all_peak_cor = do.call('rbind',all_peak_cor)
        #group_means <- sapply(colnames(beta.estimate_simu), function(g) rowMeans(all_peak_cor[, groups == g, drop = FALSE],na.rm = T))
        # Observed Spearman correlation between each peak and the gene,
        # computed separately within each group -> peaks x groups matrix
        group_means <- sapply(colnames(beta.estimate_simu), function(g) as.vector(cor(X_simu[groups==g,],Y[groups==g],method = 'spearman')))
        rownames(group_means) = colnames(X_simu)

        # Null distribution: 1000 permutations. In each one, gene expression is
        # shuffled among the cells *within each group* (preserving group
        # structure), and the group-wise Spearman correlations are recomputed.
        null_means = lapply(1:1000,function(loop){
            for(i in unique(groups)){
    Y[groups==i] = sample(Y[groups==i],size = sum(groups==i))

    }
        group_means <- sapply(colnames(beta.estimate_simu), function(g) as.vector(cor(X_simu[groups==g,],Y[groups==g],method = 'spearman')))
        rownames(group_means) = colnames(X_simu)
        return(group_means)
        })
        # Convert observed correlations to one-sided p-values against the null
        gwc_p = signac_p(null_means,group_means)

        # Per-gene output:
        #   fusion_old : fused lasso fit with fixed lambda/gamma (Step 2a)
        #   lasso      : within-group elastic-net coefficients (Step 2 baseline)
        #   marginal   : group-wise Spearman correlations (Step 3)
        #   marginal_p : permutation-based p-values for `marginal` (Step 3)
        #   fusion     : final CV-tuned fused lasso coefficients (Step 2d)
        #   fusion_cv  : selected lambda/gamma and the full CV grid results
        return(list(fusion_old=beta.estimate_simu,
                    lasso=all_group_result,
                    marginal = group_means,
                    marginal_p = gwc_p,
                    fusion = beta_estimate_final,
                    fusion_cv = list(lambda = best_lambda,  
                                     gamma = best_gamma,  
                                     cv_df = cv_df
                                    )
                   )
              )
    }, error = function(e) {
        # On failure, log the gene and the error message (appended) and return NULL for this gene
        cat(sprintf("Error in iteration %d (%s): %s\n", i, use_genes[i], e$message),
            file='fuser_region_error.log', append = TRUE)
        return(NULL)
    })
},.options = furrr_options(seed = TRUE))  # seed = TRUE: parallel-safe RNG
names(result) = use_genes

# Grouping information shared by all genes:
#   group_all   = group label of every cell (before dropping small groups)
#   groups      = integer group ID of every retained cell
group_all = group_x
group_x = group_x[is.element(group_x$group, group_count$group), ]
groups = level[group_x$group]
names(groups) = NULL

return(list(coef=result,groups=groups,group_all=group_all,group_count=group_count))
}
