################################################################################
# SAVI_tuning.R  --  Module 1 of the SAVI toolkit: choosing sp_num and time_num
#
# SAVI needs two grouping parameters: sp_num (number of spatial bins along x and y)
# and time_num (number of pseudotime bins). This module recommends them with a
# simulation on the real data:
#
#   1. Take the FIRST gene and its candidate peaks (imputed ATAC = real predictors).
#   2. Declare the first 5 peaks the true enhancers and simulate the gene's
#      expression from them, with enhancer effects that vary smoothly over space and
#      pseudotime (RBF-kernel Gaussian processes, as in the simulation study):
#           y_i = sum_{j=1..5} X_ij * coef_j(group(i)) + noise
#   3. Grid search over sp_num x time_num: group the cells exactly as SAVI does,
#      fit the fused lasso (lambda / gamma chosen by the same grouped 5-fold CV as in
#      SAVI.R), and score the estimated coefficients against the known truth.
#   4. Recommend the (sp_num, time_num) pair with the best evaluation metrics.
#
# It is independent of SAVI.R (the grouping and the CV are re-implemented here with the
# same rules) and only needs the package `myFuser` (fusedLassoProximal), dplyr and
# parallel; ggplot2 for the figure.
#
# MAIN FUNCTIONS
#   SAVI_grouping_summary()  how many groups / cells are kept for each (sp_num, time_num)
#                            (cheap, no model fitting: use it to choose the grid)
#   SAVI_simulate_gene()     the simulator (can be used on its own)
#   SAVI_tune_grid()         the grid search; returns metrics and the recommendation
#   plot_tuning_grid()       metric vs time_num, one panel per sp_num
#
# EVALUATION METRICS (per group, then averaged over groups -- "groupMean" -- and
# over replicates): AUROC and AUPR of the estimated coefficients for recovering the
# true enhancers, and precision / recall / F1 of the selected peaks (selection =
# coefficient >= elbow of the sorted coefficients, the rule used for SAVI in the
# validation). The AUROC/AUPR here are computed in base R (rank statistic; average
# precision) so no extra package is needed.
#
# NOTE: written for the authors' pipeline but not executed where it was written (no R
# available there): run SAVI_tune_grid(..., dry_run = TRUE) first, then a small
# grid with n_rep = 1 before the full search.
################################################################################

suppressPackageStartupMessages({
  library(parallel)
  library(dplyr)
})

# ------------------------------------------------------------------------------
# 1. Grouping: identical to Step 1 of fusing_bench()
#    bins = floor((value - min) / interval); x and y share the same interval
#    (x range / sp_num); groups with more than `min_cells` cells are kept;
#    G = exp(-dist^2) between groups in bin space, entries < exp(-1) set to 0.
# ------------------------------------------------------------------------------
SAVI_group_cells <- function(spatial_coord, t, sp_num, time_num, min_cells = 5) {
  grid <- data.frame(spatial_coord)
  grid$t <- t
  sp_interval <- (max(grid[, 1]) - min(grid[, 1])) / sp_num
  time_interval <- (max(grid$t) - min(grid$t)) / time_num
  cell_x_simu <- floor((grid[, 1] - min(grid[, 1])) / sp_interval)
  cell_y_simu <- floor((grid[, 2] - min(grid[, 2])) / sp_interval)
  cell_t_simu <- floor((grid$t - min(grid$t)) / time_interval)
  group_count <- data.frame(cell_x_simu, cell_y_simu, cell_t_simu) %>%
    count(cell_x_simu, cell_y_simu, cell_t_simu)
  group_count <- group_count[group_count$n > min_cells, ]
  group_count$group <- paste0(group_count$cell_x_simu, "-", group_count$cell_y_simu, "-", group_count$cell_t_simu)
  rownames(group_count) <- seq_len(nrow(group_count))
  cell_group <- paste0(cell_x_simu, "-", cell_y_simu, "-", cell_t_simu)
  keep <- cell_group %in% group_count$group
  level <- setNames(seq_len(nrow(group_count)), group_count$group)
  G <- NULL
  if (nrow(group_count) > 1) {
    # constant columns (e.g. time_num = 1) give NaN after scaling; dist() ignores them
    d <- as.matrix(dist(scale(group_count[, 1:3], center = FALSE)))
    G <- exp(-d^2)
    G[G < exp(-1)] <- 0
    diag(G) <- 0
  }
  list(group_count = group_count, keep = keep,
       groups = unname(level[cell_group[keep]]), G = G)
}

# NOTE: fusing_bench() now receives bin widths (sp_interval, time_interval) instead of
# bin counts. This file works with counts (sp_num, time_num) and converts them:
#   sp_interval = (max(x) - min(x)) / sp_num ,  time_interval = (max(t) - min(t)) / time_num
# which is exactly what SAVI_group_cells() below does, so the groupings are identical.
#
# Table of the groupings a grid would produce on the real data (no model fitting).
# Use it to drop settings with too few / too small groups before running the search.
SAVI_grouping_summary <- function(spatial_coord, t, sp_nums = c(3, 5, 7, 10),
                                  time_nums = c(1, 3, 5, 7), min_cells = 5) {
  do.call(rbind, lapply(sp_nums, function(sp) {
    do.call(rbind, lapply(time_nums, function(tn) {
      g <- SAVI_group_cells(spatial_coord, t, sp, tn, min_cells)
      data.frame(sp_num = sp, time_num = tn, n_groups = nrow(g$group_count),
                 frac_cells_kept = mean(g$keep),
                 median_cells_per_group = if (nrow(g$group_count)) median(g$group_count$n) else NA_real_)
    }))
  }))
}

# ------------------------------------------------------------------------------
# 2. Simulator
#    Real inputs : spatial coordinates, the gene's pseudotime, the imputed ATAC of the
#                  gene's candidate peaks, and the real expression scale of the gene.
#    Simulated   : enhancer effects and the expression y.
#    Effects     : coef_j(cell) = coef_mean + coef_sd * GP_j(x, y, t), GP_j with kernel
#                  exp(-d^2 / length_scale^2) on coordinates rescaled to [0,1]
#                  (sampled with random Fourier features: no n x n matrix needed);
#                  coef * (coef > coef_threshold) as in the simulation study, then
#                  averaged inside the "true" groups (true_sp_num x true_time_num).
#    Truth       : in an evaluated group, peak j is a true enhancer if the group-mean
#                  coefficient > truth_cut; non-enhancers (peaks 6..) have coefficient 0.
#    Scale       : y is rescaled so that the signal has the standard deviation of the
#                  real gene, and noise sd = noise_frac * sd(real gene).
# ------------------------------------------------------------------------------
.rff_field <- function(Z, length_scale, n_features = 300) {
  # one draw of a zero-mean, unit-variance GP with kernel exp(-d^2 / length_scale^2)
  W <- matrix(rnorm(n_features * ncol(Z), 0, sqrt(2) / length_scale), nrow = n_features)
  b <- runif(n_features, 0, 2 * pi)
  a <- rnorm(n_features)
  as.numeric(sqrt(2 / n_features) * cos(Z %*% t(W) + matrix(b, nrow(Z), n_features, byrow = TRUE)) %*% a)
}

SAVI_simulate_gene <- function(X, y_real, spatial_coord, t, n_true = 5,
                               true_sp_num = 6, true_time_num = 4,
                               length_scale = 0.5, coef_mean = 0.3, coef_sd = 0.5,
                               coef_threshold = 0.1, noise_frac = 0.5, seed = 1) {
  set.seed(seed)
  X <- as.matrix(X)
  stopifnot(nrow(X) == nrow(spatial_coord), length(t) == nrow(X), ncol(X) > n_true)
  X <- sweep(X, 2, pmax(apply(X, 2, sd), 1e-8), "/")       # peaks on a common scale (sd = 1)
  sc <- as.matrix(spatial_coord)[, 1:2]
  rs <- function(v) (v - min(v)) / max(diff(range(v)), 1e-8)
  Z <- cbind(rs(sc[, 1]), rs(sc[, 2]), rs(t))
  # cell-level effects of the n_true enhancers
  coef_cell <- sapply(seq_len(n_true), function(j) coef_mean + coef_sd * .rff_field(Z, length_scale))
  coef_cell <- coef_cell * (coef_cell > coef_threshold)
  # average inside the "true" groups (all cells, no size filter)
  tg <-paste0(floor((sc[, 1] - min(sc[, 1])) / ((diff(range(sc[, 1]))) / true_sp_num)), "-",
               floor((sc[, 2] - min(sc[, 2])) / ((diff(range(sc[, 1]))) / true_sp_num)), "-",
               floor((t - min(t)) / (diff(range(t)) / true_time_num)))
  coef_true <- apply(coef_cell, 2, function(v) ave(v, tg, FUN = mean))
  signal <- rowSums(X[, seq_len(n_true), drop = FALSE] * coef_true)
  scale_f <- sd(y_real) / max(sd(signal), 1e-8)
  y <- scale_f * signal + rnorm(nrow(X), 0, noise_frac * sd(y_real))
  list(X = X, y = y, coef_cell = coef_true, n_true = n_true, scale = scale_f)
}

# ------------------------------------------------------------------------------
# 3. Fused-lasso fit with grouped 5-fold CV (same procedure as SAVI.R, Step 2)
# ------------------------------------------------------------------------------
SAVI_default_lambda <- c(seq(0, 0.001, 0.0002), 0.002, 0.005, 0.008, 0.01)
SAVI_default_gamma  <- c(seq(0.0001, 0.001, 0.0002), 0.002, 0.005, 0.008, 0.01)

.savi_cv_fused <- function(X, y, groups, G, lambda_seq, gamma_seq, nfolds = 5,
                           intercept = FALSE, seed = 123) {
  if (!requireNamespace("myFuser", quietly = TRUE)) stop("package 'myFuser' is required (fusedLassoProximal)")
  y <- matrix(y, ncol = 1)
  set.seed(seed)
  group_ids <- unique(groups)
  group_folds <- lapply(group_ids, function(gid) {
    idx <- which(groups == gid)
    split(idx, sample(rep(seq_len(nfolds), length.out = length(idx))))
  })
  names(group_folds) <- group_ids
  grid <- expand.grid(lam = lambda_seq, gam = gamma_seq)
  cv_rmse <- vapply(seq_len(nrow(grid)), function(r) {
    tryCatch({
      errs <- numeric(nfolds)
      for (f in seq_len(nfolds)) {
        test_idx <- unlist(lapply(group_ids, function(gid) group_folds[[as.character(gid)]][[f]]))
        train_idx <- unlist(lapply(group_ids, function(gid) unlist(group_folds[[as.character(gid)]][-f])))
        b <- suppressWarnings(myFuser::fusedLassoProximal(
          X[train_idx, , drop = FALSE], y[train_idx, , drop = FALSE], groups[train_idx],
          lambda = grid$lam[r], gamma = grid$gam[r], G, intercept = intercept,
          scaling = TRUE, tol = 1e-6, num.it = 1000))
        cn <- colnames(X)
        pred <- diag(X[test_idx, cn, drop = FALSE] %*% b[cn, groups[test_idx]])
        if (intercept) pred <- pred + b[nrow(b), groups[test_idx]]
        errs[f] <- mean((y[test_idx] - pred)^2)
      }
      sqrt(mean(errs))
    }, error = function(e) NA_real_)
  }, numeric(1))
  best <- which.min(cv_rmse)
  if (!length(best)) stop("every (lambda, gamma) fit failed")
  beta <- myFuser::fusedLassoProximal(X, y, groups, lambda = grid$lam[best], gamma = grid$gam[best],
                                      G, intercept = intercept, scaling = TRUE, tol = 1e-6, num.it = 1000)
  list(beta = beta, lambda = grid$lam[best], gamma = grid$gam[best], cv_rmse = cv_rmse[best])
}

# ------------------------------------------------------------------------------
# 4. Scoring against the truth
# ------------------------------------------------------------------------------
.elbow <- function(y, threshold) {                    # same rule as get.elbow.points.indices()
  y <- sort(y)
  y[y < median(y)] <- median(y)
  d2 <- diff(diff(y))
  i <- which(abs(d2) > threshold)[1]
  if (any(abs(d2) > threshold)) y[i + 2] else max(y) + 1
}
.auroc <- function(score, label) {
  n1 <- sum(label); n0 <- sum(!label)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  (sum(rank(score)[label]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
.aupr <- function(score, label) {                      # average precision
  if (sum(label) == 0) return(NA_real_)
  o <- order(score, decreasing = TRUE)
  l <- label[o]
  mean((cumsum(l) / seq_along(l))[l])
}

# beta: peaks x groups estimate; truth: logical peaks x groups matrix (same peak order)
.score_groups <- function(beta, truth, elbow_threshold = 0.001) {
  beta <- as.matrix(beta)
  beta <- beta[!(rownames(beta) %in% c("", "(Intercept)")), , drop = FALSE]
  beta <- beta[rownames(truth), , drop = FALSE]
  per <- lapply(seq_len(ncol(beta)), function(j) {
    s <- beta[, j]; l <- truth[, j]
    sel <- s >= .elbow(s, elbow_threshold)
    tp <- sum(sel & l)
    prec <- if (sum(sel) == 0) 0 else tp / sum(sel)
    rec <- if (sum(l) == 0) NA_real_ else tp / sum(l)
    f1 <- if (tp == 0) 0 else 2 * prec * rec / (prec + rec)
    c(auroc = .auroc(s, l), aupr = .aupr(s, l), precision = if (sum(l) == 0) NA_real_ else prec,
      recall = rec, f1 = if (sum(l) == 0) NA_real_ else f1)
  })
  per <- do.call(rbind, per)
  c(colMeans(per, na.rm = TRUE), setNames(apply(per, 2, sd, na.rm = TRUE), paste0(colnames(per), "_sd")))
}

# ------------------------------------------------------------------------------
# 5. Grid search
#   imputation, spatial_coord, pseudoT, use_genes, use_peaks : as in fusing_bench()
#                pseudoT = ONE pseudotime per cell shared by all genes: a numeric vector
#                (or a data frame with column `x`; a legacy list of per-gene data frames
#                with column `x` is also accepted, the simulated gene's element is used)
#   gene         index or name of the gene used for the simulation (default: the first)
#   true_peaks   which of the gene's candidate peaks are the true enhancers (default 1:5)
#   sp_nums, time_nums   the grid
#   n_rep        simulation replicates (different noise / effect draws)
#   criterion    "composite" = mean(AUROC, AUPR, F1) | "auroc" | "aupr" | "f1"
#   min_keep_frac  settings that discard more cells than this (groups with <= 5 cells
#                  are dropped by SAVI) are not recommended
#   tol          settings within `tol` of the best score are treated as ties; the tie is
#                broken in favour of the coarser grouping (fewer groups)
#   sim_args     list of further arguments of SAVI_simulate_gene() (length_scale, noise_frac, ...)
#   dry_run      only return the grouping summary (no fitting)
# Returns list(recommended, summary, results, grouping, simulation_info)
# ------------------------------------------------------------------------------
.savi_pseudotime <- function(pseudoT, g, gi) {
  if (is.numeric(pseudoT)) return(as.numeric(pseudoT))
  if (is.data.frame(pseudoT)) return(as.numeric(pseudoT$x))
  if (is.list(pseudoT)) {                                   # legacy: one data frame per gene
    el <- if (!is.null(names(pseudoT)) && g %in% names(pseudoT)) pseudoT[[g]] else pseudoT[[gi]]
    return(as.numeric(if (is.data.frame(el)) el$x else el))
  }
  stop("pseudoT must be a numeric vector (one value per cell)")
}

SAVI_tune_grid <- function(imputation, spatial_coord, pseudoT, use_genes, use_peaks,
                           gene = 1, true_peaks = 1:5,
                           sp_nums = c(3, 5, 7, 10), time_nums = c(1, 3, 5, 7),
                           n_rep = 2, lambda_seq = SAVI_default_lambda, gamma_seq = SAVI_default_gamma,
                           criterion = c("composite", "auroc", "aupr", "f1"),
                           min_keep_frac = 0.8, tol = 0.01, truth_cut = 0.1,
                           sim_args = list(), n_cores = 1, seed = 2025, dry_run = FALSE) {
  criterion <- match.arg(criterion)
  g <- if (is.numeric(gene)) use_genes[gene] else gene
  gi <- match(g, use_genes)
  pk <- if (!is.null(names(use_peaks)) && g %in% names(use_peaks)) use_peaks[[g]] else use_peaks[[gi]]
  stopifnot(length(true_peaks) <= length(pk))
  t_real <- .savi_pseudotime(pseudoT, g, gi)
  coord <- as.matrix(spatial_coord)[, 1:2]
  X_real <- as.matrix(imputation$atac[, pk, drop = FALSE])
  # the simulator treats the first n_true columns as the true enhancers: reorder if needed
  ord <- c(true_peaks, setdiff(seq_along(pk), true_peaks))
  X_real <- X_real[, ord, drop = FALSE]
  y_real <- as.numeric(imputation$rna[, g])

  grouping <- SAVI_grouping_summary(coord, t_real, sp_nums, time_nums)
  message("Gene used for the simulation: ", g, " | candidate peaks: ", length(pk),
          " | true enhancers: ", length(true_peaks), " | cells: ", nrow(coord))
  if (dry_run) return(list(grouping = grouping))

  combos <- grouping[grouping$frac_cells_kept > 0 & grouping$n_groups > 1, c("sp_num", "time_num")]
  res_all <- list()
  for (rep_i in seq_len(n_rep)) {
    sim <- do.call(SAVI_simulate_gene,
                   c(list(X = X_real, y_real = y_real, spatial_coord = coord, t = t_real,
                          n_true = length(true_peaks), seed = seed + rep_i), sim_args))
    rows <- parallel::mclapply(seq_len(nrow(combos)), function(r) {
      sp <- combos$sp_num[r]; tn <- combos$time_num[r]
      tryCatch({
        gr <- SAVI_group_cells(coord, t_real, sp, tn)
        Xk <- sim$X[gr$keep, , drop = FALSE]
        colnames(Xk) <- paste0("peak", seq_len(ncol(Xk)))
        rownames(Xk) <- seq_len(nrow(Xk))
        yk <- sim$y[gr$keep] + rnorm(sum(gr$keep), 0, 1e-3)   # same tie-breaking noise as SAVI
        fit <- .savi_cv_fused(Xk, yk, gr$groups, gr$G, lambda_seq, gamma_seq)
        # truth of every evaluated group: mean true coefficient of its cells > truth_cut
        cc <- sim$coef_cell[gr$keep, , drop = FALSE]
        tm <- matrix(FALSE, ncol(Xk), max(gr$groups), dimnames = list(colnames(Xk), NULL))
        for (k in seq_len(max(gr$groups))) {
          tm[seq_len(sim$n_true), k] <- colMeans(cc[gr$groups == k, , drop = FALSE]) > truth_cut
        }
        sc <- .score_groups(fit$beta, tm)
        data.frame(rep = rep_i, sp_num = sp, time_num = tn, n_groups = nrow(gr$group_count),
                   frac_cells_kept = mean(gr$keep), lambda = fit$lambda, gamma = fit$gamma,
                   t(sc))
      }, error = function(e) {
        message("sp_num=", sp, " time_num=", tn, " failed: ", conditionMessage(e)); NULL
      })
    }, mc.cores = n_cores)
    res_all[[rep_i]] <- do.call(rbind, rows)
    message("replicate ", rep_i, "/", n_rep, " done")
  }
  results <- do.call(rbind, res_all)
  if (is.null(results)) stop("no setting could be fitted")
  rownames(results) <- NULL

  summ <- results %>%
    group_by(sp_num, time_num) %>%
    summarise(n_groups = first(n_groups), frac_cells_kept = first(frac_cells_kept),
              auroc = mean(auroc, na.rm = TRUE), aupr = mean(aupr, na.rm = TRUE),
              precision = mean(precision, na.rm = TRUE), recall = mean(recall, na.rm = TRUE),
              f1 = mean(f1, na.rm = TRUE),
              auroc_rep_sd = sd(auroc, na.rm = TRUE), aupr_rep_sd = sd(aupr, na.rm = TRUE),
              f1_rep_sd = sd(f1, na.rm = TRUE), .groups = "drop") %>%
    as.data.frame()
  summ$score <- switch(criterion,
                       composite = rowMeans(summ[, c("auroc", "aupr", "f1")]),
                       auroc = summ$auroc, aupr = summ$aupr, f1 = summ$f1)
  summ <- summ[order(-summ$score), ]
  rownames(summ) <- NULL

  ok <- summ[summ$frac_cells_kept >= min_keep_frac & !is.na(summ$score), ]
  if (!nrow(ok)) {
    warning("no setting keeps >= ", min_keep_frac, " of the cells; choosing among all settings")
    ok <- summ[!is.na(summ$score), ]
  }
  near <- ok[ok$score >= max(ok$score) - tol, ]
  rec <- near[order(near$n_groups, -near$score), ][1, ]
  message(sprintf("Recommended: sp_num = %d, time_num = %d  (score %.3f, %d groups, %.0f%% of cells kept)",
                  rec$sp_num, rec$time_num, rec$score, rec$n_groups, 100 * rec$frac_cells_kept))
  # fusing_bench() takes bin WIDTHS: sp_interval = range(x) / sp_num, time_interval = range(t) / time_num
  list(recommended = list(sp_num = rec$sp_num, time_num = rec$time_num,
                          sp_interval = diff(range(coord[, 1])) / rec$sp_num,
                          time_interval = diff(range(t_real)) / rec$time_num),
       recommended_row = rec, summary = summ, results = results, grouping = grouping,
       simulation_info = list(gene = g, n_peaks = length(pk), true_peaks = pk[true_peaks],
                              criterion = criterion, n_rep = n_rep, truth_cut = truth_cut))
}

# metric vs time_num, one panel per sp_num (as in the simulation study); the
# recommended setting is circled
plot_tuning_grid <- function(tuning, metrics = c("auroc", "aupr", "f1")) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required")
  s <- tuning$summary
  long <- do.call(rbind, lapply(metrics, function(m)
    data.frame(sp_num = s$sp_num, time_num = s$time_num, metric = toupper(m), value = s[[m]])))
  rec <- tuning$recommended
  long$rec <- long$sp_num == rec$sp_num & long$time_num == rec$time_num
  ggplot2::ggplot(long, ggplot2::aes(time_num, value, color = metric)) +
    ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::geom_point(data = long[long$rec, ], shape = 1, size = 5, color = "black") +
    ggplot2::facet_wrap(~sp_num, labeller = ggplot2::label_both) +
    ggplot2::labs(x = "time_num (pseudotime bins)", y = "mean over groups and replicates",
                  title = sprintf("Recommended: sp_num = %d, time_num = %d", rec$sp_num, rec$time_num))
}
