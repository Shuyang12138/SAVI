################################################################################
# SAVI_evaluation.R  --  Module 2 of the SAVI toolkit: baselines + evaluation
#
# Independent of SAVI.R: it only reads the object returned by
# fusing_bench(), so it can be sourced next to it.
#
# A. BASELINES  (every baseline is a list named by gene)
#     baseline_global_lasso()   elastic net fitted on all cells            -> "Lasso(g)"
#     baseline_moran()          Moran's R between each peak and the gene   -> "Moran's R"
#     baseline_eSpatial()       eSpatial domain-based peak-gene FDRs       -> "eSpatial"
#     build_link_object() + baseline_signac_links()   correlation p-values  -> "Signac"
#     baseline_cicero()         Cicero co-accessibility with the promoter  -> "Cicero"
#     SAVI_baselines()          wrapper that runs the requested baselines
#
# B. EVALUATION
#     eval_reproducibility()       cross-modality reproducibility (Spearman, Jaccard)
#                                  needs a SECOND MODALITY (e.g. H3K27ac): optional
#     eval_link_reproducibility()  same for the global link baselines (Signac, Cicero)
#     eval_validation_pairs()      known gene-CRE pairs vs. every method's peak scores
#     eval_validation_metrics()    per-gene AUROC / AUPR
#     eval_topk_precision()        precision / recall of the top-k peaks
#     plot_*()                     figures in the style of the paper
#     SAVI_evaluate()              wrapper: skips everything that needs missing inputs
#
# INPUT CONVENTIONS (same as the SAVI example)
#     res        output of fusing_bench(): res$coef[[gene]]$fusion / $lasso /
#                $marginal_p / $group_count ...
#     imputation list(rna = cells x genes, atac = cells x peaks)  (prep$imputed)
#     use_genes  character vector; use_peaks  list of candidate peaks per gene
#     spatial_coord  cells x 2
#
# Required packages: Matrix, glmnet, parallel (base).  Optional, only for the part
# that needs them: pROC + PRROC (validation), ggplot2 (plots), Seurat + Signac +
# GenomicRanges (Signac / Cicero baselines, gene-CRE overlap), cicero + monocle3.
################################################################################

library(Matrix)
library(glmnet)
library(parallel)

`%||%` <- function(a, b) if (is.null(a)) b else a

.need <- function(pkg, why = "") {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required", if (nzchar(why)) paste0(" ", why), ". Please install it.",
         call. = FALSE)
  }
}

# candidate peaks of gene g (list named by gene, or aligned with use_genes)
.peaks_of <- function(use_peaks, g, i) {
  if (!is.null(names(use_peaks)) && g %in% names(use_peaks)) use_peaks[[g]] else use_peaks[[i]]
}
.as_num <- function(x) as.numeric(as.matrix(x))
.rowmax <- function(m) {
  r <- suppressWarnings(apply(m, 1, max, na.rm = TRUE))
  r[!is.finite(r)] <- NA_real_
  r
}

################################################################################
# A. BASELINES
################################################################################

# ------------------------------------------------------------------------------
# A1. Global lasso: elastic net (alpha = 0.5), lambda by 5-fold CV, all cells.
#     Returns, per gene, the sparse coefficient column (first row = intercept).
# ------------------------------------------------------------------------------
baseline_global_lasso <- function(imputation, use_genes, use_peaks, alpha = 0.5,
                                  nfolds = 5, n_cores = 1, seed = 2025) {
  set.seed(seed, kind = "L'Ecuyer-CMRG")
  out <- parallel::mclapply(seq_along(use_genes), function(i) {
    g <- use_genes[i]
    tryCatch({
      y <- imputation$rna[, g]
      X <- as.matrix(imputation$atac[, .peaks_of(use_peaks, g, i), drop = FALSE])
      cv <- glmnet::cv.glmnet(X, y, alpha = alpha, nfolds = nfolds)
      coef(glmnet::glmnet(X, y, alpha = alpha, lambda = cv$lambda.min, intercept = TRUE))
    }, error = function(e) {
      message("global lasso: gene ", g, " skipped (", conditionMessage(e), ")")
      NULL
    })
  }, mc.cores = n_cores)
  names(out) <- use_genes
  out
}

# ------------------------------------------------------------------------------
# A2. Moran's R
#     W = Gaussian kernel of the spatial distances (zero diagonal). The default
#     bandwidth is (range of the first coordinate) / 20.
#     NOTE: the original analysis script used (max(coord[,1]) - min(coord[,2])) / 20,
#     which mixes the two coordinates; pass bandwidth = that value to reproduce it.
#     The variance term uses sum(W %*% t(W)) = sum(colSums(W)^2) (identical value,
#     O(n^2) instead of O(n^3)) and is computed once for all genes.
# ------------------------------------------------------------------------------
moran_weights <- function(spatial_coord, bandwidth = NULL) {
  spatial_coord <- as.matrix(spatial_coord)[, 1:2, drop = FALSE]
  if (is.null(bandwidth)) bandwidth <- (max(spatial_coord[, 1]) - min(spatial_coord[, 1])) / 20
  message("Moran's R: Gaussian kernel bandwidth = ", signif(bandwidth, 4))
  d <- as.matrix(dist(spatial_coord))
  W <- exp(-(d^2) / (2 * bandwidth^2))
  diag(W) <- 0
  W
}

.moran_setup <- function(W) {
  n <- nrow(W)
  W <- W * n / sum(W)
  rs <- rowSums(W)
  cs <- colSums(W)
  varR <- (n^2 * sum(cs^2) - 2 * n * sum(rs * cs) + sum(W)^2) / (n^2 * (n - 1)^2)
  list(W = W, varR = varR)
}

# one-sided p-values of Moran's R between every column of X_mat and Y
moransR <- function(X_mat, Y, W = NULL, setup = NULL) {
  if (is.null(setup)) setup <- .moran_setup(W)
  Yc <- Y - mean(Y)
  dY <- sqrt(sum(Yc^2))
  Wy <- setup$W %*% Yc
  Xc <- sweep(X_mat, 2, colMeans(X_mat))
  dX <- sqrt(colSums(Xc^2))
  R <- colSums(Xc * matrix(Wy, nrow(Xc), ncol(Xc))) / (dX * dY)
  z <- R / sqrt(setup$varR)
  p <- pnorm(z, lower.tail = FALSE)
  p[is.na(p)] <- 1          # constant peaks carry no evidence
  p
}

baseline_moran <- function(imputation, use_genes, use_peaks, spatial_coord,
                           bandwidth = NULL, n_cores = 1) {
  setup <- .moran_setup(moran_weights(spatial_coord, bandwidth))
  out <- parallel::mclapply(seq_along(use_genes), function(i) {
    g <- use_genes[i]
    tryCatch({
      X <- as.matrix(imputation$atac[, .peaks_of(use_peaks, g, i), drop = FALSE])
      p <- moransR(X, imputation$rna[, g], setup = setup)
      names(p) <- colnames(X)
      p
    }, error = function(e) {
      message("Moran's R: gene ", g, " skipped (", conditionMessage(e), ")")
      NULL
    })
  }, mc.cores = n_cores)
  names(out) <- use_genes
  out
}

# ------------------------------------------------------------------------------
# A3. eSpatial (computed beforehand with the eSpatial code; this only reshapes it)
#       GPTab        data frame with Gene, Peak, FDR, class ("corr" rows are used)
#       enh_pattern  data frame with original_type ("domainA_domainB") and peaks
#     the *_h arguments (second modality) are optional. A gene-peak pair is kept
#     only if it has an FDR in every modality supplied.
#     Returns list(atac = ..., h = ...): per gene a peaks x domains matrix with the
#     FDR where the peak belongs to the domain's enhancer set and 1 elsewhere.
# ------------------------------------------------------------------------------
.espatial_domain_peaks <- function(domains, pattern) {
  sel <- setNames(vector("list", length(domains)), domains)
  types <- strsplit(as.character(pattern$original_type), "_")
  peaks <- as.character(pattern$peaks)
  for (p in seq_along(types)) {
    for (d in types[[p]]) if (d %in% domains) sel[[d]] <- c(sel[[d]], peaks[p])
  }
  sel
}

baseline_eSpatial <- function(GPTab, enh_pattern, GPTab_h = NULL, enh_pattern_h = NULL,
                              exclude_domains = "0-radom") {
  # "0-radom" is spelled as in the eSpatial output of the original analysis
  GPTab <- as.data.frame(GPTab)
  GPTab <- GPTab[GPTab$class == "corr", ]
  has_h <- !is.null(GPTab_h) && !is.null(enh_pattern_h)
  domains <- unique(unlist(strsplit(as.character(enh_pattern$original_type), "_")))
  if (has_h) {
    GPTab_h <- as.data.frame(GPTab_h)
    GPTab_h <- GPTab_h[GPTab_h$class == "corr", ]
    domains <- unique(c(domains, unlist(strsplit(as.character(enh_pattern_h$original_type), "_"))))
  }
  domains <- setdiff(domains, exclude_domains)
  sel_a <- .espatial_domain_peaks(domains, enh_pattern)
  sel_h <- if (has_h) .espatial_domain_peaks(domains, enh_pattern_h) else NULL

  if (has_h) {
    joined <- merge(GPTab[, c("Gene", "Peak", "FDR")], GPTab_h[, c("Gene", "Peak", "FDR")],
                    by = c("Gene", "Peak"), suffixes = c(".x", ".y"))
    joined <- joined[!is.na(joined$FDR.x + joined$FDR.y), ]
  } else {
    joined <- GPTab[!is.na(GPTab$FDR), c("Gene", "Peak", "FDR")]
    names(joined)[3] <- "FDR.x"
  }
  idx <- split(seq_len(nrow(joined)), joined$Gene)
  build <- function(fdr_col, sel) {
    lapply(idx, function(ii) {
      a <- joined[[fdr_col]][ii]
      pk <- as.character(joined$Peak[ii])
      m <- matrix(rep(a, length(domains)), nrow = length(a))
      rownames(m) <- pk
      colnames(m) <- domains
      for (d in domains) m[!(pk %in% sel[[d]]), d] <- 1
      m
    })
  }
  list(atac = build("FDR.x", sel_a), h = if (has_h) build("FDR.y", sel_h) else NULL)
}

# ------------------------------------------------------------------------------
# A4. Signac-style links: Pearson correlation p-value between every candidate
#     peak (TF-IDF) and the gene (SCT data), computed on a Seurat object.
#     build_link_object() creates that object from the spatial Seurat object:
#       sp_data     Seurat object with an SCT assay (counts and data)
#       fragments   path to the fragment file, or a Signac Fragment object / list
#       annotation  GRanges gene annotation with gene_name (e.g. Annotation(<chromatin assay>))
#     It also returns the promoter peak of each gene (-2000/+500 bp), needed by Cicero.
# ------------------------------------------------------------------------------
build_link_object <- function(sp_data, use_genes, use_peaks, fragments, annotation,
                              rna_assay = NULL, atac_assay = "ATAC",
                              promoter_upstream = 2000, promoter_downstream = 500,
                              lsi_dims = 2:30) {
  .need("Seurat"); .need("Signac"); .need("GenomicRanges")
  library(Seurat)
  library(Signac)
  prom <- GenomicRanges::promoters(annotation, upstream = promoter_upstream,
                                   downstream = promoter_downstream)
  prom <- prom[prom$gene_name %in% use_genes]
  prom <- prom[!duplicated(prom$gene_name)]
  gene_promoter_peak <- paste0(GenomicRanges::seqnames(prom), "-",
                               GenomicRanges::start(prom), "-", GenomicRanges::end(prom))
  names(gene_promoter_peak) <- prom$gene_name

  all_peaks <- unique(c(unlist(use_peaks), gene_promoter_peak))
  if (is.character(fragments)) {
    fragments <- Signac::CreateFragmentObject(fragments, cells = colnames(sp_data))
  }
  counts <- Signac::FeatureMatrix(fragments = fragments,
                                  features = Signac::StringToGRanges(all_peaks, sep = c("-", "-")),
                                  cells = colnames(sp_data))
  if (nrow(counts) != length(all_peaks)) stop("FeatureMatrix returned an unexpected number of peaks.")
  rownames(counts) <- all_peaks

  # expression assay: SCT when the object has one, otherwise the spatial RNA assay
  # (e.g. "Spatial"), on which SCTransform() is run here to obtain an SCT assay
  if (is.null(rna_assay)) {
    rna_assay <- if ("SCT" %in% Seurat::Assays(sp_data)) "SCT" else
      if ("Spatial" %in% Seurat::Assays(sp_data)) "Spatial" else Seurat::DefaultAssay(sp_data)
  }
  if (!(rna_assay %in% Seurat::Assays(sp_data)))
    stop("assay '", rna_assay, "' not in sp_data; available: ", paste(Seurat::Assays(sp_data), collapse = ", "))
  message("build_link_object: expression assay = ", rna_assay)
  obj <- CreateSeuratObject(counts = GetAssayData(sp_data, assay = rna_assay, layer = "counts"))
  if (rna_assay == "SCT") {
    obj[["SCT"]] <- sp_data[["SCT"]]
  } else {
    obj <- SCTransform(obj, verbose = FALSE)        # builds the "SCT" assay (counts + data) used for the correlations
    rna_assay <- "SCT"
  }
  obj[[atac_assay]] <- CreateChromatinAssay(counts = counts, fragments = fragments)
  DefaultAssay(obj) <- atac_assay
  obj <- RunTFIDF(obj)
  obj <- FindTopFeatures(obj, min.cutoff = "q0")
  obj <- RunSVD(obj, reduction.name = "lsi")
  obj <- RunUMAP(obj, reduction = "lsi", dims = lsi_dims, reduction.name = "umap.atac",
                 reduction.key = "ATACUMAP_")
  DefaultAssay(obj) <- rna_assay
  list(obj = obj, gene_promoter_peak = gene_promoter_peak,
       peak_assay = atac_assay, rna_assay = rna_assay)
}

baseline_signac_links <- function(link_object, use_genes, use_peaks,
                                  peak_layer = "data", rna_layer = "data") {
  .need("Seurat")
  library(Seurat)
  obj <- link_object$obj
  expr_mat <- GetAssayData(obj, assay = link_object$rna_assay, layer = rna_layer)
  peak_mat <- GetAssayData(obj, assay = link_object$peak_assay, layer = peak_layer)
  out <- setNames(vector("list", length(use_genes)), use_genes)
  for (i in seq_along(use_genes)) {
    g <- use_genes[i]
    if (!g %in% rownames(expr_mat)) next
    peaks <- intersect(.peaks_of(use_peaks, g, i), rownames(peak_mat))
    if (length(peaks) == 0) next
    expr <- as.numeric(expr_mat[g, ])
    pm <- peak_mat[peaks, , drop = FALSE]
    pv <- vapply(seq_along(peaks), function(j) {
      tryCatch(suppressWarnings(stats::cor.test(as.numeric(pm[j, ]), expr)$p.value),
               error = function(e) NA_real_)
    }, numeric(1))
    names(pv) <- peaks
    out[[g]] <- pv
  }
  out
}

# ------------------------------------------------------------------------------
# A5. Cicero: co-accessibility between each candidate peak and the gene promoter peak
#     (needs the packages cicero and monocle3). Peaks without a connection get 0.
# ------------------------------------------------------------------------------
.seurat_to_cicero_input <- function(mat) {
  mat <- as(mat, "dgCMatrix")
  idx <- Matrix::summary(mat)
  data.frame(Peak = rownames(mat)[idx$i], Cell = colnames(mat)[idx$j], Count = idx$x,
             stringsAsFactors = FALSE)
}

.make_genomic_coords <- function(peak_names) {
  .need("Signac")
  gr <- Signac::StringToGRanges(peak_names, sep = c("-", "-"))
  gc <- data.frame(chromosome = as.character(GenomicRanges::seqnames(gr)),
                   length = GenomicRanges::end(gr), stringsAsFactors = FALSE)
  gc <- aggregate(length ~ chromosome, gc, max)
  gc <- gc[order(gc$chromosome), ]
  rownames(gc) <- NULL
  gc
}

baseline_cicero <- function(link_object, use_genes, use_peaks, genomic_coords = NULL,
                            umap_name = "umap.atac") {
  .need("cicero"); .need("monocle3"); .need("Seurat")
  library(Seurat)
  obj <- link_object$obj
  mat <- GetAssayData(obj, assay = link_object$peak_assay, layer = "counts")
  input <- .seurat_to_cicero_input(mat)
  cds <- cicero::make_atac_cds(input)
  coords <- Embeddings(obj, umap_name)
  coords <- coords[intersect(unique(input$Cell), rownames(coords)), , drop = FALSE]
  cicero_cds <- cicero::make_cicero_cds(cds, reduced_coordinates = coords)
  if (is.null(genomic_coords)) genomic_coords <- .make_genomic_coords(rownames(mat))
  conns <- cicero::run_cicero(cicero_cds, genomic_coords = genomic_coords)
  gpp <- link_object$gene_promoter_peak
  out <- lapply(seq_along(use_genes), function(i) {
    g <- use_genes[i]
    peaks <- .peaks_of(use_peaks, g, i)
    vals <- setNames(rep(0, length(peaks)), peaks)
    pp <- gpp[g]
    if (!is.na(pp)) {
      idx <- (conns$Peak1 == pp & conns$Peak2 %in% peaks) | (conns$Peak2 == pp & conns$Peak1 %in% peaks)
      sub <- conns[which(idx), ]
      if (nrow(sub) > 0) {
        nm <- ifelse(sub$Peak1 == pp, sub$Peak2, sub$Peak1)
        v <- sub$coaccess
        v[is.na(v)] <- 0
        vals[nm] <- v
      }
    }
    vals
  })
  names(out) <- use_genes
  out
}

# ------------------------------------------------------------------------------
# A6. Wrapper
#       methods         any of "global_lasso", "moran", "eSpatial", "signac", "cicero"
#       eSpatial_input  list of the arguments of baseline_eSpatial()
#       link_object     output of build_link_object() (for "signac" and "cicero")
#     eSpatial for the second modality is stored as $eSpatial_h (see SAVI_evaluate).
# ------------------------------------------------------------------------------
SAVI_baselines <- function(imputation, use_genes, use_peaks, spatial_coord,
                           methods = c("global_lasso", "moran"),
                           eSpatial_input = NULL, link_object = NULL,
                           moran_bandwidth = NULL, n_cores = 1) {
  out <- list()
  if ("global_lasso" %in% methods) {
    message("baseline: global lasso")
    out$global_lasso <- baseline_global_lasso(imputation, use_genes, use_peaks, n_cores = n_cores)
  }
  if ("moran" %in% methods) {
    message("baseline: Moran's R")
    out$moran <- baseline_moran(imputation, use_genes, use_peaks, spatial_coord,
                                bandwidth = moran_bandwidth, n_cores = n_cores)
  }
  if ("eSpatial" %in% methods) {
    if (is.null(eSpatial_input)) stop("methods includes 'eSpatial' but eSpatial_input is NULL.")
    message("baseline: eSpatial")
    esp <- do.call(baseline_eSpatial, eSpatial_input)
    out$eSpatial <- esp$atac
    if (!is.null(esp$h)) out$eSpatial_h <- esp$h
  }
  if (any(c("signac", "cicero") %in% methods) && is.null(link_object)) {
    stop("methods includes 'signac' / 'cicero' but link_object is NULL (see build_link_object()).")
  }
  if ("signac" %in% methods) {
    message("baseline: Signac links")
    out$signac <- baseline_signac_links(link_object, use_genes, use_peaks)
  }
  if ("cicero" %in% methods) {
    message("baseline: Cicero")
    out$cicero <- baseline_cicero(link_object, use_genes, use_peaks)
  }
  out
}

################################################################################
# B. EVALUATION
################################################################################

# ---- helpers copied from the analysis scripts --------------------------------
get.elbow.points.indices <- function(y, threshold) {
  y <- sort(y)
  y[y < median(y)] <- median(y)
  d1 <- diff(y)       # first derivative
  d2 <- diff(d1)      # second derivative
  indices <- which(abs(d2) > threshold)[1]
  return(ifelse(sum(abs(d2) > threshold) > 0, y[indices + 2], max(y) + 1))
}

# Jaccard index of the selected peaks of two vectors
#   elbow = c(1,1)      -> SAVI fusion coefficients: elbow thresholds (see get.elbow.points.indices)
#   elbow = c(1e-5,1e-5)-> coefficients: selected if >= 1e-5
#   elbow = NULL        -> p-values: BH-adjusted p < 0.05 (eSpatial = TRUE: raw FDR < 0.05)
jaccard_similarity <- function(vec1, vec2, elbow = NULL, eSpatial = FALSE) {
  if (!is.null(elbow)) {
    if (elbow[1] > 1e-4) {
      elbow <- c(get.elbow.points.indices(sort(vec1), 0.001),
                 get.elbow.points.indices(sort(vec2), 0.001))
    }
    vec1 <- vec1 >= elbow[1]
    vec2 <- vec2 >= elbow[2]
  } else if (eSpatial) {
    vec1 <- vec1 < 0.05
    vec2 <- vec2 < 0.05
  } else {
    vec1 <- p.adjust(vec1, method = "BH") < 0.05
    vec2 <- p.adjust(vec2, method = "BH") < 0.05
  }
  if (length(vec1) != length(vec2)) stop("Vectors must be of the same length")
  intersection <- sum(vec1 & vec2, na.rm = TRUE)
  union <- sum(vec1 | vec2, na.rm = TRUE)
  if (union == 0) return(0)
  intersection / union
}

# Spearman correlation and Jaccard index between the results of two modalities,
# averaged over the groups they share (columns picked by atac_use_level / h_use_level)
#   setting "fusion" drops the last row (intercept), "lasso" drops the first row
Eva_reprod_func <- function(atac_use_level, h3k27ac_use_level, atac_result, h3k27ac_result,
                            setting, elbow = NULL, eSpatial = FALSE) {
  use_cell <- if (setting == "fusion") -nrow(atac_result)
              else if (setting == "lasso") -1
              else seq_len(nrow(atac_result))
  atac_result <- atac_result[use_cell, , drop = FALSE]
  h3k27ac_result <- h3k27ac_result[use_cell, , drop = FALSE]
  n <- length(atac_use_level)
  spear <- suppressWarnings(vapply(seq_len(n), function(i) {
    stats::cor(atac_result[, atac_use_level[i]], h3k27ac_result[, h3k27ac_use_level[i]],
               method = "spearman")
  }, numeric(1)))
  jac <- suppressWarnings(vapply(seq_len(n), function(i) {
    jaccard_similarity(atac_result[, atac_use_level[i]], h3k27ac_result[, h3k27ac_use_level[i]],
                       elbow, eSpatial = eSpatial)
  }, numeric(1)))
  data.frame(corre = mean(spear, na.rm = TRUE),
             corre_sd = sqrt(var(spear, na.rm = TRUE) / length(spear)),
             corre_pos = mean(jac, na.rm = TRUE),
             corre_pos_sd = sqrt(var(jac, na.rm = TRUE) / length(jac)))
}

.get_group_count <- function(result, g) {
  gc <- result$coef[[g]]$group_count
  if (is.null(gc)) gc <- result$group_count     # fusing_bench(): one grouping for all genes
  gc
}

# matrix with the same vector repeated in every group column (global methods)
.const_matrix <- function(v, ncol) matrix(rep(.as_num(v), ncol), nrow = length(v))

# ---- B1. cross-modality reproducibility (needs a second modality) ------------
.reprod_one_gene <- function(g, atac_result, h_result, bl, bl_h, shift_x, shift_y, intercept) {
  a <- atac_result$coef[[g]]
  h <- h_result$coef[[g]]
  if (length(a) < 5 || length(h) < 5) return(NULL)        # SAVI failed for this gene
  ga <- .get_group_count(atac_result, g)
  gh <- .get_group_count(h_result, g)
  gh$shifted_group <- paste0(gh$cell_x_simu + shift_x, "-", gh$cell_y_simu + shift_y, "-", gh$cell_t_simu)
  use_group <- intersect(ga$group, gh$shifted_group)
  a_level <- setNames(seq_len(nrow(ga)), ga$group)
  h_level <- setNames(seq_len(nrow(gh)), gh$shifted_group)
  a_use <- a_level[use_group]
  h_use <- h_level[use_group]
  if (nrow(a$fusion) != nrow(h$fusion)) {
    stop("the two modalities have a different number of candidate peaks for this gene ",
         "(rows are matched by position)")
  }
  rows <- list()
  # SAVI fusion: elbow thresholds; the last row is the intercept only if intercept = TRUE
  r <- Eva_reprod_func(a_use, h_use, a$fusion, h$fusion, if (intercept) "fusion" else "fusion_", c(1, 1))
  rows$fusion <- cbind(r, method = "fusion")
  r <- Eva_reprod_func(a_use, h_use, a$lasso, h$lasso, "lasso", c(1e-5, 1e-5))
  rows$lasso <- cbind(r, method = "lasso")
  r <- Eva_reprod_func(a_use, h_use, a$marginal_p, h$marginal_p, "gwc")
  rows$marginal <- cbind(r, method = "marginal")

  if (!is.null(bl$eSpatial[[g]]) && !is.null(bl_h$eSpatial[[g]])) {
    ea <- bl$eSpatial[[g]]
    eh <- bl_h$eSpatial[[g]]
    r <- Eva_reprod_func(colnames(ea), colnames(eh), ea, eh, "gwc", eSpatial = TRUE)
    rows$eSpatial <- cbind(r, method = "eSpatial")
  }
  if (!is.null(bl$global_lasso[[g]]) && !is.null(bl_h$global_lasso[[g]])) {
    r <- Eva_reprod_func(a_use, h_use, .const_matrix(bl$global_lasso[[g]], nrow(ga)),
                         .const_matrix(bl_h$global_lasso[[g]], nrow(gh)), "lasso", c(1e-5, 1e-5))
    rows$global_lasso <- cbind(r, method = "global_lasso")
  }
  if (!is.null(bl$moran[[g]]) && !is.null(bl_h$moran[[g]])) {
    r <- Eva_reprod_func(a_use, h_use, .const_matrix(bl$moran[[g]], nrow(ga)),
                         .const_matrix(bl_h$moran[[g]], nrow(gh)), "gwc")
    rows$moran <- cbind(r, method = "Moran's I")
  }
  out <- do.call(rbind, rows)
  cbind(gene = g, out)
}

# atac_result / h_result : SAVI results of the two modalities (same genes)
# baselines / baselines_h: lists made by SAVI_baselines() for each modality
#                          (eSpatial of the 2nd modality is taken from baselines$eSpatial_h)
# shift_x, shift_y       : offset between the spatial grids of the two modalities
# intercept              : TRUE only if SAVI was run with intercept = TRUE
# Returns a data frame: one row per gene and method with corre (Spearman) and
# corre_pos (Jaccard index), averaged over the shared groups.
eval_reproducibility <- function(atac_result, h_result, baselines = list(), baselines_h = list(),
                                 shift_x = 0, shift_y = 0, intercept = FALSE,
                                 n_cores = 1, genes = NULL) {
  genes <- intersect(genes %||% names(atac_result$coef), names(h_result$coef))
  baselines_h$eSpatial <- baselines_h$eSpatial %||% baselines$eSpatial_h
  rows <- parallel::mclapply(genes, function(g) {
    tryCatch(.reprod_one_gene(g, atac_result, h_result, baselines, baselines_h,
                              shift_x, shift_y, intercept),
             error = function(e) {
               message("reproducibility: gene ", g, " skipped (", conditionMessage(e), ")")
               NULL
             })
  }, mc.cores = n_cores)
  out <- do.call(rbind, rows)
  if (is.null(out)) {
    warning("no gene could be evaluated")
    return(NULL)
  }
  rownames(out) <- NULL
  out
}

# Reproducibility of the global link baselines (Signac p-values, Cicero co-accessibility):
# one vector per gene and modality, peaks matched by position in the candidate lists.
#   type "pvalue": Jaccard of BH-adjusted p < 0.05;  "score": Jaccard of score >= score_cutoff
# (written to replace Eva_link_reprod(), which was not part of the shared scripts)
.align_scores <- function(s, peaks, fill) {
  v <- rep(fill, length(peaks))
  names(v) <- peaks
  if (length(s)) {
    k <- intersect(names(s), peaks)
    v[k] <- s[k]
  }
  v
}

eval_link_reproducibility <- function(scores_a, scores_h, peaks_a, peaks_h, method,
                                      type = c("pvalue", "score"), score_cutoff = 0.25) {
  type <- match.arg(type)
  fill <- if (type == "pvalue") 1 else 0
  genes <- intersect(names(scores_a), names(scores_h))
  rows <- lapply(genes, function(g) {
    if (is.null(peaks_a[[g]]) || is.null(peaks_h[[g]])) return(NULL)
    a <- .align_scores(scores_a[[g]], peaks_a[[g]], fill)
    b <- .align_scores(scores_h[[g]], peaks_h[[g]], fill)
    if (length(a) != length(b) || length(a) < 3) return(NULL)
    corre <- suppressWarnings(stats::cor(a, b, method = "spearman"))
    if (type == "pvalue") {
      ja <- p.adjust(a, "BH") < 0.05
      jb <- p.adjust(b, "BH") < 0.05
    } else {
      ja <- a >= score_cutoff
      jb <- b >= score_cutoff
    }
    u <- sum(ja | jb, na.rm = TRUE)
    data.frame(gene = g, corre = corre, corre_sd = NA_real_,
               corre_pos = if (u == 0) 0 else sum(ja & jb, na.rm = TRUE) / u,
               corre_pos_sd = NA_real_, method = method)
  })
  out <- do.call(rbind, rows)
  if (!is.null(out)) rownames(out) <- NULL
  out
}

# ---- B2. validation against known gene-CRE pairs ------------------------------
# kind "coef"   : higher coefficient = stronger link; selected if >= cutoff in some group
# kind "pvalue" : BH-adjusted p-value per group; selected if <= cutoff in some group
# kind "fusion" : SAVI coefficients; per group the elbow of the sorted coefficients
#                 is the selection threshold (cutoff = threshold of the 2nd derivative)
# The ranking statistic (used for AUROC / AUPR / top-k) is the best value over groups:
# max coefficient, or -(smallest adjusted p-value).
.method_cfg <- list(
  fusion       = list(kind = "fusion", cutoff = 0.001, label = "SAVI"),
  lasso        = list(kind = "coef",   cutoff = 1e-4,  label = "Lasso(w)"),
  marginal     = list(kind = "pvalue", cutoff = 1e-6,  label = "Marginal"),
  global_lasso = list(kind = "coef",   cutoff = 1e-4,  label = "Lasso(g)"),
  eSpatial     = list(kind = "pvalue", cutoff = 0.05,  label = "eSpatial"),
  moran        = list(kind = "pvalue", cutoff = 1e-6,  label = "Moran's R"),
  signac       = list(kind = "pvalue", cutoff = 1e-6,  label = "Signac"),
  cicero       = list(kind = "coef",   cutoff = 0.25,  label = "Cicero")
)
.method_colors <- c("SAVI" = "#ec6bff", "Marginal" = "#04adf0", "Lasso(w)" = "#00c179",
                    "Lasso(g)" = "#a7a501", "eSpatial" = "#ed7469", "Moran's R" = "#ff9f1c",
                    "Cicero" = "#7b2cbf", "Signac" = "#1d4ed8")
.method_label <- function(m) {
  map <- c(fusion = "SAVI", lasso = "Lasso(w)", marginal = "Marginal", global_lasso = "Lasso(g)",
           eSpatial = "eSpatial", moran = "Moran's R", signac = "Signac", cicero = "Cicero",
           "Moran's I" = "Moran's R")
  out <- unname(map[m])
  out[is.na(out)] <- m[is.na(out)]
  out
}

.select_stat <- function(mat, kind, cutoff) {
  mat <- as.matrix(mat)
  if (length(mat) == 0 || is.null(rownames(mat))) return(NULL)
  keep <- !is.na(rownames(mat)) & rownames(mat) != "" & rownames(mat) != "(Intercept)"
  mat <- mat[keep, , drop = FALSE]
  peaks <- rownames(mat)
  if (kind == "pvalue") {
    for (j in seq_len(ncol(mat))) mat[, j] <- -p.adjust(mat[, j], method = "BH")
    sel <- rowSums(mat >= -cutoff, na.rm = TRUE) > 0
  } else if (kind == "fusion") {
    sel <- rep(FALSE, nrow(mat))
    for (j in seq_len(ncol(mat))) {
      thr <- get.elbow.points.indices(mat[, j], cutoff)
      sel <- sel | (mat[, j] >= thr) %in% TRUE
    }
  } else {
    sel <- rowSums(mat >= cutoff, na.rm = TRUE) > 0
  }
  data.frame(select = sel, stat = .rowmax(mat), row.names = peaks)
}

.method_matrix <- function(method, g, res, baselines) {
  switch(method,
         fusion   = res$coef[[g]]$fusion,
         lasso    = res$coef[[g]]$lasso,
         marginal = res$coef[[g]]$marginal_p,
         baselines[[method]][[g]])
}

# One row per gene-peak pair: the selection flag and the ranking statistic of every
# method (columns <method>_select and <method>), and geneCRE = the peak overlaps a
# known CRE of that gene.
#   gene_cre  data frame with columns Coordinate1, Gene, Coordinate2 (chr_start_end)
#   baselines list from SAVI_baselines(); every method present is evaluated
# Rows with a missing statistic in any method are dropped so that all methods are
# compared on the same gene-peak pairs.
eval_validation_pairs <- function(res, use_genes, use_peaks, gene_cre, baselines = list(),
                                  genes = NULL, cfg = .method_cfg) {
  .need("Signac"); .need("GenomicRanges"); .need("IRanges")
  gene_cre <- unique(as.data.frame(gene_cre)[, c("Coordinate1", "Gene", "Coordinate2")])
  methods <- c("fusion", "lasso", "marginal", intersect(names(baselines), names(cfg)))
  cre1 <- Signac::StringToGRanges(gene_cre$Coordinate1, sep = c("_", "_"))
  cre2 <- Signac::StringToGRanges(gene_cre$Coordinate2, sep = c("_", "_"))
  genes_all <- intersect(gene_cre$Gene, names(res$coef))
  genes_all <- genes_all[vapply(genes_all, function(g) length(res$coef[[g]]) >= 5, logical(1))]
  if (!is.null(genes)) genes_all <- intersect(genes_all, genes)
  out <- lapply(genes_all, function(g) {
    peaks <- .peaks_of(use_peaks, g, match(g, use_genes))
    df <- data.frame(gene = g, peak = peaks, stringsAsFactors = FALSE)
    for (m in methods) {
      mat <- .method_matrix(m, g, res, baselines)
      ss <- if (is.null(mat)) NULL else .select_stat(mat, cfg[[m]]$kind, cfg[[m]]$cutoff)
      if (is.null(ss)) {
        df[[paste0(m, "_select")]] <- NA
        df[[m]] <- NA_real_
      } else {
        df[[paste0(m, "_select")]] <- ss[df$peak, "select"]
        df[[m]] <- ss[df$peak, "stat"]
      }
    }
    sub <- gene_cre$Gene == g
    peak_gr <- Signac::StringToGRanges(df$peak, sep = c("-", "-"))
    df$geneCRE <- IRanges::overlapsAny(peak_gr, cre1[sub]) |
                  IRanges::overlapsAny(peak_gr, cre2[sub])
    df
  })
  pairs <- do.call(rbind, out)
  if (is.null(pairs)) stop("no gene of the SAVI result appears in gene_cre")
  pairs <- pairs[stats::complete.cases(pairs[, methods, drop = FALSE]), ]
  rownames(pairs) <- NULL
  pairs
}

.pair_methods <- function(pairs, cfg = .method_cfg) intersect(names(cfg), colnames(pairs))

# genes with at least one known CRE peak and one other peak
.valid_genes <- function(pairs) {
  ok <- tapply(pairs$geneCRE, pairs$gene, function(x) any(x) && !all(x))
  names(ok)[ok %in% TRUE]
}

# AUROC and AUPR per gene and method (pROC / PRROC, as in the analysis scripts)
eval_validation_metrics <- function(pairs, genes = NULL, cfg = .method_cfg) {
  .need("pROC"); .need("PRROC")
  genes <- genes %||% .valid_genes(pairs)
  genes <- intersect(genes, .valid_genes(pairs))
  methods <- .pair_methods(pairs, cfg)
  rows <- lapply(genes, function(g) {
    sub <- pairs$gene == g
    y <- pairs$geneCRE[sub]
    do.call(rbind, lapply(methods, function(m) {
      s <- pairs[[m]][sub]
      auroc <- tryCatch(as.numeric(pROC::auc(pROC::roc(y, s, direction = "<", quiet = TRUE))),
                        error = function(e) NA_real_)
      aupr <- tryCatch(PRROC::pr.curve(scores.class0 = s, weights.class0 = as.numeric(y))$auc.integral,
                       error = function(e) NA_real_)
      data.frame(gene = g, method = m, Method = .method_label(m), auroc = auroc, aupr = aupr)
    }))
  })
  out <- do.call(rbind, rows)
  if (!is.null(out)) rownames(out) <- NULL
  out
}

# precision and recall of the k best-ranked peaks of one gene
precision_recall_curve <- function(scores, labels) {
  ord <- order(scores, decreasing = TRUE)
  labels_sorted <- labels[ord]
  tp <- cumsum(labels_sorted == 1)
  k_vals <- seq_along(labels_sorted)
  data.frame(k = k_vals, precision = tp / k_vals, recall = tp / sum(labels == 1))
}

# precision / recall at k = 1.. averaged over genes (a gene stops contributing beyond its number of peaks)
eval_topk_precision <- function(pairs, genes = NULL, cfg = .method_cfg) {
  genes <- genes %||% .valid_genes(pairs)
  genes <- intersect(genes, .valid_genes(pairs))
  methods <- .pair_methods(pairs, cfg)
  rows <- lapply(genes, function(g) {
    sub <- pairs$gene == g
    do.call(rbind, lapply(methods, function(m) {
      r <- precision_recall_curve(pairs[[m]][sub], pairs$geneCRE[sub])
      r$method <- m
      r$Method <- .method_label(m)
      r
    }))
  })
  all <- do.call(rbind, rows)
  if (is.null(all)) return(NULL)
  aggregate(cbind(precision, recall) ~ k + method + Method, data = all, FUN = mean)
}

# ---- B3. figures ---------------------------------------------------------------
plot_reproducibility <- function(repro) {
  .need("ggplot2")
  lab <- .method_label(repro$method)
  lv <- names(.method_colors)[names(.method_colors) %in% lab]
  d <- data.frame(Method = factor(rep(lab, 2), levels = lv),
                  Metric = rep(c("Spearman Correlation", "Jaccard Index"), each = nrow(repro)),
                  value = c(repro$corre, repro$corre_pos))
  ggplot2::ggplot(d, ggplot2::aes(x = Method, y = value, fill = Method)) +
    ggplot2::geom_boxplot() + ggplot2::facet_wrap(~Metric) +
    ggplot2::scale_fill_manual(values = .method_colors) +
    ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::ylab("Reproducibility Metrics")
}

plot_validation_metrics <- function(metrics) {
  .need("ggplot2")
  lv <- names(.method_colors)[names(.method_colors) %in% metrics$Method]
  d <- data.frame(Method = factor(rep(metrics$Method, 2), levels = lv),
                  Metric = rep(c("AUROC", "AUPR"), each = nrow(metrics)),
                  value = c(metrics$auroc, metrics$aupr))
  ggplot2::ggplot(d, ggplot2::aes(x = Method, y = value, fill = Method)) +
    ggplot2::geom_boxplot() + ggplot2::facet_wrap(~Metric, scales = "free_y") +
    ggplot2::scale_fill_manual(values = .method_colors) +
    ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::ylab("Agreement with known gene-CRE pairs")
}

plot_topk_precision <- function(topk, min_k = 2, max_k = NULL) {
  .need("ggplot2")
  d <- topk[topk$k >= min_k & (is.null(max_k) | topk$k <= (max_k %||% Inf)), ]
  ggplot2::ggplot(d, ggplot2::aes(x = k, y = precision, color = Method)) +
    ggplot2::geom_line() + ggplot2::geom_point(size = 0.5) +
    ggplot2::scale_color_manual(values = .method_colors) +
    ggplot2::labs(x = "Number of Selected Peaks", y = "Precision")
}

# ---- B4. wrapper ------------------------------------------------------------------
# Runs every evaluation whose inputs are available and says what it skipped.
#   gene_cre        known gene-CRE pairs (Coordinate1, Gene, Coordinate2)   -> validation
#   res_h           SAVI result of a second modality (e.g. H3K27ac)         -> reproducibility
#   baselines_h     baselines computed on that modality
#   use_peaks_h     candidate peaks of that modality (for Signac / Cicero reproducibility)
# Returns a list: validation_pairs, validation_metrics, topk_precision,
# reproducibility, link_reproducibility, plots.
SAVI_evaluate <- function(res, use_genes, use_peaks, baselines = list(),
                          gene_cre = NULL, validation_genes = NULL,
                          res_h = NULL, baselines_h = list(), use_peaks_h = NULL,
                          h_shift = c(0, 0), intercept = FALSE, n_cores = 1, make_plots = TRUE) {
  out <- list()
  if (!is.null(gene_cre)) {
    out$validation_pairs <- eval_validation_pairs(res, use_genes, use_peaks, gene_cre, baselines)
    out$validation_metrics <- eval_validation_metrics(out$validation_pairs, genes = validation_genes)
    out$topk_precision <- eval_topk_precision(out$validation_pairs, genes = validation_genes)
  } else {
    message("No gene_cre: skipping validation against known gene-CRE pairs.")
  }
  if (!is.null(res_h)) {
    out$reproducibility <- eval_reproducibility(res, res_h, baselines, baselines_h,
                                                shift_x = h_shift[1], shift_y = h_shift[2],
                                                intercept = intercept, n_cores = n_cores)
    if (!is.null(use_peaks_h)) {
      link <- list()
      if (!is.null(baselines$signac) && !is.null(baselines_h$signac)) {
        link$signac <- eval_link_reproducibility(baselines$signac, baselines_h$signac,
                                                 use_peaks, use_peaks_h, "Signac", "pvalue")
      }
      if (!is.null(baselines$cicero) && !is.null(baselines_h$cicero)) {
        link$cicero <- eval_link_reproducibility(baselines$cicero, baselines_h$cicero,
                                                 use_peaks, use_peaks_h, "Cicero", "score")
      }
      if (length(link)) out$link_reproducibility <- do.call(rbind, link)
    }
  } else {
    message("No second modality (res_h): skipping cross-modality reproducibility.")
  }
  if (make_plots && requireNamespace("ggplot2", quietly = TRUE)) {
    out$plots <- list()
    rep_all <- rbind(out$reproducibility[, c("gene", "method", "corre", "corre_pos")],
                     out$link_reproducibility[, c("gene", "method", "corre", "corre_pos")])
    if (!is.null(rep_all)) out$plots$reproducibility <- plot_reproducibility(rep_all)
    if (!is.null(out$validation_metrics)) out$plots$validation <- plot_validation_metrics(out$validation_metrics)
    if (!is.null(out$topk_precision)) out$plots$topk <- plot_topk_precision(out$topk_precision)
  }
  out
}

################################################################################
# C. SPLIT-HALF REPRODUCIBILITY (needs no second modality)
#
# The cells are split at random into two equal halves, the SAME analysis (SAVI plus
# the baselines) is run on each half, and the two results are compared with the
# cross-modality reproducibility code (eval_reproducibility(): Spearman correlation and
# Jaccard index of the coefficients / p-values over the groups the halves share;
# eval_link_reproducibility() for Signac / Cicero).
#
#   imputation, spatial_coord, pseudoT, use_genes, use_peaks, cell_types : as for fusing_bench()
#   sp_interval, time_interval : bin widths, the same for both halves
#   methods      baselines to compute on each half: "global_lasso", "moran" and, when
#                link_object is given, "signac" / "cicero". (eSpatial has to be run
#                with the eSpatial code on each half and is not computed here.)
#   link_object  OPTIONAL output of build_link_object() on all cells; it is subset to each
#                half (spatial_coord must then have the cell names as rownames)
#   n_splits     number of independent random splits (results are stacked, column `split`)
#   fit_fn       the SAVI function (default: fusing_bench); fit_args: extra arguments for it
#                (e.g. list(intercept = TRUE))
# Use future::plan(multisession) beforehand, as for the normal SAVI run.
#
# Note: bins are anchored at the minimum x / y / pseudotime of the data given to
# fusing_bench(), which differs slightly between the two halves, so group labels can be
# shifted by one bin at the edges; groups are matched by label as in the original code.
# Returns list(reproducibility, summary, splits, plot).
################################################################################
SAVI_split_reproducibility <- function(imputation, spatial_coord, pseudoT, use_genes, use_peaks,
                                       cell_types = NULL, sp_interval, time_interval,
                                       methods = c("global_lasso", "moran"),
                                       link_object = NULL, n_splits = 1, seed = 2025,
                                       moran_bandwidth = NULL, fit_fn = NULL, fit_args = list(),
                                       intercept = FALSE, n_cores = 1) {
  fit_fn <- fit_fn %||% get("fusing_bench", mode = "function")
  coord <- as.matrix(spatial_coord)[, 1:2, drop = FALSE]
  n <- nrow(coord)
  stopifnot(length(pseudoT) == n, nrow(imputation$rna) == n, nrow(imputation$atac) == n)
  if (is.null(names(use_peaks))) names(use_peaks) <- use_genes
  if ("eSpatial" %in% methods) {
    message("split-half: eSpatial is not computed here (run the eSpatial code on each half); skipping it")
    methods <- setdiff(methods, "eSpatial")
  }
  use_link <- !is.null(link_object) && any(c("signac", "cicero") %in% methods)
  if (!use_link) methods <- setdiff(methods, c("signac", "cicero"))
  if (use_link && is.null(rownames(coord)))
    stop("spatial_coord needs the cell names as rownames to subset link_object")
  # same kernel bandwidth in both halves (taken from all cells)
  if (is.null(moran_bandwidth)) moran_bandwidth <- diff(range(coord[, 1])) / 20

  run_half <- function(idx) {
    imp <- list(rna = imputation$rna[idx, , drop = FALSE], atac = imputation$atac[idx, , drop = FALSE])
    res <- do.call(fit_fn, c(list(sp_interval = sp_interval, time_interval = time_interval,
                                  spatial_coord = coord[idx, , drop = FALSE], pseudoT = pseudoT[idx],
                                  use_peaks = use_peaks, use_genes = use_genes, imputation = imp,
                                  cell_types = if (is.null(cell_types)) NULL else cell_types[idx]),
                             fit_args))
    lo <- NULL
    if (use_link) {
      lo <- link_object
      lo$obj <- base::subset(lo$obj, cells = rownames(coord)[idx])
    }
    bl <- SAVI_baselines(imp, use_genes, use_peaks, coord[idx, , drop = FALSE], methods = methods,
                         link_object = lo, moran_bandwidth = moran_bandwidth, n_cores = n_cores)
    list(res = res, baselines = bl)
  }

  all_rep <- list(); splits <- list()
  for (s in seq_len(n_splits)) {
    set.seed(seed + s)
    m <- floor(n / 2)
    perm <- sample.int(n)
    h1 <- sort(perm[seq_len(m)]); h2 <- sort(perm[m + seq_len(m)])
    splits[[s]] <- list(half1 = h1, half2 = h2)
    message("split ", s, "/", n_splits, ": ", m, " cells per half")
    a <- run_half(h1); b <- run_half(h2)
    rep_s <- eval_reproducibility(a$res, b$res, a$baselines, b$baselines,
                                  intercept = intercept, n_cores = n_cores)
    if (use_link) {
      for (mm in intersect(c("signac", "cicero"), methods)) {
        lk <- eval_link_reproducibility(a$baselines[[mm]], b$baselines[[mm]], use_peaks, use_peaks,
                                        method = if (mm == "signac") "Signac" else "Cicero",
                                        type = if (mm == "signac") "pvalue" else "score")
        rep_s <- rbind(rep_s, lk)
      }
    }
    if (!is.null(rep_s)) { rep_s$split <- s; all_rep[[s]] <- rep_s }
  }
  rep_all <- do.call(rbind, all_rep)
  if (is.null(rep_all)) { warning("no gene could be evaluated"); return(NULL) }
  rownames(rep_all) <- NULL
  summ <- aggregate(cbind(corre, corre_pos) ~ method, data = rep_all, FUN = mean)
  list(reproducibility = rep_all, summary = summ, splits = splits,
       plot = if (requireNamespace("ggplot2", quietly = TRUE)) plot_reproducibility(rep_all) else NULL)
}
