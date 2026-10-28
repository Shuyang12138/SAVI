################################################################################
# SAVI_data_prep.R
#
# End-to-end data preparation for SAVI, in three stages:
#
#   (1) R       SAVI_prepare_gleam_inputs(sc_data, sp_data, ...)
#               two Seurat objects  ->  files for GLEAM in `work_dir`
#   (2) Python  python GLEAM_impute.py --work_dir <work_dir>
#               run separately (own conda environment, not called from R);
#               trains GLEAM and writes the imputation files to `work_dir`
#   (3) R       SAVI_load_imputation(work_dir)
#               GLEAM output  ->  raw + imputed data and the other SAVI inputs
#
# INPUT of stage (1)
#   sc_data : single-cell multiome Seurat object (e.g. e18_brain from
#             process_multiome(): RNA + ATAC ChromatinAssay with annotation,
#             TF-IDF already run on ATAC).
#   sp_data : spatial Seurat object (e.g. P22mouse: assays "Spatial" (RNA
#             counts) and "peaks", reduction "spatial", metadata "RNA_clusters").
#   fragment_path : ATAC fragment file (.tsv.gz) of the spatial data.
#
# OUTPUT of stage (3): a list
#   $raw      list(rna, atac)  observed data of the spatial cells
#   $imputed  list(rna, atac)  GLEAM-imputed data (cells x features)
#             --> pass as `imputation =` to fusing_bench() / fusing_bench_multiveloT()
#   $spatial_coord, $cell_types, $use_genes, $use_peaks (peaks per gene),
#   $sc_peaks, $sp_peaks, $sp_cells
#   All matrices are cells x features, rows in the order of colnames(sp_data).
#
# REQUIREMENTS
#   Seurat >= 5, Signac, GenomicRanges, Matrix, data.table. This file is
#   self-contained: it does not need SAVI.R (the two helper functions it uses,
#   find_peak() and CollapseToLongestTranscript(), are defined below).
################################################################################

library(Seurat)
library(Signac)
library(GenomicRanges)
library(Matrix)
library(data.table)

################################################################################
# Helper functions (same as in SAVI.R)
################################################################################

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

################################################################################
# Internal helpers
################################################################################

# kNN object -> data.frame(orig, neighbor, neg_dist) with neg_dist = 10 / distance
# (the format of sc_cca_knn.csv / sp_cca_knn.csv; same as merge_knn() before)
.savi_merge_knn <- function(knn_neighbor) {
  idx <- knn_neighbor@nn.idx
  dst <- knn_neighbor@nn.dist
  data.frame(orig = rep(seq_len(nrow(idx)), ncol(idx)),
             neighbor = as.vector(idx),
             neg_dist = 10 / as.vector(dst))
}

# Within-dataset RNA + ATAC weighted-nearest-neighbour graph (as multi_emb()
# before): SCT + PCA for RNA, TF-IDF + SVD (LSI) for ATAC, then WNN.
# Returns data.frame(orig, neighbor, neg_dist) with neg_dist = 1 - distance.
.savi_wnn_knn <- function(counts, atac_assay) {
  obj <- CreateSeuratObject(counts)
  obj[["ATAC"]] <- atac_assay
  DefaultAssay(obj) <- "RNA"
  obj <- SCTransform(obj, verbose = FALSE)
  obj <- RunPCA(obj, verbose = FALSE)
  DefaultAssay(obj) <- "ATAC"
  obj <- RunTFIDF(obj, verbose = FALSE)
  obj <- FindTopFeatures(obj, min.cutoff = "q0", verbose = FALSE)
  obj <- RunSVD(obj, verbose = FALSE)
  obj <- FindMultiModalNeighbors(obj, reduction.list = list("pca", "lsi"),
                                 dims.list = list(1:50, 2:50), verbose = FALSE)
  idx <- obj@neighbors$weighted.nn@nn.idx
  dst <- obj@neighbors$weighted.nn@nn.dist
  data.frame(orig = rep(seq_len(nrow(idx)), ncol(idx)),
             neighbor = as.vector(idx),
             neg_dist = 1 - as.vector(dst))
}

# Scale each feature by its root mean square, without centering, clipped at 10
# (ScaleData(do.center = FALSE)); non-finite values are set to 0.
.savi_scale <- function(m) {
  out <- as.matrix(ScaleData(m, do.center = FALSE, verbose = FALSE))
  dimnames(out) <- dimnames(m)
  out[!is.finite(out)] <- 0
  out
}

# Block of -1 (= "not measured") with named rows/cols, used for padding
.savi_pad <- function(row_names, col_names) {
  matrix(-1, nrow = length(row_names), ncol = length(col_names),
         dimnames = list(row_names, col_names))
}

.savi_write <- function(x, work_dir, file) {
  fwrite(as.data.frame(x), file.path(work_dir, file))
}

################################################################################
# STAGE 1: SAVI_prepare_gleam_inputs()
################################################################################
SAVI_prepare_gleam_inputs <- function(sc_data, sp_data, fragment_path, work_dir,
                                      sc_rna_assay = "RNA", sc_atac_assay = "ATAC",
                                      sp_rna_assay = "Spatial", sp_atac_assay = "peaks",
                                      sp_reduction = "spatial", cell_type_col = "RNA_clusters",
                                      n_sc_var_peaks = 60000, n_sc_var_genes = 1000,
                                      n_sp_var_peaks = 1000,
                                      window = 500000, min_peaks_per_gene = 50,
                                      additional_gene = c(), k = 20, seed = 2024) {
  # Arguments (besides the inputs above)
  #   *_assay / sp_reduction / cell_type_col : names of assays, spatial reduction
  #       and cell-label metadata column in the two objects
  #   n_sc_var_peaks, n_sc_var_genes, n_sp_var_peaks : variable features used
  #       for selection
  #   window : half-width (bp) of the window around each gene's start
  #   min_peaks_per_gene : a gene is kept if more than this many sc variable
  #       peaks fall in its window
  #   additional_gene : genes to add to the sc variable genes given to GLEAM
  #   k : number of sc <-> spatial neighbours in CCA space
  # Writes (to work_dir): the csv files read by GLEAM_impute.py,
  #   generated_data.Rdata (feature names needed to read GLEAM's output) and
  #   SAVI_selection.rds (selected genes / peaks).

  set.seed(seed)
  old_ver <- getOption("Seurat.object.assay.version")
  on.exit(options(Seurat.object.assay.version = old_ver), add = TRUE)
  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
  sep <- c("-", "-")
  sp_cells <- colnames(sp_data)
  sc_cells <- colnames(sc_data)

  ##############################################################################
  # STEP 1: select genes and peaks
  ##############################################################################
  message("Step 1: selecting genes and peaks")

  # 1a. Variable features: sc ATAC peaks, sc RNA genes, spatial ATAC peaks
  sc_data <- FindVariableFeatures(sc_data, assay = sc_atac_assay,
                                  nfeatures = n_sc_var_peaks, verbose = FALSE)
  sc_data <- FindVariableFeatures(sc_data, assay = sc_rna_assay,
                                  nfeatures = n_sc_var_genes, verbose = FALSE)
  sp_data <- FindVariableFeatures(sp_data, assay = sp_atac_assay,
                                  nfeatures = n_sp_var_peaks, verbose = FALSE)
  sc_var_genes <- VariableFeatures(sc_data, assay = sc_rna_assay)
  sc_var_peaks <- VariableFeatures(sc_data, assay = sc_atac_assay)
  sp_var_peaks <- VariableFeatures(sp_data, assay = sp_atac_assay)

  # 1b. One window per sc variable gene: [start - window, start + window]
  #     (windows are defined from the gene start, as in the original script)
  gene_coords <- as.data.frame(
    CollapseToLongestTranscript(ranges = Annotation(sc_data[[sc_atac_assay]]))
  )
  gene_win <- gene_coords[, c("seqnames", "start", "end", "gene_name")]
  gene_win$begin <- ifelse(gene_win$start < window, 0, gene_win$start - window)
  gene_win$end <- gene_win$start + window
  gene_win$win <- paste0(gene_win$seqnames, "-", gene_win$begin, "-", gene_win$end)
  gene_win <- gene_win[gene_win$gene_name %in% c(sc_var_genes, additional_gene), ]
  if (nrow(gene_win) == 0) stop("No sc variable genes found in the ATAC annotation.")

  # peak x gene matrix: which sc variable peaks lie in which gene window
  peak_gene <- find_peak(
    peaks = StringToGRanges(sc_var_peaks, sep = sep),
    genes = StringToGRanges(gene_win$win, sep = sep)
  )
  rownames(peak_gene) <- sc_var_peaks
  colnames(peak_gene) <- gene_win$gene_name

  # keep genes with enough nearby peaks; sc_peaks = sc variable peaks near them
  keep <- Matrix::colSums(peak_gene) > min_peaks_per_gene
  if (!any(keep)) stop("No gene has more than min_peaks_per_gene peaks in its window.")
  sc_peaks <- rownames(peak_gene)[Matrix::rowSums(peak_gene[, keep, drop = FALSE]) > 0]

  # spatial-specific peaks: spatial variable peaks that do not overlap any sc_peaks
  sp_vs_sc <- find_peak(
    peaks = StringToGRanges(sp_var_peaks, sep = sep),
    genes = StringToGRanges(sc_peaks, sep = sep)
  )
  sp_peaks <- sp_var_peaks[Matrix::rowSums(sp_vs_sc) == 0]
  all_peaks <- c(sc_peaks, sp_peaks)

  # 1c. Genes given to GLEAM: sc variable genes (+ additional_gene). Genes
  #     measured in the spatial data come first; the others are padded with -1
  #     in the spatial matrix and imputed by GLEAM.
  sc_counts <- GetAssayData(sc_data, assay = sc_rna_assay, layer = "counts")
  sp_counts <- GetAssayData(sp_data, assay = sp_rna_assay, layer = "counts")
  rna_genes <- intersect(unique(c(sc_var_genes, additional_gene)), rownames(sc_counts))
  sp_genes <- intersect(rna_genes, rownames(sp_counts))
  gleam_genes <- c(sp_genes, setdiff(rna_genes, sp_genes))
  sc_only_genes <- setdiff(rna_genes, sp_genes)

  # 1d. Genes analysed by SAVI and the candidate peaks of each gene
  use_genes <- intersect(colnames(peak_gene)[keep], sp_genes)
  if (length(use_genes) == 0) stop("No genes left after filtering.")
  cand <- find_peak(
    peaks = StringToGRanges(all_peaks, sep = sep),
    genes = StringToGRanges(gene_win$win[match(use_genes, gene_win$gene_name)], sep = sep)
  )
  use_peaks <- lapply(seq_along(use_genes), function(j) all_peaks[cand[, j] > 0])
  names(use_peaks) <- use_genes
  message("  genes: ", length(use_genes), " | sc peaks: ", length(sc_peaks),
          " | spatial-specific peaks: ", length(sp_peaks))

  ##############################################################################
  # STEP 2: raw data matrices for GLEAM (features x cells; -1 = not measured)
  ##############################################################################
  message("Step 2: building raw data matrices")

  # RNA: scale the (unnormalised) counts of the selected genes, as before
  sc_rna_raw <- .savi_scale(sc_counts[gleam_genes, , drop = FALSE])
  sp_rna_raw <- .savi_scale(sp_counts[sp_genes, , drop = FALSE])
  sp_rna_raw <- rbind(sp_rna_raw, .savi_pad(sc_only_genes, sp_cells))

  # ATAC (sc): TF-IDF values of the sc peaks; spatial-specific peaks are not
  # measured in sc and are padded with -1
  sc_atac <- .savi_scale(
    GetAssayData(sc_data, assay = sc_atac_assay, layer = "data")[sc_peaks, , drop = FALSE]
  )
  sc_atac <- rbind(sc_atac, .savi_pad(sp_peaks, sc_cells))

  # ATAC (spatial): count fragments in all selected peaks, then scale
  sp_frag <- CreateFragmentObject(path = fragment_path, cells = sp_cells, verbose = FALSE)
  sp_atac_counts <- FeatureMatrix(
    fragments = sp_frag,
    features = StringToGRanges(all_peaks, sep = sep),
    verbose = FALSE
  )
  if (!all(all_peaks %in% rownames(sp_atac_counts))) {
    stop("FeatureMatrix did not return all selected peaks; check peak name format.")
  }
  sp_atac <- .savi_scale(sp_atac_counts[all_peaks, sp_cells])

  .savi_write(sc_rna_raw, work_dir, "sc_rna_raw.csv")
  .savi_write(sp_rna_raw, work_dir, "sp_rna_raw.csv")
  .savi_write(sc_atac, work_dir, "sc_atac_scaled.csv")
  .savi_write(sp_atac, work_dir, "sp_atac_scaled.csv")

  ##############################################################################
  # STEP 3: graphs for GLEAM
  ##############################################################################
  message("Step 3: joint embedding and kNN graphs")

  # 3a. CCA integration of sc and spatial RNA. Cell names get a prefix so that
  #     identical barcodes in the two datasets cannot collide.
  options(Seurat.object.assay.version = "v5")
  ifnb <- merge(CreateSeuratObject(sc_counts), y = CreateSeuratObject(sp_counts),
                add.cell.ids = c("sc", "sp"))
  ifnb <- NormalizeData(ifnb, verbose = FALSE)
  ifnb <- FindVariableFeatures(ifnb, verbose = FALSE)
  ifnb <- ScaleData(ifnb, verbose = FALSE)
  ifnb <- RunPCA(ifnb, verbose = FALSE)
  ifnb$source <- ifelse(startsWith(colnames(ifnb), "sc_"), "sc", "sp")
  ifnb <- IntegrateLayers(object = ifnb, method = CCAIntegration,
                          orig.reduction = "pca", new.reduction = "integrated.cca",
                          verbose = FALSE)
  emb <- Embeddings(ifnb, reduction = "integrated.cca")
  sc_emb <- emb[paste0("sc_", sc_cells), , drop = FALSE]
  sp_emb <- emb[paste0("sp_", sp_cells), , drop = FALSE]
  .savi_write(sc_emb, work_dir, "sc_cca_emb.csv")
  .savi_write(sp_emb, work_dir, "sp_cca_emb.csv")

  # 3b. sc <-> spatial kNN in CCA space
  sp_cca_knn <- .savi_merge_knn(FindNeighbors(object = sc_emb, query = sp_emb,
                                              k.param = k, return.neighbor = TRUE, verbose = FALSE))
  sc_cca_knn <- .savi_merge_knn(FindNeighbors(object = sp_emb, query = sc_emb,
                                              k.param = k, return.neighbor = TRUE, verbose = FALSE))
  .savi_write(sc_cca_knn, work_dir, "sc_cca_knn.csv")
  .savi_write(sp_cca_knn, work_dir, "sp_cca_knn.csv")

  # 3c. within-dataset RNA + ATAC (WNN) kNN graphs; the original (not the
  #     selected-peak) ATAC assays are used here, as before
  options(Seurat.object.assay.version = "v3")
  sc_rna_atac_knn <- .savi_wnn_knn(sc_counts, sc_data[[sc_atac_assay]])
  sp_rna_atac_knn <- .savi_wnn_knn(sp_counts, sp_data[[sp_atac_assay]])
  .savi_write(sc_rna_atac_knn, work_dir, "sc_rna_atac_knn.csv")
  .savi_write(sp_rna_atac_knn, work_dir, "sp_rna_atac_knn.csv")

  # 3d. spatial kNN of the spatial cells (self-links, distance 0, are dropped)
  sp_coord <- Embeddings(sp_data, reduction = sp_reduction)[sp_cells, , drop = FALSE]
  sp_spatial_knn <- FindNeighbors(scale(sp_coord), return.neighbor = TRUE,
                                  k.param = 21, verbose = FALSE)
  sp_spatial_knn_result <- data.frame(
    orig = rep(seq_len(nrow(sp_spatial_knn@nn.idx)), ncol(sp_spatial_knn@nn.idx)),
    neighbor = as.vector(sp_spatial_knn@nn.idx),
    neg_dist = 1 - as.vector(sp_spatial_knn@nn.dist)
  )
  sp_spatial_knn_result <- sp_spatial_knn_result[sp_spatial_knn_result$neg_dist < 1, ]
  .savi_write(sp_spatial_knn_result, work_dir, "sp_spatial_knn.csv")

  # spatial coordinates and cell labels (used later by SAVI)
  .savi_write(sp_coord, work_dir, "sp_spatial_coord.csv")
  if (cell_type_col %in% colnames(sp_data@meta.data)) {
    .savi_write(data.frame(x = as.character(sp_data@meta.data[sp_cells, cell_type_col])),
                work_dir, "sp_RNA_clusters.csv")
  } else {
    warning("cell_type_col '", cell_type_col, "' not found in sp_data metadata; no cell labels written")
  }

  ##############################################################################
  # Feature names and selection, needed to read GLEAM's output
  ##############################################################################
  save(list = c("sc_rna_raw", "sp_rna_raw", "sc_atac", "sp_atac"),
       file = file.path(work_dir, "generated_data.Rdata"))
  saveRDS(list(use_genes = use_genes, use_peaks = use_peaks,
               sc_peaks = sc_peaks, sp_peaks = sp_peaks),
          file = file.path(work_dir, "SAVI_selection.rds"))
  message("Done. Next: run GLEAM in its own conda environment: python GLEAM_impute.py --work_dir ", work_dir)
  invisible(list(work_dir = work_dir, use_genes = use_genes, use_peaks = use_peaks,
                 sc_peaks = sc_peaks, sp_peaks = sp_peaks))
}

################################################################################
# STAGE 3: SAVI_load_imputation()
#   Reads the files written by GLEAM_impute.py and assembles the SAVI inputs.
#   (read_imputation() before, plus the other SAVI inputs.)
################################################################################
SAVI_load_imputation <- function(work_dir, restrict_genes = NULL) {
  # restrict_genes : optional vector of genes to keep (e.g. names(P22_multipseudoT),
  #                  genes that have a pseudotime estimate)

  # feature names / cell order used by GLEAM, and the selection
  e <- new.env()
  load(file.path(work_dir, "generated_data.Rdata"), envir = e)
  sel <- readRDS(file.path(work_dir, "SAVI_selection.rds"))
  gene <- rownames(e$sc_rna_raw)   # gene order in GLEAM's RNA output
  peak <- rownames(e$sp_atac)      # peak order in GLEAM's ATAC output
  sp_cells <- colnames(e$sp_rna_raw)

  # GLEAM csv: first column = node index, then one column per feature, one row
  # per spatial cell (in the order of sp_cells)
  read_gleam <- function(file, col_names) {
    m <- fread(file.path(work_dir, file), data.table = FALSE)
    m <- as.matrix(m[, -1, drop = FALSE])
    if (nrow(m) != length(sp_cells) || ncol(m) != length(col_names)) {
      stop(file, ": expected ", length(sp_cells), " x ", length(col_names),
           " but found ", nrow(m), " x ", ncol(m))
    }
    dimnames(m) <- list(sp_cells, col_names)
    m
  }
  rna_exist <- read_gleam("sp_rna_exist_withatac_all.csv", gene)
  rna_pred <- read_gleam("sp_rna_pred_withatac_all.csv", gene)
  atac_exist <- read_gleam("sp_atac_exist_withatac.csv", peak)

  # imputed RNA = existence probability x predicted level (kept where > 0);
  # imputed ATAC = predicted accessibility
  rna_aggregate <- rna_exist * rna_pred * (rna_pred > 0)

  # genes / peaks analysed by SAVI
  use_genes <- intersect(sel$use_genes, gene)
  if (!is.null(restrict_genes)) use_genes <- intersect(use_genes, restrict_genes)
  if (length(use_genes) == 0) stop("No genes left after filtering.")
  use_peaks <- sel$use_peaks[use_genes]
  needed_peaks <- peak[peak %in% unique(unlist(use_peaks))]

  # spatial coordinates and cell labels
  spatial_coord <- as.matrix(fread(file.path(work_dir, "sp_spatial_coord.csv"), data.table = FALSE))
  rownames(spatial_coord) <- sp_cells
  ct_file <- file.path(work_dir, "sp_RNA_clusters.csv")
  cell_types <- if (file.exists(ct_file)) fread(ct_file, data.table = FALSE)$x else NULL

  list(
    raw = list(rna = t(e$sp_rna_raw[use_genes, , drop = FALSE]),
               atac = t(e$sp_atac[needed_peaks, , drop = FALSE])),
    imputed = list(rna = rna_aggregate[, use_genes, drop = FALSE],
                   atac = atac_exist[, needed_peaks, drop = FALSE]),
    spatial_coord = spatial_coord, cell_types = cell_types,
    use_genes = use_genes, use_peaks = use_peaks,
    sc_peaks = sel$sc_peaks, sp_peaks = sel$sp_peaks, sp_cells = sp_cells
  )
}

################################################################################
# EXAMPLE (not run); see SAVI_data_prep_example.ipynb for a complete run
################################################################################
# # 1. R: Seurat objects -> GLEAM inputs
# SAVI_prepare_gleam_inputs(
#   sc_data = e18_brain, sp_data = P22mouse,
#   fragment_path = "/storage10/shuyang/GLEAM_real/P22_GSM6758285_MouseBrain_20um_100barcodes_ATAC_fragments.tsv.gz",
#   work_dir = "/storage10/shuyang/GLEAM_real/P22_E18"
# )
#
# # 2. shell: train GLEAM and write the imputation
# #    python GLEAM_impute.py --work_dir /storage10/shuyang/GLEAM_real/P22_E18 --gpu 1
#
# # 3. R: load raw + imputed data and the other SAVI inputs, then save them
# prep <- SAVI_load_imputation("/storage10/shuyang/GLEAM_real/P22_E18")
# saveRDS(prep, "/storage10/shuyang/GLEAM_real/P22_E18/SAVI_input.rds")
#
# # SAVI.R then reads SAVI_input.rds (see SAVI_example.ipynb)
