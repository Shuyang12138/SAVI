library(Seurat)
library(Signac)
library(ggplot2)
library(dplyr)
library(data.table)
library(GenomicRanges)
library(Matrix)
find_peak <- function(peaks, genes, distance = 0, sep = c("-", "-")) {
  overlaps <- findOverlaps(
    query = peaks,
    subject = genes,
    type = 'any',
    select = 'all'
  )
  hit_matrix <- sparseMatrix(
    i = queryHits(x = overlaps),
    j = subjectHits(x = overlaps),
    x = 1,
    dims = c(length(x = peaks), length(x = genes))
  )
  rownames(x = hit_matrix) <- GRangesToString(grange = peaks, sep = sep)
  colnames(x = hit_matrix) <- GRangesToString(grange = genes, sep = sep)
  return(hit_matrix)
}
CollapseToLongestTranscript <- function(ranges) {
  seqnames=start=end=strand=gene_biotype=gene_name=NULL
  range.df <- as.data.table(x = ranges)
  range.df$strand <- as.character(x = range.df$strand)
  range.df$strand <- ifelse(
    test = range.df$strand == "*",
    yes = "+",
    no = range.df$strand
  )
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
  collapsed$gene_name <- make.unique(names = collapsed$gene_name)
  gene.ranges <- makeGRangesFromDataFrame(
    df = collapsed,
    keep.extra.columns = TRUE
  )
  return(gene.ranges)
}
library(glmnet)
gaussian_kernel <- function(d, bandwidth) {
          exp(- (d^2) / (2 * bandwidth^2))
        }
gw_correlation <- function(i, X, Y, W) {
          w <- W[i, ]  # Extract weights for location i
          w <- w / sum(w)  # Normalize weights
          
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
signac_p = function(background_corr,corre){
    background_array <- simplify2array(background_corr)

# Compute mean and standard deviation across the third dimension (i.e., over the list)
    mean_mat <- apply(background_array, c(1, 2), mean)
    sd_mat <- apply(background_array, c(1, 2), sd)

# Compute Z-score
    z_score <- (corre - mean_mat) / sd_mat
    pval_mat <- pnorm(z_score, lower.tail = FALSE)
    return(pval_mat)
}
fusing_bench_multiveloT = function(sp_num,time_num,spatial_coord,pseudoT,use_peaks,use_genes,imputation,cell_types,intercept=F){# need to make group_all = group_x much ealier!!!!!!!!!!!!!!!!
set.seed(2025, kind = "L'Ecuyer-CMRG")
tryCatch({
    result = future_map(1:length(use_genes),function(i){
grid = data.frame(spatial_coord)
grid$t = pseudoT[[i]]$x
        sp_interval <- (max(grid[,1]) - min(grid[,1])) / sp_num
      time_interval <- (max(grid$t) - min(grid$t)) / time_num
cell_x_simu = floor((grid[,1]-min(grid[,1]))/sp_interval)
cell_y_simu = floor((grid[,2]-min(grid[,2]))/sp_interval)
cell_t_simu = floor((grid$t-min(grid$t))/time_interval)
library("dplyr")
group_count = data.frame(cell_x_simu,cell_y_simu,cell_t_simu)%>%count(cell_x_simu,cell_y_simu,cell_t_simu)
group_count = group_count[group_count$n>5,]
group_count$group = paste0(group_count$cell_x_simu,'-',group_count$cell_y_simu,'-',group_count$cell_t_simu)
group_count$region = paste0(group_count$cell_x_simu,'-',group_count$cell_y_simu)
rownames(group_count) = 1:nrow(group_count)
group_dist = as.matrix(dist(scale(group_count[,1:3],center = F)))
level = 1:nrow(group_count)
names(level) = group_count$group
group_x = data.frame(cell_x_simu,cell_y_simu,cell_t_simu)
group_x$group = paste0(group_x$cell_x_simu,'-',group_x$cell_y_simu,'-',group_x$cell_t_simu)
group_x$region = paste0(group_x$cell_x_simu,'-',group_x$cell_y_simu)
G_simu = exp(-group_dist^2)
G_simu[G_simu<exp(-1)] = 0
diag(G_simu) = 0


        # Your existing function logic
        y_simu = imputation$rna[, use_genes[i]]
        X_simu = imputation$atac[, use_peaks[[i]]]
        rownames(X_simu) = 1:nrow(X_simu)

        gene_exp = y_simu
        X_simu = X_simu[is.element(group_x$group, group_count$group), ]
        X_simu = as.matrix(X_simu)
        y_simu = y_simu[is.element(group_x$group, group_count$group)]
        y_simu = y_simu+rnorm(n = length(y_simu),mean=0,sd = 1e-3)
    cell_types =  cell_types[is.element(group_x$group, group_count$group)]
        group_x = group_x[is.element(group_x$group, group_count$group), ]

        groups = level[group_x$group]
        names(groups) = NULL
    #######################################################
        beta.estimate_simu = myFuser::fusedLassoProximal(
            X_simu, matrix(y_simu, ncol=1), groups, lambda=0.001, 
            tol=1e-6, gamma=0.001, G_simu, intercept=intercept, scaling=T, num.it=5000
        )
        lambda_seq <- c(seq(0, 0.001, 0.0002),0.002,0.005,0.008,0.01)
        gamma_seq <- c(seq(0.0001, 0.001, 0.0002),0.002,0.005,0.008,0.01)
        param_grid <- expand.grid(lam = lambda_seq, gam = gamma_seq)
        set.seed(123)
        group_ids <- unique(groups)
        group_folds <- lapply(group_ids, function(gid) {
              idx <- which(groups == gid)
              fold_assignments <- sample(rep(1:5, length.out = length(idx)))
              split(idx, fold_assignments)
        })
        names(group_folds) <- group_ids

        process_combination_cv_grouped <- function(params) {
            tryCatch({
              lam <- params$lam
              gam <- params$gam
              all_fold_errors <- numeric(5)

        for (fold in 1:5) {
            train_idx <- c()
            test_idx <- c()

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
    col_names <- colnames(X_test)
    y_pred <- diag(X_test[, col_names, drop = FALSE] %*% beta_estimate[col_names, groups_test])
    if(intercept) y_pred = y_pred+beta_estimate[nrow(beta_estimate),groups_test]

    # Compute squared error
    fold_error <- mean((y_test - y_pred)^2)
    all_fold_errors[fold] <- fold_error
  }

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
param_grid <- expand.grid(lam = lambda_seq, gam = gamma_seq)


param_list <- split(param_grid, seq(nrow(param_grid)))

# Run in parallel using future_map
cv_results <- lapply(
  param_list,
  function(params) suppressWarnings(process_combination_cv_grouped(params))
)

cv_df <- do.call(rbind, lapply(cv_results, as.data.frame))
best_params <- cv_df[which.min(cv_df$cv_rmse), ]
best_lambda <- best_params$lambda
best_gamma <- best_params$gamma
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
        all_group_result = c()
        for(g in colnames(beta.estimate_simu)){
                x_group = X_simu[groups==g,]
                y_group = y_simu[groups==g]
                if(length(y_group)>3){
                cv_model = cv.glmnet(x_group, y_group, alpha = 0.5,nfolds = 5)
                best_lambda <- cv_model$lambda.min
                    }
            else{
                best_lambda = 0.002
            }
                lasso_withingroup = glmnet(x_group, y_group, alpha = 0.5,lambda = best_lambda,intercept = intercept)
                all_group_result = cbind(all_group_result,coef(lasso_withingroup))
        }

        
        
        Y = y_simu[1:length(y_simu)]
        #all_peak_cor = mclapply(1:ncol(X_simu),function(j){
        #local_correlations <- sapply(1:length(Y), FUN=function(k) gw_correlation(k, X_simu[,j], Y, W))
        
        #return(local_correlations)
#},mc.cores=2)
        #all_peak_cor = do.call('rbind',all_peak_cor)
        #group_means <- sapply(colnames(beta.estimate_simu), function(g) rowMeans(all_peak_cor[, groups == g, drop = FALSE],na.rm = T))
        group_means <- sapply(colnames(beta.estimate_simu), function(g) as.vector(cor(X_simu[groups==g,],Y[groups==g],method = 'spearman')))
        rownames(group_means) = colnames(X_simu)

        null_means = lapply(1:1000,function(loop){
            for(i in unique(groups)){
    Y[groups==i] = sample(Y[groups==i],size = sum(groups==i))

    }
            #all_peak_cor = lapply(1:ncol(X_simu),function(j){
        #local_correlations <- sapply(1:length(Y), FUN=function(k) gw_correlation(k, X_simu[,j], Y, W))
        
        #return(local_correlations)
#})
        #all_peak_cor = do.call('rbind',all_peak_cor)
        #group_means <- sapply(colnames(beta.estimate_simu), function(g) rowMeans(all_peak_cor[, groups == g, drop = FALSE],na.rm = T))
        group_means <- sapply(colnames(beta.estimate_simu), function(g) as.vector(cor(X_simu[groups==g,],Y[groups==g],method = 'spearman')))
        rownames(group_means) = colnames(X_simu)
        return(group_means)
        })
        gwc_p = signac_p(null_means,group_means)
        group_all = group_x
        group_x = group_x[is.element(group_x$group, group_count$group), ]

        groups = level[group_x$group]
        names(groups) = NULL
        # cell type-specific GWC and spearman correlation 
            
        return(list(fusion_old=beta.estimate_simu,
                    lasso=all_group_result,
                    marginal = group_means,
                    marginal_p = gwc_p,
                    fusion = beta_estimate_final,
                    fusion_cv = list(lambda = best_lambda,  
                                     gamma = best_gamma,  
                                     cv_df = cv_df
                                    ),
                    groups = groups,
                    group_all = group_all,
                    group_count = group_count
                   )
              )
    }, .options = furrr_options(seed = TRUE))
                              
},error = function(e) {
        cat(sprintf(e$message), file='fuser_multivelo_error.log')
        return(NULL)  # Return NULL if an error occurs
    })
names(result) = use_genes
    
names(result) = use_genes

return(list(coef=result))
#return(list(initial_setting = list(coef=result,groups=groups,group_all=group_all,group_count=group_count)))
}
folder_path <- "code/MultiVeloT/"

# List all files ending with _t.csv
file_list <- list.files(path = folder_path, pattern = "_t\\.csv$", full.names = TRUE)

# Read all files into a list
P22_multipseudoT <- lapply(file_list, read.csv)
names(P22_multipseudoT) <- tools::file_path_sans_ext(basename(file_list))
names(P22_multipseudoT) <- sub("_t$", "", names(P22_multipseudoT))