perform_kfold_cv_wiser <- function(geno_df, raw_pheno_df, trait_,
                                   whitening_method,
                                   reg_method, alpha_frob,
                                   pred_method, k_folds) {
  # set seed for reproducibility and get set of indices
  set.seed(123)
  idx <- 1:nrow(geno_df)

  # create folds for k-folds cv
  folds <- cvFolds(nrow(geno_df), K = k_folds, type = "consecutive")

  # compute wiser phenotypes
  pheno_obj <- estimate_wiser_phenotype(
    geno_df, raw_pheno_df, trait_,
    whitening_method = whitening_method,
    regularization_method = reg_method,
    alpha_frob_ = alpha_frob,
    reduce_raw_dataset_size_ = F
  )
  v_hat <- pheno_obj$v_hat

  # use future_lapply for folds
  results <- future_lapply(1:k_folds,
    future.seed = T,
    function(fold) {
      idx_train <- idx[folds$which != fold]
      idx_val <- idx[folds$which == fold]

      # train and predict with random forest (using ranger package)
      if (pred_method == "rf") {
        rf_model <- ranger(
          y = v_hat[idx_train],
          x = geno_df[idx_train, ],
          mtry = ncol(geno_df) / 3,
          num.trees = 1000
        )
        f_hat_val_rf <- predict(rf_model, geno_df[idx_val, ])
        pa_ <- cor(
          f_hat_val_rf$predictions,
          v_hat[idx_val]
        )
        # train and predict with non linear svr (using kernlab package)
      } else if (pred_method == "svr") {
        c_par <- max(
          abs(mean(v_hat[idx_train])
          + 3 * sd(v_hat[idx_train])),
          abs(mean(v_hat[idx_train])
          - 3 * sd(v_hat[idx_train]))
        )
        gaussian_svr_model <- ksvm(
          x = as.matrix(geno_df[idx_train, ]),
          y = v_hat[idx_train],
          scaled = FALSE, type = "eps-svr",
          kernel = "rbfdot",
          kpar = "automatic", C = c_par, epsilon = 0.1
        )
        f_hat_val_gaussian_svr <- predict(
          gaussian_svr_model,
          as.matrix(geno_df[idx_val, ])
        )
        pa_ <- cor(
          f_hat_val_gaussian_svr,
          v_hat[idx_val]
        )
        # train and predict with gblup (using KRMM package)
      } else if (pred_method == "gblup") {
        linear_krmm_model <- krmm(
          Y = v_hat[idx_train],
          Matrix_covariates = geno_df[idx_train, ],
          method = "GBLUP"
        )
        f_hat_val_linear_krmm <- predict_krmm(linear_krmm_model,
          Matrix_covariates = geno_df[idx_val, ],
          add_fixed_effects = T
        )
        pa_ <- cor(
          f_hat_val_linear_krmm,
          v_hat[idx_val]
        )
      } else if (pred_method == "rkhs") {
        gaussian_krmm_model <- krmm(
          Y = v_hat[idx_train],
          Matrix_covariates = geno_df[idx_train, ],
          method = "RKHS", kernel = "Gaussian",
          rate_decay_kernel = 0.1
        )
        f_hat_val_gaussian_krmm <- predict_krmm(gaussian_krmm_model,
          Matrix_covariates = geno_df[idx_val, ],
          add_fixed_effects = T
        )
        pa_ <- cor(
          f_hat_val_gaussian_krmm,
          v_hat[idx_val]
        )
      } else {
        cv_fit_lasso_model <- cv.glmnet(
          intercept = TRUE, y = v_hat[idx_train],
          x = as.matrix(geno_df[idx_train, ]),
          type.measure = "mse", alpha = 1.0, nfold = 10,
          parallel = TRUE
        )
        f_hat_val_lasso <- predict(cv_fit_lasso_model,
          newx = as.matrix(geno_df[idx_val, ]),
          s = "lambda.min"
        )
        pa_ <- suppressWarnings(cor(
          f_hat_val_lasso,
          v_hat[idx_val]
        ))
      }

      data.frame(pa = pa_)
    }, future.packages = c(
      "ranger", "KRMM", "kernlab",
      "glmnet", "cvTools", "dplyr",
      "stringr", "matrixcalc",
      "Matrix", "whitening", "mixOmics"
    )
  )

  df_results <- do.call(rbind, results)
  mean_pa <- mean(df_results$pa)
  return(mean_pa)
}
