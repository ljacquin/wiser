optimize_whitening_and_regularization <- function(
    geno_df, raw_pheno_df, trait_,
    fixed_effects_vars = c(
      "Envir", "Country", "Year",
      "Row", "Position", "Management"
    ),
    random_effects_vars = "Genotype",
    prediction_method = c("rf", "svr", "gblup", "rkhs", "lasso"),
    whitening_method_grid = c("ZCA-cor", "PCA-cor", "Cholesky"),
    regularization_method_ = "frobenius_norm",
    alpha_frob_grid = c(0.01, 0.1),
    reduce_raw_dataset_size_ = T,
    nrow_lim_raw_dataset_ = 5e3,
    parallelized_cholesky_ = T,
    k_folds_ = 5) {
  # remove all rows with na w.r.t to trait_
  raw_pheno_df <- raw_pheno_df %>% drop_na(all_of(trait_))

  # define variables of interest
  sel_vars_ <- c(fixed_effects_vars, random_effects_vars, trait_)

  # get only variables of interest from raw_pheno_df
  raw_pheno_df <- raw_pheno_df[, sel_vars_]
  raw_pheno_df <- na.omit(raw_pheno_df)

  # downsize dataset for computation time optimization purpose
  if (reduce_raw_dataset_size_) {
    set.seed(123)
    raw_pheno_df <- as.data.frame(
      reduce_dataset_based_on_genotypes(
        df_ = raw_pheno_df,
        nrow_lim = nrow_lim_raw_dataset_
      )
    )
  }

  # create a grid combining whitening methods, regularization parameter and
  # prediction methods
  grid_ <- expand.grid(
    whitening_method = whitening_method_grid,
    alpha_frob = alpha_frob_grid,
    pred_method = prediction_method
  )

  # configure parallelization
  plan(multisession, workers = parallel::detectCores())

  df_results <- future_lapply(1:nrow(grid_),
    future.seed = F,
    function(i) {
      mean_pa <- perform_kfold_cv_wiser(
        geno_df, raw_pheno_df, trait_,
        whitening_method = grid_$whitening_method[i],
        reg_method = regularization_method_,
        alpha_frob = grid_$alpha_frob[i],
        pred_method = grid_$pred_method[i],
        k_folds = k_folds_
      )
      data.frame(
        "whitening_method" = grid_$whitening_method[i],
        "alpha_frob" = grid_$alpha_frob[i],
        "prediction_method" = grid_$pred_method[i],
        "mean_pa" = mean_pa
      )
    }, future.packages = c(
      "ranger", "KRMM", "kernlab",
      "glmnet", "cvTools", "dplyr",
      "stringr", "matrixcalc",
      "Matrix", "whitening", "mixOmics"
    )
  )
  df_results <- na.omit(do.call(rbind, df_results))

  # get optimal whitening method based on mean pa for each prediction method
  df_opt_ <- data.frame()
  for (method_ in prediction_method) {
    df_res_method_ <- df_results[
      df_results$prediction_method == method_,
    ]
    df_res_method_$white_reg_combination <- paste0(
      df_res_method_$whitening_method,
      "/", df_res_method_$alpha_frob
    )
    df_opt_ <- rbind(
      df_opt_,
      df_res_method_[
        which.max(df_res_method_$mean_pa),
      ]
    )
  }
  opt_mode_ <- compute_vect_mode(df_opt_$white_reg_combination)
  opt_mode_ <- unlist(str_split(opt_mode_, pattern = "/"))
  opt_whitening_method <- opt_mode_[1]
  opt_alpha_frob <- as.numeric(opt_mode_[2])

  # stop parallelization
  plan(sequential)

  return(list(
    "opt_results" = df_opt_,
    "opt_alpha_frob" = opt_alpha_frob,
    "opt_whitening_method" = opt_whitening_method
  ))
}
