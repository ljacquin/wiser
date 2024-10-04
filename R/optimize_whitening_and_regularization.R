# function which finds the optimal whitening method and regularization parameter
optimize_whitening_and_regularization <- function(
    omic_df, raw_pheno_df, trait_,
    fixed_effects_vars = c(
      "Envir", "Country", "Year",
      "Row", "Position", "Management"
    ),
    fixed_effects_vars_computed_as_factor = c(
      "Envir", "Country", "Year",
      "Row", "Position", "Management"
    ),
    site_var = "Country",
    fixed_effects_vars_computed_as_factor_by_site = c("Row", "Position"),
    random_effects_vars = "Genotype",
    prediction_method = c("rf", "svr", "gblup", "rkhs", "lasso"),
    whitening_method_grid = c("ZCA-cor", "PCA-cor", "Cholesky"),
    regularization_method_ = "frobenius_norm_regularization",
    alpha_grid = c(0.01, 0.1),
    reduce_raw_dataset_size_ = T,
    nrow_approx_lim_raw_dataset_ = 5e3,
    parallelized_cholesky = T,
    k_folds_ = 5,
    nb_cores_ = 12) {
  # remove all rows with na w.r.t to trait_
  raw_pheno_df <- raw_pheno_df %>% drop_na(all_of(trait_))

  # get raw phenotypes and marker data based on common genotypes
  raw_pheno_df <- match_indices(raw_pheno_df, omic_df)
  omic_df <- omic_df[rownames(omic_df) %in% raw_pheno_df$Genotype, ]

  # define variables of interest
  sel_vars_ <- c(fixed_effects_vars, random_effects_vars, trait_)

  # get only variables of interest from raw_pheno_df
  raw_pheno_df <- raw_pheno_df[, sel_vars_]
  raw_pheno_df <- na.omit(raw_pheno_df)

  # downsize dataset for computation time optimization purpose
  if (reduce_raw_dataset_size_ && 
      (nrow(raw_pheno_df) > nrow_approx_lim_raw_dataset_)) {
    set.seed(123)
    raw_pheno_df <- as.data.frame(
      reduce_dataset_based_on_genotypes(
        df_ = raw_pheno_df,
        nrow_approx_lim = nrow_approx_lim_raw_dataset_
      )
    )
  }

  # create a grid combining whitening methods, regularization parameter and
  # prediction methods
  grid_ <- expand.grid(
    whitening_method = whitening_method_grid,
    alpha_ = alpha_grid,
    pred_method = prediction_method
  )

  # pre-compute unique wiser object for unique combinations
  wiser_cache <- list()
  unique_combinations <- unique(grid_[, c("whitening_method", "alpha_")])

  for (j in 1:nrow(unique_combinations)) {
    method <- unique_combinations$whitening_method[j]
    alpha <- unique_combinations$alpha_[j]

    wiser_obj <- estimate_wiser_phenotype(
      omic_df, raw_pheno_df, trait_,
      fixed_effects_vars,
      fixed_effects_vars_computed_as_factor,
      site_var,
      fixed_effects_vars_computed_as_factor_by_site,
      random_effects_vars,
      whitening_method = method,
      regularization_method = regularization_method_,
      alpha_ = alpha,
      reduce_raw_dataset_size_ = FALSE
    )

    cache_key <- paste(method, alpha, sep = "_")
    wiser_cache[[cache_key]] <- wiser_obj
  }

  # configure parallelization
  plan(multisession, workers = nb_cores_)

  df_results <- future_lapply(
    1:nrow(grid_),
    future.seed = T,
    function(i) {
      # retrieve the precomputed wiser_obj for this combination
      cache_key <- paste(grid_$whitening_method[i], grid_$alpha_[i], sep = "_")
      wiser_obj_local <- wiser_cache[[cache_key]]

      mean_pa <- tryCatch(
        {
          perform_kfold_cv_wiser(
            omic_df, raw_pheno_df, trait_,
            fixed_effects_vars,
            fixed_effects_vars_computed_as_factor,
            site_var,
            fixed_effects_vars_computed_as_factor_by_site,
            random_effects_vars,
            whitening_method = grid_$whitening_method[i],
            reg_method = regularization_method_,
            alpha_ = grid_$alpha_[i],
            pred_method = grid_$pred_method[i],
            k_folds = k_folds_,
            wiser_obj_local = wiser_obj_local
          )
        },
        error = function(e) {
          cat("Error during iteration", i, ":", e$message, "\n")
          return(NA)
        }
      )
      data.frame(
        "whitening_method" = grid_$whitening_method[i],
        "alpha_" = grid_$alpha_[i],
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
  for (method_ in df_results$prediction_method) {
    df_res_method_ <- df_results[
      df_results$prediction_method == method_,
    ]
    df_res_method_$white_reg_combination <- paste0(
      df_res_method_$whitening_method,
      "/", df_res_method_$alpha_
    )
    df_opt_ <- rbind(
      df_opt_,
      unique(df_res_method_[
        which.max(df_res_method_$mean_pa)[1],
      ])
    )
  }
  opt_mode_ <- compute_vect_mode(df_opt_$white_reg_combination)
  opt_mode_ <- unlist(str_split(opt_mode_, pattern = "/"))
  opt_whitening_method <- opt_mode_[1]
  opt_alpha_ <- as.numeric(opt_mode_[2])

  # stop parallelization
  plan(sequential)

  return(list(
    "opt_results" = unique(df_opt_),
    "opt_alpha_" = opt_alpha_,
    "opt_whitening_method" = opt_whitening_method
  ))
}
