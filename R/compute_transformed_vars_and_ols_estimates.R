# function which computes transformed fixed variables and least squares
compute_transformed_vars_and_ols_estimates <- function(
    omic_df, raw_pheno_df, trait_,
    fixed_effect_vars,
    fixed_effect_vars_computed_as_factor,
    envir_var,
    fixed_effect_vars_computed_as_factor_by_envir,
    random_effect_vars,
    sigma2_u, sigma2_e, kernel_type,
    whitening_method,
    regularization_method,
    alpha_,
    parallelized_cholesky,
    reduce_raw_dataset_size_,
    nrow_approx_lim_raw_dataset_zca_cor,
    nrow_approx_lim_raw_dataset_pca_cor,
    nrow_approx_lim_raw_dataset_chol) {
  tryCatch(
    {
      # get raw phenotypes and omic data based on common genotypes
      raw_data_obj <-
        raw_pheno_and_marker_based_on_trait_common_genotypes(
          raw_pheno_df,
          omic_df,
          trait_,
          fixed_effect_vars,
          random_effect_vars
        )
      raw_pheno_df <- raw_data_obj$raw_pheno_df
      omic_df <- raw_data_obj$omic_df

      # should raw dataset size be reduced wrt to selected whitening method ?
      if (reduce_raw_dataset_size_ && (
        (whitening_method == "ZCA-cor" &&
          nrow(raw_pheno_df) > nrow_approx_lim_raw_dataset_zca_cor) ||
          (whitening_method == "PCA-cor" &&
            nrow(raw_pheno_df) > nrow_approx_lim_raw_dataset_pca_cor) ||
          (whitening_method == "Cholesky" &&
            nrow(raw_pheno_df) > nrow_approx_lim_raw_dataset_chol)
      )
      ) {
        raw_pheno_df <- reduce_dataset_based_on_selected_whitening(
          whitening_method,
          raw_pheno_df,
          nrow_approx_lim_raw_dataset_zca_cor,
          nrow_approx_lim_raw_dataset_pca_cor,
          nrow_approx_lim_raw_dataset_chol
        )
      }

      # computes fixed effect vars as factors for those declared as
      raw_pheno_df <- compute_fixed_effect_vars_declared_as_factors(
        raw_pheno_df,
        fixed_effect_vars,
        fixed_effect_vars_computed_as_factor,
        envir_var,
        fixed_effect_vars_computed_as_factor_by_envir
      )

      # get omic data associated to common genotypes
      omic_df <- omic_df[rownames(omic_df) %in% unique(raw_pheno_df$Genotype), ]

      # compute Gram matrix (e.g. genomic covariance matrix)
      k_mat <- compute_gram_matrix(omic_df, kernel_type)

      # remove fixed effects with no variance or unique level for factors
      if (!is.null(ncol(raw_pheno_df[, fixed_effect_vars])) &&
        ncol(raw_pheno_df[, fixed_effect_vars]) > 1) {
        fixed_effect_vars <- find_columns_with_multiple_unique_values(
          raw_pheno_df[, fixed_effect_vars]
        )
      }

      # get incidence matrices for fixed and random effects
      # NB. column of ones is added for intercept associated to fixed effects
      incid_obj <- compute_incidence_matrices_fixed_and_random_effects(
        fixed_effect_vars,
        fixed_effect_vars_computed_as_factor,
        random_effect_vars,
        raw_pheno_df
      )
      x_mat <- incid_obj$x_mat
      z_mat <- incid_obj$z_mat

      # compute Σu, i.e. sig_mat_ here
      sig_mat_ <- sigma2_u * crossprod(t(z_mat), tcrossprod(k_mat, z_mat))

      # compute the whitening matrix for Σu based on the selected
      # whitening method
      white_obj <- compute_whitening_matrix_for_sig_mat_(
        whitening_method,
        regularization_method,
        parallelized_cholesky,
        sig_mat_, alpha_
      )
      w_mat <- white_obj$w_mat
      sig_mat_ <- white_obj$sig_mat_

      # whiten x_mat using w_mat
      # NB. intercept is already present in x_mat and x_mat_tilde
      x_mat_tilde <- w_mat %*% x_mat

      # get raw phenotypes associated to common genotypes
      y <- as.numeric(raw_pheno_df[, trait_])

      # get ols estimates for fixed effects and xi
      beta_hat <- ginv(t(x_mat_tilde) %*% x_mat_tilde) %*% t(x_mat_tilde) %*% y
      y_hat <- x_mat_tilde %*% beta_hat
      xi_hat <- y - y_hat
      
      # add the individual estimated phenotype with fixed effects eliminated
      # and corrected for the genetic covariance structure.
      raw_pheno_df$xi_hat <- xi_hat

      return(list(
        "omic_df" = omic_df,
        "sig_mat_u" = sig_mat_,
        "w_mat" = w_mat,
        "x_mat" = x_mat,
        "x_mat_tilde" = x_mat_tilde,
        "z_mat" = z_mat,
        "k_mat" = k_mat,
        "beta_hat" = beta_hat,
        "y_hat" = y_hat,
        "xi_hat" = xi_hat,
        "y" = y,
        "xi_phenotypes" = raw_pheno_df
      ))
    },
    error = function(e) {
      cat(
        "Error with : ", conditionMessage(e), "\n"
      )
    }
  )
}
