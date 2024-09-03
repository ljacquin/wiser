compute_transformed_vars_and_ols_estimates <- function(
    omic_df, raw_pheno_df, trait_,
    fixed_effects_vars,
    fixed_effects_vars_computed_as_factor,
    site_var,
    fixed_effects_vars_computed_as_factor_by_site,
    random_effects_vars,
    sigma2_u, sigma2_e, kernel_type,
    whitening_method,
    regularization_method,
    alpha_frob_,
    percent_eig_,
    non_zero_precision_eig_,
    parallelized_cholesky,
    reduce_raw_dataset_size_,
    nrow_lim_raw_dataset_zca_cor,
    nrow_lim_raw_dataset_pca_cor,
    nrow_lim_raw_dataset_chol) {
  tryCatch(
    {
      # get raw phenotypes and omic data based on common genotypes
      raw_data_obj <-
        raw_pheno_and_marker_based_on_trait_common_genotypes(
          raw_pheno_df,
          omic_df,
          trait_,
          fixed_effects_vars,
          random_effects_vars
        )
      raw_pheno_df <- raw_data_obj$raw_pheno_df
      omic_df <- raw_data_obj$omic_df

      # should raw dataset size be reduced wrt to selected whitening method ?
      if (reduce_raw_dataset_size_) {
        raw_pheno_df <- reduce_dataset_based_on_selected_whitening(
          whitening_method,
          raw_pheno_df,
          nrow_lim_raw_dataset_zca_cor,
          nrow_lim_raw_dataset_pca_cor,
          nrow_lim_raw_dataset_chol
        )
      }

      # computes fixed effect vars as factors for those declared as
      raw_pheno_df <- compute_fixed_effect_vars_declared_as_factors(
        raw_pheno_df,
        fixed_effects_vars,
        fixed_effects_vars_computed_as_factor,
        site_var,
        fixed_effects_vars_computed_as_factor_by_site
      )

      # get omic data associated to common genotypes
      omic_df <- omic_df[rownames(omic_df) %in% unique(raw_pheno_df$Genotype), ]

      # compute Gram matrix (e.g. genomic covariance matrix)
      k_mat <- compute_gram_matrix(omic_df, kernel_type)

      # remove fixed effects with no variance or unique level for factors
      if (!is.null(ncol(raw_pheno_df[, fixed_effects_vars])) &&
        ncol(raw_pheno_df[, fixed_effects_vars]) > 1) {
        fixed_effects_vars <- find_columns_with_multiple_unique_values(
          raw_pheno_df[, fixed_effects_vars]
        )
      }

      # get incidence matrices for fixed and random effects
      # NB. column of ones is added for intercept associated to fixed effects
      incid_obj <- compute_incidence_matrices_fixed_and_random_effects(
        fixed_effects_vars,
        fixed_effects_vars_computed_as_factor,
        random_effects_vars,
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
        sig_mat_, k_mat, sigma2_u,
        percent_eig_,
        non_zero_precision_eig_,
        alpha_frob_
      )
      w_mat <- white_obj$w_mat
      sig_mat_ <- white_obj$sig_mat_

      # whiten x_mat using w_mat
      # NB. intercept is already present in x_mat and x_mat_tilde
      x_mat_tilde <- w_mat %*% x_mat

      # get raw phenotypes associated to common genotypes
      y <- raw_pheno_df[, trait_]

      # get ols estimates for fixed effects and xi
      beta_hat <- ginv(t(x_mat_tilde) %*% x_mat_tilde) %*% t(x_mat_tilde) %*% y
      y_hat <- x_mat_tilde %*% beta_hat
      xi_hat <- y - y_hat

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
        "y" = y
      ))
    },
    error = function(e) {
      cat(
        "Error with : ", conditionMessage(e), "\n"
      )
    }
  )
}
