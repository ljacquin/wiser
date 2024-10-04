# function which computes phenotypes approximating genetic values using whitening
estimate_wiser_phenotype <- function(omic_df, raw_pheno_df, trait_,
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
                                     init_sigma2_u = 1,
                                     init_sigma2_e = 1,
                                     prior_scale_factor = 1e3,
                                     n_sim_abc = 100,
                                     seed_abc = 123,
                                     quantile_threshold_abc = 0.05,
                                     nb_iter_abc = 1,
                                     kernel_type = "linear",
                                     whitening_method = "ZCA-cor",
                                     regularization_method = "frobenius_norm_regularization",
                                     alpha_ = 0.01,
                                     parallelized_cholesky = T,
                                     reduce_raw_dataset_size_ = T,
                                     nrow_approx_lim_raw_dataset_zca_cor = 20e3,
                                     nrow_approx_lim_raw_dataset_pca_cor = 20e3,
                                     nrow_approx_lim_raw_dataset_chol = 40e3) {
  tryCatch(
    {
      # compute transformed variables associated to fixed effects and least-squares
      # to estimate these
      transform_and_ls_obj <- compute_transformed_vars_and_ols_estimates(
        omic_df, raw_pheno_df, trait_,
        fixed_effects_vars,
        fixed_effects_vars_computed_as_factor,
        site_var,
        fixed_effects_vars_computed_as_factor_by_site,
        random_effects_vars,
        sigma2_u = init_sigma2_u,
        sigma2_e = init_sigma2_e,
        kernel_type,
        whitening_method,
        regularization_method,
        alpha_,
        parallelized_cholesky,
        reduce_raw_dataset_size_,
        nrow_approx_lim_raw_dataset_zca_cor,
        nrow_approx_lim_raw_dataset_pca_cor,
        nrow_approx_lim_raw_dataset_chol
      )

      # get an upper bound for sigma2_u et sigma2_e priors
      prior_sigma2_upper_bound <- var(
        transform_and_ls_obj$y
      )

      for (iter_ in 1:nb_iter_abc) {
        # compute variance components using abc
        var_comp_abc_obj <- abc_variance_component_estimation(
          y = transform_and_ls_obj$y,
          x_mat = transform_and_ls_obj$x_mat,
          z_mat = transform_and_ls_obj$z_mat,
          k_mat = transform_and_ls_obj$k_mat,
          beta_hat = transform_and_ls_obj$beta_hat,
          prior_sigma2_u = c(
            ceiling(prior_sigma2_upper_bound / prior_scale_factor),
            prior_sigma2_upper_bound
          ),
          prior_sigma2_e = c(
            ceiling(prior_sigma2_upper_bound / prior_scale_factor),
            prior_sigma2_upper_bound
          ),
          n_sim_abc, seed_abc,
          quantile_threshold_abc
        )
        # compute variance components again with abc using new estimates
        transform_and_ls_obj <- compute_transformed_vars_and_ols_estimates(
          omic_df, raw_pheno_df, trait_,
          fixed_effects_vars,
          fixed_effects_vars_computed_as_factor,
          site_var,
          fixed_effects_vars_computed_as_factor_by_site,
          random_effects_vars,
          sigma2_u = var_comp_abc_obj$sigma2_u_hat_mean,
          sigma2_e = var_comp_abc_obj$sigma2_e_hat_mean,
          kernel_type,
          whitening_method,
          regularization_method,
          alpha_,
          parallelized_cholesky,
          reduce_raw_dataset_size_,
          nrow_approx_lim_raw_dataset_zca_cor,
          nrow_approx_lim_raw_dataset_pca_cor,
          nrow_approx_lim_raw_dataset_chol
        )
      }

      # get estimated and/or modified components after abc

      # extract marker data in case of modification
      omic_df <- transform_and_ls_obj$omic_df

      # extract estimated fixed effects
      beta_hat <- transform_and_ls_obj$beta_hat

      # compute phenotypic values using ols
      v_hat <- ginv(t(transform_and_ls_obj$z_mat) %*% transform_and_ls_obj$z_mat) %*%
        t(transform_and_ls_obj$z_mat) %*% transform_and_ls_obj$xi_hat

      # save wiser fixed effects estimates in a data frame
      wiser_pheno_df <- data.frame(
        "Genotype" = str_replace_all(
          colnames(transform_and_ls_obj$z_mat),
          pattern = "Genotype_",
          replacement = ""
        ),
        "v_hat" = v_hat
      )

      # save wiser phenotypes in a data frame
      wiser_fix_eff_df <- data.frame(
        "fixed_effect_var" = colnames(transform_and_ls_obj$x_mat),
        "beta_hat_var" = beta_hat
      )

      return(list(
        "wiser_omic_data" = omic_df,
        "sig_mat_u" = transform_and_ls_obj$sig_mat_,
        "w_mat" = transform_and_ls_obj$w_mat,
        "wiser_fixed_effect_estimates" = wiser_fix_eff_df,
        "wiser_abc_variance_component_estimates" = var_comp_abc_obj,
        "wiser_phenotypes" = wiser_pheno_df,
        "wiser_z_mat" = transform_and_ls_obj$z_mat,
        "wiser_x_mat" = transform_and_ls_obj$x_mat,
        "wiser_x_mat_tilde" = transform_and_ls_obj$x_mat_tilde,
        "wiser_xi_hat" = transform_and_ls_obj$xi_hat,
        "wiser_y_hat" = transform_and_ls_obj$y_hat,
        "wiser_y" = transform_and_ls_obj$y
      ))
    },
    error = function(e) {
      cat(
        "Error with : ", conditionMessage(e), "\n"
      )
    }
  )
}
