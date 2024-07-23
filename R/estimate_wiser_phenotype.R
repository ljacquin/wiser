# function which computes phenotypes approximating genetic values using whitening
estimate_wiser_phenotype <- function(geno_df, raw_pheno_df, trait_,
                                     fixed_effects_vars = c(
                                       "Envir", "Country", "Year",
                                       "Row", "Position", "Management"
                                     ),
                                     random_effects_vars = "Genotype",
                                     init_sigma2_u = 1,
                                     init_sigma2_e = 1,
                                     n_sim_abc = 100,
                                     seed_abc = 123,
                                     quantile_threshold_abc = 0.05,
                                     nb_iter_abc = 1,
                                     kernel_type = "linear",
                                     rate_decay_kernel = 0.1,
                                     whitening_method = "Cholesky") {
  tryCatch(
    {
      # compute transformed variables associated to fixed effects and least-squares
      # to estimate these
      transform_and_ls_obj <- compute_transformed_vars_and_ols_estimates(
        geno_df, raw_pheno_df, fixed_effects_vars, random_effects_vars, trait_,
        sigma2_u = init_sigma2_u,
        sigma2_e = init_sigma2_e,
        kernel_type, rate_decay_kernel,
        whitening_method
      )

      # get an upper bound for sigma2_u et sigma2_e priors
      prior_sigma2_upper_bound <- var(
        transform_and_ls_obj$y
      )

      for (iter_ in 1:nb_iter_abc) {
        # print(paste0('iter : ', iter_))
        # compute variance components using abc
        var_comp_abc_obj <- abc_variance_component_estimation(
          y = transform_and_ls_obj$y,
          x_mat = transform_and_ls_obj$x_mat,
          z_mat = transform_and_ls_obj$z_mat,
          k_mat = transform_and_ls_obj$k_mat,
          beta_hat = transform_and_ls_obj$beta_hat,
          prior_sigma2_u = c(1e-2, prior_sigma2_upper_bound),
          prior_sigma2_e = c(1e-2, prior_sigma2_upper_bound),
          n_sim_abc, seed_abc,
          quantile_threshold_abc
        )
        # compute variance components again with abc using new estimates
        transform_and_ls_obj <- compute_transformed_vars_and_ols_estimates(
          geno_df, raw_pheno_df, fixed_effects_vars, random_effects_vars, trait_,
          sigma2_u = var_comp_abc_obj$sigma2_u_hat_mean,
          sigma2_e = var_comp_abc_obj$sigma2_e_hat_mean,
          kernel_type, rate_decay_kernel,
          whitening_method
        )
      }

      # get estimated components after abc

      # extract estimated fixed effects
      beta_hat <- transform_and_ls_obj$beta_hat

      # compute phenotypic values using ols
      v_hat <- Matrix::solve(t(transform_and_ls_obj$z_mat) %*% transform_and_ls_obj$z_mat) %*%
        t(transform_and_ls_obj$z_mat) %*% transform_and_ls_obj$xi_hat

      return(list(
        "var_comp_abc_obj" = var_comp_abc_obj,
        "beta_hat" = beta_hat,
        "v_hat" = v_hat
      ))
    },
    error = function(e) {
      cat(
        "Error with : ", conditionMessage(e), "\n"
      )
    }
  )
}
