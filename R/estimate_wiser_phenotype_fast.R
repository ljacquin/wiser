# function which accelarates computation of phenotypes  approximating
# genetic values using whitening
estimate_wiser_phenotype_fast <- function(geno_df, raw_pheno_df, trait_,
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
                                          whitening_method = "Cholesky",
                                          n_geno_subsets_ = 5) {
  tryCatch(
    {
      # compute several sets of randomized indices for parallelization
      set.seed(123)
      idx_geno_shuff_ <- sample(1:nrow(geno_df), replace = F)
      n_idx_geno_subset <- floor(nrow(geno_df) / n_geno_subsets_)
      
      # Configure parallel environment
      plan(multisession, workers = floor(detectCores() / 2))
      
      # apply parallel function
      v_hat_and_var_comp_list <- future_lapply(1:n_geno_subsets_, function(i) {
        estimate_wiser_phenotype_subset(
          i, idx_geno_shuff_, n_idx_geno_subset,
          geno_df, raw_pheno_df, trait_,
          fixed_effects_vars, random_effects_vars,
          init_sigma2_u, init_sigma2_e,
          n_sim_abc, seed_abc, quantile_threshold_abc, nb_iter_abc,
          kernel_type, rate_decay_kernel, whitening_method,
          n_geno_subsets_
        )
      })
      v_hat_and_var_comp <- do.call(rbind, v_hat_and_var_comp_list)
      v_hat_and_var_comp <- v_hat_and_var_comp[
        match(rownames(geno_df), rownames(v_hat_and_var_comp)),
      ]
      
      sigma2_u_hat_mean <- colMeans(
        v_hat_and_var_comp["sigma2_u_hat_mean_subset_number"]
      )
      names(sigma2_u_hat_mean) <- NULL
      
      sigma2_e_hat_mean <- colMeans(
        v_hat_and_var_comp["sigma2_e_hat_mean_subset_number"]
      )
      names(sigma2_e_hat_mean) <- NULL
      
      return(list(
        "v_hat_and_var_comp" = v_hat_and_var_comp,
        "sigma2_u_hat_mean" = sigma2_u_hat_mean,
        "sigma2_e_hat_mean" = sigma2_e_hat_mean
      ))
    },
    error = function(e) {
      cat(
        "Error with : ", conditionMessage(e), "\n"
      )
    }
  )
}
