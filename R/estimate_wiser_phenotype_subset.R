# function which estimate wiser phenotype for a subset of genotypes
estimate_wiser_phenotype_subset <- function(
    i, idx_geno_shuff_, n_idx_geno_subset,
    geno_df, raw_pheno_df, trait_,
    fixed_effects_vars, random_effects_vars,
    init_sigma2_u, init_sigma2_e,
    n_sim_abc, seed_abc, quantile_threshold_abc, nb_iter_abc,
    kernel_type, rate_decay_kernel, whitening_method, n_geno_subsets_) {
  tryCatch(
    {
      if (i < n_geno_subsets_) {
        idx_geno_subset_ <- idx_geno_shuff_[
          ((i - 1) * n_idx_geno_subset + 1):(i * n_idx_geno_subset)
        ]
      } else {
        idx_geno_subset_ <- idx_geno_shuff_[
          ((i - 1) * n_idx_geno_subset + 1):length(idx_geno_shuff_)
        ]
      }
      
      sub_geno_df <- geno_df[idx_geno_subset_, ]
      
      sub_raw_pheno_df <- match_indices(raw_pheno_df, sub_geno_df)
      
      sub_wiser_obj <- estimate_wiser_phenotype(
        geno_df = sub_geno_df,
        raw_pheno_df = sub_raw_pheno_df,
        trait_ = trait_,
        fixed_effects_vars = fixed_effects_vars,
        random_effects_vars = random_effects_vars,
        init_sigma2_u = init_sigma2_u,
        init_sigma2_e = init_sigma2_e,
        n_sim_abc = n_sim_abc,
        seed_abc = seed_abc,
        quantile_threshold_abc = quantile_threshold_abc,
        nb_iter_abc = nb_iter_abc,
        kernel_type = kernel_type,
        rate_decay_kernel = rate_decay_kernel,
        whitening_method = whitening_method
      )
      sub_v_hat_and_var_comp <- data.frame(
        "v_hat" = sub_wiser_obj$v_hat,
        "sigma2_u_hat_mean_subset_number" =
          sub_wiser_obj$var_comp_abc_obj$sigma2_u_hat_mean,
        "sigma2_e_hat_mean_subset_number" =
          sub_wiser_obj$var_comp_abc_obj$sigma2_e_hat_mean,
        "parallelized_subset_number" = i
      )
      rownames(sub_v_hat_and_var_comp) <- str_replace_all(
        rownames(sub_v_hat_and_var_comp),
        pattern = "Genotype_", replacement = ""
      )
      return(sub_v_hat_and_var_comp)
    },
    error = function(e) {
      cat("Error in subset", i, ":", conditionMessage(e), "\n")
      return(NULL)
    }
  )
}
