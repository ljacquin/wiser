# function which computes transformed fixed variables and least squares
compute_transformed_vars_and_ols_estimates <- function(
    geno_df, raw_pheno_df, fixed_effects_vars, random_effects_vars, trait_,
    compute_row_and_position_as_factors,
    sigma2_u, sigma2_e, kernel_type,
    rate_decay_kernel,
    whitening_method,
    regularization_method,
    alpha_frob_,
    percent_eig_,
    non_zero_precision_eig_) {
  tryCatch(
    {
      # remove all rows with na w.r.t to trait_
      raw_pheno_df <- raw_pheno_df %>% drop_na(all_of(trait_))

      # define variables of interest
      sel_vars_ <- c(fixed_effects_vars, random_effects_vars, trait_)

      # get only variables of interest from raw_pheno_df
      raw_pheno_df <- raw_pheno_df[, sel_vars_]
      raw_pheno_df <- na.omit(raw_pheno_df)

      # compute Gram matrix (i.e. genomic covariance matrix)
      geno_names <- rownames(geno_df)
      geno_df <- apply(geno_df, 2, as.numeric)

      if (kernel_type == "linear") {
        k_mat <- tcrossprod(scale(apply(geno_df, 2, as.numeric),
          center = T, scale = F
        ))
      } else {
        # kernel identity is not recommended due to constrained hypothesis about
        # genotypes independence which may lead to low precision
        k_mat <- as.matrix(diag(nrow(geno_df)))
      }

      # test positive definiteness and force it if necessary
      if (!is.positive.definite(k_mat, tol = 1e-8)) {
        k_mat <- as.matrix(nearPD(k_mat)$mat)
      }

      # assign genotype rownames and colnames to k_mat
      colnames(k_mat) <- rownames(k_mat) <- geno_names

      # get common geontype between raw_pheno_df and geno_df
      raw_pheno_df <- match_indices(raw_pheno_df, k_mat)

      # convert fixed effects variables to factors, and remove
      # buffer for management if exists
      for (fix_eff_var_ in fixed_effects_vars) {
        if ((fix_eff_var_ == "Row" || fix_eff_var_ == "Position") &&
          !compute_row_and_position_as_factors) {
          raw_pheno_df[, fix_eff_var_] <- raw_pheno_df[, fix_eff_var_]
        } else {
          raw_pheno_df[, fix_eff_var_] <- as.factor(
            raw_pheno_df[, fix_eff_var_]
          )
        }
        if ("BUFFER" %in% raw_pheno_df[, fix_eff_var_]) {
          raw_pheno_df <- raw_pheno_df[
            raw_pheno_df[, fix_eff_var_] != "BUFFER",
          ]
        }
      }
      # droplevels in order to remove levels which don't exist anymore
      raw_pheno_df <- droplevels(raw_pheno_df)

      # get raw phenotypes associated to common genotypes
      y <- raw_pheno_df[, trait_]

      # get incidence matrices for fixed and random effects
      # NB. column of ones is added for intercept associated to fixed effects

      # define list of incidence matrices for fixed effects
      list_x_mat <- vector("list", length(fixed_effects_vars))
      names(list_x_mat) <- fixed_effects_vars

      # add incidence matrix for first fixed effect to list of matrices
      fix_eff_var_ <- fixed_effects_vars[1]
      list_x_mat[[fix_eff_var_]] <- model.matrix(
        as.formula(paste0("~", fix_eff_var_)),
        data = raw_pheno_df
      )
      colnames(list_x_mat[[fix_eff_var_]]) <- str_replace_all(
        colnames(list_x_mat[[fix_eff_var_]]),
        pattern = fix_eff_var_, replacement = paste0(fix_eff_var_, "_")
      )

      # add incidence matrices (without intercept) for other fixed effects to list
      for (fix_eff_var_ in fixed_effects_vars[-1]) {
        if ((fix_eff_var_ == "Row" || fix_eff_var_ == "Position") &&
          !compute_row_and_position_as_factors) {
          list_x_mat[[fix_eff_var_]] <- raw_pheno_df[, fix_eff_var_]
          names(list_x_mat[[fix_eff_var_]]) <- fix_eff_var_
        } else {
          list_x_mat[[fix_eff_var_]] <- model.matrix(
            as.formula(paste0("~", fix_eff_var_, " - 1")),
            data = raw_pheno_df
          )
          colnames(list_x_mat[[fix_eff_var_]]) <- str_replace_all(
            colnames(list_x_mat[[fix_eff_var_]]),
            pattern = fix_eff_var_, replacement = paste0(fix_eff_var_, "_")
          )
        }
      }
      x_mat <- do.call(cbind, list_x_mat)
      x_mat <- apply(x_mat, 2, as.numeric)

      # define list of incidence matrices for random effects
      list_z_mat <- vector("list", length(random_effects_vars))
      names(list_z_mat) <- random_effects_vars

      # add incidence matrices for random effects to list
      for (rand_eff_var in random_effects_vars) {
        list_z_mat[[rand_eff_var]] <- model.matrix(
          as.formula(paste0("~", rand_eff_var, " - 1")),
          data = raw_pheno_df
        )
        colnames(list_z_mat[[rand_eff_var]]) <- str_replace_all(
          colnames(list_z_mat[[rand_eff_var]]),
          pattern = rand_eff_var, replacement = paste0(rand_eff_var, "_")
        )
      }
      z_mat <- do.call(cbind, list_z_mat)
      z_mat <- apply(z_mat, 2, as.numeric)

      # compute Σu, i.e. sig_mat_ here
      sig_mat_ <- sigma2_u * crossprod(t(z_mat), tcrossprod(k_mat, z_mat))

      # regularize covariance matrix, by adding a strictly positive value to the
      # diagonal of Σ, to ensure its positive definiteness
      if (regularization_method == "mean_small_eigenvalues") {
        sig_mat_ <- regularize_covariance_mean_small_eigenvalues(
          sig_mat_, k_mat, sigma2_u, percent_eig_, non_zero_precision_eig_
        )
      } else if (regularization_method == "mean_eigenvalues") {
        sig_mat_ <- regularize_covariance_mean_eigenvalues(
          sig_mat_
        )
      } else if (regularization_method == "frobenius_norm") {
        sig_mat_ <- regularize_covariance_frobenius_norm(
          sig_mat_, alpha_frob_
        )
      }

      # compute whitening matrix, either from cholesky decomposition,
      # i.e. Σ = LL', or ZCA-cor
      if (whitening_method == "Cholesky") {
        # compute w_mat = L^−1 from Cholesky decomposition
        L <- t(cholesky(sig_mat_, parallel = T))
        w_mat <- forwardsolve(L, diag(nrow(L)))
      } else {
        # compute w_mat from ZCA-cor
        w_mat <- whiteningMatrix(sig_mat_, method = "ZCA-cor")
      }

      # NB. intercept is already present in x_mat and x_mat_tilde
      x_mat_tilde <- w_mat %*% x_mat

      # get ols estimates for fixed effects and xi
      beta_hat <- ginv(t(x_mat_tilde) %*% x_mat_tilde) %*% t(x_mat_tilde) %*% y
      xi_hat <- y - x_mat_tilde %*% beta_hat

      return(list(
        "x_mat" = x_mat,
        "z_mat" = z_mat,
        "k_mat" = k_mat,
        "beta_hat" = beta_hat,
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
