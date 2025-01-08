# function which computes incidence matrices for fixed and random effects
# NB. column of ones is added for intercept, and is associated to first
# fixed effect (which must be factor or numeric) during construction
compute_incidence_matrices_fixed_and_random_effects <- function(
    fixed_effect_vars,
    fixed_effect_vars_computed_as_factor,
    random_effects_vars,
    raw_pheno_df) {
  # define list of incidence matrices for fixed effects
  list_x_mat <- vector("list", length(fixed_effect_vars))
  names(list_x_mat) <- fixed_effect_vars
  
  # add incidence matrix for first fixed effect to list of matrices
  fix_eff_var_ <- fixed_effect_vars[1]
  if ((length(fixed_effect_vars_computed_as_factor) > 0) &&
      (fix_eff_var_ %in% fixed_effect_vars_computed_as_factor)
  ) {
    list_x_mat[[fix_eff_var_]] <- model.matrix(
      as.formula(paste0("~", fix_eff_var_)),
      data = raw_pheno_df
    )
    colnames(list_x_mat[[fix_eff_var_]]) <- str_replace_all(
      colnames(list_x_mat[[fix_eff_var_]]),
      pattern = fix_eff_var_, replacement = paste0(fix_eff_var_, "_")
    )
  } else {
    list_x_mat[[fix_eff_var_]] <- cbind(
      rep(1, nrow(raw_pheno_df)),
      raw_pheno_df[, fix_eff_var_]
    )
    colnames(list_x_mat[[fix_eff_var_]]) <- c("Intercept", fix_eff_var_)
  }
  # add incidence matrices (without intercept) for other fixed effects to list
  for (fix_eff_var_ in fixed_effect_vars[-1]) {
    if ((length(fixed_effect_vars_computed_as_factor) > 0) &&
        (fix_eff_var_ %in% fixed_effect_vars_computed_as_factor)
    ) {
      list_x_mat[[fix_eff_var_]] <- model.matrix(
        as.formula(paste0("~", fix_eff_var_, " - 1")),
        data = raw_pheno_df
      )
      colnames(list_x_mat[[fix_eff_var_]]) <- str_replace_all(
        colnames(list_x_mat[[fix_eff_var_]]),
        pattern = fix_eff_var_, replacement = paste0(fix_eff_var_, "_")
      )
    } else {
      list_x_mat[[fix_eff_var_]] <- raw_pheno_df[, fix_eff_var_]
      names(list_x_mat[[fix_eff_var_]]) <- fix_eff_var_
    }
  }
  x_mat <- do.call(cbind, list_x_mat)
  x_mat <- apply(x_mat, 2, as.numeric)
  
  # define list of incidence matrices for random effects
  list_z_mat <- vector("list", length(random_effects_vars))
  names(list_z_mat) <- random_effects_vars
  
  # add incidence matrices for random effects to list
  for (rand_eff_var in random_effects_vars) {
    # make sure effect is indeed a factor
    raw_pheno_df[, rand_eff_var] <- as.factor(raw_pheno_df[, rand_eff_var])
    # build incidence matrix for random effect rand_eff_var
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
  return(
    list(
      "x_mat" = x_mat,
      "z_mat" = z_mat
    )
  )
}