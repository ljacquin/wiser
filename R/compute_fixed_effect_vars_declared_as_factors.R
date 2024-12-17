# function which computes fixed effect vars as factors for those declared as
compute_fixed_effect_vars_declared_as_factors <- function(
    raw_pheno_df,
    fixed_effects_vars,
    fixed_effects_vars_computed_as_factor,
    envir_var,
    fixed_effects_vars_computed_as_factor_by_envir) {
  # convert fixed effects variables to factors for those declared
  for (fix_eff_var_ in fixed_effects_vars) {
    if ((length(fixed_effects_vars_computed_as_factor) > 0) &&
      (fix_eff_var_ %in% fixed_effects_vars_computed_as_factor)
    ) {
      # if fix_eff_var_ is equal to management, remove buffer if exists
      if ("buffer" %in% tolower(raw_pheno_df[, fix_eff_var_])) {
        raw_pheno_df <- raw_pheno_df[
          tolower(raw_pheno_df[, fix_eff_var_]) != "buffer",
        ]
      }
      # compute fix_eff_var_ as factor by envir, for those declared as,
      # otherwise compute fix_eff_var_ as a factor only
      if ((length(envir_var) > 0 &&
        length(fixed_effects_vars_computed_as_factor_by_envir) > 0) &&
        (fix_eff_var_ %in% fixed_effects_vars_computed_as_factor_by_envir)
      ) {
        envir_var_fix_eff_var_ <- paste0(
          raw_pheno_df[, envir_var], "_", fix_eff_var_, "_",
          raw_pheno_df[, fix_eff_var_]
        )
        raw_pheno_df[, fix_eff_var_] <- as.factor(envir_var_fix_eff_var_)
      } else {
        raw_pheno_df[, fix_eff_var_] <- as.factor(
          raw_pheno_df[, fix_eff_var_]
        )
      }
    }
  }
  return(droplevels(raw_pheno_df))
}
