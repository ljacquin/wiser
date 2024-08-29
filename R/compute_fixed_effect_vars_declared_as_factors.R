compute_fixed_effect_vars_declared_as_factors <- function(
    raw_pheno_df,
    fixed_effects_vars,
    fixed_effects_vars_computed_as_factor,
    site_var,
    fixed_effects_vars_computed_as_factor_by_site) {
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
      # compute fix_eff_var_ as factor by site, for those declared as,
      # otherwise compute fix_eff_var_ as a factor only
      if ((length(site_var) > 0 &&
           length(fixed_effects_vars_computed_as_factor_by_site) > 0) &&
          (fix_eff_var_ %in% fixed_effects_vars_computed_as_factor_by_site)
      ) {
        site_var_fix_eff_var_ <- paste0(
          raw_pheno_df[, site_var], "_", fix_eff_var_, "_",
          raw_pheno_df[, fix_eff_var_]
        )
        raw_pheno_df[, fix_eff_var_] <- as.factor(site_var_fix_eff_var_)
      } else {
        raw_pheno_df[, fix_eff_var_] <- as.factor(
          raw_pheno_df[, fix_eff_var_]
        )
      }
    }
  }
  return(droplevels(raw_pheno_df))
}