raw_pheno_and_marker_based_on_trait_common_genotypes <- function(
    raw_pheno_df,
    omic_df,
    trait_,
    fixed_effects_vars,
    random_effects_vars) {
  # remove all rows with na w.r.t to trait_
  raw_pheno_df <- raw_pheno_df %>% drop_na(all_of(trait_))
  
  # define variables of interest
  sel_vars_ <- c(fixed_effects_vars, random_effects_vars, trait_)
  
  # get only variables of interest from raw_pheno_df
  raw_pheno_df <- raw_pheno_df[, sel_vars_]
  raw_pheno_df <- na.omit(raw_pheno_df)
  
  # droplevels in order to remove levels which don't exist anymore
  raw_pheno_df <- droplevels(raw_pheno_df)
  
  # get phenotypes and marker data based on common genotypes
  raw_pheno_df <- match_indices(raw_pheno_df, omic_df)
  omic_df <- omic_df[rownames(omic_df) %in% unique(raw_pheno_df$Genotype), ]
  
  return(list(
    "raw_pheno_df" = raw_pheno_df,
    "omic_df" = omic_df
  ))
}