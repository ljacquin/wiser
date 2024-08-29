match_indices <- function(raw_pheno_df, omic_df) {
  common_indices <- unlist(lapply(rownames(omic_df), function(g) {
    which(raw_pheno_df$Genotype == g)
  }))
  return(raw_pheno_df[common_indices, ])
}
