# function which match indices (not only the first one)
match_indices <- function(raw_pheno_df, mat_) {
  common_indices <- unlist(lapply(rownames(mat_), function(g) {
    which(raw_pheno_df$Genotype == g)
  }))
  return(raw_pheno_df[common_indices, ])
}
