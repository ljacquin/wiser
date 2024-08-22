reduce_dataset_based_on_genotypes <- function(df_, nrow_lim = 10e3,
                                              min_samples = 3) {
  K <- nrow(df_) / nrow_lim
  
  df_reduced <- df_ %>%
    group_by(Genotype) %>%
    group_modify(~ {
      n_genotype <- nrow(.x)
      samples_to_take <- max(min_samples, ceiling(n_genotype / K))
      sample_n(.x, min(samples_to_take, n_genotype))
    }) %>%
    ungroup()
  
  return(df_reduced)
}