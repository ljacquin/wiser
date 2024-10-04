# reduce dataset based on genotype counts
reduce_dataset_based_on_genotypes <- function(df_, nrow_approx_lim = 10e3,
                                              min_samples = 30) {
  K <- nrow(df_) / nrow_approx_lim
  df_reduced <- df_ %>%
    group_by(Genotype) %>%
    group_modify(~ {
      n_genotype <- nrow(.x)
      if (n_genotype < min_samples) {
        samples_to_take <- n_genotype
      } else {
        samples_to_take <- ceiling(n_genotype / K)
      }
      sample_n(.x, samples_to_take)
    }) %>%
    ungroup()

  return(df_reduced)
}
