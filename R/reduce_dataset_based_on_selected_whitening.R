reduce_dataset_based_on_selected_whitening <- function(
    whitening_method,
    raw_pheno_df,
    nrow_lim_raw_dataset_zca_cor,
    nrow_lim_raw_dataset_pca_cor,
    nrow_lim_raw_dataset_chol) {
  set.seed(123)
  if (whitening_method == "ZCA-cor") {
    raw_pheno_df <- as.data.frame(
      reduce_dataset_based_on_genotypes(
        df_ = raw_pheno_df,
        nrow_lim = nrow_lim_raw_dataset_zca_cor
      )
    )
  } else if (whitening_method == "PCA-cor") {
    raw_pheno_df <- as.data.frame(
      reduce_dataset_based_on_genotypes(
        df_ = raw_pheno_df,
        nrow_lim = nrow_lim_raw_dataset_pca_cor
      )
    )
  } else {
    raw_pheno_df <- as.data.frame(
      reduce_dataset_based_on_genotypes(
        df_ = raw_pheno_df,
        nrow_lim = nrow_lim_raw_dataset_chol
      )
    )
  }
  return(raw_pheno_df)
}