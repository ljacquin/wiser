compute_gram_matrix <- function(omic_df, kernel_type) {
  geno_names <- rownames(omic_df)
  omic_df <- apply(omic_df, 2, as.numeric)
  if (kernel_type == "linear") {
    k_mat <- tcrossprod(scale(apply(omic_df, 2, as.numeric),
                              center = T, scale = F
    ))
  } else {
    # kernel identity is not recommended due to constrained hypothesis about
    # genotypes independence which may lead to low precision
    k_mat <- as.matrix(diag(nrow(omic_df)))
  }
  # test positive definiteness and force it if necessary
  if (!is.positive.definite(k_mat, tol = 1e-8)) {
    k_mat <- as.matrix(nearPD(k_mat)$mat)
  }
  # assign genotype rownames and colnames to k_mat
  colnames(k_mat) <- rownames(k_mat) <- geno_names
  return(k_mat)
}
