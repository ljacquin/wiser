compute_whitening_matrix_for_sig_mat_ <- function(whitening_method,
                                                  regularization_method,
                                                  parallelized_cholesky,
                                                  sig_mat_, k_mat, sigma2_u,
                                                  percent_eig_,
                                                  non_zero_precision_eig_,
                                                  alpha_frob_) {
  # regularize covariance matrix, by adding a strictly positive value to the
  # diagonal of Σ, to ensure its positive definiteness
  if (regularization_method == "mean_small_eigenvalues") {
    sig_mat_ <- regularize_covariance_mean_small_eigenvalues(
      sig_mat_, k_mat, sigma2_u, percent_eig_, non_zero_precision_eig_
    )
  } else if (regularization_method == "mean_eigenvalues") {
    sig_mat_ <- regularize_covariance_mean_eigenvalues(
      sig_mat_
    )
  } else if (regularization_method == "frobenius_norm") {
    sig_mat_ <- regularize_covariance_frobenius_norm(
      sig_mat_, alpha_frob_
    )
  }
  # compute whitening matrix, either from ZCA-cor or cholesky
  # decomposition (i.e. Σ = LL' )
  if (whitening_method == "ZCA-cor") {
    # compute w_mat from ZCA-cor
    w_mat <- whiteningMatrix(sig_mat_, method = "ZCA-cor")
  } else if (whitening_method == "PCA-cor") {
    # compute w_mat from ZCA-cor
    w_mat <- whiteningMatrix(sig_mat_, method = "PCA-cor")
  } else {
    # compute w_mat = L^−1 from Cholesky decomposition
    if (parallelized_cholesky) {
      L <- t(cholesky(sig_mat_, parallel = T))
    } else {
      L <- t(cholesky(sig_mat_, parallel = F))
    }
    w_mat <- forwardsolve(L, diag(nrow(L)))
  }
  return(list(
    "sig_mat_" = sig_mat_,
    "w_mat" = w_mat
  ))
}
