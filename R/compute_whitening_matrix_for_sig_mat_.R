# function which computes the whitening matrix for sig_mat_ based on the
# selected whitening method
compute_whitening_matrix_for_sig_mat_ <- function(whitening_method,
                                                  regularization_method,
                                                  parallelized_cholesky,
                                                  sig_mat_, alpha_) {
  # regularize covariance matrix, by adding a strictly positive value to the
  # diagonal of Σ, to ensure its positive definiteness
  if (regularization_method == "frobenius_norm_regularization") {
    sig_mat_ <- frobenius_norm_regularization(
      sig_mat_, alpha_
    )
  } else if (regularization_method == "frobenius_norm_shrinkage") {
    sig_mat_ <- frobenius_norm_shrinkage(
      sig_mat_, alpha_
    )
  } else if (regularization_method == "frobenius_norm_partial_shrinkage") {
    sig_mat_ <- frobenius_norm_partial_shrinkage(
      sig_mat_, alpha_
    )
  } else if (regularization_method == "trace_sample_variance_regularization") {
    sig_mat_ <- trace_sample_variance_regularization(
      sig_mat_, alpha_
    )
  } else if (regularization_method == "trace_sample_variance_shrinkage") {
    sig_mat_ <- trace_sample_variance_shrinkage(
      sig_mat_, alpha_
    )
  } else if (regularization_method == "trace_sample_variance_partial_shrinkage") {
    sig_mat_ <- trace_sample_variance_partial_shrinkage(
      sig_mat_, alpha_
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
