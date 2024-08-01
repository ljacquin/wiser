# function which computes KL divergence in a more stable way numerically
kl_divergence_ <- function(Sigma_1, Sigma_2, alpha_frob_) {
  Sigma_1 <- regularize_covariance_frobenius_norm(Sigma_1, alpha_frob_)
  Sigma_2 <- regularize_covariance_frobenius_norm(Sigma_2, alpha_frob_)

  inv_Sigma_2 <- ginv(Sigma_2)
  term_1 <- sum(diag(inv_Sigma_2 %*% Sigma_1))
  term_2 <- log_det(Sigma_2) - log_det(Sigma_1)

  k <- nrow(Sigma_1)
  kl_div <- 0.5 * (term_1 - k + term_2)

  return(kl_div)
}
