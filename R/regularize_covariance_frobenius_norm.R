# function which makes a covariance matrix positive definite by adding a
# positive value delta to the diagonal, based on l2 norm
regularize_covariance_frobenius_norm <- function(cov_mat_, alpha_frob_) {
  n <- nrow(cov_mat_)
  delta_ <- alpha_frob_ * frobenius_norm(cov_mat_)
  cov_mat_ <- cov_mat_ + delta_ * diag(n)
  return(cov_mat_)
}