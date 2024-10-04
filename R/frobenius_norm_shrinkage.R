# function which makes a covariance matrix positive definite by using
# a convex shrinkage estimator which adds a strictly positive diagonal matrix
# based on the Frobenius norm
frobenius_norm_shrinkage <- function(cov_mat_, alpha_) {
  n <- nrow(cov_mat_)
  cov_mat_ <- (1 - alpha_) * cov_mat_ + alpha_ * frobenius_norm(cov_mat_) * diag(n)
  return(cov_mat_)
}
