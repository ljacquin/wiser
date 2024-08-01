# function which makes a covariance matrix positive definite by adding a
# positive value delta to the diagonal, based on the trace
regularize_covariance_mean_eigenvalues <- function(cov_mat_) {
  n <- nrow(cov_mat_)
  delta_ <- trace_mat(cov_mat_) / n
  cov_mat_ <- cov_mat_ + delta_ * diag(n)
  return(cov_mat_)
}