# function which regularizes a covariance matrix by adding a small
# positive value delta to the diagonal
regularize_covariance <- function(cov_mat_, alpha_ = 1e-2) {
  n <- nrow(cov_mat_)
  delta_ <- alpha_ * trace_mat(cov_mat_)
  cov_mat_ <- cov_mat_ + delta_ * diag(n)
  return(cov_mat_)
}
