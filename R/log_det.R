# function which computes the log of determinant
log_det <- function(Sigma) {
  return(2 * sum(log(diag(cholesky(Sigma, parallel = T)))))
}
