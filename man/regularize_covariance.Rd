\name{regularize_covariance}
\alias{regularize_covariance}
\title{Regularize a covariance matrix}
\description{
  Regularizes a covariance matrix by adding a small positive value to the diagonal to ensure positive semi-definiteness.
}
\usage{
  regularize_covariance(cov_mat_, alpha_ = 1e-2)
}
\arguments{
  \item{cov_mat_}{A numeric matrix representing the covariance matrix to be regularized.}
  \item{alpha_}{A numeric value specifying the regularization parameter.}
}
\value{
  A numeric matrix representing the regularized covariance matrix.
}
\examples{
  # Example usage
  # reg_cov_mat <- regularize_covariance(cov_mat_)
}
