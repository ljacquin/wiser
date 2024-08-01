\name{regularize_covariance_frobenius_norm}
\alias{regularize_covariance_frobenius_norm}
\title{Regularize a covariance matrix based on Frobenius norm to ensure positive definiteness}
\description{
  Regularizes a covariance matrix to ensure it is positive definite by adding a positive value \code{delta} to the diagonal. The value of \code{delta} is computed based on the Frobenius norm of the matrix.
}
\usage{
  regularize_covariance_frobenius_norm(cov_mat_, alpha_frob_)
}
\arguments{
  \item{cov_mat_}{A numeric matrix representing the covariance matrix to be regularized.}
  \item{alpha_frob_}{A numeric value used for regularization, based on the Frobenius norm of \code{cov_mat_}.}
}
\value{
  A numeric matrix representing the regularized covariance matrix, which is positive definite.
}
\details{
  To ensure that the covariance matrix \code{cov_mat_} is positive definite, a regularization term \code{delta} is added to the diagonal elements of the matrix. The value of \code{delta} is computed as follows:

  \deqn{\delta = \alpha_{\text{frob}} \times \| \Sigma \|_F}

  where \eqn{\| \Sigma \|_F} is the Frobenius norm of the covariance matrix \eqn{\Sigma}.

  The regularized matrix \eqn{\Sigma_{\text{reg}}} is then computed by adding \eqn{\delta} to the diagonal elements of \eqn{\Sigma}:

  \deqn{\Sigma_{\text{reg}} = \Sigma + \delta \mathbf{I}}

  Here, \eqn{\mathbf{I}} denotes the identity matrix of the same dimensions as \eqn{\Sigma}. This regularization step ensures that the resulting matrix \eqn{\Sigma_{\text{reg}}} is positive definite.
}
\author{
  Laval Jacquin, \email{jacquin.julien@gmail.com}
}
\examples{
  # Example usage
  cov_mat_ <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
  alpha_frob_ <- 0.01
  regularized_cov_mat <- regularize_covariance_frobenius_norm(cov_mat_, alpha_frob_)
  print(regularized_cov_mat)
}
