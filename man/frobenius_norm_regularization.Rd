\name{frobenius_norm_regularization}
\alias{frobenius_norm_regularization}
\title{Regularization of covariance matrix using Frobenius norm}
\description{
  This function regularizes a covariance matrix by adding a scaled identity matrix, where the scaling factor is determined by the Frobenius norm of the original covariance matrix. The regularization intensity is governed by \code{alpha_}.
  
  The regularization is performed using the formula:
    \deqn{\Sigma_{\text{regularized}} = \Sigma + \alpha \|\Sigma\|_F I_n}
  where \eqn{\|\Sigma\|_F} is the Frobenius norm of the covariance matrix, and \eqn{I_n} is the identity matrix of size \eqn{n}.
}
\usage{
  frobenius_norm_regularization(cov_mat_, alpha_)
}
\arguments{
  \item{cov_mat_}{A square numeric matrix representing the sample covariance matrix.}
  \item{alpha_}{A numeric value controlling the regularization intensity. A higher value of \code{alpha_} increases the influence of the Frobenius norm on the regularization.}
}
\value{
  A regularized covariance matrix with diagonal elements scaled by the Frobenius norm.
}
\author{
  Laval Jacquin, \email{jacquin.julien@gmail.com}
}
\examples{
  sample_cov <- matrix(c(1, 0.8, 0.8, 1), nrow = 2)
  reg_cov <- frobenius_norm_regularization(sample_cov, alpha_ = 0.5)
  print(reg_cov)
}
