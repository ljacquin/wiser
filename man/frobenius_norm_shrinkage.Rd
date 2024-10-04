\name{frobenius_norm_shrinkage}
\alias{frobenius_norm_shrinkage}
\title{Shrinkage of covariance matrix towards a scaled identity matrix using Frobenius norm}
\description{
  This function shrinks a covariance matrix towards a scaled identity matrix, where the scaling factor is determined by the Frobenius norm of the matrix. The shrinkage intensity is controlled by \code{alpha_}.
  
  The shrinkage is performed according to the formula:
    \deqn{\Sigma_{\text{shrink}} = (1 - \alpha) \Sigma + \alpha \|\Sigma\|_F I_n}
  where \eqn{\|\Sigma\|_F} is the Frobenius norm of the covariance matrix, defined as:
    \deqn{\|\Sigma\|_F = \sqrt{\sum_{i,j} \Sigma_{ij}^2}}
  and \eqn{I_n} is the identity matrix of size \eqn{n}.
}
\usage{
  frobenius_norm_shrinkage(cov_mat_, alpha_)
}
\arguments{
  \item{cov_mat_}{A square numeric matrix representing the sample covariance matrix.}
  \item{alpha_}{A numeric value between 0 and 1 that controls the shrinkage intensity. A value of 0 corresponds to no shrinkage, and a value of 1 shrinks fully towards a scaled identity matrix using the Frobenius norm.}
}
\value{
  A shrunken covariance matrix with off-diagonal elements reduced, and diagonal elements scaled by the Frobenius norm of the matrix.
}
\author{
  Laval Jacquin, \email{jacquin.julien@gmail.com}
}
\examples{
  sample_cov <- matrix(c(1, 0.8, 0.8, 1), nrow = 2)
  shrunk_cov <- frobenius_norm_shrinkage(sample_cov, alpha_ = 0.5)
  print(shrunk_cov)
}
