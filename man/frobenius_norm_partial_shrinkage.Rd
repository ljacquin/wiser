\name{frobenius_norm_partial_shrinkage}
\alias{frobenius_norm_partial_shrinkage}
\title{Partial shrinkage of diagonal elements of a covariance matrix using Frobenius norm}
\description{
  This function applies partial shrinkage to the diagonal elements of a covariance matrix by adjusting them based on the Frobenius norm of the matrix. The off-diagonal elements remain unchanged. The shrinkage intensity is controlled by \code{alpha_}.
  
  The shrinkage is performed according to the formula:
    \deqn{\Sigma_{\text{shrink},ii} = (1 - \alpha) \Sigma_{ii} + \alpha \|\Sigma\|_F}
  where \eqn{\|\Sigma\|_F} is the Frobenius norm of the covariance matrix, defined as:
    \deqn{\|\Sigma\|_F = \sqrt{\sum_{i,j} \Sigma_{ij}^2}}
  and the off-diagonal elements of \eqn{\Sigma} remain unchanged.
}
\usage{
  frobenius_norm_partial_shrinkage(cov_mat_, alpha_)
}
\arguments{
  \item{cov_mat_}{A square numeric matrix representing the sample covariance matrix.}
  \item{alpha_}{A numeric value between 0 and 1 that controls the shrinkage intensity for the diagonal elements. A value of 0 corresponds to no shrinkage, and a value of 1 fully replaces the diagonal elements with the Frobenius norm of the matrix.}
}
\value{
  A covariance matrix with its diagonal elements partially shrunk towards the Frobenius norm, while the off-diagonal elements remain unchanged.
}
\author{
  Laval Jacquin, \email{jacquin.julien@gmail.com}
}
\examples{
  # example usage of frobenius_norm_partial_shrinkage
  
  sample_cov <- matrix(c(1, 0.8, 0.8, 1), nrow = 2)
  
  # apply partial shrinkage with alpha = 0.5
  shrunk_cov <- frobenius_norm_partial_shrinkage(sample_cov, alpha_ = 0.5)
  
  # print the shrunken covariance matrix
  print(shrunk_cov)
}
