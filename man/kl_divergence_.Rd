\name{kl_divergence_}
\alias{kl_divergence_}
\title{Compute the Kullback-Leibler divergence between two covariance matrices for gaussian distributions}
\description{
  Computes the Kullback-Leibler (KL) divergence between two covariance matrices in a numerically stable way. The KL divergence measures the difference between two probability distributions. This implementation assumes that both distributions are multivariate Gaussian with zero means.
}
\usage{
  kl_divergence_(Sigma_1, Sigma_2, alpha_frob_)
}
\arguments{
  \item{Sigma_1}{A numeric matrix representing the first covariance matrix.}
  \item{Sigma_2}{A numeric matrix representing the second covariance matrix.}
  \item{alpha_frob_}{A numeric value used for regularization of the covariance matrices to ensure numerical stability.}
}
\value{
  A numeric value representing the Kullback-Leibler divergence between the two covariance matrices.
}
\details{
  The Kullback-Leibler (KL) divergence between two covariance matrices \eqn{\Sigma_1} and \eqn{\Sigma_2} is given by:

  \deqn{D_{KL}(\Sigma_1 \| \Sigma_2) = \frac{1}{2} \left( \text{tr}(\Sigma_2^{-1} \Sigma_1) - k + \log \frac{\text{det}(\Sigma_2)}{\text{det}(\Sigma_1)} \right)}{D_{KL}(\Sigma_1 || \Sigma_2) = \frac{1}{2} \left( tr(\Sigma_2^{-1} \Sigma_1) - k + \log \frac{det(\Sigma_2)}{det(\Sigma_1)} \right)}

  where:
  
  - \eqn{tr()} denotes the trace of a matrix,
  
  - \eqn{det()} denotes the determinant of a matrix,
  
  - \eqn{k} is the dimensionality of the covariance matrices.

  This formula assumes that both multivariate Gaussian distributions have zero means. 

 }
\author{
  Laval Jacquin, \email{jacquin.julien@gmail.com}
}
\examples{
  # Example usage
  Sigma_1 <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
  Sigma_2 <- matrix(c(1, 0.3, 0.3, 1), nrow = 2)
  alpha_frob_ <- 0.01
  kl_divergence_(Sigma_1, Sigma_2, alpha_frob_)
}
