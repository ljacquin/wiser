\name{kl_divergence_}
\alias{kl_divergence_}
\title{Compute Kullback-Leibler (KL) divergence between two covariance matrices}
\description{
  Computes the KL divergence between two covariance matrices in a numerically stable manner.
}
\usage{
  kl_divergence_(Sigma_1, Sigma_2)
}
\arguments{
  \item{Sigma_1}{A numeric matrix representing the first covariance matrix.}
  \item{Sigma_2}{A numeric matrix representing the second covariance matrix.}
}
\value{
  A numeric value representing the KL divergence between the two covariance matrices.
}
\examples{
  # Example usage
  # kl_div <- kl_divergence_(Sigma_1, Sigma_2)
}
