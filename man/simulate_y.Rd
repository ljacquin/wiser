\name{simulate_y}
\alias{simulate_y}
\title{Simulate phenotypic data}
\description{
  Simulates phenotypic values based on fixed effects, random effects, and residual variance components.
}
\usage{
  simulate_y(x_mat, z_mat, beta_hat, sigma2_u, sigma2_e, k_mat)
}
\arguments{
  \item{x_mat}{A numeric matrix of incidence variables for fixed effects.}
  \item{z_mat}{A numeric matrix of incidence variables for random effects.}
  \item{beta_hat}{A numeric vector of estimated fixed effects.}
  \item{sigma2_u}{A numeric value for the variance component of random effects.}
  \item{sigma2_e}{A numeric value for the variance component of residuals.}
  \item{k_mat}{A numeric matrix representing the genomic covariance matrix.}
}
\value{
  A numeric vector of simulated phenotypic values.
}
\author{
  Laval Jacquin, \email{jacquin.julien@gmail.com}
}
\examples{
  # Example usage
  # y_sim <- simulate_y(x_mat, z_mat, beta_hat, sigma2_u, sigma2_e, k_mat)
}
