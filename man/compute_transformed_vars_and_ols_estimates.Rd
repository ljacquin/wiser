\name{compute_transformed_vars_and_ols_estimates}
\alias{compute_transformed_vars_and_ols_estimates}
\title{Compute transformed variables and OLS estimates for fixed effects}
\description{
  Computes transformed variables for the design matrix of fixed effects, and ordinary least squares (OLS) estimates for these effects.
}
\usage{
  compute_transformed_vars_and_ols_estimates(
    geno_df, raw_pheno_df, fixed_effects_vars, random_effects_vars, trait_,
    sigma2_u, sigma2_e, kernel_type, rate_decay_kernel, whitening_method
  )
}
\arguments{
  \item{geno_df}{A data frame of genomic data with genotypes as rows and SNPs as columns.}
  \item{raw_pheno_df}{A data frame of raw phenotypic data.}
  \item{fixed_effects_vars}{A character vector of fixed effects variables.}
  \item{random_effects_vars}{A character vector of random effects variables.}
  \item{trait_}{A character string specifying the trait of interest.}
  \item{sigma2_u}{A numeric value for the variance component of random effects.}
  \item{sigma2_e}{A numeric value for the variance component of residuals.}
  \item{kernel_type}{A character string specifying the kernel type ("linear", "gaussian", or "identity").}
  \item{rate_decay_kernel}{A numeric value for the rate decay parameter used in the Gaussian kernel.}
  \item{whitening_method}{A character string specifying the whitening method ("Cholesky" or "ZCA-cor").}
}
\value{
  A list containing:
    \item{x_mat}{The transformed incidence matrix of fixed effects variables.}
  \item{z_mat}{The incidence matrix of random effects variables.}
  \item{k_mat}{The genomic covariance matrix.}
  \item{beta_hat}{The estimated fixed effects (adjusted for the genetic covariance structure between individuals in the experimental design).}
  \item{xi_hat}{The estimated residuals.}
  \item{y}{The vector of raw phenotypic measurements for the trait.}
}
\examples{
  # Example usage
  # result <- compute_transformed_vars_and_ols_estimates(
  #   geno_df, raw_pheno_df, fixed_effects_vars = c("Envir", "Country"),
  #   random_effects_vars = "Genotype", trait_ = "Trait1",
  #   sigma2_u = 1, sigma2_e = 1, kernel_type = "linear",
  #   rate_decay_kernel = 0.1, whitening_method = "Cholesky"
  # )
}
