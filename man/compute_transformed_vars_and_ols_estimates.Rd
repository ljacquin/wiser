\name{compute_transformed_vars_and_ols_estimates}
\alias{compute_transformed_vars_and_ols_estimates}
\title{Compute transformed variables and OLS estimates for fixed effects}
\description{
  Computes transformed variables for the design matrix of fixed effects, and ordinary least squares (OLS) estimates for these effects.
}
\usage{
compute_transformed_vars_and_ols_estimates(
  geno_df, raw_pheno_df, fixed_effects_vars, random_effects_vars, trait_,
  compute_row_and_position_as_factors, sigma2_u, sigma2_e, kernel_type,
  rate_decay_kernel, whitening_method, regularization_method, alpha_frob_,
  percent_eig_, non_zero_precision_eig_
)
}
\arguments{
  \item{geno_df}{A data frame of genomic data with genotypes as rows and SNPs as columns.}
  \item{raw_pheno_df}{A data frame of raw phenotypic data.}
  \item{fixed_effects_vars}{A character vector of fixed effects variables.}
  \item{random_effects_vars}{A character vector of random effects variables.}
  \item{trait_}{A character string specifying the trait of interest.}
  \item{compute_row_and_position_as_factors}{A logical value indicating whether to treat "Row" and "Position" as factors.}
  \item{sigma2_u}{A numeric value for the variance component of random effects.}
  \item{sigma2_e}{A numeric value for the variance component of residuals.}
  \item{kernel_type}{A character string specifying the kernel type ("linear", "gaussian", or "identity").}
  \item{rate_decay_kernel}{A numeric value for the rate decay parameter used in the Gaussian kernel.}
  \item{whitening_method}{A character string specifying the whitening method ("Cholesky" or "ZCA-cor").}
  \item{regularization_method}{A character string specifying the regularization method ("mean_small_eigenvalues", "mean_eigenvalues", or "frobenius_norm").}
  \item{alpha_frob_}{A numeric value for the Frobenius norm regularization parameter.}
  \item{percent_eig_}{A numeric value between 0 and 1 representing the percentage of the smallest strictly positive eigenvalues to be considered.}
  \item{non_zero_precision_eig_}{A small positive value used to filter out near-zero eigenvalues.}
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
\author{
  Laval Jacquin, \email{jacquin.julien@gmail.com}
}
\examples{
# Example usage
if(requireNamespace("mixOmics", quietly = TRUE)) {
# -- load wiser library and data subsets from the refpop dataset
# load library
library(wiser)

# load raw individual phenotypic measurements (for subset of 80 genotypes from refpop data)
data("refpop_raw_indiv_pheno_data_subset")
head(refpop_raw_indiv_pheno_data_subset)

# load SNP marker data (for subset of 80 genotypes from refpop data)
data("refpop_geno_data_subset")
head(refpop_geno_data_subset)[,1:10]

# -- estimate wiser phenotype for a specified trait,i.e. Flowering_begin here
# define trait
trait_ <- "Flowering_begin"

result <- compute_transformed_vars_and_ols_estimates(
  refpop_geno_data_subset, refpop_raw_indiv_pheno_data_subset,
  fixed_effects_vars,
  random_effects_vars,
  trait_,
  compute_row_and_position_as_factors = T,
  sigma2_u = 1, sigma2_e = 1,
  kernel_type = 'linear',
  rate_decay_kernel = 0.1,
  whitening_method = 'Cholesky',
  regularization_method = 'frobenius_norm',
  alpha_frob_ = 0.01,
  percent_eig_ = 0.05,
  non_zero_precision_eig_ = 1e-5
  )
  print(result)
} else {
  message("wiser package is required for this example.")
 }
}
