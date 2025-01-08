\name{compute_transformed_vars_and_ols_estimates}
\alias{compute_transformed_vars_and_ols_estimates}
\title{Compute transformed variables and OLS estimates for fixed effects}
\description{
  This function computes transformed variables for the design matrix of fixed effects, based on whitening, and ordinary least squares (OLS) estimates for these effects. It also computes the genetic covariance matrix (Gram matrix) using omic data (e.g., genomic, metabolomic, or phenomic data) and derives the associated whitening matrix according to the specified whitening method.
}
\usage{
compute_transformed_vars_and_ols_estimates(
  omic_df, raw_pheno_df, trait_,
  fixed_effect_vars, fixed_effect_vars_computed_as_factor,
  envir_var, fixed_effect_vars_computed_as_factor_by_envir,
  random_effects_vars, sigma2_u, sigma2_e, kernel_type,
  whitening_method, regularization_method, alpha_,
  parallelized_cholesky, reduce_raw_dataset_size_,
  nrow_approx_lim_raw_dataset_zca_cor,
  nrow_approx_lim_raw_dataset_pca_cor,
  nrow_approx_lim_raw_dataset_chol
)
}
\arguments{
  \item{omic_df}{A data frame containing the omic data (e.g., genotypes as rows and SNPs, or metabolites, or wavelength reflectances/absorbances as columns).}
  \item{raw_pheno_df}{A data frame containing the raw phenotypic data.}
  \item{trait_}{A character string specifying the trait of interest.}
  \item{fixed_effect_vars}{A character vector specifying the fixed effects variables.}
  \item{fixed_effect_vars_computed_as_factor}{A character vector specifying which fixed effects variables should be treated as factors.}
  \item{envir_var}{A character string specifying the environment variable.}
  \item{fixed_effect_vars_computed_as_factor_by_envir}{A character vector specifying which fixed effects variables should be treated as factors by environment.}
  \item{random_effects_vars}{A character vector specifying the random effects variables.}
  \item{sigma2_u}{A numeric value specifying the variance component for random effects.}
  \item{sigma2_e}{A numeric value specifying the variance component for residuals.}
  \item{kernel_type}{A character string specifying the type of kernel used in computing the Gram matrix. Can be \code{"linear"} or \code{"identity"}.}
  \item{whitening_method}{A character string specifying the whitening method. Can be \code{"ZCA-cor"}, \code{"PCA-cor"}, or \code{"Cholesky"}.}
  \item{regularization_method}{A character string specifying the regularization or shrinkage method used to ensure positive definiteness of the genetic covariance matrix in the whitening process. Options include \code{"frobenius_norm_regularization"}, \code{"frobenius_norm_shrinkage"}, \code{"frobenius_norm_partial_shrinkage"}, \code{"trace_sample_variance_regularization"}, \code{"trace_sample_variance_shrinkage"}, or \code{"trace_sample_variance_partial_shrinkage"}.}
  \item{alpha_}{A numeric value for the regularization parameter associated with the regularization or shrinkage method. Defaults to \code{0.01}.}
  \item{parallelized_cholesky}{A logical value indicating whether to use parallelized Cholesky decomposition. Defaults to \code{TRUE}.}
  \item{reduce_raw_dataset_size_}{A logical value indicating whether to reduce the size of the raw dataset based on the selected whitening method. Defaults to \code{TRUE}.}
  \item{nrow_approx_lim_raw_dataset_zca_cor}{An integer specifying the maximum number of rows in the raw dataset when using the ZCA-cor whitening method. Defaults to \code{20e3}. This limit is an estimate, as genotypes with few repeated phenotypic measurements (less than 30) will not be sampled. The reduction factor \code{K} used is computed as \code{K = nrow(raw_pheno_df) / nrow_approx_lim_raw_dataset_zca_cor}.}
  \item{nrow_approx_lim_raw_dataset_pca_cor}{An integer specifying the maximum number of rows in the raw dataset when using the PCA-cor whitening method. Defaults to \code{20e3}. This limit is an estimate, as genotypes with few repeated phenotypic measurements (less than 30) will not be sampled. The reduction factor \code{K} used is computed as \code{K = nrow(raw_pheno_df) / nrow_approx_lim_raw_dataset_pca_cor}.}
  \item{nrow_approx_lim_raw_dataset_chol}{An integer specifying the maximum number of rows in the raw dataset when using the Cholesky whitening method. Defaults to \code{40e3}. This limit is an estimate, as genotypes with few repeated phenotypic measurements (less than 30) will not be sampled. The reduction factor \code{K} used is computed as \code{K = nrow(raw_pheno_df) / nrow_approx_lim_raw_dataset_chol}.}
}
\value{
  A list containing the following elements:
  \item{omic_df}{The omic data used in the estimation process.}
  \item{sig_mat_u}{The computed genetic covariance matrix for individuals (i.e., \eqn{\Sigma_u}).}
  \item{w_mat}{The whitening matrix derived from \code{sig_mat_u}.}
  \item{x_mat}{The incidence matrix of fixed effects variables.}
  \item{x_mat_tilde}{The transformed incidence matrix of fixed effects variables after applying the whitening matrix.}
  \item{z_mat}{The incidence matrix of random effects variables (genetic effects specifically).}
  \item{k_mat}{The genetic covariance matrix (Gram matrix) computed using the omic data.}
  \item{beta_hat}{The estimated coefficients for the fixed effects from the OLS regression.}
  \item{y_hat}{The predicted values from the OLS regression.}
  \item{xi_hat}{The residuals (difference between observed and predicted values).}
  \item{y}{The vector of raw phenotypic measurements for the specified trait.}
}
\author{
  Laval Jacquin, \email{jacquin.julien@gmail.com}
}
\examples{
# example usage
if(requireNamespace("wiser", quietly = TRUE)) {
  # load the wiser library and example data subsets from the refpop dataset
  library(wiser)

  # load raw individual phenotypic measurements (subset of 80 genotypes from refpop data)
  data("refpop_raw_indiv_pheno_data_subset")
  head(refpop_raw_indiv_pheno_data_subset)

  # load SNP marker data (subset of 80 genotypes from refpop data)
  data("refpop_geno_data_subset")
  head(refpop_geno_data_subset)[,1:10]

  # define the trait of interest
  trait_ <- "Flowering_begin"

  # compute transformed variables and OLS estimates
  result <- compute_transformed_vars_and_ols_estimates(
    omic_df = refpop_geno_data_subset,
    raw_pheno_df = refpop_raw_indiv_pheno_data_subset,
    trait_ = trait_,
    fixed_effect_vars = c("Envir", "Country", "Year", "Row", "Position", "Management"),
    fixed_effect_vars_computed_as_factor = c("Envir", "Country", "Year"),
    envir_var = "Country",
    fixed_effect_vars_computed_as_factor_by_envir = c("Row", "Position"),
    random_effects_vars = "Genotype",
    sigma2_u = 1, 
    sigma2_e = 1,
    kernel_type = "linear",
    whitening_method = "ZCA-cor",
    regularization_method = "frobenius_norm_shrinkage",
    alpha_ = 0.01,
    parallelized_cholesky = TRUE,
    reduce_raw_dataset_size_ = TRUE,
    nrow_approx_lim_raw_dataset_zca_cor = 5000,
    nrow_approx_lim_raw_dataset_pca_cor = 5000,
    nrow_approx_lim_raw_dataset_chol = 10000
  )

  # print the result
  print(result)
} else {
  message("The wiser package is required for this example.")
}
}
