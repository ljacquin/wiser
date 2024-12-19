\name{estimate_wiser_phenotype}
\alias{estimate_wiser_phenotype}
\title{Estimate phenotypes approximating genetic values using whitening, successive OLS, and ABC}
\description{
  Estimates phenotypes that approximate genetic values by applying a combination of whitening, successive ordinary least squares (OLS), and approximate Bayesian computation (ABC) for variance component estimation. The function iteratively refines variance component and phenotypic value estimates through multiple rounds of ABC.
}
\usage{
  estimate_wiser_phenotype(omic_df, 
   raw_pheno_df, trait_,
   fixed_effects_vars = c(
     "Envir", "Row", "Position"
   ),
   fixed_effects_vars_computed_as_factor = c(
     "Envir", "Row", "Position"
   ),
   envir_var = "Envir",
   fixed_effects_vars_computed_as_factor_by_envir = c("Row", "Position"),
   random_effects_vars = "Genotype",
   init_sigma2_u = 1,
   init_sigma2_e = 1,
   prior_sigma2_lower_bound = 1e-2,
   n_sim_abc = 100,
   seed_abc = 123,
   quantile_threshold_abc = 0.05,
   nb_iter_abc = 1,
   kernel_type = "linear",
   whitening_method = "ZCA-cor",
   regularization_method = "frobenius_norm_regularization",
   alpha_ = 0.01,
   parallelized_cholesky = TRUE,
   reduce_raw_dataset_size_ = TRUE,
   nrow_approx_lim_raw_dataset_zca_cor = 20e3,
   nrow_approx_lim_raw_dataset_pca_cor = 20e3,
   nrow_approx_lim_raw_dataset_chol = 40e3)
}
\arguments{
  \item{omic_df}{A data frame containing the omic data (e.g. genotypes as rows and SNPs, or metabolites, or wavelength reflectances/absorbances as columns).}
  \item{raw_pheno_df}{A data frame containing the raw phenotypic data.}
  \item{trait_}{A character string specifying the trait of interest.}
  \item{fixed_effects_vars}{A character vector specifying the fixed effects variables. Defaults to \code{c("Envir", "Row", "Position")}.}
  \item{fixed_effects_vars_computed_as_factor}{A character vector specifying which fixed effects variables should be treated as factors. Defaults to \code{c("Envir", "Row", "Position")}. Note: declared fixed effects variables which are not specified to be computed as factors will be treated as quantitative variables.}
  \item{envir_var}{A character string specifying the environment variable. Defaults to \code{"Envir"} and should be set to \code{"NULL"} if there is no environment variable.}
  \item{fixed_effects_vars_computed_as_factor_by_envir}{A character vector specifying which fixed effects variables should be treated as factors by environment. Defaults to \code{c("Row", "Position")}.}
  \item{random_effects_vars}{A character vector specifying the random effects variables. Defaults to \code{"Genotype"}.}
  \item{init_sigma2_u}{A numeric value for the initial variance component of random effects. Defaults to \code{1}.}
  \item{init_sigma2_e}{A numeric value for the initial variance component of residuals. Defaults to \code{1}.}
  \item{prior_sigma2_lower_bound}{A numeric value specifying the lower bound for the prior distribution of the variance components (\eqn{\sigma_u^2} and \eqn{\sigma_e^2}) in ABC estimation. Defaults to \code{1e-2}.}

  \item{n_sim_abc}{An integer specifying the number of simulations for ABC. Defaults to \code{100}.}
  \item{seed_abc}{An integer to set the seed for random number generation. Defaults to \code{123}.}
  \item{quantile_threshold_abc}{A numeric value for the quantile threshold in ABC. Defaults to \code{0.05}.}
  \item{nb_iter_abc}{An integer specifying the number of iterations for ABC. Defaults to \code{1}.}
  \item{kernel_type}{A character string specifying the kernel type. Can be \code{"linear"} or \code{"identity"}. Defaults to \code{"linear"}.}
  \item{whitening_method}{A character string specifying the whitening method. Can be \code{"ZCA-cor"}, \code{"PCA-cor"}, or \code{"Cholesky"}. Defaults to \code{"ZCA-cor"}.}
  \item{regularization_method}{A character string specifying the regularization or shrinkage method used for ensuring positive definiteness of the genetic covariance matrix, implicated in the whitening process. Can be \code{"frobenius_norm_regularization"}, \code{"frobenius_norm_shrinkage"}, \code{"frobenius_norm_partial_shrinkage"}, \code{"trace_sample_variance_regularization"}, \code{"trace_sample_variance_shrinkage"}, or \code{"trace_sample_variance_partial_shrinkage"}. Defaults to \code{"frobenius_norm_regularization"}.}
  \item{alpha_}{A numeric value for the regularization parameter associated with the regularization or shrinkage method. Defaults to \code{0.01}.}
  \item{parallelized_cholesky}{A logical value indicating whether to use parallelized Cholesky decomposition. Defaults to \code{TRUE}.}
  \item{reduce_raw_dataset_size_}{A logical value indicating whether to reduce the size of the raw dataset for performance optimization. Defaults to \code{TRUE}.}
  \item{nrow_approx_lim_raw_dataset_zca_cor}{An integer specifying the maximum number of rows in the raw dataset when using the ZCA-cor whitening method. Defaults to \code{20e3}. This limit is an estimate, as genotypes with few repeated phenotypic measurements (less than 30) will not be sampled. The reduction factor \code{K} used is computed as \code{K = nrow(raw_pheno_df) / nrow_approx_lim_raw_dataset_zca_cor}.}
  \item{nrow_approx_lim_raw_dataset_pca_cor}{An integer specifying the maximum number of rows in the raw dataset when using the PCA-cor whitening method. Defaults to \code{20e3}. This limit is an estimate, as genotypes with few repeated phenotypic measurements (less than 30) will not be sampled. The reduction factor \code{K} used is computed as \code{K = nrow(raw_pheno_df) / nrow_approx_lim_raw_dataset_pca_cor}.}
  \item{nrow_approx_lim_raw_dataset_chol}{An integer specifying the maximum number of rows in the raw dataset when using the Cholesky whitening method. Defaults to \code{40e3}. This limit is an estimate, as genotypes with few repeated phenotypic measurements (less than 30) will not be sampled. The reduction factor \code{K} used is computed as \code{K = nrow(raw_pheno_df) / nrow_approx_lim_raw_dataset_chol}.}
}
\value{
  A list containing:
  \item{wiser_omic_data}{The omic data used in the estimation process.}
  \item{sig_mat_u}{The computed genetic covariance matrix associated to individuals (i.e. \eqn{\Sigma_u}).}
  \item{w_mat}{The transformation matrix from the whitening process.}
  \item{wiser_fixed_effect_estimates}{A data frame of the estimated fixed effects.}
  \item{wiser_abc_variance_component_estimates}{The results of the ABC variance component estimation.}
  \item{wiser_phenotypes}{A data frame containing the estimated phenotypic values (i.e. \code{v_hat}), for the genotypes, from the second OLS regression.}
  \item{wiser_z_mat}{The design matrix for the genetic effects.}
  \item{wiser_x_mat}{The design matrix for the fixed effects.}
  \item{wiser_x_mat_tilde}{The transformed design matrix after whitening.}
  \item{wiser_xi_hat}{The estimated residuals from the first OLS regression.}
  \item{wiser_y_hat}{The predicted values from the first OLS regression.}
  \item{wiser_y}{The raw observed phenotype used by wiser.}
}
\author{
  Laval Jacquin, \email{jacquin.julien@gmail.com}
}
