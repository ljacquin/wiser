\name{perform_kfold_cv_wiser}
\alias{perform_kfold_cv_wiser}
\title{Perform K-fold cross-validation using WISER-estimated phenotypes}
\description{
  This function performs k-fold cross-validation (CV) to evaluate the mean squared error (MSE) of a model based on phenotypes estimated using the WISER method. 
  The WISER-estimated phenotypes used are computed beforehand using a specific combination of whitening method and regularization parameter.
  }
\usage{
  perform_kfold_cv_wiser(omic_df, raw_pheno_df, trait_,
                         fixed_effect_vars,
                         fixed_effect_vars_computed_as_factor,
                         envir_var,
                         fixed_effect_vars_computed_as_factor_by_envir,
                         random_effects_vars,
                         whitening_method, reg_method, alpha_,
                         pred_method, k_folds, wiser_cache)
}
\arguments{
  \item{omic_df}{A data frame containing the omic data (e.g. genotypes as rows and SNPs, or metabolites, or wavelength reflectances/absorbances as columns).}
  \item{raw_pheno_df}{A data frame containing the raw phenotypic data with observations across different environments and conditions.}
  \item{trait_}{A character string specifying the trait of interest to be analyzed.}
  \item{fixed_effect_vars}{A character vector listing the fixed effects variables to be included in the model.}
  \item{fixed_effect_vars_computed_as_factor}{A character vector specifying the fixed effects variables that should be treated as factors during modeling.}
  \item{envir_var}{A character string specifying the environment variable (e.g., combination of "Site", "Year", "Management" and "Block").}
  \item{fixed_effect_vars_computed_as_factor_by_envir}{A character vector listing the fixed effects variables that should be treated as factors within each environment.}
  \item{random_effects_vars}{A character string specifying the random effects variable (e.g., "Genotype").}
  \item{whitening_method}{A character string specifying the whitening method to be used in the WISER phenotype estimation (e.g., \code{"ZCA-cor"}, \code{"PCA-cor"}, or \code{"Cholesky"}).}
  \item{reg_method}{A character string specifying the regularization or shrinkage method used in conjunction with whitening (e.g., \code{"frobenius_norm_regularization"}).}
  \item{alpha_}{A numeric value for the regularization parameter, associated to the regularization or shrinkage method used for having positive definiteness.}
  \item{pred_method}{A character string specifying the prediction method to be used (e.g., \code{"rf"}, \code{"svr"}, \code{"gblup"}, \code{"rkhs"}, or \code{"lasso"}).}
  \item{k_folds}{An integer specifying the number of folds to use in cross-validation.}
  \item{wiser_cache}{A list used to store precomputed WISER phenotypes to avoid redundant computations and speed up the k-fold cross-validation process.}
}
\value{
  A numeric value representing the MSE across all folds, calculated as the MSE between predicted and observed phenotypic values in the validation sets.
}
\author{
  Laval Jacquin, \email{jacquin.julien@gmail.com}
}
\examples{
  # example usage

  # load required data
  data("refpop_raw_indiv_pheno_data_subset")
  data("refpop_geno_data_subset")

  # define parameters
  trait <- "Flowering_begin"
  fixed_effect_vars <- c("Envir", "Country", "Year", "Row", "Position", "Management")
  fixed_effect_vars_computed_as_factor <- c("Envir", "Country", "Year", "Row", "Position", "Management")
  envir_var <- "Country"
  fixed_effect_vars_computed_as_factor_by_envir <- c("Row", "Position")
  random_effects_vars <- "Genotype"
  whitening_method <- "ZCA-cor"
  reg_method <- "frobenius_norm_shrinkage"
  alpha_ <- 0.01
  pred_method <- "rf"
  k_folds <- 5
  wiser_cache <- list()

  # perform k-fold cross-validation using WISER phenotypes
  mse <- perform_kfold_cv_wiser(
    omic_df = refpop_geno_data_subset,
    raw_pheno_df = refpop_raw_indiv_pheno_data_subset,
    trait_ = trait,
    fixed_effect_vars = fixed_effect_vars,
    fixed_effect_vars_computed_as_factor = fixed_effect_vars_computed_as_factor,
    envir_var = envir_var,
    fixed_effect_vars_computed_as_factor_by_envir = fixed_effect_vars_computed_as_factor_by_envir,
    random_effects_vars = random_effects_vars,
    whitening_method = whitening_method,
    reg_method = reg_method,
    alpha_ = alpha_,
    pred_method = pred_method,
    k_folds = k_folds,
    wiser_cache = wiser_cache
  )

  # print the MSE
  print(mse)
}
