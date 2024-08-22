\name{perform_kfold_cv_wiser}
\alias{perform_kfold_cv_wiser}
\title{Perform K-fold cross-validation using WISER-estimated phenotypes}
\description{
  This function performs k-fold cross-validation (CV) to evaluate the predictive ability of a model for phenotypes estimated using the WISER method. The function computes WISER phenotypes for a specific whitening method and regularization parameter, and then applies cross-validation to assess model performance.
}
\usage{
  perform_kfold_cv_wiser(geno_df, raw_pheno_df, trait_,
                         whitening_method, reg_method, alpha_frob,
                         pred_method, k_folds)
}
\arguments{
  \item{geno_df}{A data frame containing the genomic data, where rows represent genotypes and columns represent SNPs.}
  \item{raw_pheno_df}{A data frame containing the raw phenotypic data.}
  \item{trait_}{A character string specifying the trait of interest to be analyzed.}
  \item{whitening_method}{A character string specifying the whitening method to be used in the WISER phenotype estimation (i.e. \code{"ZCA-cor"}, \code{"PCA-cor"} or \code{"Cholesky"}).}
  \item{reg_method}{A character string specifying the regularization method used in conjunction with whitening (i.e. \code{"frobenius_norm"}, \code{"mean_small_eigenvalues"} or \code{"mean_eigenvalues"}).}
  \item{alpha_frob}{A numeric value specifying the regularization parameter for the Frobenius norm.}
  \item{pred_method}{A character string specifying the prediction method to be used (i.e. \code{"rf"}, \code{"svr"}, \code{"gblup"}, \code{"rkhs"} or \code{"lasso"}).}
  \item{k_folds}{An integer specifying the number of folds to use in cross-validation.}
}
\value{
  A numeric value representing the mean predictive ability (PA) across all folds, calculated as the mean Pearson correlation coefficient between predicted and observed phenotypic values in the validation sets.
}
\author{
  Laval Jacquin, \email{jacquin.julien@gmail.com}
}
\examples{
  # example usage
  # load library
  library(wiser)

  # load required data
  data("refpop_raw_indiv_pheno_data_subset")
  data("refpop_geno_data_subset")

  # define parameters
  trait <- "Flowering_begin"
  whitening_method <- "ZCA-cor"
  reg_method <- "frobenius_norm"
  alpha_frob <- 0.01
  pred_method <- "rf"
  k_folds <- 5

  # perform k-fold cross-validation using WISER phenotypes
  mean_pa <- perform_kfold_cv_wiser(
    geno_df = refpop_geno_data_subset,
    raw_pheno_df = refpop_raw_indiv_pheno_data_subset,
    trait_ = trait,
    whitening_method = whitening_method,
    reg_method = reg_method,
    alpha_frob = alpha_frob,
    pred_method = pred_method,
    k_folds = k_folds
    )

  # print the mean predictive ability (note the dataset used is quite small and
  # is for example usage purposes only, descent result is not expected)
  print(mean_pa)
}
