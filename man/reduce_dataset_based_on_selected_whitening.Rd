\name{reduce_dataset_based_on_selected_whitening}
\alias{reduce_dataset_based_on_selected_whitening}
\title{Reduce dataset size based on selected whitening method}
\description{
  This function reduces the size of a dataset based on the specified whitening method. It selects a subset of the data based on predefined row limits associated with different whitening methods: ZCA-correlation, PCA-correlation, or Cholesky.
}
\usage{
  reduce_dataset_based_on_selected_whitening(
    whitening_method,
    raw_pheno_df,
    nrow_approx_lim_raw_dataset_zca_cor,
    nrow_approx_lim_raw_dataset_pca_cor,
    nrow_approx_lim_raw_dataset_chol
  )
}
\arguments{
  \item{whitening_method}{A character string specifying the whitening method used to determine the dataset reduction approach. Possible values are \code{"ZCA-cor"}, \code{"PCA-cor"}, or \code{"Cholesky"}.}
  \item{raw_pheno_df}{A data frame containing the raw phenotypic data to be reduced.}
  \item{nrow_approx_lim_raw_dataset_zca_cor}{An integer specifying the maximum number of rows in the raw dataset when using the ZCA-cor whitening method. Defaults to \code{20e3}. This limit is an estimate, as genotypes with few repeated phenotypic measurements (less than 30) will not be sampled. The reduction factor \code{K} is computed as \code{K = nrow(raw_pheno_df) / nrow_approx_lim_raw_dataset_zca_cor}.}
  \item{nrow_approx_lim_raw_dataset_pca_cor}{An integer specifying the maximum number of rows in the raw dataset when using the PCA-cor whitening method. Defaults to \code{20e3}. This limit is an estimate, as genotypes with few repeated phenotypic measurements (less than 30) will not be sampled. The reduction factor \code{K} is computed as \code{K = nrow(raw_pheno_df) / nrow_approx_lim_raw_dataset_pca_cor}.}
  \item{nrow_approx_lim_raw_dataset_chol}{An integer specifying the maximum number of rows in the raw dataset when using the Cholesky whitening method. Defaults to \code{40e3}. This limit is an estimate, as genotypes with few repeated phenotypic measurements (less than 30) will not be sampled. The reduction factor \code{K} is computed as \code{K = nrow(raw_pheno_df) / nrow_approx_lim_raw_dataset_chol}.}
}
\value{
  A data frame containing the reduced dataset based on the selected whitening method. 
}
\author{
  Laval Jacquin, \email{jacquin.julien@gmail.com}
}
\examples{
  # example usage
  
  # load example data
  data("example_raw_pheno_df")
  
  # define parameters
  whitening_method <- "PCA-cor"
  nrow_approx_lim_raw_dataset_zca_cor <- 1000
  nrow_approx_lim_raw_dataset_pca_cor <- 500
  nrow_approx_lim_raw_dataset_chol <- 750
  
  # reduce dataset based on the selected whitening method
  reduced_pheno_df <- reduce_dataset_based_on_selected_whitening(
    whitening_method = whitening_method,
    raw_pheno_df = example_raw_pheno_df,
    nrow_approx_lim_raw_dataset_zca_cor = nrow_approx_lim_raw_dataset_zca_cor,
    nrow_approx_lim_raw_dataset_pca_cor = nrow_approx_lim_raw_dataset_pca_cor,
    nrow_approx_lim_raw_dataset_chol = nrow_approx_lim_raw_dataset_chol
  )
  
  # print the dimensions of the reduced dataset
  print(dim(reduced_pheno_df))
}
