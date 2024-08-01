\name{match_indices}
\alias{match_indices}
\title{Match indices of genotypes between data frames}
\description{
  Matches indices of genotypes between a raw phenotype data frame and a matrix, ensuring that only common genotypes are retained.
}
\usage{
  match_indices(raw_pheno_df, mat_)
}
\arguments{
  \item{raw_pheno_df}{A data frame of raw phenotypic data with genotypes as rows.}
  \item{mat_}{A matrix with genotypes as row names.}
}
\value{
  A data frame of raw phenotypic data with rows matched to the genotypes in the matrix.
}
\author{
  Laval Jacquin, \email{jacquin.julien@gmail.com}
}
\examples{
  # Example usage
  # matched_df <- match_indices(raw_pheno_df, mat_)
}
