\name{match_indices}
\alias{match_indices}
\title{Match indices of genotypes between data frames}
\description{
  Matches indices of genotypes between a raw phenotype data frame and an omic data frame, ensuring that only common genotypes are retained.
}
\usage{
  match_indices(raw_pheno_df, omic_df)
}
\arguments{
  \item{raw_pheno_df}{A data frame of raw phenotypic data containing a column named "Genotype".}
  \item{omic_df}{An omic data frame with genotypes as row names.}
}
\value{
  A data frame of raw phenotypic data with rows matched to the genotypes in the omic data frame.
}
\author{
  Laval Jacquin, \email{jacquin.julien@gmail.com}
}
\examples{
  # example usage
  # matched_df <- match_indices(raw_pheno_df, omic_df)
}
