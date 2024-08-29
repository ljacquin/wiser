\name{remove_monomorphic_markers}
\alias{remove_monomorphic_markers}
\title{Remove monomorphic markers from genotypic data}
\description{
  This function removes monomorphic markers (i.e., markers with no variation) from a genotypic data frame. It also provides a list of the removed monomorphic markers.
}
\usage{
  remove_monomorphic_markers(geno_df)
}
\arguments{
  \item{geno_df}{A data frame containing genotypic data, where columns represent markers and rows represent genotypes. It may optionally include a column named \code{"Genotype"} to identify the genotypes.}
}
\value{
  A list with two elements:
    \itemize{
      \item{\code{"filtered_df"}}{A data frame with monomorphic markers removed. If the original data frame contained a \code{"Genotype"} column, it is preserved as the first column.}
      \item{\code{"monomorphic_markers"}}{A character vector containing the names of the monomorphic markers that were removed. If no monomorphic markers were found, this is \code{NULL}.}
    }
}
\details{
  The function first checks if the \code{"Genotype"} column is present in the input data frame. If so, it removes this column temporarily for the purpose of identifying monomorphic markers. It then identifies markers where all values are identical (i.e., monomorphic) and removes these markers from the data frame. If the \code{"Genotype"} column was originally present, it is added back to the filtered data frame.
}
\author{
  Laval Jacquin, \email{jacquin.julien@gmail.com}
}
\examples{
  # example usage
  
  # define a sample data frame with genotypic data
  geno_df <- data.frame(
    Genotype = c("G1", "G2", "G3"),
    Marker1 = c(1, 1, 1),
    Marker2 = c(0, 1, 0),
    Marker3 = c(2, 2, 2)
  )
  
  # demove monomorphic markers
  result <- remove_monomorphic_markers(geno_df)
  
  # print the filtered data frame and monomorphic markers
  print(result$filtered_df)
  print(result$monomorphic_markers)
}
