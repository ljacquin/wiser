\name{compute_vect_mode}
\alias{compute_vect_mode}
\title{Compute the statistical mode of a vector}
\description{
  This function calculates the statistical mode of a given vector. The mode is the value that appears most frequently in the vector. If multiple values have the same highest frequency, the function returns the first one encountered.
}
\usage{
compute_vect_mode(vect)
}
\arguments{
  \item{vect}{A vector of any type (e.g., numeric, character, factor) for which the mode is to be computed.}
}
\value{
  A single value representing the mode of the input vector. If there are multiple modes, only the first one (in order of appearance) is returned.
}
\author{
  Laval Jacquin, \email{jacquin.julien@gmail.com}
}
\examples{
  # Example usage with a numeric vector
  numeric_vect <- c(1, 2, 2, 3, 4, 4, 4, 5)
  mode_numeric <- compute_vect_mode(numeric_vect)
  print(mode_numeric) # Should return 4

  # Example usage with a character vector
  char_vect <- c("apple", "banana", "apple", "orange", "banana", "banana")
  mode_char <- compute_vect_mode(char_vect)
  print(mode_char) # Should return "banana"

  # Example usage with a vector having multiple modes
  multi_mode_vect <- c("a", "b", "b", "c", "a")
  mode_multi <- compute_vect_mode(multi_mode_vect)
  print(mode_multi) # Should return "a" or "b", depending on order of appearance
}
