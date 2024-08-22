\name{find_columns_with_multiple_unique_values}
\alias{find_columns_with_multiple_unique_values}
\title{Find columns with multiple unique Values}
\description{
  Identifies columns in a data frame that contain more than one unique value. This function is useful for detecting columns that exhibit variability across rows.
}
\usage{
  find_columns_with_multiple_unique_values(df_)
}
\arguments{
  \item{df_}{A data frame. The data frame in which to search for columns with multiple unique values.}
}
\value{
  A character vector containing the names of the columns in the data frame that have more than one unique value.
}
\author{
  Laval Jacquin, \email{jacquin.julien@gmail.com}
}
\examples{
  # example usage
  # create a sample data frame
  df <- data.frame(
    A = c(1, 1, 1),
    B = c(1, 2, 1),
    C = c(3, 3, 3),
    D = c(4, 5, 6)
  )
  
  # find columns with multiple unique values
  multiple_value_columns <- find_columns_with_multiple_unique_values(df)
  
  # print the result
  print(multiple_value_columns) # Expected output: "B" "D"
}
