\name{generate_row_column_variables_by_environment}
\alias{generate_row_column_variables_by_environment}
\title{Generate row and column variables by environment}
\description{
  This function generates new columns for rows and columns in a dataset, grouped by unique combinations of environment-related factors such as site, year, management, and block. The generated columns indicate the row and column values for each specific environment. The rows and columns are treated as quantitative variables.
}
\usage{
  generate_row_column_variables_by_environment(df)
}
\arguments{
  \item{df}{A data frame containing the input data with at least the following columns: \code{Site}, \code{year}, \code{Management}, \code{block}, \code{Row}, and \code{Column}.}
}
\value{
  A data frame with additional columns for each environment. Each new column contains the corresponding row or column value for the specific environment, with 0 elsewhere.
}
\details{
  The function iterates through all unique combinations (i.e., environments) of the \code{Site}, \code{year}, \code{Management}, and \code{block} variables. For each environment, two new columns are added to the data frame:
    \itemize{
      \item{A row variable column, named as \code{"<site>_<year>_<management>_block_<block_num>_row"}}
      \item{A column variable column, named as \code{"<site>_<year>_<management>_block_<block_num>_column"}}
    }
  These columns store the row and column values for observations matching the environment, with all other rows set to 0.
}
\author{
  Laval Jacquin, \email{jacquin.julien@gmail.com}
}
\examples{
  # example usage
  # define a sample dataset
  df <- data.frame(
    Site = c("A", "A", "B", "B"),
    year = c(2020, 2020, 2021, 2021),
    Management = c("Conv", "Conv", "Org", "Org"),
    block = c(1, 2, 1, 2),
    Row = c(10, 20, 30, 40),
    Column = c(1, 2, 3, 4)
  )
  
  # generate row and column variables by environment
  df_updated <- generate_row_column_variables_by_environment(df)
  
  # view the updated dataset
  print(df_updated)
}
