\name{generate_latitude_longitude_variables_by_environment}
\alias{generate_latitude_longitude_variables_by_environment}
\title{Generate latitude and longitude variables by environment}
\description{
  This function creates new columns for latitude and longitude in a dataset, grouped by unique combinations of environment-related factors such as site, year, and block. The new columns contain the latitude and longitude values for each specific environment, with zero values elsewhere.
  The latitudes and longitudes are treated as quantitative variables.
  }
\usage{
  generate_latitude_longitude_variables_by_environment(df)
}
\arguments{
  \item{df}{A data frame containing the input data with at least the following columns: \code{Site}, \code{Year}, \code{Block}, \code{Latitude}, and \code{Longitude}.}
}
\value{
  A data frame with additional columns for each environment. Each new column contains the corresponding latitude or longitude value for the specific environment, with 0 elsewhere.
}
\details{
  The function iterates through all unique combinations (i.e., environments) of the \code{Site}, \code{Year}, and \code{Block} variables. For each environment, two new columns are added to the data frame:
    \itemize{
      \item{A latitude variable column, named as \code{"<site>_<year>_block_<block_num>_latitude"}}
      \item{A longitude variable column, named as \code{"<site>_<year>_block_<block_num>_longitude"}}
    }
  These columns store the latitude and longitude values for observations matching the environment, with all other rows set to 0.
}
\author{
  Laval Jacquin, \email{jacquin.julien@gmail.com}
}
\examples{
  # example usage
  
  # define a sample dataset
  df <- data.frame(
    Site = c("A", "A", "B", "B"),
    Year = c(2020, 2020, 2021, 2021),
    Block = c(1, 2, 1, 2),
    Latitude = c(45.5, 45.6, 46.0, 46.1),
    Longitude = c(-73.6, -73.7, -74.0, -74.1)
  )
  
  # generate latitude and longitude variables by environment
  df_updated <- generate_latitude_longitude_variables_by_environment(df)
  
  # view the updated dataset
  print(df_updated)
}
