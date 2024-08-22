\name{squared_l2_norm}
\alias{squared_l2_norm}
\title{Compute squared L2 norm distance}
\description{
  Computes the squared L2 norm distance between two numeric vectors.
}
\usage{
  squared_l2_norm(y, y_sim)
}
\arguments{
  \item{y}{A numeric vector of observed phenotypic measurements.}
  \item{y_sim}{A numeric vector of simulated phenotypic measurements.}
}
\value{
  A numeric value representing the squared L2 norm distance between the observed and simulated vectors.
}
\author{
  Laval Jacquin, \email{jacquin.julien@gmail.com}
}
\examples{
  # example usage
  # distance <- squared_l2_norm(y, y_sim)
}
