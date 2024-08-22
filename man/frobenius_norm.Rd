\name{frobenius_norm}
\alias{frobenius_norm}
\title{Compute the Frobenius norm of a matrix}
\description{
  Computes the Frobenius norm of a given matrix. The Frobenius norm is a measure of the size of a matrix, defined as the square root of the sum of the absolute squares of its entries.
}
\usage{
  frobenius_norm(matrix_)
}
\arguments{
  \item{matrix_}{A numeric matrix for which the Frobenius norm is to be computed.}
}
\value{
  A numeric value representing the Frobenius norm of the input matrix.
}
\details{
  The Frobenius norm of a matrix \code{A} is defined as:

  \deqn{\| \mathbf{A} \|_F = \sqrt{\sum_{i,j} a_{ij}^2}}{||A||_F = sqrt(sum_{i,j} a_{ij}^2)}

  where \code{a_{ij}} denotes the entry in the \code{i}-th row and \code{j}-th column of \code{A}. It is also equivalent to the square root of the trace of the matrix \code{A A^T}, where \code{A^T} is the transpose of \code{A}. This can be demonstrated as follows:

  \deqn{\| \mathbf{A} \|_F^2 = \text{trace}(\mathbf{A} \mathbf{A}^T)}{||A||_F^2 = trace(A A^T)}

  **Proof:**

  1. Expand \code{trace(A A^T)}:

     \deqn{\text{trace}(A A^T) = \sum_{i} (A A^T)_{ii}}{trace(A A^T) = sum_{i} (A A^T)_{ii}}

  2. Compute the (i,i)-th entry of \code{A A^T}:

     \deqn{(A A^T)_{ii} = \sum_{j} a_{ij}^2}{(A A^T)_{ii} = sum_{j} a_{ij}^2}

  3. Therefore,

     \deqn{\text{trace}(A A^T) = \sum_{i} \sum_{j} a_{ij}^2}{trace(A A^T) = sum_{i} sum_{j} a_{ij}^2}

  4. Taking the square root gives:

     \deqn{\| \mathbf{A} \|_F = \sqrt{\text{trace}(A A^T)}}{||A||_F = sqrt(trace(A A^T))}
}
\author{
  Laval Jacquin, \email{jacquin.julien@gmail.com}
}
\examples{
  # example usage
  mat <- matrix(c(1, 2, 3, 4), nrow = 2)
  frobenius_norm(mat)
}
