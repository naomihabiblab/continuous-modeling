#' @title Matrix Utility Functions
#' @description Simple utility functions for matrix operations.

#' Check if matrix is square
#'
#' @param m Matrix to check
#'
#' @return Logical indicating if dimensions are equal
is_square_matrix <- function(m) {
  if (dim(m)[1] == dim(m)[2]) {
    TRUE
  }
  else {
    FALSE
  }
}

#' Sample rows from a matrix
#'
#' @param m Matrix to sample from
#' @param n Number of rows to sample (capped at nrow(m))
#'
#' @return Matrix with n randomly sampled rows
sample_rows <- function (m, n) {
  number_of_rows <- nrow(m)
  n <- min(n, number_of_rows)
  m[sample(number_of_rows, n),]
}
