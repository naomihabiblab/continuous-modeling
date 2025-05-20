
is_square_matrix <- function(m) {
  if (dim(m)[1] == dim(m)[2]) {
    TRUE
  }
  else {
    FALSE
  }
}

sample_rows <- function (m, n) {
  number_of_rows <- nrow(m)
  n <- min(n, number_of_rows)
  m[sample(number_of_rows, n),]
}
