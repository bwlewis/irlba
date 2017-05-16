# Tests for svdr
require("irlba")
loc <- ""

test <- function()
{
  on.exit(message("Error occured in: ", loc))
  # Dense matrix
  loc <<- "svdr dense"
  set.seed(1)
  A <- matrix(rnorm(16), 4)
  L <- svdr(A, 1)
  S <- svd(A, nu=1, nv=1)
  stopifnot(isTRUE(all.equal(L$d, S$d[1])))

  loc <<- "svdr dense m > n"
  A <- matrix(rnorm(50 * 40), 50)
  L <- svdr(A, 5, 10, extra=15)
  S <- svd(A, nu=5, nv=5)
  stopifnot(isTRUE(all.equal(L$d, S$d[1:5])))

  loc <<- "svdr dense m < n"
  A <- matrix(rnorm(50 * 40), 40)
  L <- svdr(A, 5, 10, extra=15)
  S <- svd(A, nu=5, nv=5)
  stopifnot(isTRUE(all.equal(L$d, S$d[1:5])))

  # Sparse but not dgCMatrix (issue #6)
  loc <<- "svdr misc sparse"
  A <- Matrix(matrix(rnorm(100), 10))
  L <- svdr(A, 1)
  S <- svd(A, nu=1, nv=1)
  stopifnot(isTRUE(all.equal(L$d, S$d[1])))

  loc <<- "svdr logical sparse"
  A <- Matrix(sample(c(FALSE, TRUE), 100, replace=TRUE), 10, 10)
  L <- svdr(A, 1)
  S <- svd(A, nu=1, nv=1)
  stopifnot(isTRUE(all.equal(L$d, S$d[1])))

  loc <<- "svdr center only, sparse"
  A <- Matrix(matrix(rnorm(100), 10))
  m <- colMeans(A)
  L <- svdr(A, 3, center=m)
  S <- svd(scale(A, center=TRUE, scale=FALSE))
  stopifnot(isTRUE(all.equal(L$d, S$d[1:3])))

  on.exit()
}
test()
