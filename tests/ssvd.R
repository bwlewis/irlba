# Tests for sparse SVD/PCA
require("irlba")
loc <- ""

test <- function()
{
  on.exit(message("Error occured in: ", loc))
  loc <<- "sparse SVD"
  set.seed(1)
  x <- matrix(rnorm(100), 10)
  s <- ssvd(x, 1, n=5)
  stopifnot(isTRUE(all.equal(sqrt(drop(crossprod(x %*% s$v - s$u %*% s$d))), 0)))

  loc <<- "sparse PCA"
  set.seed(1)
  x <- matrix(rnorm(100), 10)
  s <- ssvd(x, 1, n=5, center=TRUE)
  stopifnot(isTRUE(all.equal(sqrt(drop(crossprod(scale(x, center=TRUE, scale=FALSE) %*% s$v - s$u %*% s$d))), 0)))

  loc <<- "sparse PCA + scale"
  set.seed(1)
  x <- matrix(rnorm(100), 10)
  s <- ssvd(x, 1, n=5, center=TRUE, scale.=TRUE)
  isTRUE(all.equal(sqrt(drop(crossprod(scale(x, center=TRUE, scale=TRUE) %*% s$v - s$u %*% s$d))), 0))

  loc <<- "sparse scaled"
  set.seed(1)
  x <- matrix(rnorm(100), 10)
  s <- ssvd(x, 1, n=5, center=FALSE, scale.=TRUE)
  isTRUE(all.equal(sqrt(drop(crossprod(scale(x, center=FALSE, scale=TRUE) %*% s$v - s$u %*% s$d))), 0))

  on.exit()
}
test()
