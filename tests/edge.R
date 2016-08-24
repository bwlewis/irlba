# Tests for a few edge cases

require("irlba")
# Dense matrix
set.seed(1)
A <- matrix(rnorm(16),4)
L <- irlba(A, nu=1, nv=1, tol=1e-9, fastpath=FALSE)
L1 <- irlba(A, nu=1, nv=1, tol=1e-9, fastpath=TRUE)
S <- svd(A, nu=1, nv=1)
if(!isTRUE(all.equal(L$d, S$d[1])))
{
  stop("Failed tiny reference example ")
}
if(!isTRUE(all.equal(L1$d, S$d[1])))
{
  stop("Failed tiny fastpath example")
}
