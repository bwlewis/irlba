# Tests for a few edge cases

require("irlba")
# Dense matrix
set.seed(1)
A <- matrix(rnorm(16), 4)
L <- irlba(A, nu=1, nv=1, tol=1e-9, fastpath=FALSE)
L1 <- irlba(A, nu=1, nv=1, tol=1e-9, fastpath=TRUE)
S <- svd(A, nu=1, nv=1)
if (!isTRUE(all.equal(L$d, S$d[1])))
{
  stop("Failed tiny reference example ")
}
if (!isTRUE(all.equal(L1$d, S$d[1])))
{
  stop("Failed tiny fastpath example")
}

# Tickle misc. checks
L <- irlba(A, nv=3, tol=1e-9, fastpath=FALSE, work=2, v=rep(0, nrow(A)))
set.seed(1)
A <- matrix(rnorm(100), 10)
L <- tryCatch(irlba(A, nv=3, tol=1e-9, fastpath=FALSE, work=2, v=rep(0, nrow(A))), error=function(e) "NULLSPACE")
S <- svd(A)
L <- irlba(A, nv=3, tol=1e-9, fastpath=FALSE, work=2, v=S$v[, 1])
A <- S$u %*% diag(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1e-12)) %*% t(S$v)
L <- irlba(A, nv=3, tol=1e-9, fastpath=FALSE, work=2, reorth=FALSE)

# Convergence
A <- S$u %*% (c(1e-5, rep(1e-9, 9)) * t(S$v))
for(tol in 10^-(7:12))
{
  L <- irlba(A, 3, tol=tol)
  converged <- svd(A %*% L$v - L$u  %*% diag(L$d))$d[1] < tol * L$d[1]
  stopifnot(converged)
}
