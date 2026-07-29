# Tests for some edge cases
require("irlba")


testex <- function(message, expr)
{
  cat(message, "\t", file=stderr())
  tryCatch(expr, warning=function(w) cat("Warning: ", conditionMessage(w), "\t", file=stderr()))
  cat("OK\n", file=stderr()) 
}

set.seed(1)

testex("Convergence", {
A <- matrix(rnorm(100), 10)
L <- tryCatch(irlba(A, nv=3, tol=1e-9, work=2, v=rep(0, nrow(A))), error=function(e) "NULLSPACE")
S <- svd(A)
A <- S$u %*% (c(1e-5, rep(1e-9, 9)) * t(S$v))
for (tol in 10 ^ - (7:11))
  { 
    L <- irlba(A, 3, tol=tol, svtol=Inf)
    converged <- svd(A %*% L$v - L$u  %*% diag(L$d))$d[1] < tol * L$d[1]
    stopifnot(converged)
  }
})

testex("dgeMatrix", {
  X <- Matrix(matrix(rnorm(100), 10))
  L <- irlba(X, nv=1)
  S <- svd(X, nu=1, nv=1)
  stopifnot(all.equal(L$d, S$d[1]))
})

testex("lgCMatrix", {
  X <- Matrix(sample(c(FALSE, TRUE), 100, replace=TRUE), 10, 10)
  L <- irlba(X, nv=1)
  S <- svd(X, nu=1, nv=1)
  stopifnot(all.equal(L$d, S$d[1]))
})

testex("GitHub issue #7", {
  mx <- matrix(sample(1:10, 10 * 100, replace=TRUE), nrow=10)
  S <- withCallingHandlers(
    irlba(mx, nv = 2, center = colMeans(mx), right_only = TRUE),
    warning = function(w) {
      cat("caught atomic type coercion warning\t", file=stderr())
      invokeRestart("muffleWarning")
    }
  )
  stopifnot(all.equal(S$d, svd(sweep(mx, 2, colMeans(mx), FUN=`-`))$d[1:2]))
})

testex("GitHub issue #9", {
  set.seed(2)
  s1 <- irlba(diag(c(1, 2, 3, 4, 5, 0, 0, 0, 0)), 4)
  set.seed(2)
  s2 <- irlba(diag(c(1, 2, 3, 4, 5, 0, 0, 0, 0)), 4)
  stopifnot(all.equal(s1$d, s2$d))
  # Repeat this test with different seed
  set.seed(3)
  s2 <- irlba(diag(c(1, 2, 3, 4, 5, 0, 0, 0, 0)), 4)
  stopifnot(all.equal(s1$d, s2$d))
})

testex("GitHub issue #26", {
  set.seed(1)
  r <- 10
  n <- 1000
  X1 <- matrix(rnorm(n * r), n)
  X2 <- matrix(rnorm(n * r), n)
  X <- X1 %*% t(X2)
  l <- irlba(X, 20)$d
  stopifnot(all.equal(tail(l, 10), rep(0, 10)))
  l <- irlba(X, 20)$d
  stopifnot(all.equal(tail(l, 10), rep(0, 10)))
})
