require("irlba")

for (FAST in c(FALSE, TRUE))
{
  # Dense matrix
  set.seed(1)
  A <- matrix(rnorm(400), 20)
  L <- irlba(A, nu=2, nv=2, tol=1e-9, fastpath=FAST)
  S <- svd(A, nu=2, nv=2)
  if (!isTRUE(all.equal(L$d, S$d[1:2])))
  {
    stop("Failed simple dense singular value test", " fastpath=", FAST)
  }

  # restart
  L1 <- irlba(A, nv=3, v=L, fastpath=FAST)
  if (!isTRUE(all.equal(L1$d, S$d[1:3])))
  {
    stop("Failed restart", " fastpath=", FAST)
  }

  # Scaling and centering, dense
  s <- sqrt(apply(A, 2, crossprod))
  m <- colMeans(A)
  L <- irlba(A, 3, tol=1e-9, center=m, scale=s, fastpath=FAST)
  S <- svd(scale(A, center=TRUE, scale=s))
  if (!isTRUE(all.equal(L$d, S$d[1:3])))
  {
    stop("Failed scaling/centering test", " fastpath=", FAST)
  }
  # Scale only, non-square, dense
  A <- matrix(rnorm(200), 10)
  s <- seq(1, ncol(A))
  m <- colMeans(A)
  L <- irlba(A, 3, tol=1e-9, scale=s, fastpath=FAST)
  S <- svd(scale(A, center=FALSE, scale=s))
  if (!isTRUE(all.equal(L$d, S$d[1:3])))
  {
    stop("Failed dense scaling test", " fastpath=", FAST)
  }
  # Center only, non-square, dense
  L <- irlba(A, 3, tol=1e-9, center=m, fastpath=FAST)
  S <- svd(scale(A, center=TRUE, scale=FALSE))
  if (!isTRUE(all.equal(L$d, S$d[1:3])))
  {
    stop("Failed dense centering test", " fastpath=", FAST)
  }
  # Sparse matrix
  require("Matrix")
  K <- 400
  N <- 2000
  i <- sample(K, size=N, replace=TRUE)
  j <- sample(K, size=N, replace=TRUE)
  A <- sparseMatrix(i, j, x=rnorm(N))
  L <- irlba(A, nu=2, nv=2, tol=1e-9, fastpath=FAST)
  S <- svd(A, nu=2, nv=2)
  if (!isTRUE(all.equal(L$d, S$d[1:2])))
  {
    stop("Failed simple sparse singular value test", " fastpath=", FAST)
  }
  # Center only, sparse
  m <- colMeans(A)
  L <- irlba(A, 3, tol=1e-9, center=m, fastpath=FAST)
  S <- svd(scale(A, center=TRUE, scale=FALSE))
  if (!isTRUE(all.equal(L$d, S$d[1:3])))
  {
    stop("Failed sparse centering test", " fastpath=", FAST)
  }
  # scale only, spase
  s <- seq(1, ncol(A))
  L <- irlba(A, 3, tol=1e-9, scale=s, fastpath=FAST)
  S <- svd(scale(A, center=FALSE, scale=s))
  if (!isTRUE(all.equal(L$d, S$d[1:3])))
  {
    stop("Failed sparse scaling test", " fastpath=", FAST)
  }

  # Symmetric partial eigendecomposition
  set.seed(1)
  V <- qr.Q(qr(matrix(runif(100), nrow=10)))
  x <- V %*% diag(c(10, -9, 8, -7, 6, -5, 4, -3, 2, -1)) %*% t(V)
  if (!isTRUE(all.equal(partial_eigen(x, 3, fastpath=FAST)$values, c(10, 8, 6))))
  {
    stop("Failed partial_eigen test", " fastpath=", FAST)
  }

  # Test right-only option
  L <- irlba(A, 2, tol=1e-3, right_only=TRUE, fastpath=FAST)
  S <- svd(A, nu=2, nv=2)
  if (!isTRUE(all.equal(L$d, S$d[1:2])))
  {
    stop("Failed right_only test", " fastpath=", FAST)
  }

  # Dense complex-valued matrix
  A <- matrix(rnorm(400), 20) + 1i * matrix(rnorm(400), 20)
  L <- irlba(A, nu=2, nv=2, tol=1e-9, fastpath=FAST)
  S <- svd(A, nu=2, nv=2)
  if (!isTRUE(all.equal(L$d, S$d[1:2])))
  {
    stop("Failed complex-valued dense singular value test", " fastpath=", FAST)
  }

  # test extra reorthogonalization
  L <- irlba(A, nu=2, nv=2, tol=1e-9, reorth=TRUE, fastpath=FAST)
  if (!isTRUE(all.equal(L$d, S$d[1:2])))
  {
    stop("Failed reorthogonalization test", " fastpath=", FAST)
  }

  # prcomp convenience function
  x  <- matrix(rnorm(200), nrow=20)
  p1 <- prcomp_irlba(x, n=3, fastpath=FAST)
  p2 <- prcomp(x, tol=0.7)
  if (!isTRUE(all.equal(p1$sdev[1:2], p2$sdev[1:2])))
  {
    stop("Failed prcomp test", " fastpath=", FAST)
  }

  # very non-square dense matrices
  set.seed(1)
  A <- matrix(rnorm(2000), 20)
  L1 <- irlba(A, nu=2, nv=2, tol=1e-9, fastpath=FAST)
  L2 <- irlba(t(A), nu=2, nv=2, tol=1e-9, fastpath=FAST)
  if (!isTRUE(all.equal(L1$d, L2$d)))
  {
    stop("Failed nonsquare test", " fastpath=", FAST)
  }

# This pathological example was provided by Giuseppe Rodriguez, http://bugs.unica.it/~gppe/
# The singular values cluster at 1 and 0, making it hard to converge to a truncated
# subspace containing the largest few singular values (they are all very close).
# Or, for that matter, the smallest.
#
#   Reference:
#   J. M. Varah. The Prolate matrix. Linear Algebra and Appl.,
#   187:269-278, 1993.
#  Michela Redivo-Zaglia, University of Padova, Italy
#       Email: Michela.RedivoZaglia@unipd.it
#  Giuseppe Rodriguez, University of Cagliari, Italy
#       Email: rodriguez@unica.it
  tprolate <- function(n, w=0.25)
  {
    a <- rep(0, n)
    a[1] <- 2 * w
    a[2:n] <- sin(2 * pi * w * (1:(n-1))) / (pi * (1:(n-1)))
    toeplitz(a)
  }

  x <- tprolate(512)
  set.seed(1)
  l <- irlba(x, nv=20, fastpath=FAST)
  if (isTRUE(max(abs(l$d - 1)) > 1e-3))
  {
    stop("Failed tprolate test fastpath=", FAST)
  }
}
