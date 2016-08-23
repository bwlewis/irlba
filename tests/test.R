require("irlba")
# Dense matrix
set.seed(1)
A <- matrix(rnorm(400), 20)
L <- irlba(A, nu=2, nv=2, tol=1e-9)
S <- svd(A, nu=2, nv=2)
if(!isTRUE(all.equal(L$d, S$d[1:2])))
{
  stop("Failed simple dense signular value test")
}
L_slow <- irlba(A, nu=2, nv=2, tol=1e-9, fastpath=FALSE)
if(!isTRUE(all.equal(L_slow$d, S$d[1:2])))
{
  stop("Failed simple dense signular value test (fastpath=FALSE reference implementation)")
}

# restart
L1 <- irlba(A, nv=3, v=L)
if(!isTRUE(all.equal(L1$d, S$d[1:3])))
{
  stop("Failed restart")
}

# Scaling and centering
s <- sqrt(apply(A, 2, crossprod))
m <- colMeans(A)
L <- irlba(A, 3, tol=1e-9, center=m, scale=s)
S <- svd(scale(A, center=TRUE, scale=s))
if(!isTRUE(all.equal(L$d, S$d[1:3])))
{
  stop("Failed scaling/centering test")
}

# Sparse matrix
require("Matrix")
K <- 400
N <- 2000
i <- sample(K, size=N, replace=TRUE)
j <- sample(K, size=N, replace=TRUE)
A <- sparseMatrix(i, j, x=rnorm(N))
L <- irlba(A, nu=2, nv=2, tol=1e-9)
S <- svd(A, nu=2, nv=2)
if(!isTRUE(all.equal(L$d, S$d[1:2])))
{
  stop("Failed simple sparse signular value test")
}

# Symmetric partial eigendecomposition
set.seed(1)
V <- qr.Q(qr(matrix(runif(100), nrow=10)))
x <- V %*% diag(c(10, -9, 8, -7, 6, -5, 4, -3, 2, -1)) %*% t(V)
if(!isTRUE(all.equal(partial_eigen(x, 3)$values, c(10,8,6))))
{
  stop("Failed partial_eigen test")
}

# Test right-only option
L <- irlba(A, 2, tol=1e-9, right_only=TRUE)
S <- svd(A, nu=2, nv=2)
if(!isTRUE(all.equal(L$d, S$d[1:2])))
{
  stop("Failed right_only test")
}

# Dense complex-valued matrix
A <- matrix(rnorm(400), 20) + 1i * matrix(rnorm(400), 20)
L <- irlba(A, nu=2, nv=2, tol=1e-9)
S <- svd(A, nu=2, nv=2)
if(!isTRUE(all.equal(L$d, S$d[1:2])))
{
  stop("Failed complex-valued dense signular value test")
}

# test extra reorthogonalization
L <- irlba(A, nu=2, nv=2, tol=1e-9, reorth=TRUE)
if(!isTRUE(all.equal(L$d, S$d[1:2])))
{
  stop("Failed reorthogonalization test")
}

# prcomp convenience function
x  <- matrix(rnorm(200), nrow=20)
p1 <- prcomp_irlba(x, n=3)
p2 <- prcomp(x, tol=0.7)
if(!isTRUE(all.equal(p1$sdev[1:2], p2$sdev[1:2])))
{
  stop("Failed prcomp test")
}
