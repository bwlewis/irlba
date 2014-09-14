# Basic package tests. To run these tests during package
# build, set the system environment variable
# IRLBA_TEST=TRUE
# (TRUE can be any string) before building the package.
# To disable tests, unset the system environment variable
# IRLBA_TEST.

GO <- nchar(Sys.getenv("IRLBA_TEST")) > 0

if(GO)
{
  require("irlba")
# Dense matrix
  set.seed(1)
  A <- matrix(rnorm(400),20)
  L <- irlba(A,nu=2,nv=2,tol=1e-8)
  S <- svd(A,nu=2,nv=2)
  if(!isTRUE(all.equal(L$d, S$d[1:2])))
  {
    stop("Failed simple dense signular value test")
  }

# Sparse matrix
  require("Matrix")
  K <- 400
  N <- 2000
  i <- sample(K, size=(N), replace=TRUE)
  j <- sample(K, size=(N), replace=TRUE)
  A <- sparseMatrix(i,j,x=rnorm(N))
  L <- irlba(A,nu=2,nv=2,tol=1e-8)
  S <- svd(A,nu=2,nv=2)
  if(!isTRUE(all.equal(L$d, S$d[1:2])))
  {
    stop("Failed simple sparse signular value test")
  }
}
