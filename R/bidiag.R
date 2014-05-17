# The Golub-Kahan Lanczos bidiagonalization, Paige and Saunders variation.

# Input
# A: m*n real matrix, m>=n (otherwise call with t(A) for efficiency).
# k: Dimension of subspace
# u: Optional starting vector of length m.
# Output
# A list with entries U, V, B matrices such that U and V have orthonoral columns
# and B is lower bidiagonal and A %*% V = U %*% B

bidiag = function(A, k, u)
{
  m = nrow(A)
  n = ncol(A)
  if(n>m) stop("# cols exceeds #rows, use t(A) instead.")
  if(missing(u)) u = rnorm(m)
  eps = .Machine$double.eps
  U = matrix(0,m,k+1)
  V = matrix(0,n,k)
  B = matrix(0,k+1,k)
  u = u/sqrt(crossprod(u)[1])
  U[,1] = u
  v = 0
  beta = 0
  for(j in 1:k)
  {
    r = t(crossprod(u,A)) - beta * v
    alpha = sqrt(crossprod(r)[1])
# XXX add check
    v = r/alpha
    p = A %*% v - alpha * u
    beta = sqrt(crossprod(p)[1])
# XXX add check
    u = p/beta
    B[j,j] = alpha
    B[j+1,j] = beta
    V[,j] = as.vector(v)
    U[,j+1] = as.vector(u)
  }
  list(V=V,U=U,B=B)
}
