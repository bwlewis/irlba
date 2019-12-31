# ---------------------------------------------------------------------
# Internal supporting functions
# ---------------------------------------------------------------------
# General real/complex crossprod
cross <- function(x, y)
{
  if (missing(y))
  {
    if (is.complex(x)) return(abs(Conj(t(x)) %*% x))
    return(crossprod(x))
  }
  if (!is.complex(x) && !is.complex(y)) return(crossprod(x, y))
  Conj(t(x)) %*% y
}

# Euclidean norm
norm2 <- function(x)
{
  drop(sqrt(cross(x)))
}

# Orthogonalize vectors Y against vectors X.
orthog <- function(Y, X)
{
  dx2 <- dim(X)[2]
  if (is.null(dx2)) dx2 <- 1
  dy2 <- dim(Y)[2]
  if (is.null(dy2)) dy2 <- 1
  if (dx2 < dy2) doty <- cross(X, Y)
  else doty <- Conj(t(cross(Y, X)))
  Y - X %*% doty
}

# Convergence tests
# Input parameters
# Bsz            Number of rows of the bidiagonal matrix B
# tol
# k_org
# Bsvd           svd list of small matrix B
# residuals
# k
# Smax
# lastsv, svtol, work, S
#
# Output parameter list
# converged      TRUE/FALSE
# k              Number of singular vectors returned
convtests <- function(Bsz, tol, k_org, Bsvd, residuals, k, Smax, lastsv, svtol, maxritz, work, S)
{
# Converged singular triplets
  subspace_converged <- residuals[1:k_org] < tol * Smax
# Converged fixed point triplets
  if (is.null(lastsv)) lastsv <- 0 * Bsvd$d
  delta_converged <- (abs(Bsvd$d[1:k_org] - lastsv[1:k_org]) / Bsvd$d[1:k_org])  < svtol
  len_res <- sum(subspace_converged & delta_converged) # both
  if (is.na(len_res)) len_res <- 0
  if (len_res >= k_org) return(list(converged=TRUE, k=k))
  if (S == 0) return(list(converged=TRUE, k=k))
# Not converged yet...
# Adjust k to include more vectors as the number of vectors converge, but not
# too many (maxritz):
  augment <- min(sum(subspace_converged), maxritz)
  k <- min(max(k, k_org + augment), work - 1)
  list(converged=FALSE, k=k)
}

message_once <- function(..., flag)
{
  if (flag$flag) return()
  flag$flag <- TRUE
  message(...)
}

irlba_mult_single <- function(X, V) {
    drop(X %*% V)
}

irlba_cross_single <- function(X, V) {
    drop(crossprod(X, V))
}
