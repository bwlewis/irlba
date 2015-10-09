# ---------------------------------------------------------------------
# Internal supporting functions
# ---------------------------------------------------------------------
# General real/complex crossprod
cross <- function(x,y)
{
  if(missing(y))
  {
    if(is.complex(x)) return(abs(Conj(t(x)) %*% x))
    return(crossprod(x))
  }
  if(!is.complex(x) && !is.complex(y)) return(crossprod(x,y))
  Conj(t(x)) %*% y
}

# Euclidean norm
norm2 <- function (x)
{
  as.numeric(sqrt(cross(x)))
}
# Orthogonalize vectors Y against vectors X. Y and X must be R matrix
# objects (they must have a dim attribute).
# Note: this function unnecessarily copies the contents of Y
orthog <- function (Y, X)
 {
  dx2 <- dim(X)[2]
  if(is.null(dx2)) dx2 <- 1
  dy2 <- dim(Y)[2]
  if(is.null(dy2)) dy2 <- 1
  if (dx2 < dy2) doty <- cross(X, Y)
  else doty <- Conj(t(cross(Y, X)))
  return (Y - X %*% doty)
 }

# Convergence tests
# Input parameters
# Bsz            Number of rows of the bidiagonal matrix B
# tol
# k_org
# U_B            Left singular vectors of small matrix B
# S_B            Singular values of B
# V_B            Right singular vectors of B
# residuals
# k
# SVTol
# Smax
#
# Output parameter list
# converged      TRUE/FALSE
# U_B            Left singular vectors of small matrix B
# S_B            Singular values of B
# V_B            Right singular vectors of B
# k              Number of singular vectors returned
convtests <- function (Bsz, tol, k_org, U_B, S_B, V_B,
                       residuals, k, SVTol, Smax)
 {
  len_res <- sum(residuals[1:k_org] < tol * Smax)
  if(is.na(len_res)) len_res <- 0
  if (len_res == k_org) {
    return (list(converged=TRUE, U_B=U_B[,1:k_org, drop=FALSE],
                  S_B=S_B[1:k_org, drop=FALSE],
                  V_B=V_B[,1:k_org, drop=FALSE], k=k))
  }
# Not converged yet...
# Adjust k to include more vectors as the number of vectors converge.
  len_res <- sum(residuals[1:k_org] < SVTol * Smax)
  k <- max(k, k_org + len_res)
  if (k > Bsz - 3) k <- max(Bsz - 3,1)
  return (list(converged=FALSE, U_B=U_B, S_B=S_B, V_B=V_B, k=k) )
 }
