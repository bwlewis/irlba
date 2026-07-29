# ---------------------------------------------------------------------
# Internal functions
# ---------------------------------------------------------------------

which_matrix <- function(x)
{
   is.matrix(x) && !isS4(x)
}

# General real/complex crossproduct
cross <- function(x, y)
{
  if(missing(y))
  { 
    if(is.null(dim(x))) {  # dense vector norm-squared case (real or complex)
      return(abs(drop(gemm(x,x,TRUE))))
    }
    if(which_matrix(x)) { # dense matrix case
      return(gemm(x, x, TRUE))
    }
    verbose("cross-product wrapper generic case")
    return(Conj(t(x)) %*% x)  # generic matrix object case (slow, memory is copied)
  }
  if((which_matrix(x) || is.null(dim(x))) && which_matrix(y)) {  # dense matrix case
    return(gemm(x, y, TRUE))
  }
  verbose("cross-product wrapper generic case")
  Conj(t(x)) %*% y    # generic case
}

# Euclidean norm of vectors
norm2 <- function(x)
{
 sqrt(abs(drop(cross(x))))
}

# Orthogonalize vectors Y against vectors X.
orthog <- function(Y, X)
{
  dx2 <- dim(X)[2]
  if (is.null(dx2)) dx2 <- 1
  dy2 <- dim(Y)[2]
  if (is.null(dy2)) dy2 <- 1
  if (dx2 < dy2) doty <-  cross(X, Y)
  else doty <-  Conj(t(cross(Y, X)))
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

oknum <- function(x) {
  if(is.atomic(x)) {
    return(.Call("okatomic", x))
  }
  is_mat <- inherits(x, "Matrix")
  if(!is_mat) return(FALSE)     # rules out any non-atomic matrix or Matrix objects.
  .Call("okatomic", as.numeric(x@x))
}

verbose <- function(...) {
  if(isTRUE(getOption("irlba.verbose"))) message(...)
}

gemm <- function(A, B, transA = FALSE, transB = FALSE) {
  transA <- as.logical(transA)
  transB <- as.logical(transB)
  n <- if(is.null(dim(A))) {
    length(A)
  } else {
    dim(A)[(!transA) + 1]
  }
  p <- if(is.null(dim(B))) {
    length(B)
  } else {
    dim(B)[transB + 1]
  }
  stopifnot(n==p)
  if(is.complex(A) && !is.complex(B)) B <- as.complex(B)
  if(!is.complex(A) && is.complex(B)) A <- as.complex(A)
  if (is.complex(A) && is.complex(B)) {
    .Call("direct_zgemm_c", A, B, transA, transB)
  } else {
    .Call("direct_dgemm_c", A, B, transA, transB)
  }
}


# limited matrix multiplication wrapper, only for use internally This is highly
# idiosyncratic to irlba, avoid using more generally. See "gemm" for general
# use.
# x, y: matrices to multiply
# interchange: boolean, if TRUE compute y %*% x, otherwise x %*% y
# U: optional deflation matrix, replace either x or y with, e.g., x = x- UU'x
# i: which matrix to deflate, i=1 deflate x, otherwise deflate y
# complex: boolean TRUE/FALSE
# tx: x <- (conjugate) transpose (x)
# ty: y <- (conjugate) transpose (y)
mult <- function(x,y,interchange=FALSE,u=NULL,i=1,tx=FALSE,ty=FALSE) {
  
  mtype <- if( (which_matrix(x) || is.null(dim(x))) && (which_matrix(y) || is.null(dim(y))) ) {
    "matrix"
  } else {
    "Matrix"
  }

  verbose("mult interchange, deflate, mtype, tx, ty: ", interchange, ",", !is.null(u), ",", mtype, ",", tx, ",", ty)
  DEFLATE = !is.null(u)
  if(isFALSE(DEFLATE)) {
    if(interchange) {
      ty <- is.null(dim(y))
      tx <- FALSE
      if(mtype == "matrix") {
        return(gemm(y, x, ty, tx))
      } else return(y %*% x)
    }
    if(mtype == "matrix") {
      tx <- is.null(dim(x))
      return(gemm(x,y,tx,ty))
    } else return(x %*% y)
  }
# Deflation...
  if(isFALSE(mtype == "matrix")) {  # Not a dense matrix, so "Matrix" or something else...
    if(i==1) { # Case i=1: deflate x
      if(interchange) {
        return(y %*% x - (y %*% u) %*% cross(u, x))
      }
      a <- x %*% y
      return(a - u %*% cross(u, a))
    }
    if(interchange) { # otherwise: deflate y
      a <- y %*% x
      return(a - u %*% cross(u, a))
    }
    return(x %*% y - (x %*% u) %*% cross(u, y))
  }
# Dense deflation case...
  if(i==1) { # Case i=1: deflate x
    if(interchange) {
      return(gemm(y,x,is.null(dim(y))) - gemm(gemm(y, u, is.null(dim(y))), cross(u, x)))
    }
    a <- gemm(x, y, is.null(dim(x)))
    return(a - u %*% cross(u, a))
  }
  if(interchange) { # otherwise: deflate y
    a <- gemm(y, x, is.null(dim(y)))
    return(a - gemm(u, cross(u, a)))
  }
  return(gemm(x,y,is.null(dim(x))) - gemm(gemm(x, u, is.null(dim(x))), cross(u, y)))
}
