#' Find a few approximate largest singular values and corresponding
#' singular vectors of a matrix.
#'
#' The augmented implicitly restarted Lanczos bidiagonalization algorithm
#' (IRLBA) finds a few approximate largest singular values and corresponding
#' singular vectors of a sparse or dense matrix using a method of Baglama and
#' Reichel.  It is a fast and memory-efficient way to compute a partial SVD.
#'
#' @param A numeric real- or complex-valued matrix or real-valued sparse matrix.
#' @param nv number of right singular vectors to estimate.
#' @param nu number of left singular vectors to estimate (defaults to \code{nv}).
#' @param maxit maximum number of iterations.
#' @param work working subspace dimension, larger values can speed convergence at the cost of more memory use.
#' @param reorth if \code{TRUE}, apply full reorthogonalization to both SVD bases, otherwise
#'   only apply reorthogonalization to the right SVD basis vectors; the latter case is cheaper per
#'   iteration but, overall, may require more iterations for convergence.
#' @param tol convergence is determined when \eqn{\|AV - US\| < tol\|A\|}{||AV - US|| < tol*||A||},
#'   where the spectral norm ||A|| is approximated by the
#'   largest estimated singular value, and U, V, S are the matrices corresponding
#'   to the estimated left and right singular vectors, and diagonal matrix of
#'   estimated singular values, respectively.
#' @param v optional starting vector or output from a previous run of \code{irlba} used
#'   to restart the algorithm from where it left off (see the notes).
#' @param right_only logical value indicating return only the right singular vectors
#'  (\code{TRUE}) or both sets of vectors (\code{FALSE}). The right_only option can be
#'  cheaper to compute and use much less memory when \code{nrow(A) >> ncol(A)}.
#' @param verbose logical value that when \code{TRUE} prints status messages during the computation.
#' @param scale optional column scaling vector whose values divide each column of \code{A};
#'   must be as long as the number of columns of \code{A} (see notes).
#' @param center optional column centering vector whose values are subtracted from each
#'   column of \code{A}; must be as long as the number of columns of \code{A} and may
#"   not be used together with the deflation options below (see notes).
#' @param du DEPRECATED optional subspace deflation vector (see notes).
#' @param ds DEPRECATED optional subspace deflation scalar (see notes).
#' @param dv DEPRECATED optional subspace deflation vector (see notes).
#' @param shift optional shift value (square matrices only, see notes).
#' @param mult optional custom matrix multiplication function (default is \code{\%*\%}, see notes).
#' @param fastpath try a fast C algorithm implementation if possible; set \code{fastpath=FALSE} to use the reference R implementation. See notes.
#'
#' @return
#' Returns a list with entries:
#' \describe{
#'   \item{d:}{ max(nu, nv) approximate singular values}
#'   \item{u:}{ nu approximate left singular vectors (only when right_only=FALSE)}
#'   \item{v:}{ nv approximate right singular vectors}
#'   \item{iter:}{ The number of Lanczos iterations carried out}
#'   \item{mprod:}{ The total number of matrix vector products carried out}
#' }
#'
#' @note
#' The syntax of \code{irlba} partially follows \code{svd}, with an important
#' exception. The usual R \code{svd} function always returns a complete set of
#' singular values, even if the number of singular vectors \code{nu} or \code{nv}
#' is set less than the maximum. The \code{irlba} function returns a number of
#' estimated singular values equal to the maximum of the number of specified
#' singular vectors \code{nu} and \code{nv}.
#'
#' Use the optional \code{scale} parameter to implicitly scale each column of
#' the matrix \code{A} by the values in the \code{scale} vector, computing the
#' truncated SVD of the column-scaled \code{sweep(A, 2, scale, FUN=`/`)}, or
#' equivalently, \code{A \%*\% diag(1 / scale)}, without explicitly forming the
#' scaled matrix. \code{scale} must be a non-zero vector of length equal
#' to the number of columns of \code{A}.
#'
#' Use the optional \code{center} parameter to implicitly subtract the values
#' in the \code{center} vector from each column of \code{A}, computing the
#' truncated SVD of \code{sweep(A, 2, center, FUN=`-`)},
#' without explicitly forming the centered matrix. This option may not be
#' used together with the general rank 1 deflation options. \code{center}
#' must be a vector of length equal to the number of columns of \code{A}.
#' This option may be used to efficiently compute principal components without
#' explicitly forming the centered matrix (which can, importantly, preserve
#' sparsity in the matrix). See the examples.
#'
#' The optional deflation parameters are deprecated and will be removed in
#' a future version. They could be used to compute the rank-one deflated
#' SVD of \eqn{A - ds \cdot du dv^T}{A - ds*du \%*\% t(dv)}, where
#' \eqn{du^T A - ds\cdot dv^T = 0}{t(du) \%*\% A - ds * t(dv) == 0}. For
#' example, the triple \code{ds, du, dv} may be a known singular value
#' and corresponding singular vectors. Or \code{ds=m} and \code{dv}
#' and \code{du} represent a vector of column means of \code{A} and of ones,
#' respectively, where \code{m} is the number of rows of \code{A}.
#' This functionality can be effectively replaced with custom matrix
#' product functions.
#'
#' Specify an optional alternative matrix multiplication operator in the
#' \code{mult} parameter. \code{mult} must be a function of two arguments,
#' and must handle both cases where one argument is a vector and the other
#' a matrix. See the examples.
#'
#' Use the \code{v} option to supply a starting vector for the iterative
#' method. A random vector is used by default. Optionally set \code{v} to
#' the output of a previous run of \code{irlba} to restart the method, adding
#' additional singular values/vectors without recomputing the solution
#' subspace. See the examples.
#'
#' The function may generate the following warnings:
#' \itemize{
#'   \item{"did not converge--results might be invalid!; try increasing maxit or fastpath=FALSE" means that the algorithm didn't
#'   converge -- this is potentially a serious problem and the returned results may not be valid. \code{irlba}
#'   reports a warning here instead of an error so that you can inspect whatever is returned. If this
#'   happens, carefully heed the warning and inspect the result.}
#'   \item{"You're computing a large percentage of total singular values, standard svd might work better!"
#'     \code{irlba} is designed to efficiently compute a few of the largest singular values and associated
#'      singular vectors of a matrix. The standard \code{svd} function will be more efficient for computing
#'      large numbers of singular values than \code{irlba}.}
#'    \item{"convergence criterion below machine epsilon" means that the product of \code{tol} and the
#'      largest estimated singular value is really small and the normal convergence criterion is only
#'      met up to round off error.}
#' }
#' The function might return an error for several reasons including a situation when the starting
#' vector \code{v} is near the null space of the matrix. In that case, try a different \code{v}.
#'
#' If your matrix is sparse and \code{fastpath=TRUE} it will try to be coerced to the \code{dgCMatrix} class. Use \code{fastpath=FALSE} to bypass this coercion.
#'
#' @references
#' Augmented Implicitly Restarted Lanczos Bidiagonalization Methods, J. Baglama and L. Reichel, SIAM J. Sci. Comput. 2005.
#'
#' @examples
#' set.seed(1)
#'
#' A <- matrix(runif(400), nrow=20)
#' S <- irlba(A, 3)
#' S$d
#'
#' # Compare with svd
#' svd(A)$d[1:3]
#'
#' # Restart the algorithm to compute more singular values
#' # (starting with an existing solution S)
#' S1 <- irlba(A, 5, v=S)
#'
#' # Principal components (see also prcomp_irlba)
#' P <- irlba(A, nv=1, center=colMeans(A))
#'
#' # Compare with prcomp and prcomp_irlba (might vary up to sign)
#' cbind(P$v,
#'       prcomp(A)$rotation[, 1],
#'       prcomp_irlba(A)$rotation[, 1])
#'
#' # A custom matrix multiplication function that scales the columns of A
#' # (cf the scale option). This function scales the columns of A to unit norm.
#' col_scale <- sqrt(apply(A, 2, crossprod))
#' mult <- function(x, y)
#'         {
#'           # check if x is a  vector
#'           if (is.vector(x))
#'           {
#'             return((x %*% y) / col_scale)
#'           }
#'           # else x is the matrix
#'           x %*% (y / col_scale)
#'         }
#' irlba(A, 3, mult=mult)$d
#'
#' # Compare with:
#' irlba(A, 3, scale=col_scale)$d
#'
#' # Compare with:
#' svd(sweep(A, 2, col_scale, FUN=`/`))$d[1:3]
#'
#' @seealso \code{\link{svd}}, \code{\link{prcomp}}, \code{\link{partial_eigen}}
#' @import Matrix
#' @importFrom stats rnorm
#' @importFrom methods slotNames
#' @useDynLib irlba
#' @export
irlba <-
function (A,                     # data matrix
          nv=5, nu,              # number of singular vectors to estimate
          maxit=1000,            # maximum number of iterations
          work=nv + 7,           # working subspace size
          reorth=TRUE,           # TRUE=full reorthogonalization
          tol=1e-5,              # stopping tolerance
          v=NULL,                # optional starting vector or restart
          right_only=FALSE,      # TRUE=only return V
          verbose=FALSE,         # display status messages
          scale,                 # optional column scaling
          center,                # optional column centering
          du, ds, dv,            # optional general rank 1 deflation
          shift,                 # optional shift for square matrices
          mult,                  # optional custom matrix multiplication func.
          fastpath=TRUE)         # use the faster C implementation if possible
{
# ---------------------------------------------------------------------
# Check input parameters
# ---------------------------------------------------------------------
  ropts <- options(warn=1) # immediately show warnings
  on.exit(options(ropts))  # reset on exit
  eps <- .Machine$double.eps
  deflate <- missing(du) + missing(ds) + missing(dv)
  if (deflate == 3)
  {
    deflate <- FALSE
  } else if (deflate == 0)
  {
    deflate <- TRUE
    warning("The deflation options are deprecated and will be removed in a future version.")
    if (length(ds) > 1) stop("deflation limited to one dimension")
    if (!is.null(dim(du))) du <- du[, 1]
    if (!is.null(dim(dv))) dv <- dv[, 1]
  } else stop("all three du ds dv parameters must be specified for deflation")
  if (!missing(center))
  {
    if (deflate) stop("the center parameter can't be specified together with deflation parameters")
    if (length(center) != ncol(A)) stop("center must be a vector of length ncol(A)")
    if (fastpath && ! right_only) du <- NULL
    else du <- 1
    ds <- 1
    dv <- center
    deflate <- TRUE
  }
  iscomplex <- is.complex(A)
  m <- nrow(A)
  n <- ncol(A)
  if (missing(nu)) nu <- nv
  if (!missing(mult) && deflate) stop("the mult parameter can't be specified together with deflation parameters")
  missingmult <- FALSE
  if (missing(mult))
  {
    missingmult <- TRUE
    mult <- `%*%`
  }
  k <- max(nu, nv)
  if (k <= 0)  stop("max(nu, nv) must be positive")
  if (k > min(m - 1, n - 1)) stop("max(nu, nv) must be strictly less than min(nrow(A), ncol(A))")
  if (k >= 0.5 * min(m, n))
  {
    warning("You're computing a large percentage of total singular values, standard svd might work better!")
  }
  if (work <= 1) stop("work must be greater than 1")
  if (tol < 0) stop("tol must be non-negative")
  if (maxit <= 0) stop("maxit must be positive")
  if (work <= k) work <- k + 1 # work must be strictly larger than requested subspace dimension
  if (work >= min(n, m))
  {
    work <- min(n, m)
    if (work <= k)
    {
      k <- work - 1  # the best we can do! Need to reduce output subspace dimension
      warning("Requested subspace dimension too large! Reduced to ", k)
    }
  }
  k_org <- k
  w_dim <- work
  if (right_only)
  {
    w_dim <- 1
    work <- min(min(m, n), work + 10 ) # typically need this to help convergence
  }

  if (verbose)
  {
    message ("Working dimension size ", work)
  }
# Check for tiny problem, use standard SVD in that case. Make definition of 'tiny' larger?
  if (min(m, n) < 6)
  {
    if (verbose) warning("Tiny problem detected, using standard `svd` function.")
    algorithm = "svd"
    if (!missing(scale)) A <- A / scale
    if (!missing(shift)) A <- A + diag(shift)
    if (deflate)
    {
      if (is.null(du)) du <- rep(1, nrow(A))
      A <- A - (ds * du) %*% t(dv)
    }
    s <- svd(A)
    return(list(d=s$d[1:k], u=s$u[, 1:nu, drop=FALSE],
              v=s$v[, 1:nv, drop=FALSE], iter=0, mprod=0))
  }

# Try to use the fast C-language code path
  if (deflate) fastpath <- fastpath && is.null(du)
# Only dgCMatrix supported by fastpath for now
  if ("Matrix" %in% attributes(class(A)) && !("dgCMatrix" %in% class(A)))
  {
    fastpath <- FALSE
  }
  if (fastpath && missingmult && !iscomplex && !right_only)
  {
    RESTART <- 0
    RV <- RW <- RS <- NULL
    if (is.null(v))
      v <- rnorm(n)
    else if (is.list(v))  # restarted case
    {
      if (is.null(v$v) || is.null(v$d) || is.null(v$u)) stop("restart requires left and right singular vectors")
      if (max(nu, nv) <= min(ncol(v$u), ncol(v$v))) return(v) # Nothing to do!
      RESTART <- as.integer(length(v$d))
      RND <- rnorm(n)
      RND <- orthog(RND, v$v)
      RV <- cbind(v$v, RND / norm2(RND))
      RW <- v$u
      RS <- v$d
      v <- NULL
    }

    SP <- ifelse(is.matrix(A), 0L, 1L)
    if (verbose) message("irlba: using fast C implementation")
    SCALE <- NULL
    SHIFT <- NULL
    CENTER <- NULL
    if (!missing(scale))
    {
      if (length(scale) != ncol(A)) stop("scale length must mactch number of matrix columns")
      SCALE <- as.double(scale)
    }
    if (!missing(shift))
    {
      if (length(shift) != 1) stop("shift length must be 1")
      SHIFT <- as.double(shift)
    }
    if (deflate)
    {
      if (length(center) != ncol(A)) stop("the centering vector length must match the number of matrix columns")
      CENTER <- as.double(center)
    }
    ans <- .Call("IRLB", A, as.integer(k), as.double(v), as.integer(work),
                 as.integer(maxit), as.double(tol), .Machine$double.eps, as.integer(SP),
                 RESTART, RV, RW, RS, SCALE, SHIFT, CENTER, PACKAGE="irlba")
    if (ans[[6]] == 0 || ans[[6]] == -2)
    {
      names(ans) <- c("d", "u", "v", "iter", "mprod", "err")
      ans$u <- matrix(head(ans$u, m * nu), nrow=m, ncol=nu)
      ans$v <- matrix(head(ans$v, n * nv), nrow=n, ncol=nv)
      if(tol * ans$d[1] < eps) warning("convergence criterion below machine epsilon")
      if (ans[[6]] == -2) warning("did not converge--results might be invlaid!; try increasing maxit or fastpath=FALSE")
      return(ans[-6])
    }
    errors <- c("invalid dimensions",
                "didn't converge",
                "out of memory",
                "starting vector near the null space",
                "linear dependency encountered")
    erridx <- abs(ans[[6]])
    if (erridx > 1)
      warning("fast code path error ", errors[erridx], "; re-trying with fastpath=FALSE.")
  }

# Allocate memory for W and F:
  W <- matrix(0.0, m, w_dim)
  F <- matrix(0.0, n, 1)
  restart <- FALSE
  if (is.list(v))
  {
    if (is.null(v$v) || is.null(v$d) || is.null(v$u)) stop("restart requires left and right singular vectors")
    if (max(nu, nv) <= min(ncol(v$u), ncol(v$v))) return(v) # Nothing to do!
    right_only <- FALSE
    W[, 1:ncol(v$u)] <- v$u
    d <- v$d
    V <- matrix(0.0, n, work)
    V[, 1:ncol(v$v)] <- v$v
    restart <- TRUE
  } else if (is.null(v))
  {
# If starting matrix v is not given then set V to be an (n x 1) matrix of
# normally distributed random numbers.  In any case, allocate V appropriate to
# problem size:
    V <- matrix(0.0, n, work)
    V[, 1] <- rnorm(n)
  } else
  {
# user-supplied starting subspace
    V <- matrix(0.0, n, work)
    V[1:length(v)] <- v
  }

# ---------------------------------------------------------------------
# Initialize local variables
# ---------------------------------------------------------------------
  B <- NULL                  # Bidiagonal matrix
  Bsz <- NULL                # Size of B
  eps23 <- eps ^ (2 / 3)     # Used for Smax/avoids using zero
  eps2 <- 2 * eps
  iter <- 1                  # Man loop iteration count
  mprod <- 0                 # Number of matrix-vector products
  R_F <- NULL                # 2-norm of residual vector F
  sqrteps <- sqrt(eps)       #
  Smax <- 1                  # Max value of all computed singular values of
                             # B est. ||A||_2
  Smin <- NULL               # Min value of all computed singular values of
                             # B est. cond(A)
  SVTol <- tol  # Tolerance for singular vector convergence

# Check for user-supplied restart condition
  if (restart)
  {
    B <- cbind(diag(d), 0)
    k <- length(d)

    F <- rnorm(n)
    F <- orthog(F, V[, 1:k])
    V[, k + 1] <- F / norm2(F)
  }

# ---------------------------------------------------------------------
# Main iteration
# ---------------------------------------------------------------------
  while (iter <= maxit)
  {
# ---------------------------------------------------------------------
# Compute the Lanczos bidiagonal decomposition:
# such that AV  = WB
# and       t(A)W = VB + Ft(E)
# This routine updates W, V, F, B, mprod
#
# Note on scale and center: These options are applied implicitly below
# for maximum computational efficiency. This complicates their application
# somewhat, but saves a few flops.
# ---------------------------------------------------------------------
    j <- 1
#   Normalize starting vector:
    if (iter == 1 && !restart)
    {
      V[, 1] <- V[, 1] / norm2(V[, 1])
    }
    else j <- k + 1
#   j_w is used here to support the right_only=TRUE case.
    j_w <- ifelse(w_dim > 1, j, 1)

#   Compute W=AV
#   Optionally apply scale
    VJ <- V[, j]
    if (!missing(scale))
    {
      VJ <- VJ / scale
    }
#   Handle sparse products.
    avj <- mult(A, VJ)
    if ("Matrix" %in% attributes(class(avj)) && "x" %in% slotNames(avj))
    {
      if (length(avj@x) == nrow(W)) avj <- slot(avj, "x")
      else avj <- as.vector(avj)
    }
    W[, j_w] <- avj
    mprod <- mprod + 1

#   Optionally apply shift
    if (!missing(shift))
    {
      W[, j_w] <- W[, j_w] + shift * VJ
    }

#   Optionally apply deflation
    if (deflate)
    {
      W[, j_w] <- W[, j_w] - ds * drop(cross(dv, VJ)) * du
    }

#   Orthogonalize W
    if (iter != 1 && w_dim > 1 && reorth)
    {
      W[, j] <- orthog (W[, j, drop=FALSE], W[, 1:(j - 1), drop=FALSE])
    }

    S <- norm2(W[, j_w, drop=FALSE])
#   Check for linearly dependent vectors
    if (S < eps2 && j == 1) stop("starting vector near the null space")
    if (S < eps2)
    {
      W[, j_w] <- rnorm(nrow(W))
      if (w_dim > 1) W[, j] <- orthog(W[, j], W[, 1:(j - 1)])
      W[, j_w] <- W[, j_w] / norm2(W[, j_w])
      S <- 0
    }
    else W[, j_w] <- W[, j_w] / S

#   Lanczos process
    while (j <= work)
    {
      j_w <- ifelse(w_dim > 1, j, 1)
      if (iscomplex)
      {
        F <- Conj(t(drop(mult(Conj(drop(W[, j_w])), A))))
      }
      else F <- t(drop(mult(drop(W[, j_w]), A)))
#     Optionally apply shift and scale
      if (!missing(shift)) F <- F + shift * W[, j_w]
      if (!missing(scale)) F <- F / scale
      mprod <- mprod + 1
      F <- drop(F - S * V[, j])
#     Orthogonalize
      F <- orthog(F, V[, 1:j, drop=FALSE])
      if (j + 1 <= work)
      {
        R <- norm2(F)
#       Check for linear dependence
        if (R < eps2)
        {
          F <- matrix(rnorm(dim(V)[1]), dim(V)[1], 1)
          F <- orthog(F, V[, 1:j, drop=FALSE])
          V[, j + 1] <- F / norm2(F)
          R <- 0
        }
        else V[, j + 1] <- F / R

#       Compute block diagonal matrix
        if (is.null(B)) B <- cbind(S, R)
        else            B <- rbind(cbind(B, 0), c(rep(0, ncol(B) - 1), S, R))

        jp1_w <- ifelse(w_dim > 1, j + 1, 1)
        w_old <- W[, j_w]

#       Optionally apply scale
        VJP1 <- V[, j + 1]
        if (!missing(scale))
        {
          VJP1 <- VJP1 / scale
        }
        W[, jp1_w] <- drop(mult(A, drop(VJP1)))
        mprod <- mprod + 1

#       Optionally apply shift
        if (!missing(shift))
        {
          W[, jp1_w] <- W[, jp1_w] + shift * VJP1
        }

#       Optionally apply deflation
        if (deflate)
        {
          W[, jp1_w] <- W[, jp1_w] - ds * drop(cross(dv, VJP1)) * du
        }

#       One step of the classical Gram-Schmidt process
        W[, jp1_w] <- W[, jp1_w] - R * w_old

#       Full reorthogonalization of W
        if (reorth && w_dim > 1) W[, j + 1] <- orthog(W[, j + 1], W[, 1:j])
        S <- norm2(W[, jp1_w])
#       Check for linear dependence
        if (S < eps2)
        {
          W[, jp1_w] <- rnorm(nrow(W))
          if (w_dim > 1) W[, j + 1] <- orthog(W[, j + 1], W[, 1:j])
          W[, jp1_w] <- W[, jp1_w] / norm2(W[, jp1_w])
          S <- 0
        }
        else W[, jp1_w] <- W[, jp1_w] / S
      }
      else
      {
#       Add a last block to matrix B
        B <- rbind(B, c(rep(0, j - 1), S))
      }
      j <- j + 1
    }
    if (verbose)
    {
      message("Lanczos iter = ", iter, ", dim = ", j - 1, ", mprod = ", mprod)
    }
# ---------------------------------------------------------------------
# (End of the Lanczos bidiagonalization part)
# ---------------------------------------------------------------------
    Bsz <- nrow(B)
    R_F <- norm2(F)
    F <- F / R_F
#   Compute singular triplets of B. Expect svd to return s.v.s in order
#   from largest to smallest.
    Bsvd <- svd(B)

#   Estimate ||A|| using the largest singular value over all iterations
#   and estimate the cond(A) using approximations to the largest and
#   smallest singular values. If a small singular value is less than sqrteps
#   require two-sided reorthogonalization.
    if (iter == 1)
    {
      Smax <- Bsvd$d[1]
      Smin <- Bsvd$d[Bsz]
    }
    else
    {
      Smax <- max(Smax, Bsvd$d[1])
      Smin <- min(Smin, Bsvd$d[Bsz])
    }
    Smax <- max(eps23, Smax)
    if (Smin / Smax < sqrteps && !reorth)
    {
      warning("The matrix is ill-conditioned. Basis will be reorthogonalized.")
      reorth <- TRUE
    }

#   Compute the residuals
    R <- R_F * Bsvd$u[Bsz, , drop=FALSE]
#   Check for convergence
    ct <- convtests(Bsz, tol, k_org, Bsvd$u,
                    Bsvd$d, Bsvd$v, abs(R), k, SVTol, Smax)
    k <- ct$k

#   If all desired singular values converged, then exit main loop
    if (ct$converged) break
    if (iter >= maxit) break

#   Compute the starting vectors and first block of B[1:k, 1:(k+1), drop=FALSE]
#   using the Ritz vectors
      V[, 1:(k + dim(F)[2])] <- cbind(V[, 1:(dim(Bsvd$v)[1]),
                                     drop=FALSE] %*% Bsvd$v[, 1:k], F)
      B <- cbind( diag(Bsvd$d[1:k], nrow=k), R[1:k])

#   Update the left approximate singular vectors
    if (w_dim > 1)
    {
      W[, 1:k] <- W[, 1:(dim(Bsvd$u)[1]), drop=FALSE] %*% Bsvd$u[, 1:k]
    }

    iter <- iter + 1
  }
# ---------------------------------------------------------------------
# End of the main iteration loop
# Output results
# ---------------------------------------------------------------------
  if (iter > maxit) warning("did not converge--results might be invalid!; try increasing maxit")
  d <- Bsvd$d[1:k_org]
  if (!right_only)
  {
    u <- W[, 1:(dim(Bsvd$u)[1]), drop=FALSE] %*% Bsvd$u[, 1:k_org, drop=FALSE]
  }
  v <- V[, 1:(dim(Bsvd$v)[1]), drop=FALSE] %*% Bsvd$v[, 1:k_org, drop=FALSE]
  if(tol * d[1] < eps) warning("convergence criterion below machine epsilon")
  if (right_only)
    return(list(d=d, v=v[, 1:nv, drop=FALSE], iter=iter, mprod=mprod))
  return(list(d=d, u=u[, 1:nu, drop=FALSE],
              v=v[, 1:nv, drop=FALSE], iter=iter, mprod=mprod))
}
