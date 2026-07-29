#' Find a few approximate singular values and corresponding
#' singular vectors of a matrix.
#'
#' The augmented implicitly restarted Lanczos bidiagonalization algorithm
#' (IRLBA) finds a few approximate largest (or, optionally, smallest) singular
#' values and corresponding
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
#' @param tol convergence is determined when \eqn{\|A^TU - VS\| < tol\|A\|}{||A^T U - VS|| < tol*||A||},
#'   and when the maximum relative change in estimated singular values from one iteration to the
#'   next is less than \code{svtol = tol} (see \code{svtol} below),
#'   where the spectral norm ||A|| is approximated by the
#'   largest estimated singular value, and U, V, S are the matrices corresponding
#'   to the estimated left and right singular vectors, and diagonal matrix of
#'   estimated singular values, respectively.
#' @param v optional starting vector or output from a previous run of \code{irlba} used
#'   to restart the algorithm from where it left off (see the notes).
#' @param right_only logical value indicating return only the right singular vectors
#'  (\code{TRUE}) or both sets of vectors (\code{FALSE}). The right_only option can be
#'  cheaper to compute and use much less memory when \code{nrow(A) >> ncol(A)} but note
#'  that obtained solutions typically lose accuracy due to lack of re-orthogonalization in the
#'  algorithm and that \code{right_only = TRUE} (only use this option
#'  for really large problems that run out of memory and when \code{nrow(A) >> ncol(A)}).
#'  Consider increasing the \code{work} option to improve accuracy with \code{right_only=TRUE}.
#' @param verbose logical value that when \code{TRUE} prints status messages during the computation.
#' @param scale optional column scaling vector whose values divide each column of \code{A};
#'   must be as long as the number of columns of \code{A} (see notes).
#' @param center optional column centering vector whose values are subtracted from each
#'   column of \code{A}; must be as long as the number of columns of \code{A} and may
#'   not be used together with the deflation options below (see notes).
#' @param shift optional shift value (square matrices only, see notes).
#' @param svtol additional stopping tolerance on maximum allowed absolute relative change across each
#' estimated singular value between iterations.
#' The default value of this parameter is to set it to \code{tol}. You can set \code{svtol=Inf} to
#' effectively disable this stopping criterion. Setting \code{svtol=Inf} allows the method to
#' terminate on the first Lanczos iteration if it finds an invariant subspace, but with less certainty
#' that the converged subspace is the desired one. (It may, for instance, miss some of the largest
#' singular values in difficult problems.)
#' @param smallest set \code{smallest=TRUE} to estimate the smallest singular values and associated
#' singular vectors. WARNING: this option is somewhat experimental, and may produce poor
#' estimates for ill-conditioned matrices.
#' @param ... optional additional arguments used to support experimental and deprecated features.
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
#' without explicitly forming the centered matrix. \code{center}
#' must be a vector of length equal to the number of columns of \code{A}.
#' This option may be used to efficiently compute principal components without
#' explicitly forming the centered matrix (importantly, preserving
#' sparsity in the matrix). See the examples.
#'
#' The optional \code{shift} scalar valued argument applies only to square matrices; use it
#' to estimate the partial svd of \code{A + diag(shift, nrow(A), nrow(A))}
#' (without explicitly forming the shifted matrix).
#'
#' Use the \code{v} option to supply a starting vector for the iterative
#' method. A random vector is used by default (precede with \code{set.seed()}
#' for reproducibility). Optionally set \code{v} to
#' the output of a previous run of \code{irlba} to restart the method, adding
#' additional singular values/vectors without recomputing the solution
#' subspace. See the examples.
#'
#' The function may generate the following warnings:
#' \itemize{
#'   \item{"did not converge--results might be invalid!; try increasing work or maxit"
#'   means that the algorithm didn't
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
#' @references
#' Baglama, James, and Lothar Reichel. "Augmented implicitly restarted Lanczos bidiagonalization methods." SIAM Journal on Scientific Computing 27.1 (2005): 19-42.
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
#' # Restart the algorithm to compute 5 more singular values and vectors
#' # (starting with a previous solution S). S1 will contain a solution
#' # of dimension 8.
#' S1 <- irlba(A, 5, v=S)
#'
#' # Estimate smallest singular values
#' irlba(A, 3, smallest=TRUE)$d
#'
#' #Compare with
#' tail(svd(A)$d, 3)
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
#' setClass("scaled_matrix", contains="matrix", slots=c(scale="numeric"))
#' setMethod("%*%", signature(x="scaled_matrix", y="numeric"),
#'    function(x ,y) x@.Data %*% (y / x@scale))
#' setMethod("%*%", signature(x="numeric", y="scaled_matrix"),
#'    function(x ,y) (x %*% y@.Data) / y@scale)
#' a <- new("scaled_matrix", A, scale=col_scale)
#' irlba(a, 3)$d
#'
#' # Compare with:
#' svd(sweep(A, 2, col_scale, FUN=`/`))$d[1:3]
#'
#'
#' @seealso \code{\link{svd}}, \code{\link{prcomp}}, \code{\link{partial_eigen}}, \code{\link{svdr}}
#' @import Matrix
#' @importFrom stats rnorm
#' @importFrom methods slotNames
#' @useDynLib irlba, .registration=TRUE, .fixes="C_"
#' @export
irlba <-
function(A,                     # data matrix
         nv=5, nu=nv,           # number of singular vectors to estimate
         maxit=1000,            # maximum number of iterations
         work=nv + 7,           # working subspace size
         reorth=TRUE,           # TRUE=full reorthogonalization
         tol=1e-5,              # stopping tolerance
         v=NULL,                # optional starting vector or previous run output for deflated restart
         right_only=FALSE,      # TRUE=only return V (may not work with `smallest`)
         verbose=FALSE,         # display status messages
         scale=NULL,            # optional column scaling
         center=NULL,           # optional column centering, mutually exclusive with deflation
         shift=NULL,            # optional shift for square matrices
         svtol=tol,             # stopping tolerance percent change in estimated svs
         smallest=FALSE,        # set to TRUE to estimate subspaces associated w/smallest singular values
         ...)                   # optional experimental and deprecated arguments
{
# ---------------------------------------------------------------------
# Check input parameters
# ---------------------------------------------------------------------
  if(is.atomic(A) && !isTRUE(typeof(A) %in% c("double", "complex"))) {
    warning("A is neither real nor complex-valued. Attempting to coerce A to real values (double-precision numbers). This might not be what you intend, beware!")
    storage.mode(A) <- "double"
  }
  if(!oknum(A)) {
    stop("Data matrix A must contain finite, non-missing numeric (real or complex) values (atomic matrix- or Matrix-classed).")
  }
  if(is.null(dim(A)) || min(dim(A), na.rm=TRUE) < 1) {
    stop("A is not matrix-like.")
  }
  COMPLEX <- is.complex(A)
  if(!is.logical(smallest)) stop("smallest must be a valid logical value")
  ropts <- options(warn=1, irlba.verbose=verbose) # immediately show warnings, set message log level
  on.exit(options(ropts))
  mflag <- new.env()
  mflag$flag <- FALSE
  INTERCHANGE <- FALSE
  eps <- .Machine$double.eps
  random <- rnorm # default RNG

  mcall <- as.list(match.call())
  warn_deprecated <- if(isTRUE(mcall[["warn_on_deprecated"]])) {
    warning
  } else {
    invisible
  }
  if(!is.null(mcall[["rng"]])) random <- eval(mcall[["rng"]])
  # Maximum number of Ritz vectors to use in augmentation, may be less depending on workspace size.
  maxritz <- eval(mcall[["maxritz"]])
  if (is.null(maxritz)) maxritz <- 3
  eps2 <- eval(mcall[["invariant_subspace_tolerance"]])
  if (is.null(eps2)) eps2 <- eps ^ (4 / 5)
  if(!is.null(mcall[["fastpath"]])) warn_deprecated("The `fastpath` option is deprecated.")

  if(isFALSE(scale) || !oknum(scale)) scale <- NULL
  if(isFALSE(shift) || !oknum(shift)) shift <- NULL
  if(isFALSE(center)) center <- NULL
  CENTER <- if(!is.null(center))
  {
    if(isTRUE(center)) center <- colMeans(A)
    if(length(center) != ncol(A)) stop("center must be a vector of length ncol(A)")
    if(!oknum(center)) stop("center must contain valid numeric values")
    TRUE
  } else {
    FALSE
  }
  m <- nrow(A)
  n <- ncol(A)
  if(is.null(nu)) nu <- nv
  k <- max(nu, nv)
  if(!oknum(k) || isTRUE(k<=0)) stop("nv and/or nu must be finite positive integer values")
  if(k > min(m - 1, n - 1)) stop("max(nu, nv) must be strictly less than min(nrow(A), ncol(A))")
  if(k >= 0.5 * min(m, n))
  {
    warning("You're computing too large a percentage of total singular values, use a standard svd instead.")
  }
  work = round(work)
  if(work <= 1 || !oknum(work)) stop("work must be greater than 1")
  if(tol < 0 || !oknum(tol)) stop("tol must be a non-negative real number")
  if(maxit <= 0 || !oknum(maxit)) stop("maxit must be positive")
# work must be strictly larger than requested subspace dimension, except see right_only below
  if(work <= k && ! right_only) work <- k + 1
  if(work >= min(n, m))
  {
    work <- min(n, m)
    if(work <= k)
    {
      k <- work - 1  # the best we can do! Need to reduce output subspace dimension
      warning("Requested subspace dimension too large! Reduced to ", k)
    }
  }
  k_org <- k
  w_dim <- work
  if(right_only)
  {
    w_dim <- 1
  }
  if(n > m && smallest)
  {
    if(right_only) stop("`right_only` and `smallest` are mutually exclusive here, remove one.")
    verbose("working on transpose problem for smallest singular values...")
# Interchange dimensions m,n so that dim(A'A) = min(m,n) when seeking the
# smallest singular values; avoids finding zero-valued smallest singular values.
    INTERCHANGE <- TRUE
    temp <- m
    m <- n
    n <- temp
  }

  verbose("Working dimension size ", work)
# Check for tiny problem, use standard SVD in that case. Make definition of 'tiny' larger?
  if(min(m, n) < 6)
  {
    A <- as.matrix(A) # avoid need to define "+" and "/" for arbitrary matrix types.
    verbose("Tiny problem detected, using standard `svd` function.")
    if(!is.null(scale)) {
      A <- sweep(A, 2, scale, "/")
      dv <- dv / scale # scale the centering vector.
    }
    if(!is.null(shift)) A <- A + diag(shift, nrow(A), ncol(A))
    if(CENTER)
    {
      A <- A - rep(1, NROW(A)) %*% t(center)
    }
    s <- svd(A)
    if(smallest)
    {
      return(list(d=tail(s$d, k), u=s$u[, tail(seq(ncol(s$u)), k), drop=FALSE],
              v=s$v[, tail(seq(ncol(s$v), k)), drop=FALSE], iter=0, mprod=0))
    }
    return(list(d=s$d[1:k], u=s$u[, 1:nu, drop=FALSE],
              v=s$v[, 1:nv, drop=FALSE], iter=0, mprod=0))
  }

  W <- matrix(0.0, m, w_dim)
  F <- matrix(0.0, n, 1)
  DEFLATE <- if(is.list(v))
  {
    if(CENTER) stop("Sorry, centering is not compatible with deflation restarting in this implementation.")
    if(right_only) stop("Deflation (restarting) is not compatible with the right_only=TRUE option.")
    verbose("double-sided deflation")
    V <- matrix(0.0, n, work)
    V[, 1] <- random(n)
    if(m == n) {
      V[, 1] <- orthog(V[,1], v$u)
      V[, 1] <- orthog(V[,1], v$v)
    } else V[, 1] <- if(INTERCHANGE) {
      orthog(V[,1], v$u)
    } else {
      orthog(V[,1], v$v)
    }
    TRUE
  } else if(is.null(v))
  {
# If starting matrix v is not given then set V to be an (n x 1) matrix of
# normally distributed random numbers.  In any case, allocate V appropriate to
# problem size:
    V <- matrix(0.0, n, work)
    V[, 1] <- random(n)
    FALSE
  } else
  {
# user-supplied starting subspace
    V <- matrix(0.0, n, work)
    V[1:length(v)] <- v
    v <- NULL
    FALSE
  }

  if(COMPLEX) {
    W <- W + 0i
    F <- F + 0i
    V <- V + 0i
  }

# ---------------------------------------------------------------------
# Initialize local variables
# ---------------------------------------------------------------------
  B <- NULL                  # Bidiagonal matrix
  Bsz <- NULL                # Size of B
  eps23 <- eps ^ (2 / 3)     # Used for Smax/avoids using zero
  iter <- 1                  # Man loop iteration count
  mprod <- 0                 # Number of matrix-vector products
  R_F <- NULL                # 2-norm of residual vector F
  sqrteps <- sqrt(eps)       #
  Smax <- 1                  # Max value of all computed singular values of
                             # B est. ||A||_2
  Smin <- NULL               # Min value of all computed singular values of
                             # B est. cond(A)
  lastsv <- c()              # estimated sv in last iteration

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
# Note on scale and center: These options are applied implicitly below for
# computational efficiency. This complicates their application somewhat, but
# saves a few flops. This also applies to deflation (algorithm restarting).
# ---------------------------------------------------------------------
    j <- 1
#   Normalize starting vector:
    if(iter == 1)
    {
      V[, 1] <- V[, 1] / norm2(V[, 1])
    }
    else j <- k + 1
#   j_w is used here to support the right_only=TRUE case.
    j_w <- ifelse(w_dim > 1, j, 1)

#   Compute W=AV
#   Optionally apply scale
    VJ <- V[, j]
    if(!is.null(scale))
    {
      VJ <- VJ / scale
    }
    W[, j_w] <- drop(mult(A, VJ, INTERCHANGE, u=v$u))
    mprod <- mprod + 1

#   Optionally apply shift
    if(!is.null(shift))
    {
      W[, j_w] <- W[, j_w] + shift * VJ
    }

#   Optionally apply centering
    if(CENTER)
    {
      W[, j_w] <- W[, j_w] - drop(cross(center, VJ))
    }

#   Orthogonalize W
    if(iter != 1 && w_dim > 1 && reorth)
    {
      W[, j] <- orthog(W[, j, drop=FALSE], W[, 1:(j - 1), drop=FALSE])
    }

    S <- norm2(W[, j_w, drop=FALSE])
#   Check for linearly dependent vectors
    if(is.na(S) || S < eps2 && j == 1) stop("starting vector near the null space")
    if(is.na(S) || S < eps2)
    {
      if(isTRUE(getOption("irlba.verbose"))) message_once("invariant subspace found", flag=mflag)
      W[, j_w] <- random(nrow(W))
      if(w_dim > 1) W[, j] <- orthog(W[, j], W[, 1:(j - 1)])
      W[, j_w] <- W[, j_w] / norm2(W[, j_w])
      S <- 0
    }
    else W[, j_w] <- W[, j_w] / S

#   Lanczos process
    while (j <= work)
    {
      j_w <- ifelse(w_dim > 1, j, 1)
      if(COMPLEX)
      {
        F <- Conj(t(mult(W[, j_w], A, INTERCHANGE, u=v$u, i=2, tx=TRUE)))
      }
      else
      {
        F <- t(mult(W[, j_w], A, INTERCHANGE, u=v$u, i=2, tx=TRUE))
      }
#     Optionally apply shift, scale, center
      if(!is.null(shift)) F <- F + shift * W[, j_w]
      if(!is.null(scale)) F <- F / scale
      if(CENTER) {
        sub <- sum(W[, j_w]) * center
        if(!is.null(scale)) sub <- sub / scale
        F <- F - sub
      }
      mprod <- mprod + 1
      F <- drop(F - S * V[, j])
#     Orthogonalize
      F <- orthog(F, V[, 1:j, drop=FALSE])
      if(j + 1 <= work)
      {
        R <- norm2(F)
#       Check for linear dependence
        if(R < eps2)
        {
          if(isTRUE(getOption("irlba.verbose"))) message_once("invariant subspace found", flag=mflag)
          F <- matrix(random(dim(V)[1]), dim(V)[1], 1)
          F <- orthog(F, V[, 1:j, drop=FALSE])
          V[, j + 1] <- F / norm2(F)
          R <- 0
        }
        else V[, j + 1] <- F / R

#       Compute block diagonal matrix
        if(is.null(B)) B <- cbind(S, R)
        else           B <- rbind(cbind(B, 0), c(rep(0, ncol(B) - 1), S, R))

        jp1_w <- ifelse(w_dim > 1, j + 1, 1)
        w_old <- W[, j_w]

#       Optionally apply scale
        VJP1 <- V[, j + 1]
        if(!is.null(scale))
        {
          VJP1 <- VJP1 / scale
        }
        W[, jp1_w] <- drop(mult(A, VJP1, INTERCHANGE, u=v$u))
        mprod <- mprod + 1

#       Optionally apply shift
        if(!is.null(shift))
        {
          W[, jp1_w] <- W[, jp1_w] + shift * VJP1
        }

#       Optionally apply centering
        if(CENTER)
        {
          W[, jp1_w] <- W[, jp1_w] - drop(cross(center, VJP1))
        }

#       One step of the classical Gram-Schmidt process
        W[, jp1_w] <- W[, jp1_w] - R * w_old

#       Full reorthogonalization of W
        if(reorth && w_dim > 1) W[, j + 1] <- orthog(W[, j + 1], W[, 1:j])
        S <- norm2(W[, jp1_w])
#       Check for linear dependence
        if(S < eps2)
        {
          if(isTRUE(getOption("irlba.verbose"))) message_once("invariant subspace found", flag=mflag)
          W[, jp1_w] <- random(nrow(W))
          if(w_dim > 1) W[, j + 1] <- orthog(W[, j + 1], W[, 1:j])
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
# ---------------------------------------------------------------------
# (End of the Lanczos bidiagonalization part)
# ---------------------------------------------------------------------
    Bsz <- nrow(B)
    R_F <- norm2(F)
    F <- F / R_F
#   Compute singular triplets of B, svd must return ordered singular
#   values from largest to smallest.
    Bsvd <- svd(B)

#   Estimate ||A|| using the largest singular value over all iterations
#   and estimate the cond(A) using approximations to the largest and
#   smallest singular values. If a small singular value is less than sqrteps
#   require two-sided reorthogonalization.
    if(iter == 1)
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
    if(! reorth && Smin / Smax < sqrteps)
    {
      warning("The matrix is ill-conditioned. Basis will be reorthogonalized.")
      reorth <- TRUE
    }
    if(smallest)
    {
      jj <- seq(ncol(Bsvd$u), 1, by = -1)
      Bsvd$u <- Bsvd$u[, jj]
      Bsvd$d <- Bsvd$d[jj]
      Bsvd$v <- Bsvd$v[, jj]
    }

#   Compute the residuals
    R <- R_F * Bsvd$u[Bsz, , drop=FALSE]
#   Check for convergence
    ct <- convtests(Bsz, tol, k_org, Bsvd, abs(R), k, Smax, lastsv, svtol, maxritz, work, S)
    verbose("iter= ", iter,
              ", mprod= ", mprod,
              ", sv[", k_org, "]=", sprintf("%.2e", Bsvd$d[k_org]),
              ", %change=", sprintf("%.2e", (Bsvd$d[k_org] - lastsv[k_org])/Bsvd$d[k_org]),
              ", k=", ct$k)
    lastsv <- Bsvd$d
    k <- ct$k

#   If all desired singular values converged, then exit main loop
    if(ct$converged) break
    if(iter >= maxit) break

#   Compute the starting vectors and first block of B
    if(smallest && (Smin / Smax > sqrteps))
    {
#     Update the SVD of B to be the svd of [B ||F||E_m]
      Bsvd2.d <- Bsvd$d
      Bsvd2.d <- diag(Bsvd2.d, nrow=length(Bsvd2.d))
      Bsvd2 <- svd(cbind(Bsvd2.d, t(R)))
      jj <- seq(ncol(Bsvd2$u), 1, by=-1)
      Bsvd2$u <- Bsvd2$u[, jj]
      Bsvd2$d <- Bsvd2$d[jj]
      Bsvd2$v <- Bsvd2$v[, jj]

      Bsvd$d <- Bsvd2$d
      Bsvd$u <- Bsvd$u %*% Bsvd2$u
      Bsvd$v <- cbind(rbind(Bsvd$v, rep(0, Bsz)), c(rep(0, Bsz), 1)) %*% Bsvd2$v
      V_B_last <- Bsvd$v[Bsz + 1, 1:k]
      s <- R_F * solve(B, cbind(c(rep(0, Bsz - 1), 1)))
      Bsvd$v <- Bsvd$v[1:Bsz, , drop=FALSE] + s %*% Bsvd$v[Bsz + 1, ]

      qrv <- qr(cbind(rbind(Bsvd$v[, 1:k], 0), rbind(-s, 1)))
      Bsvd$v <- qr.Q(qrv)
      R <- qr.R(qrv)
      V[, 1:(k + 1)] <- cbind(V, F) %*% Bsvd$v

#  Update and compute the k by k+1 part of B
      UT <- t(R[1:(k + 1), 1:k] + R[, k + 1] %*% rbind(V_B_last))
      B <- diag(Bsvd$d[1:k], nrow=k) %*% (UT * upper.tri(UT, diag=TRUE))[1:k, 1:(k+1)]
    } else
    {
#   using the Ritz vectors
      V[, 1:(k + dim(F)[2])] <- cbind(V[, 1:(dim(Bsvd$v)[1]), drop=FALSE] %*% Bsvd$v[, 1:k], F)
      B <- cbind(diag(Bsvd$d[1:k], nrow=k), R[1:k])
    }

#   Update the left approximate singular vectors
    if(w_dim > 1)
    {
      W[, 1:k] <- W[, 1:(dim(Bsvd$u)[1]), drop=FALSE] %*% Bsvd$u[, 1:k]
    }

    iter <- iter + 1
  }
# ---------------------------------------------------------------------
# End of the main iteration loop
# Output results
# ---------------------------------------------------------------------
  if(!ct$converged) warning("did not converge--results might be invalid!; try increasing maxit or work")
  d <- Bsvd$d[1:k_org]
  if(!right_only)
  {
    U <- W[, 1:(dim(Bsvd$u)[1]), drop=FALSE] %*% Bsvd$u[, 1:k_org, drop=FALSE]
  }
  V <- V[, 1:(dim(Bsvd$v)[1]), drop=FALSE] %*% Bsvd$v[, 1:k_org, drop=FALSE]
  if(smallest)
  {
    reverse <- seq(length(d), 1)
    d <- d[reverse]
    if(!right_only) U <- U[, reverse]
    V <- V[, reverse]
  }
  if(tol * d[1] < eps) warning("convergence criterion below machine epsilon")
  if(right_only) {
    return(list(d=d, v=V[, 1:nv, drop=FALSE], iter=iter, mprod=mprod))
  }
  U = U[, 1:nu, drop=FALSE]
  V = V[, 1:nv, drop=FALSE]
  if(DEFLATE) {
    if(smallest) {
      U = if(INTERCHANGE) cbind(U, v$v) else cbind(U, v$u)
      V = if(INTERCHANGE) cbind(V, v$u) else cbind(V, v$v)
      d = c(d, v$d)
    } else {
      U = cbind(v$u, U)
      V = cbind(v$v, V)
      d = c(v$d, d)
    }
  }
  if(INTERCHANGE)  {
    return(list(d=d, u=V, v=U, iter=iter, mprod=mprod))
  }
  list(d=d, u=U, v=V, iter=iter, mprod=mprod)
}
