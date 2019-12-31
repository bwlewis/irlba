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
#'   iteration but, overall, may require more iterations for convergence. Automatically \code{TRUE}
#'   when \code{fastpath=TRUE} (see below).
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
#'  algorithm and that \code{right_only = TRUE} sets \code{fastpath = FALSE} (only use this option
#'  for really large problems that run out of memory and when \code{nrow(A) >> ncol(A)}).
#'  Consider increasing the \code{work} option to improve accuracy with \code{right_only=TRUE}.
#' @param verbose logical value that when \code{TRUE} prints status messages during the computation.
#' @param scale optional column scaling vector whose values divide each column of \code{A};
#'   must be as long as the number of columns of \code{A} (see notes).
#' @param center optional column centering vector whose values are subtracted from each
#'   column of \code{A}; must be as long as the number of columns of \code{A} and may
#'   not be used together with the deflation options below (see notes).
#' @param shift optional shift value (square matrices only, see notes).
#' @param mult DEPRECATED optional custom matrix multiplication function (default is \code{\%*\%}, see notes).
#' @param fastpath try a fast C algorithm implementation if possible; set \code{fastpath=FALSE} to use the
#'     reference R implementation. See the notes for more details.
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
#' explicitly forming the centered matrix (which can, importantly, preserve
#' sparsity in the matrix). See the examples.
#'
#' The optional \code{shift} scalar valued argument applies only to square matrices; use it
#' to estimate the partial svd of \code{A + diag(shift, nrow(A), nrow(A))}
#' (without explicitly forming the shifted matrix).
#'
#' (Deprecated) Specify an optional alternative matrix multiplication operator in the
#' \code{mult} parameter. \code{mult} must be a function of two arguments,
#' and must handle both cases where one argument is a vector and the other
#' a matrix. This option is deprecated and will be removed in a future version.
#' The new preferred method simply uses R itself to define a custom matrix class
#' with your user-defined matrix multiplication operator. See the examples.
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
#'   happens, carefully heed the warning and inspect the result. You may also try setting \code{fastpath=FALSE}.}
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
#' The \code{fastpath=TRUE} option only supports real-valued matrices and sparse matrices
#' of type \code{dgCMatrix} (for now). Other problems fall back to the reference
#' R implementation.
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
#' # Restart the algorithm to compute more singular values
#' # (starting with an existing solution S)
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
         v=NULL,                # optional starting vector or restart
         right_only=FALSE,      # TRUE=only return V
         verbose=FALSE,         # display status messages
         scale=NULL,            # optional column scaling
         center=NULL,           # optional column centering
         shift=NULL,            # optional shift for square matrices
         mult=NULL,             # optional custom matrix multiplication func.
         fastpath=TRUE,         # use the faster C implementation if possible
         svtol=tol,             # stopping tolerance percent change in estimated svs
         smallest=FALSE,        # set to TRUE to estimate subspaces associated w/smallest singular values
         ...)                   # optional experimental or deprecated arguments
{
# ---------------------------------------------------------------------
# Check input parameters
# ---------------------------------------------------------------------
  ropts <- options(warn=1) # immediately show warnings
  mflag <- new.env()
  mflag$flag <- FALSE
  on.exit(options(ropts))
  interchange <- FALSE
  eps <- .Machine$double.eps
  # hidden support for old, removed (previously deprecated) parameters
  # this is here as a convenience to keep old code working without change
  # also supports experimental features not yet promoted to the api
  mcall <- as.list(match.call())
  random <- eval(mcall[["rng"]])
  if (is.null(random)) random <- rnorm # default RNG
  # Maximum number of Ritz vectors to use in augmentation, may be less
  # depending on workspace size.
  maxritz <- eval(mcall[["maxritz"]]) # experimental
  if (is.null(maxritz)) maxritz <- 3
  eps2 <- eval(mcall[["invariant_subspace_tolerance"]])
  if (is.null(eps2)) eps2 <- eps ^ (4 / 5)
  du <- eval(mcall[["du"]]) # deprecated
  dv <- eval(mcall[["dv"]]) # deprecated
  ds <- eval(mcall[["ds"]]) # deprecated
  deflate <- is.null(du) + is.null(ds) + is.null(dv)
  if (is.logical(scale) && ! scale) scale <- NULL
  if (is.logical(shift) && ! shift) shift <- NULL
  if (is.logical(center) && ! center) center <- NULL
  if (smallest) fastpath <- FALSE  # for now anyway
  if (any(dim(A) > 2 ^ 32 - 1)) fastpath <- FALSE # for now
  if (deflate == 3)
  {
    deflate <- FALSE
  } else if (deflate == 0)
  {
    deflate <- TRUE
    warning("The deflation options have been deprecated. Please modify your code to not use them.")
    if (length(ds) > 1) stop("deflation limited to one dimension")
    if (!is.null(dim(du))) du <- du[, 1]
    if (!is.null(dim(dv))) dv <- dv[, 1]
  } else stop("all three du ds dv parameters must be specified for deflation")
  if (!is.null(center))
  {
    if (is.logical(center) && center) center <- colMeans(A)
    if (deflate) stop("the center parameter can't be specified together with deflation parameters")
    if (length(center) != ncol(A)) stop("center must be a vector of length ncol(A)")
    if (fastpath && ! right_only) du <- NULL
    else du <- 1
    ds <- 1
    dv <- center
    deflate <- TRUE
  }
  if ("integer" == typeof(A)) A <- A + 0.0
  iscomplex <- is.complex(A)
  m <- nrow(A)
  n <- ncol(A)
  if (is.null(nu)) nu <- nv
  if (!is.null(mult) && deflate) stop("the mult parameter can't be specified together with deflation parameters")
  missingmult <- FALSE
  if (is.null(mult))
  {
    missingmult <- TRUE
    mult <- `%*%`
  }
  k <- max(nu, nv)
  if (k <= 0)  stop("max(nu, nv) must be positive")
  if (k > min(m - 1, n - 1)) stop("max(nu, nv) must be strictly less than min(nrow(A), ncol(A))")
  if (k >= 0.5 * min(m, n))
  {
    warning("You're computing too large a percentage of total singular values, use a standard svd instead.")
  }
  if (work <= 1) stop("work must be greater than 1")
  if (tol < 0) stop("tol must be non-negative")
  if (maxit <= 0) stop("maxit must be positive")
  # work must be strictly larger than requested subspace dimension, except see right_only below
  if (work <= k && ! right_only) work <- k + 1
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
    fastpath <- FALSE
  }
  if (n > m && smallest)
  {
    # Interchange dimensions m,n so that dim(A'A) = min(m,n) when seeking the
    # smallest singular values; avoids finding zero-valued smallest singular values.
    interchange <- TRUE
    temp <- m
    m <- n
    n <- temp
  }

  if (verbose)
  {
    message("Working dimension size ", work)
  }
# Check for tiny problem, use standard SVD in that case. Make definition of 'tiny' larger?
  if (min(m, n) < 6)
  {
    A <- as.matrix(A) # avoid need to define "+" and "/" for arbitrary matrix types.
    if (verbose) message("Tiny problem detected, using standard `svd` function.")
    if (!is.null(scale)) {
      A <- sweep(A, 2, scale, "/")
      dv <- dv / scale # scale the centering vector.
    }
    if (!is.null(shift)) A <- A + diag(shift, nrow(A), ncol(A))
    if (deflate)
    {
      if (is.null(du)) du <- rep(1, nrow(A))
      A <- A - (ds * du) %*% t(dv)
    }
    s <- svd(A)
    if (smallest)
    {
      return(list(d=tail(s$d, k), u=s$u[, tail(seq(ncol(s$u)), k), drop=FALSE],
              v=s$v[, tail(seq(ncol(s$v), k)), drop=FALSE], iter=0, mprod=0))
    }
    return(list(d=s$d[1:k], u=s$u[, 1:nu, drop=FALSE],
              v=s$v[, 1:nv, drop=FALSE], iter=0, mprod=0))
  }

# Try to use the fast C-language code path
  if (deflate) fastpath <- fastpath && is.null(du)
  if (fastpath && missingmult && !iscomplex && !right_only)
  {
    RESTART <- 0L
    RV <- RW <- RS <- NULL
    if (is.null(v))
    {
      v <- random(n)
      if (verbose) message("Initializing starting vector v with samples from standard normal distribution.
Use `set.seed` first for reproducibility.")
    } else if (is.list(v))  # restarted case
    {
      if (is.null(v$v) || is.null(v$d) || is.null(v$u)) stop("restart requires left and right singular vectors")
      if (max(nu, nv) <= min(ncol(v$u), ncol(v$v))) return(v) # Nothing to do!
      RESTART <- as.integer(length(v$d))
      RND <- random(n)
      RND <- orthog(RND, v$v)
      RV <- cbind(v$v, RND / norm2(RND))
      RW <- v$u
      RS <- v$d
      v <- NULL
    }

    SP <- if (is.matrix(A)) {
        0L
    } else if (is(A, "dgCMatrix")) {
        1L
    } else {
        2L
    }

    if (verbose) message("irlba: using fast C implementation")
    SCALE <- NULL
    SHIFT <- NULL
    CENTER <- NULL
    if (!is.null(scale))
    {
      if (length(scale) != ncol(A)) stop("scale length must match number of matrix columns")
      SCALE <- as.double(scale)
    }
    if (!is.null(shift))
    {
      if (length(shift) != 1) stop("shift length must be 1")
      SHIFT <- as.double(shift)
    }
    if (deflate)
    {
      if (length(center) != ncol(A)) stop("the centering vector length must match the number of matrix columns")
      CENTER <- as.double(center)
    }
    ans <- .Call(C_IRLB, A, as.integer(k), as.double(v), as.integer(work),
                 as.integer(maxit), as.double(tol), as.double(eps2), as.integer(SP),
                 as.integer(RESTART), RV, RW, RS, SCALE, SHIFT, CENTER, as.double(svtol),
                 nrow(A), ncol(A), environment()) 
    if (ans[[6]] == 0 || ans[[6]] == -2)
    {
      names(ans) <- c("d", "u", "v", "iter", "mprod", "err")
      ans$u <- matrix(head(ans$u, m * nu), nrow=m, ncol=nu)
      ans$v <- matrix(head(ans$v, n * nv), nrow=n, ncol=nv)
      if (tol * ans$d[1] < eps) warning("convergence criterion below machine epsilon")
      if (ans[[6]] == -2) warning("did not converge--results might be invalid!; try increasing work or maxit")
      return(ans[-6])
    }
    errors <- c("invalid dimensions",
                "didn't converge",
                "out of memory",
                "starting vector near the null space",
                "linear dependency encountered")
    erridx <- abs(ans[[6]])
    if (erridx > 1)
      warning("fast code path error ", errors[erridx], "; re-trying with fastpath=FALSE.", immediate.=TRUE)
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
    V[, 1] <- random(n)
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
  iter <- 1                  # Man loop iteration count
  mprod <- 0                 # Number of matrix-vector products
  R_F <- NULL                # 2-norm of residual vector F
  sqrteps <- sqrt(eps)       #
  Smax <- 1                  # Max value of all computed singular values of
                             # B est. ||A||_2
  Smin <- NULL               # Min value of all computed singular values of
                             # B est. cond(A)
  lastsv <- c()              # estimated sv in last iteration

# Check for user-supplied restart condition
  if (restart)
  {
    B <- cbind(diag(d), 0)
    k <- length(d)

    F <- random(n)
    F <- orthog(F, V[, 1:k])
    V[, k + 1] <- F / norm2(F)
  }

# Change du to be non-NULL, for non-fastpath'able matrices with non-NULL scale.
  if (deflate && is.null(du)) du <- 1

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
    if (!is.null(scale))
    {
      VJ <- VJ / scale
    }
    if (interchange) avj <- mult(VJ, A)
    else avj <- mult(A, VJ)

#   Handle non-ordinary arrays as products.
    W[, j_w] <- as.vector(avj)
    mprod <- mprod + 1

#   Optionally apply shift
    if (!is.null(shift))
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
      W[, j] <- orthog(W[, j, drop=FALSE], W[, 1:(j - 1), drop=FALSE])
    }

    S <- norm2(W[, j_w, drop=FALSE])
#   Check for linearly dependent vectors
    if (is.na(S) || S < eps2 && j == 1) stop("starting vector near the null space")
    if (is.na(S) || S < eps2)
    {
      if (verbose) message_once("invariant subspace found", flag=mflag)
      W[, j_w] <- random(nrow(W))
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
        if (interchange) F <- Conj(t(drop(mult(A, Conj(drop(W[, j_w]))))))
        else F <- Conj(t(drop(mult(Conj(drop(W[, j_w])), A))))
      }
      else
      {
        if (interchange) F <- t(drop(mult(A, drop(W[, j_w]))))
        else F <- t(drop(mult(drop(W[, j_w]), A)))
      }
#     Optionally apply shift, scale, deflate
      if (!is.null(shift)) F <- F + shift * W[, j_w]
      if (!is.null(scale)) F <- F / scale
      if (deflate) {
        sub <- sum(W[, j_w]) * dv
        if (!is.null(scale)) sub <- sub / scale
        F <- F - sub
      }
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
          if (verbose) message_once("invariant subspace found", flag=mflag)
          F <- matrix(random(dim(V)[1]), dim(V)[1], 1)
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
        if (!is.null(scale))
        {
          VJP1 <- VJP1 / scale
        }
        if (interchange) W[, jp1_w] <- drop(mult(drop(VJP1), A))
        else W[, jp1_w] <- drop(mult(A, drop(VJP1)))
        mprod <- mprod + 1

#       Optionally apply shift
        if (!is.null(shift))
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
          if (verbose) message_once("invariant subspace found", flag=mflag)
          W[, jp1_w] <- random(nrow(W))
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
    if (! reorth && Smin / Smax < sqrteps)
    {
      warning("The matrix is ill-conditioned. Basis will be reorthogonalized.")
      reorth <- TRUE
    }
    if (smallest)
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
    if (verbose)
    {
      message("iter= ", iter,
              ", mprod= ", mprod,
              ", sv[", k_org, "]=", sprintf("%.2e", Bsvd$d[k_org]),
              ", %change=", sprintf("%.2e", (Bsvd$d[k_org] - lastsv[k_org])/Bsvd$d[k_org]),
              ", k=", ct$k)
    }
    lastsv <- Bsvd$d
    k <- ct$k

#   If all desired singular values converged, then exit main loop
    if (ct$converged) break
    if (iter >= maxit) break

#   Compute the starting vectors and first block of B
    if (smallest && (Smin / Smax > sqrteps))
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
  if (!ct$converged) warning("did not converge--results might be invalid!; try increasing maxit or work")
  d <- Bsvd$d[1:k_org]
  if (!right_only)
  {
    u <- W[, 1:(dim(Bsvd$u)[1]), drop=FALSE] %*% Bsvd$u[, 1:k_org, drop=FALSE]
  }
  v <- V[, 1:(dim(Bsvd$v)[1]), drop=FALSE] %*% Bsvd$v[, 1:k_org, drop=FALSE]
  if (smallest)
  {
    reverse <- seq(length(d), 1)
    d <- d[reverse]
    if (!right_only) u <- u[, reverse]
    v <- v[, reverse]
  }
  if (tol * d[1] < eps) warning("convergence criterion below machine epsilon")
  if (right_only)
    return(list(d=d, v=v[, 1:nv, drop=FALSE], iter=iter, mprod=mprod))
  return(list(d=d, u=u[, 1:nu, drop=FALSE],
              v=v[, 1:nv, drop=FALSE], iter=iter, mprod=mprod))
}
