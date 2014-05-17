# This is an early prototype. The f option is not working yet, don't use it
# (that is, this version is limited to the exponential function).

topm <-
function(A,
         directed=TRUE,
         type=c("centralities","communicability",
                "starting convenience","ending convenience"),
         f=function(sigma,mu) (exp(sigma-mu) + exp(-sigma-mu))/2,
         m=5,
         q=2,
         Nmax=ncol(A),
         tol=1e-3,
         maxit=min(100,ncol(A)),
         verbose=FALSE,
         ...) # additional arguments passed on to irlba or irblb
{
# Note, we will dispatch to an appropriate undirected algorithm here...
# For now, we only support IRLBA-based methods for directed problems.
  m <- m + 1  # XXX investigate this, not converging in  last value always...
  if(!directed) stop("Not yet supported...")
  n <- nrow(A)
  if(ncol(A)!=n) stop("A must be square")
  type <- match.arg(type)
  if(type!="centralities") stop("Only centralities supported just now...")

  lh <- la <- rep(0,n)
  sh <- sa <- rep(0,n)
  N <- 0
  L <- irlba(A,nv=q,tol=tol,...)
  mprod <- L$mprod 
  if(verbose) cat("Initial irlba dim",q,"mprod",mprod,"\n")
  while(TRUE)
  {
    N <- N + 1
    if(N>maxit) stop("Failed to converge")
    fs <- f(L$d[N],L$d[1])
    t <- L$u[,N] * L$u[,N]
    sh <- sh + t
    lh <- lh + fs*t
    zh <- lh + fs*(1 - sh)
    p <- order(zh,decreasing=TRUE)
    SH <- which(zh >= lh[p[m]])
    if(verbose) cat("|SH|",length(SH),"\n")
    if(length(SH)<=m) break
    if(N %% q == 0)
    {
browser()
# We need more singular vectors for strong convergence. Restart IRLBA.
      L <- irlba(A, nv=N+q, tol=1e-3, V=L$v, d=L$d, U=L$u)
      mprod <- mprod + L$mprod 
    }
  }
  return(list(call=match.call(), hubs=SH[1:(m-1)], authorities=NA,
              iter=L$iter, mprod=mprod, N=N))
}
