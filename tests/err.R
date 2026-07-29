# Tests for expected error conditions
require("irlba")

testerr <- function(message, expr)
{
  cat(message, "\t", file=stderr())
  OK <- FALSE
  tryCatch(expr, warning=function(w) {
      OK <<- TRUE
      cat("Expected warning: ", conditionMessage(w), "\t", file=stderr())
    }, error=function(e) {
      OK <<- TRUE
      cat("Expected error: ", conditionMessage(e), "\t", file=stderr())
    })
  if(!OK) stop("Didn't throw expected error!")
  cat("OK\n", file=stderr()) 
  invisible(TRUE)
}

testok <- function (message, expr)
{
  cat(message, "\t", file=stderr())
  force(expr)
  cat("OK\n", file=stderr()) 
}

set.seed(1)

testerr("Non-numeric", {
  X <- matrix(Inf, nrow=10, ncol=10)
  irlba(X, nv=2)
})

testerr("Non-matrix", {
  irlba(runif(100), nv=2)
})

testerr("Atomic matrix or Matrix only", {
  setClass("SomeClass", slots = list(x = "matrix"))
  setMethod(
  "%*%",
  signature(x = "SomeClass", y = "ANY"),
    function(x, y) {
      x.mat <- slot(x, "x")
      x.mat %*% y
    }
  )
  X <- new("SomeClass", x = diag(100))
  irlba(X, nv=2)
})

testerr("smallest argument", {
  irlba(diag(100), nv=2, smallest=pi, warn_on_deprecated=TRUE)
})

testerr("rng argument", {
  irlba(diag(100), nv=2, smallest=pi, rng=rnorm)
})

X <- matrix(rnorm(90), ncol=10)
testerr("too many (out of bounds)", {
  irlba(X, nv=20)
})
testerr("too many", {
  irlba(X, nv=15)
})
testerr("work 0", {
  irlba(X, nv=1, work=0)
})
testerr("work 2", {
  irlba(X, nv=5, work=2)
})
testerr("tol", {
  irlba(X, nv=5, tol=-1)
})
testerr("maxit", {
  irlba(X, nv=5, maxit=-1)
})
testerr("right_only", {
  irlba(X, nv=1, smallest=TRUE, right_only=TRUE)
})
testok("tiny smallest", {
  irlba(X[1:5,1:5], nv=1, smallest=TRUE)
})
testok("tiny center", {
  irlba(X[1:5,1:5], nv=1, center=TRUE)
})
testok("tiny", {
  irlba(X[1:5,1:5], nv=1)
})
X <- matrix(rnorm(100), ncol=10)
L <-  irlba(X, nv=1) 
testok("deflate", {
  irlba(X, nv=1, v=L)
})
testerr("deflate center", {
  irlba(X, nv=1, v=L, center=TRUE)
})
testerr("deflate right_only", {
  irlba(X, nv=1, v=L, right_only=TRUE)
})






