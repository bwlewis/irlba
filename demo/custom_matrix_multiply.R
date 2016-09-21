library(irlba)
set.seed(1)
A <- matrix(runif(400), nrow=20)

# A custom matrix multiplication function that scales the columns of A
# (cf the scale option). This function scales the columns of A to unit norm.
# This approach uses the `mult` argument to irlba and is also found in the
# examples.
col_scale <- sqrt(apply(A, 2, crossprod))
mult <- function(x, y)
        {
          # check if x is a  vector
          if (is.vector(x))
          {
            return((x %*% y) / col_scale)
          }
          # else x is the matrix
          x %*% (y / col_scale)
        }
print(irlba(A, 3, mult=mult)$d)
     
# Compare with:
print(svd(sweep(A, 2, col_scale, FUN=`/`))$d[1:3])


# An alternative and maybe better approach simply uses R's native operator
# overloading. You simply need to define a matrix x vector product method and a
# vector x matrix product method.
setClass("mymat", contains="matrix", S3methods=TRUE, slots=c(scale="numeric"))
setMethod("%*%", signature(x="mymat", y="numeric"), function(x ,y)
  {
    x@.Data %*% (y / x@scale)
  })
setMethod("%*%", signature(x="numeric", y="mymat"), function(x ,y)
  {
    (x %*% y@.Data) / y@scale
  })
X <- new("mymat", A, scale=col_scale)

# Compare with other approaches
print(irlba(X, 3, fastpath=FALSE)$d)

# See https://bwlewis.github.io/1000_genomes_examples/PCA_whole_genome.html
# for a much more complex example using a custom class.
