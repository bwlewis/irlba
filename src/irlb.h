/* Compute Y = Y - X * t(X) * Y */
void
orthog( double *X,  // Input data matrix
        double *Y,  // Input data matrix
        double *T,  // work matrix size xn * yn
        int xm,     // number of columns of X
        int xn,     // number of columns of X
        int yn);    // number of columns of Y

void
convtests (int Bsz,           // Number of rows of bidiagonal matrix B
           int n,             // requested number of singular values
           double tol,        // convergence tolerance
           double svtol,      // max change each singular value tolerane
           double Smax,       // largest singular value of B
           double *svratio,   // vector of relative singular value ratio compared to last iteration
           double *residuals, // vector of residual values
           int *k,            // number of estimated singular values (INPUT)
                              // adjusted subspace size (OUTPUT)
           int *converged,    // 0 = FALSE, 1 = TRUE
           double S);         // If S == 0 then invariant subspace found.

/*
 * Simple cholmod double-precision sparse matrix times dense vector multiplication interface
 * Compute c = op(a) %*% b, c changed on output a and b unchanged.
 * where, if transpose = 't' then op(a) = t(a) and length(b) = m, length(c) = n
 *        else,  then op(a) = a and length(b) = n, length(c) = m
 */
void
dsdmult(char transpose, // 't' -> op(a) = t(a), non-transposed a otherwise 
        int m,          // number of rows of a
        int n,          // number of columns of a
        void *a,        // double precision valued sparse matrix
        double *b,      // double precision dense vector
        double *c);     // output

void 
Rmult(char trans, 
      int m, 
      int n, 
      SEXP X, 
      double * v, 
      double * out, 
      SEXP rho);

/* IRLB function for sparse or dense double-precision valued matrices */
int
irlb(double *A,     // input data matrix (dense case)
     void *AS,      // input data matrix (sparse case)
     SEXP MAT,      // input data matrix (other case)
     int mult,      // 0 -> A is double *, 1-> A is sparse double *
     int m,         // data matrix number of rows
     int n,         // data matrix number of columns
     int nu,        // dimension of solution
     int m_b,       // working dimension
     int maxit,     // maximum number of Lanzcos iterations
     int restart,   // 0 -> no restart, 1 -> restarted form
     double tol,    // convergence tolerance
     double *scale, // optional scale (NULL for no scale) length n * 2 (1st n scale values 2nd n work)
     double *shift, // optional shift (NULL for no shift) length 1
     double *center,// optional center (NULL for no center) length n
     double *s,     // output singular vectors at least length nu
     double *U,     // output left singular vectors length >= m x m_b
     double *V,     // output right singular vectors length >= n x m_b
     int *ITER,     // output number of iterations performed
     int *MPROD,    // output number of matrix vector products
     double eps,    // machine epsilon
     int lwork,     // length for some intermediate values below
     double *V1,    // working storage  n * work
     double *U1,    // working storage  m * work
     double *W,     // working storage  m * work
     double *F,     // working storage  n
     double *B,     // working storage  work * work
     double *BU,    // working storage  work * work
     double *BV,    // working storage  work * work
     double *BS,    // working storage  work
     double *BW,    // working storage  lwork * lwork
     double *res,   // working storage  work
     double *T,     // working storage lwork
     double svtol,  // svtol tolerance on maximum ratio change per singular value per iteration
     double *SVRATIO, // working storage nu
     SEXP RHO);     // R environment
