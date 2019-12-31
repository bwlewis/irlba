/*
 * irlb: Implicitly restarted Lanczos bidiagonalization partial SVD.
 * Copyright (c) 2016 by Bryan W. Lewis
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <assert.h>
#include <math.h>

#include "Rinternals.h"
#include "R_ext/BLAS.h"
#include "R_ext/Lapack.h"

#include "irlb.h"


/* orthog(X,Y,...)
 * compute Y = Y - X * t(X) * Y
 * xm,xn: nrow, ncol X
 * yn: ncol Y (ASSUMED TO BE 1)
 * On entry, number of rows of Y must be xm to compute t(X) * Y and
 * T must be allocated of at least size xn * yn.
 * Modifies contents of Y.
 */
void
orthog (double *X, double *Y, double *T, int xm, int xn, int yn)
{
  double a = 1, b = 1;
  int inc = 1;
  memset (T, 0, xn * yn * sizeof (double));
  // T = t(X) * Y
  F77_NAME (dgemv) ("t", &xm, &xn, &a, X, &xm, Y, &inc, &b, T, &inc);
  // Y = Y - X * T
  a = -1.0;
  b = 1.0;
  F77_NAME (dgemv) ("n", &xm, &xn, &a, X, &xm, T, &inc, &b, Y, &inc);
}


/*
 * Convergence tests
 * Input parameters
 * Bsz            number of rows of the bidiagonal matrix B (scalar)
 * tol            convergence tolerance (scalar)
 * svtol          max change in each singular value tolerance (scalar)
 * n              requested number of singular values
 * Smax           largest singular value of B
 * svratio        vector of abs(current - previous) / current singular value ratios
 * residuals      vector of residual values
 * k              number of estimated signular values (scalar)
 * S              check for invariant subspace when S == 0
 *
 * Output
 * converged      0 = FALSE, 1 = TRUE (all converged)
 * k              adjusted subspace size.
 */
void
convtests (int Bsz, int n, double tol, double svtol, double Smax,
           double *svratio, double *residuals, int *k, int *converged, double S)
{
  int j, Len_res = 0;
  for (j = 0; j < Bsz; j++)
    {
      if ((fabs (residuals[j]) < tol * Smax) && (svratio[j] < svtol))
        Len_res++;
    }
  if (Len_res >= n || S == 0)
    {
      *converged = 1;
      return;
    }
  if (*k < n + Len_res)
    *k = n + Len_res;
  if (*k > Bsz - 3)
    *k = Bsz - 3;
  if (*k < 1)
    *k = 1;
  *converged = 0;
  return;
}

void Rmult(char trans, int m, int n, SEXP X, double * v, double * out, SEXP rho) {
  static SEXP mult_symb = NULL;
  static SEXP cross_symb = NULL;
  
  SEXP s, t;
  t = s = PROTECT(allocList(3));
  SET_TYPEOF(s, LANGSXP);

  SEXP V, result, coerced;

  if (trans=='n') {
    if (mult_symb == NULL) {
      mult_symb = installChar(asChar(mkString("irlba_mult_single")));
    }

    V = PROTECT(allocVector(REALSXP, n));
    memcpy(REAL(V), v, n * sizeof(double));

    SETCAR(t, mult_symb); t = CDR(t);
    SETCAR(t, X); t = CDR(t);
    SETCAR(t, V);

    result = eval(s, rho);
    coerced = PROTECT(coerceVector(result, REALSXP));
    memcpy(out, REAL(coerced), m * sizeof(double));
    
  } else {
    if (cross_symb == NULL) {
      cross_symb = installChar(asChar(mkString("irlba_cross_single")));
    }

    V = PROTECT(allocVector(REALSXP, m));
    memcpy(REAL(V), v, m * sizeof(double));

    SETCAR(t, cross_symb); t = CDR(t);
    SETCAR(t, X); t = CDR(t);
    SETCAR(t, V);

    result = eval(s, rho);
    coerced = PROTECT(coerceVector(result, REALSXP));
    memcpy(out, REAL(coerced), n * sizeof(double));
  }

  UNPROTECT(3);
  return;
}
