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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <assert.h>
#include <math.h>

#include "Rinternals.h"
#include "R_ext/BLAS.h"
#include "R_ext/Lapack.h"

#include "irlb.h"


/* orthog(X,Y,...)
 * compute Y = Y - X * t(X) * Y
 * xm,xn: nrow, ncol X
 * yn: ncol Y
 * On entry, number of rows of Y must be xm to compute t(X) * Y and
 * T must be allocated of at least size xn * yn.
 * Modifies contents of Y.
 */
void
orthog (double *X, double *Y, double *T, int xm, int xn, int yn)
{
  double a = 1, b = 0;
  memset(T, 0, xn * yn * sizeof(double));
  // T = t(X) * Y
  F77_NAME(dgemm) ("t", "n", &xn, &yn, &xm, &a, X, &xm, Y, &xm, &b, T, &xm);
  // Y = Y - X * T
  a = -1.0;
  b = 1.0;
  F77_NAME(dgemm) ("n", "n", &xm, &yn, &xn, &a, X, &xm, T, &xn, &b, Y, &xm);
}


/*
 * Convergence tests
 * Input parameters
 * Bsz            Number of rows of the bidiagonal matrix B (scalar)
 * tol            convergence tolerance (scalar)
 * n              Requested number of singular values
 * Smax           Largest singular value of B
 * residuals      vector of residual values
 * k              number of estimated signular values (scalar)
 *
 * Output
 * converged      0 = FALSE, 1 = TRUE (all converged)
 * k              Adjusted subspace size.
 */
void
convtests (int Bsz, int n, double tol, double Smax,
           double *residuals, int *k, int *converged)
{
  int j, Len_res = 0;
  for (j = 0; j < Bsz; ++j)
    {
      if (fabs(residuals[j]) < tol * Smax)
        Len_res++;
    }
  if (Len_res >= n)
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
