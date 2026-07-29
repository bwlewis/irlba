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

#define USE_FC_LEN_T
#include <Rinternals.h>
#include <Rconfig.h>
#include <R_ext/BLAS.h>
#include <R_ext/Complex.h>
#include <R_ext/Utils.h>
#ifndef FCONE
# define FCONE
#endif


static inline void get_dims(SEXP x, int *rows, int *cols) {
  SEXP dim = Rf_getAttrib(x, R_DimSymbol);
  if (dim != R_NilValue && Rf_length(dim) == 2) {
      *rows = INTEGER(dim)[0];
      *cols = INTEGER(dim)[1];
  } else {
      *rows = Rf_length(x);
      *cols = 1;
  }
}

SEXP direct_dgemm_c(SEXP A, SEXP B, SEXP transa, SEXP transb) {
  int ta_flag = LOGICAL(transa)[0];
  int tb_flag = LOGICAL(transb)[0];
  char ta = ta_flag ? 'T' : 'N';
  char tb = tb_flag ? 'T' : 'N';
  int a_rows, a_cols, b_rows, b_cols;
  get_dims(A, &a_rows, &a_cols);
  get_dims(B, &b_rows, &b_cols);
  int m = ta_flag ? a_cols : a_rows;
  int k = ta_flag ? a_rows : a_cols;
  int n = tb_flag ? b_rows : b_cols;
  int lda = (ta == 'N') ? m : k;
  int ldb = (tb == 'N') ? k : n;
  SEXP C;
  PROTECT(C = Rf_allocMatrix(REALSXP, m, n));
  double alpha = 1.0;
  double beta  = 0.0;
  F77_CALL(dgemm)(&ta, &tb, &m, &n, &k,
                 &alpha, REAL(A), &lda,
                 REAL(B), &ldb,
                 &beta, REAL(C), &m FCONE FCONE);
  UNPROTECT(1);
  return C;
}

SEXP direct_zgemm_c(SEXP A, SEXP B, SEXP transa, SEXP transb) {
  int ta_flag = LOGICAL(transa)[0];
  int tb_flag = LOGICAL(transb)[0];
  char ta = ta_flag ? 'C' : 'N';
  char tb = tb_flag ? 'C' : 'N';
  int a_rows, a_cols, b_rows, b_cols;
  get_dims(A, &a_rows, &a_cols);
  get_dims(B, &b_rows, &b_cols);
  int m = ta_flag ? a_cols : a_rows;
  int k = ta_flag ? a_rows : a_cols;
  int n = tb_flag ? b_rows : b_cols;
  int lda = (ta == 'N') ? m : k;
  int ldb = (tb == 'N') ? k : n;
  SEXP C;
  PROTECT(C = Rf_allocMatrix(CPLXSXP, m, n));
  Rcomplex alpha = {{1.0, 0.0}};
  Rcomplex beta  = {{0.0, 0.0}};
  F77_CALL(zgemm)(&ta, &tb, &m, &n, &k,
                 &alpha, COMPLEX(A), &lda,
                 COMPLEX(B), &ldb,
                 &beta, COMPLEX(C), &m FCONE FCONE);
  UNPROTECT(1);
  return C;
}


SEXP okatomic(SEXP x) {
  SEXPTYPE type = TYPEOF(x);
  if (type != REALSXP && type != CPLXSXP) {
    error("Input must be a numeric (real) or complex matrix.");
  }
  R_xlen_t n = XLENGTH(x);
  if (type == REALSXP) {
    double *ptr = REAL(x);
    for (R_xlen_t i = 0; i < n; i++) {
      if (!R_finite(ptr[i])) {
        return ScalarLogical(0);
      }
    }
  } else {
    Rcomplex *ptr = COMPLEX(x);
    for (R_xlen_t i = 0; i < n; i++) {
      if (!R_finite(ptr[i].r) || !R_finite(ptr[i].i)) {
        return ScalarLogical(0);
      }
    }
  }
  return ScalarLogical(1);
}


static const R_CallMethodDef CallEntries[] = {
  {"direct_dgemm_c", (DL_FUNC) &direct_dgemm_c, 4},
  {"direct_zgemm_c", (DL_FUNC) &direct_zgemm_c, 4},
  {"okatomic", (DL_FUNC) &okatomic, 1},
  {NULL, NULL, 0}
};

#ifdef HAVE_VISIBILITY_ATTRIBUTE
__attribute__ ((visibility ("default")))
#endif
void
R_init_irlba (DllInfo * dll)
{

  R_registerRoutines (dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols (dll, 0);
}
