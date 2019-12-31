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

#include <R.h>
#define USE_RINTERNALS
#include <Rinternals.h>
#include <Rdefines.h>

#include "R_ext/BLAS.h"
#include "R_ext/Lapack.h"
#include "R_ext/Rdynload.h"
#include "R_ext/Utils.h"
#include "R_ext/Parse.h"

#include "Matrix.h"
#include "Matrix_stubs.c"

#include "irlb.h"

/* helper function for calling rnorm below */
SEXP
RNORM (int n)
{
  char buf[4096];
  SEXP cmdSexp, cmdexpr, ans = R_NilValue;
  ParseStatus status;
  cmdSexp = PROTECT (allocVector (STRSXP, 1));
  snprintf (buf, 4095, "rnorm(%d)", n);
  SET_STRING_ELT (cmdSexp, 0, mkChar (buf));
  cmdexpr = PROTECT (R_ParseVector (cmdSexp, -1, &status, R_NilValue));
  if (status != PARSE_OK)
    {
      UNPROTECT (2);
      error ("invalid call");
    }
  for (int i = 0; i < length (cmdexpr); i++)
    {
      ans = PROTECT (eval (VECTOR_ELT (cmdexpr, i), R_GlobalEnv));
      UNPROTECT (1);
    }
  UNPROTECT (2);
  return ans;
}

/* irlb C implementation wrapper for R
 *
 * X double precision input matrix
 * NU integer number of singular values/vectors to compute must be > 3
 * INIT double precision starting vector length(INIT) must equal ncol(X)
 * WORK integer working subspace dimension must be > NU
 * MAXIT integer maximum number of iterations
 * TOL double tolerance
 * EPS double invariant subspace detection tolerance
 * MULT integer 0 X is a dense matrix (dgemm), 1 sparse (cholmod)
 * RESTART integer 0 no or > 0 indicates restart of dimension n
 * RV, RW, RS optional restart V W and S values of dimension RESTART
 *    (only used when RESTART > 0)
 * SCALE either NULL (no scaling) or a vector of length ncol(X)
 * SHIFT either NULL (no shift) or a single double-precision number
 * CENTER either NULL (no centering) or a vector of length ncol(X)
 * SVTOL double tolerance max allowed per cent change in each estimated singular value
 *
 * Returns a list with 6 elements:
 * 1. vector of estimated singular values
 * 2. matrix of estimated left singular vectors
 * 3. matrix of estimated right singular vectors
 * 4. number of algorithm iterations
 * 5. number of matrix vector products
 * 6. irlb C algorithm return error code (see irlb below)
 */
SEXP
IRLB (SEXP X, SEXP NU, SEXP INIT, SEXP WORK, SEXP MAXIT, SEXP TOL, SEXP EPS,
      SEXP MULT, SEXP RESTART, SEXP RV, SEXP RW, SEXP RS, SEXP SCALE,
      SEXP SHIFT, SEXP CENTER, SEXP SVTOL,
      SEXP NROW, SEXP NCOL, SEXP RHO)
{
  SEXP ANS, S, U, V;
  double *V1, *U1, *W, *F, *B, *BU, *BV, *BS, *BW, *res, *T, *scale, *shift,
    *center, *SVRATIO;
  int i, iter, mprod, ret;
  
  int m = INTEGER(NROW)[0];
  int n = INTEGER(NCOL)[0];

  int mult = INTEGER (MULT)[0];
  void *AS = NULL;
  double *A = NULL;
  switch (mult)
    {
    case 2:
      break;
    case 1:
      AS = (void *) AS_CHM_SP (X);
      break;
    default:
      A = REAL (X);
    }
  int nu = INTEGER (NU)[0];
  int work = INTEGER (WORK)[0];
  int maxit = INTEGER (MAXIT)[0];
  double tol = REAL (TOL)[0];
  double svtol = REAL (SVTOL)[0];
  int lwork = 7 * work * (1 + work);
  int restart = INTEGER (RESTART)[0];
  double eps = REAL (EPS)[0];

  PROTECT (ANS = NEW_LIST (6));
  PROTECT (S = allocVector (REALSXP, nu));
  PROTECT (U = allocVector (REALSXP, m * work));
  PROTECT (V = allocVector (REALSXP, n * work));
  if (restart == 0)
    for (i = 0; i < n; ++i)
      (REAL (V))[i] = (REAL (INIT))[i];

  /* set up intermediate working storage */
  scale = NULL;
  shift = NULL;
  center = NULL;
  if (TYPEOF (SCALE) == REALSXP)
    {
      scale = (double *) R_alloc (n * 2, sizeof (double));
      memcpy (scale, REAL (SCALE), n * sizeof (double));
    }
  if (TYPEOF (SHIFT) == REALSXP)
    {
      shift = REAL (SHIFT);
    }
  if (TYPEOF (CENTER) == REALSXP)
    {
      center = REAL (CENTER);
    }
  SVRATIO = (double *) R_alloc (work, sizeof (double));
  V1 = (double *) R_alloc (n * work, sizeof (double));
  U1 = (double *) R_alloc (m * work, sizeof (double));
  W = (double *) R_alloc (m * work, sizeof (double));
  F = (double *) R_alloc (n, sizeof (double));
  B = (double *) R_alloc (work * work, sizeof (double));
  BU = (double *) R_alloc (work * work, sizeof (double));
  BV = (double *) R_alloc (work * work, sizeof (double));
  BS = (double *) R_alloc (work, sizeof (double));
  BW = (double *) R_alloc (lwork, sizeof (double));
  res = (double *) R_alloc (work, sizeof (double));
  T = (double *) R_alloc (lwork, sizeof (double));
  if (restart > 0)
    {
      memcpy (REAL (V), REAL (RV), n * (restart + 1) * sizeof (double));
      memcpy (W, REAL (RW), m * restart * sizeof (double));
      memset (B, 0, work * work * sizeof (double));
      for (i = 0; i < restart; ++i)
        B[i + work * i] = REAL (RS)[i];
    }
  ret =
    irlb (A, AS, X, mult, m, n, nu, work, maxit, restart, tol, scale, shift, center,
          REAL (S), REAL (U), REAL (V), &iter, &mprod, eps, lwork, V1, U1, W,
          F, B, BU, BV, BS, BW, res, T, svtol, SVRATIO, RHO);
  SET_VECTOR_ELT (ANS, 0, S);
  SET_VECTOR_ELT (ANS, 1, U);
  SET_VECTOR_ELT (ANS, 2, V);
  SET_VECTOR_ELT (ANS, 3, ScalarInteger (iter));
  SET_VECTOR_ELT (ANS, 4, ScalarInteger (mprod));
  SET_VECTOR_ELT (ANS, 5, ScalarInteger (ret));
  UNPROTECT (4);
  return ANS;
}

/* irlb: main computation function.
 * returns:
 *  0 on success,
 * -1 invalid dimensions,
 * -2 not converged
 * -3 out of memory
 * -4 starting vector near the null space of A
 *
 * all data must be allocated by caller, required sizes listed below
 */
int
irlb (double *A,                // input data matrix (double case)
      void *AS,                 // input data matrix (sparse case)
      SEXP MAT,                 // input data matrix (other case) 
      int mult,                 // 0 -> use double *A, 1 -> use AS, 2 -> use R.
      int m,                    // data matrix number of rows, must be > 3.
      int n,                    // data matrix number of columns, must be > 3.
      int nu,                   // dimension of solution
      int work,                 // working dimension, must be > 3.
      int maxit,                // maximum number of main iterations
      int restart,              // 0->no, n>0 -> restarted algorithm of dimension n
      double tol,               // convergence tolerance
      double *scale,            // optional scale (NULL for no scale) size n * 2
      double *shift,            // optional shift (NULL for no shift)
      double *center,           // optional center (NULL for no center)
      // output values
      double *s,                // output singular values at least length nu
      double *U,                // output left singular vectors  m x work
      double *V,                // output right singular vectors n x work
      int *ITER,                // ouput number of Lanczos iterations
      int *MPROD,               // output number of matrix vector products
      double eps,               // tolerance for invariant subspace detection
      // working intermediate storage, sizes shown
      int lwork, double *V1,    // n x work
      double *U1,               // m x work
      double *W,                // m x work  input when restart > 0
      double *F,                // n
      double *B,                // work x work  input when restart > 0
      double *BU,               // work x work
      double *BV,               // work x work
      double *BS,               // work
      double *BW,               // lwork
      double *res,              // work
      double *T,                // lwork
      double svtol,             // svtol limit
      double *svratio,          // convtest extra storage vector of length work
      SEXP RHO)                 // environment for R-level evaluation
{
  double d, S, R, alpha, beta, R_F, SS;
  double *x;
  int jj, kk;
  int converged;
  int info, j, k = restart;
  int inc = 1;
  int mprod = 0;
  int iter = 0;
  double Smax = 0;
  SEXP FOO;

/* Check for valid input dimensions */
  if (work < 4 || n < 4 || m < 4)
    return -1;

  if (restart == 0)
    memset (B, 0, work * work * sizeof (double));
  memset(svratio, 0, work * sizeof(double));

/* Main iteration */
  while (iter < maxit)
    {
      j = 0;
/*  Normalize starting vector */
      if (iter == 0 && restart == 0)
        {
          d = F77_NAME (dnrm2) (&n, V, &inc);
          if (d < eps)
            return -1;
          d = 1 / d;
          F77_NAME (dscal) (&n, &d, V, &inc);
        }
      else
        j = k;

/* optionally apply scale */
      x = V + j * n;
      if (scale)
        {
          x = scale + n;
          memcpy (scale + n, V + j * n, n * sizeof (double));
          for (kk = 0; kk < n; ++kk)
            x[kk] = x[kk] / scale[kk];
        }

      switch (mult)
        {
        case 2:
          Rmult('n', m, n, MAT, x, W + j * m, RHO);
          break;
        case 1:
          dsdmult ('n', m, n, AS, x, W + j * m);
          break;
        default:
          alpha = 1;
          beta = 0;
          F77_NAME (dgemv) ("n", &m, &n, &alpha, (double *) A, &m, x,
                            &inc, &beta, W + j * m, &inc);
        }
      mprod++;
      R_CheckUserInterrupt ();
/* optionally apply shift in square cases m = n */
      if (shift)
        {
          jj = j * m;
          for (kk = 0; kk < m; ++kk)
            W[jj + kk] = W[jj + kk] + shift[0] * x[kk];
        }
/* optionally apply centering */
      if (center)
        {
          jj = j * m;
          beta = F77_CALL (ddot) (&n, x, &inc, center, &inc);
          for (kk = 0; kk < m; ++kk)
            W[jj + kk] = W[jj + kk] - beta;
        }
      if (iter > 0)
        orthog (W, W + j * m, T, m, j, 1);
      S = F77_NAME (dnrm2) (&m, W + j * m, &inc);
      if (S < eps && j == 0)
        return -4;
      SS = 1.0 / S;
      F77_NAME (dscal) (&m, &SS, W + j * m, &inc);

/* The Lanczos process */
      while (j < work)
        {
          switch (mult)
            {
            case 2:
              Rmult('t', m, n, MAT, W + j * m, F, RHO);
              break;
            case 1:
              dsdmult ('t', m, n, AS, W + j * m, F);
              break;
            default:
              alpha = 1.0;
              beta = 0.0;
              F77_NAME (dgemv) ("t", &m, &n, &alpha, (double *) A, &m,
                                W + j * m, &inc, &beta, F, &inc);
            }
          mprod++;
          R_CheckUserInterrupt ();
/* optionally apply shift, scale, center */
          if (shift)
            {
              // Note, not a bug because shift only applies to square matrices
              for (kk = 0; kk < m; ++kk)
                F[kk] = F[kk] + shift[0] * W[j * m + kk];
            }
          if (scale)
            {
              for (kk = 0; kk < n; ++kk)
                F[kk] = F[kk] / scale[kk];
            }
          if (center)
            {
              beta = 0;
              for (kk = 0; kk < m; ++kk) beta += W[j *m + kk];
              if (scale)
                for (kk = 0; kk < n; ++kk)
                  F[kk] = F[kk] - beta * center[kk] / scale[kk];
              else
                for (kk = 0; kk < n; ++kk)
                  F[kk] = F[kk] - beta * center[kk];
            }
          SS = -S;
          F77_NAME (daxpy) (&n, &SS, V + j * n, &inc, F, &inc);
          orthog (V, F, T, n, j + 1, 1);

          if (j + 1 < work)
            {
              R_F = F77_NAME (dnrm2) (&n, F, &inc);
              R = 1.0 / R_F;
              if (R_F < eps)        // near invariant subspace
                {
                  FOO = RNORM (n);
                  for (kk = 0; kk < n; ++kk)
                    F[kk] = REAL (FOO)[kk];
                  orthog (V, F, T, n, j + 1, 1);
                  R_F = F77_NAME (dnrm2) (&n, F, &inc);
                  R = 1.0 / R_F;
                  R_F = 0;
                }
              memmove (V + (j + 1) * n, F, n * sizeof (double));
              F77_NAME (dscal) (&n, &R, V + (j + 1) * n, &inc);
              B[j * work + j] = S;
              B[(j + 1) * work + j] = R_F;
/* optionally apply scale */
              x = V + (j + 1) * n;
              if (scale)
                {
                  x = scale + n;
                  memcpy (x, V + (j + 1) * n, n * sizeof (double));
                  for (kk = 0; kk < n; ++kk)
                    x[kk] = x[kk] / scale[kk];
                }
              switch (mult)
                {
                case 2:
                  Rmult('n', m, n, MAT, x, W + (j + 1) * m, RHO);
                  break;
                case 1:
                  dsdmult ('n', m, n, AS, x, W + (j + 1) * m);
                  break;
                default:
                  alpha = 1.0;
                  beta = 0.0;
                  F77_NAME (dgemv) ("n", &m, &n, &alpha, (double *) A, &m,
                                    x, &inc, &beta, W + (j + 1) * m, &inc);
                }
              mprod++;
              R_CheckUserInterrupt ();
/* optionally apply shift */
              if (shift)
                {
                  jj = j + 1;
                  for (kk = 0; kk < m; ++kk)
                    W[jj * m + kk] = W[jj * m + kk] + shift[0] * x[kk];
                }
/* optionally apply centering */
              if (center)
                {
                  jj = (j + 1) * m;
                  beta = F77_CALL (ddot) (&n, x, &inc, center, &inc);
                  for (kk = 0; kk < m; ++kk)
                    W[jj + kk] = W[jj + kk] - beta;
                }
/* One step of classical Gram-Schmidt */
              R = -R_F;
              F77_NAME (daxpy) (&m, &R, W + j * m, &inc, W + (j + 1) * m,
                                &inc);
/* full re-orthogonalization of W_{j+1} */
              orthog (W, W + (j + 1) * m, T, m, j + 1, 1);
              S = F77_NAME (dnrm2) (&m, W + (j + 1) * m, &inc);
              SS = 1.0 / S;
              if (S < eps)
                {
                  FOO = RNORM (m);
                  jj = (j + 1) * m;
                  for (kk = 0; kk < m; ++kk)
                    W[jj + kk] = REAL (FOO)[kk];
                  orthog (W, W + (j + 1) * m, T, m, j + 1, 1);
                  S = F77_NAME (dnrm2) (&m, W + (j + 1) * m, &inc);
                  SS = 1.0 / S;
                  F77_NAME (dscal) (&m, &SS, W + (j + 1) * m, &inc);
                  S = 0;
                }
              else
                F77_NAME (dscal) (&m, &SS, W + (j + 1) * m, &inc);
            }
          else
            {
              B[j * work + j] = S;
            }
          j++;
        }

      memmove (BU, B, work * work * sizeof (double));   // Make a working copy of B
      int *BI = (int *) T;
      F77_NAME (dgesdd) ("O", &work, &work, BU, &work, BS, BU, &work, BV,
                         &work, BW, &lwork, BI, &info);
      R_F = F77_NAME (dnrm2) (&n, F, &inc);
      R = 1.0 / R_F;
      F77_NAME (dscal) (&n, &R, F, &inc);
/* Force termination after encountering linear dependence */
      if (R_F < eps)
        R_F = 0;

      Smax = 0;
      for (jj = 0; jj < j; ++jj)
        {
          if (BS[jj] > Smax)
            Smax = BS[jj];
          svratio[jj] = fabs (svratio[jj] - BS[jj]) / BS[jj];
        }
      for (kk = 0; kk < j; ++kk)
        res[kk] = R_F * BU[kk * work + (j - 1)];
/* Update k to be the number of converged singular values. */
      convtests (j, nu, tol, svtol, Smax, svratio, res, &k, &converged, S);

      if (converged == 1)
        {
          iter++;
          break;
        }
      for (jj = 0; jj < j; ++jj)
        svratio[jj] = BS[jj];

      alpha = 1;
      beta = 0;
      F77_NAME (dgemm) ("n", "t", &n, &k, &j, &alpha, V, &n, BV, &work, &beta,
                        V1, &n);
      memmove (V, V1, n * k * sizeof (double));
      memmove (V + n * k, F, n * sizeof (double));

      memset (B, 0, work * work * sizeof (double));
      for (jj = 0; jj < k; ++jj)
        {
          B[jj * work + jj] = BS[jj];
          B[k * work + jj] = res[jj];
        }

/*   Update the left approximate singular vectors */
      alpha = 1;
      beta = 0;
      F77_NAME (dgemm) ("n", "n", &m, &k, &j, &alpha, W, &m, BU, &work, &beta,
                        U1, &m);
      memmove (W, U1, m * k * sizeof (double));
      iter++;
    }

/* Results */
  memmove (s, BS, nu * sizeof (double));        /* Singular values */
  alpha = 1;
  beta = 0;
  F77_NAME (dgemm) ("n", "n", &m, &nu, &work, &alpha, W, &m, BU, &work, &beta,
                    U, &m);
  F77_NAME (dgemm) ("n", "t", &n, &nu, &work, &alpha, V, &n, BV, &work, &beta,
                    V1, &n);
  memmove (V, V1, n * nu * sizeof (double));

  *ITER = iter;
  *MPROD = mprod;
  return (converged == 1 ? 0 : -2);
}


cholmod_common chol_c;
/* Need our own CHOLMOD error handler */
void attribute_hidden
irlba_R_cholmod_error (int status, const char *file, int line,
                       const char *message)
{
  if (status < 0)
    error ("Cholmod error '%s' at file:%s, line %d", message, file, line);
  else
    warning ("Cholmod warning '%s' at file:%s, line %d", message, file, line);
}

static const R_CallMethodDef CallEntries[] = {
  {"IRLB", (DL_FUNC) & IRLB, 19},
  {NULL, NULL, 0}
};

#ifdef HAVE_VISIBILITY_ATTRIBUTE
__attribute__ ((visibility ("default")))
#endif
void
R_init_irlba (DllInfo * dll)
{

  R_RegisterCCallable("irlba", "orthog",
                      (DL_FUNC) &orthog);
  R_RegisterCCallable("irlba", "irlb",
                      (DL_FUNC) &irlb);


  R_registerRoutines (dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols (dll, 0);
  M_R_cholmod_start (&chol_c);
  chol_c.final_ll = 1;          /* LL' form of simplicial factorization */
  /* need own error handler, that resets  final_ll (after *_defaults()) : */
  chol_c.error_handler = irlba_R_cholmod_error;
}

void
R_unload_irlba (DllInfo * dll)
{
  M_cholmod_finish (&chol_c);
}


void
dsdmult (char transpose, int m, int n, void * a, double *b, double *c)
{
  DL_FUNC sdmult = R_GetCCallable ("Matrix", "cholmod_sdmult");
  int t = transpose == 't' ? 1 : 0;
  CHM_SP cha = (CHM_SP) a;

  cholmod_dense chb;
  chb.nrow = transpose == 't' ? m : n;
  chb.d = chb.nrow;
  chb.ncol = 1;
  chb.nzmax = chb.nrow;
  chb.xtype = cha->xtype;
  chb.dtype = 0;
  chb.x = (void *) b;
  chb.z = (void *) NULL;

  cholmod_dense chc;
  chc.nrow = transpose == 't' ? n : m;
  chc.d = chc.nrow;
  chc.ncol = 1;
  chc.nzmax = chc.nrow;
  chc.xtype = cha->xtype;
  chc.dtype = 0;
  chc.x = (void *) c;
  chc.z = (void *) NULL;

  double one[] = { 1, 0 }, zero[] = { 0, 0};
  sdmult (cha, t, one, zero, &chb, &chc, &chol_c);
}
