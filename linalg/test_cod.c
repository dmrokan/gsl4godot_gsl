/* linalg/test_cod.c
 *
 * Copyright (C) 2017 Patrick Alken
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <config.h>
#include <stdlib.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_permute_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_permutation.h>

static int test_COD_decomp_dim(const gsl_matrix * m, const double eps, const char *desc);
static int test_COD_decomp(gsl_rng *r);
static int test_COD_lssolve_dim(const gsl_matrix * m, const double * actual, const double eps, const char *desc);
static int test_COD_lssolve(void);

/* create a matrix of a given rank */
static int
create_rank_matrix(const size_t rank, gsl_matrix * m, gsl_rng * r)
{
  const size_t M = m->size1;
  const size_t N = m->size2;
  size_t i;
  gsl_vector *u = gsl_vector_alloc(M);
  gsl_vector *v = gsl_vector_alloc(N);

  gsl_matrix_set_zero(m);

  /* add several rank-1 matrices together */
  for (i = 0; i < rank; ++i)
    {
      create_random_vector(u, r);
      create_random_vector(v, r);
      gsl_blas_dger(1.0, u, v, m);
    }

  gsl_vector_free(u);
  gsl_vector_free(v);

  return GSL_SUCCESS;
}

static int
test_COD_decomp_dim(const gsl_matrix * m, const double eps, const char *desc)
{
  int s = 0;
  size_t i, j, M = m->size1, N = m->size2;
  size_t rank;

  gsl_matrix * QRZT = gsl_matrix_alloc(M, N);
  gsl_matrix * Q = gsl_matrix_alloc(M, M);
  gsl_matrix * R = gsl_matrix_alloc(M, N);
  gsl_matrix * QR = gsl_matrix_alloc(M, N);
  gsl_matrix * Z = gsl_matrix_alloc(N, N);
  gsl_vector * tau_Q = gsl_vector_alloc(GSL_MIN(M, N));
  gsl_vector * tau_Z = gsl_vector_alloc(GSL_MIN(M, N));
  gsl_vector * work = gsl_vector_alloc(N);
  gsl_matrix * lhs = gsl_matrix_alloc(M, N);
  gsl_matrix * rhs = gsl_matrix_alloc(M, N);

  gsl_permutation * perm = gsl_permutation_alloc(N);

  gsl_matrix_memcpy(QRZT, m);

  s += gsl_linalg_COD_decomp(QRZT, tau_Q, tau_Z, perm, &rank, work);
  s += gsl_linalg_COD_unpack(QRZT, tau_Q, tau_Z, rank, Q, R, Z);

  /*XXX*/
  if (rank != GSL_MIN(M, N))
    {
      fprintf(stderr, "%s (%zu,%zu) rank=%zu\n", desc, M, N, rank);
    }

  /* compute lhs = m P */
  gsl_matrix_memcpy(lhs, m);
  gsl_permute_matrix(perm, lhs);

  /* compute rhs = Q R Z */
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Q, R, 0.0, QR);
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, QR, Z, 0.0, rhs);

  for (i = 0; i < M; i++)
    {
      for (j = 0; j < N; j++)
        {
          double aij = gsl_matrix_get(rhs, i, j);
          double bij = gsl_matrix_get(lhs, i, j);

          gsl_test_rel(aij, bij, eps, "%s (%3lu,%3lu)[%lu,%lu]: %22.18g   %22.18g\n",
                       desc, M, N, i, j, aij, bij);
        }
    }

  gsl_permutation_free (perm);
  gsl_vector_free(work);
  gsl_vector_free(tau_Q);
  gsl_vector_free(tau_Z);
  gsl_matrix_free(QRZT);
  gsl_matrix_free(lhs);
  gsl_matrix_free(rhs);
  gsl_matrix_free(QR);
  gsl_matrix_free(Q);
  gsl_matrix_free(R);
  gsl_matrix_free(Z);

  return s;
}

static int
test_COD_decomp(gsl_rng *r)
{
  int s = 0;
  size_t N;
  gsl_matrix *m;

  /* test COD decomposition on Hilbert matrices */
  for (N = 1; N <= 12; ++N)
    {
      m = gsl_matrix_alloc(N, N);
      create_hilbert_matrix2(m);

      test_COD_decomp_dim(m, 256.0 * N * GSL_DBL_EPSILON, "COD_decomp hilbert");

      gsl_matrix_free(m);
    }

  /* build some matrices of a given rank and test */

  m = gsl_matrix_alloc(100, 50);
  create_rank_matrix(26, m, r);
  test_COD_decomp_dim(m, 1.0e2 * GSL_DBL_EPSILON, "COD_decomp rank 26");
  gsl_matrix_free(m);

  m = gsl_matrix_alloc(550, 200);
  create_rank_matrix(176, m, r);
  test_COD_decomp_dim(m, 1.0e3 * GSL_DBL_EPSILON, "COD_decomp rank 176");
  gsl_matrix_free(m);

  test_COD_decomp_dim(m35, 1.0e1 * GSL_DBL_EPSILON, "COD_decomp m(3,5)");
  test_COD_decomp_dim(m53, 1.0e1 * GSL_DBL_EPSILON, "COD_decomp m(5,3)");

  test_COD_decomp_dim(s35, 1.0e1 * GSL_DBL_EPSILON, "COD_decomp s(3,5)");
  test_COD_decomp_dim(s53, 1.0e1 * GSL_DBL_EPSILON, "COD_decomp s(5,3)");

  test_COD_decomp_dim(vander2, 1.0e1 * GSL_DBL_EPSILON, "COD_decomp vander(2)");
  test_COD_decomp_dim(vander3, 1.0e1 * GSL_DBL_EPSILON, "COD_decomp vander(3)");
  test_COD_decomp_dim(vander4, 1.0e1 * GSL_DBL_EPSILON, "COD_decomp vander(4)");
  test_COD_decomp_dim(vander12, 1e-3, "COD_decomp vander(12)"); /* FIXME: large tolerance needed */

  return s;
}

static int
test_COD_lssolve_dim(const gsl_matrix * m, const double * actual, const double eps, const char *desc)
{
  int s = 0;
  size_t i, M = m->size1, N = m->size2;

  gsl_vector * lhs = gsl_vector_alloc(M);
  gsl_vector * rhs = gsl_vector_alloc(M);
  gsl_matrix * QRZT  = gsl_matrix_alloc(M, N);
  gsl_vector * tau_Q = gsl_vector_alloc(GSL_MIN(M, N));
  gsl_vector * tau_Z = gsl_vector_alloc(GSL_MIN(M, N));
  gsl_vector * work = gsl_vector_alloc(N);
  gsl_vector * x = gsl_vector_alloc(N);
  gsl_vector * r = gsl_vector_alloc(M);
  gsl_vector * res = gsl_vector_alloc(M);
  gsl_permutation * perm = gsl_permutation_alloc(N);
  size_t rank;

  gsl_matrix_memcpy(QRZT, m);

  for (i = 0; i < M; i++)
    gsl_vector_set(rhs, i, i + 1.0);

  s += gsl_linalg_COD_decomp(QRZT, tau_Q, tau_Z, perm, &rank, work);
  s += gsl_linalg_COD_lssolve(QRZT, tau_Q, tau_Z, perm, rank, rhs, x, res);

  for (i = 0; i < N; i++)
    {
      double xi = gsl_vector_get(x, i);
      gsl_test_rel(xi, actual[i], eps,
                   "%s (%3lu,%3lu)[%lu]: %22.18g   %22.18g\n",
                   desc, M, N, i, xi, actual[i]);
    }

  /* compute residual r = b - m x */
  if (M == N)
    {
      gsl_vector_set_zero(r);
    }
  else
    {
      gsl_vector_memcpy(r, rhs);
      gsl_blas_dgemv(CblasNoTrans, -1.0, m, x, 1.0, r);
    }

  for (i = 0; i < N; i++)
    {
      gsl_test_rel(gsl_vector_get(res, i), gsl_vector_get(r, i), sqrt(eps),
                   "%s res (%3lu,%3lu)[%lu]: %22.18g   %22.18g\n",
                   desc, M, N, i, gsl_vector_get(res, i), gsl_vector_get(r,i));
    }

  gsl_vector_free(r);
  gsl_vector_free(res);
  gsl_vector_free(x);
  gsl_vector_free(tau_Q);
  gsl_vector_free(tau_Z);
  gsl_matrix_free(QRZT);
  gsl_vector_free(rhs);
  gsl_vector_free(lhs);
  gsl_vector_free(work);
  gsl_permutation_free(perm);

  return s;
}

static int
test_COD_lssolve(void)
{
  int s = 0;

  test_COD_lssolve_dim(m53, m53_lssolution, 1.0e2 * GSL_DBL_EPSILON, "COD_lssolve m(5,3)");

  test_COD_lssolve_dim(hilb2, hilb2_solution, 1.0e2 * GSL_DBL_EPSILON, "COD_lssolve hilbert(2)");
  test_COD_lssolve_dim(hilb3, hilb3_solution, 1.0e2 * GSL_DBL_EPSILON, "COD_lssolve hilbert(3)");
  test_COD_lssolve_dim(hilb4, hilb4_solution, 1.0e4 * GSL_DBL_EPSILON, "COD_lssolve hilbert(4)");

  test_COD_lssolve_dim(vander2, vander2_solution, 1.0e1 * GSL_DBL_EPSILON, "COD_lssolve vander(2)");
  test_COD_lssolve_dim(vander3, vander3_solution, 1.0e1 * GSL_DBL_EPSILON, "COD_lssolve vander(3)");
  test_COD_lssolve_dim(vander4, vander4_solution, 1.0e2 * GSL_DBL_EPSILON, "COD_lssolve vander(4)");

  /* rank-1 least squares problem from the 'lin2' test dataset for multifit_nlinear */
  {
    const size_t M = 20;
    const size_t N = 5;

    /* unique minimum norm solution from "Factorize" matlab package */
    const double x_sol[] = { 1.818181818181817e-02, 3.636363636363636e-02, 5.454545454545454e-02,
                             7.272727272727272e-02, 9.090909090909088e-02 };

    gsl_matrix *m = gsl_matrix_alloc(M, N);
    size_t i, j;

    for (i = 0; i < M; ++i)
      {
        for (j = 0; j < N; ++j)
          {
            gsl_matrix_set(m, i, j, (i + 1.0) * (j + 1.0));
          }
      }

    test_COD_lssolve_dim(m, x_sol, 1.0e2 * GSL_DBL_EPSILON, "COD_lssolve lin2");

    gsl_matrix_free(m);
  }

  return s;
}
