/* linalg/test_ldlt.c
 *
 * Copyright (C) 2018 Patrick Alken
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
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>

static int
test_ldlt_decomp_eps(const gsl_matrix * m, const double eps, const char * desc)
{
  int s = 0;
  size_t i, j, N = m->size2;

  gsl_matrix * V  = gsl_matrix_alloc(N, N);
  gsl_matrix * A  = gsl_matrix_alloc(N, N);
  gsl_matrix * L  = gsl_matrix_calloc(N, N);
  gsl_matrix * LT = gsl_matrix_calloc(N, N);
  gsl_vector_view d;

  gsl_matrix_memcpy(V, m);

  s += gsl_linalg_ldlt_decomp(V);

  /* compute L and LT */
  gsl_matrix_tricpy('L', 0, L, V);
  d = gsl_matrix_diagonal(L);
  gsl_vector_set_all(&d.vector, 1.0);
  gsl_matrix_transpose_tricpy('L', 1, LT, L);

  /* compute L <- L D */
  d = gsl_matrix_diagonal(V);
  for (i = 0; i < N; ++i)
    {
      gsl_vector_view c = gsl_matrix_column(L, i);
      double di = gsl_vector_get(&d.vector, i);
      gsl_vector_scale(&c.vector, di);
    }

  /* compute A = L D LT */
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, L, LT, 0.0, A);

  for (i = 0; i < N; i++)
    {
      for (j = 0; j < N; j++)
        {
          double Aij = gsl_matrix_get(A, i, j);
          double mij = gsl_matrix_get(m, i, j);

          gsl_test_rel(Aij, mij, eps,
                       "%s: (%3lu,%3lu)[%lu,%lu]: %22.18g   %22.18g\n",
                       desc, N, N, i, j, Aij, mij);
        }
    }

  gsl_matrix_free(V);
  gsl_matrix_free(A);
  gsl_matrix_free(L);
  gsl_matrix_free(LT);

  return s;
}

static int
test_ldlt_decomp(gsl_rng * r)
{
  int s = 0;
  const size_t N_max = 50;
  size_t N;

  for (N = 1; N <= N_max; ++N)
    {
      gsl_matrix * m = gsl_matrix_alloc(N, N);

      create_posdef_matrix(m, r);

      test_ldlt_decomp_eps(m, 1.0e2 * N * GSL_DBL_EPSILON, "ldlt_decomp random");

      if (N <= 12)
        {
          create_hilbert_matrix2(m);
          test_ldlt_decomp_eps(m, N * GSL_DBL_EPSILON, "ldlt_decomp hilbert");
        }

      gsl_matrix_free(m);
    }

  return s;
}

int
test_ldlt_solve_eps(const gsl_matrix * m, const gsl_vector * rhs,
                    const gsl_vector * sol, const double eps,
                    const char * desc)
{
  int s = 0;
  size_t i, N = m->size1;
  gsl_matrix * u  = gsl_matrix_alloc(N, N);
  gsl_vector * x = gsl_vector_calloc(N);

  gsl_matrix_memcpy(u, m);

  s += gsl_linalg_ldlt_decomp(u);
  s += gsl_linalg_ldlt_solve(u, rhs, x);

  for (i = 0; i < N; i++)
    {
      double xi = gsl_vector_get(x, i);
      double yi = gsl_vector_get(sol, i);

      gsl_test_rel(xi, yi, eps,
                   "%s: %3lu[%lu]: %22.18g   %22.18g\n",
                   desc, N, i, xi, yi);
    }

  gsl_vector_free(x);
  gsl_matrix_free(u);

  return s;
}

static int
test_ldlt_solve(gsl_rng * r)
{
  int s = 0;
  const size_t N_max = 50;
  size_t N;

  for (N = 1; N <= N_max; ++N)
    {
      gsl_matrix * m = gsl_matrix_alloc(N, N);
      gsl_vector * rhs = gsl_vector_alloc(N);
      gsl_vector * sol = gsl_vector_alloc(N);

      create_posdef_matrix(m, r);
      create_random_vector(sol, r);
      gsl_blas_dsymv(CblasLower, 1.0, m, sol, 0.0, rhs);

      test_ldlt_solve_eps(m, rhs, sol, 64.0 * N * GSL_DBL_EPSILON, "ldlt_solve random");

      if (N <= 3)
        {
          create_hilbert_matrix2(m);
          gsl_blas_dsymv(CblasLower, 1.0, m, sol, 0.0, rhs);
          test_ldlt_solve_eps(m, rhs, sol, 1024.0 * N * GSL_DBL_EPSILON, "ldlt_solve hilbert");
        }

      gsl_matrix_free(m);
      gsl_vector_free(rhs);
      gsl_vector_free(sol);
    }

  return s;
}

static int
test_ldlt_band_decomp_eps(const size_t p, const gsl_matrix * m, const double eps, const char * desc)
{
  int s = 0;
  size_t i, j, N = m->size2;

  gsl_matrix * V  = gsl_matrix_alloc(N, p + 1);
  gsl_matrix * A  = gsl_matrix_alloc(N, N);
  gsl_matrix * L  = gsl_matrix_calloc(N, N);
  gsl_matrix * LT = gsl_matrix_calloc(N, N);
  gsl_vector * D = gsl_vector_alloc(N);
  gsl_vector * work = gsl_vector_alloc(N);

  /* convert m to packed banded format */
  symm2band_matrix(p, m, V);

  s += gsl_linalg_ldlt_band_decomp(V, work);

  /* compute L and LT */
  gsl_linalg_ldlt_band_unpack(V, L, D);
  gsl_matrix_transpose_tricpy('L', 1, LT, L);

  /* compute L <- L D */
  for (i = 0; i < N; ++i)
    {
      gsl_vector_view c = gsl_matrix_column(L, i);
      double di = gsl_vector_get(D, i);
      gsl_vector_scale(&c.vector, di);
    }

  /* compute A = L D LT */
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, L, LT, 0.0, A);

  for (i = 0; i < N; i++)
    {
      for (j = 0; j < N; j++)
        {
          double Aij = gsl_matrix_get(A, i, j);
          double mij = gsl_matrix_get(m, i, j);

          gsl_test_rel(Aij, mij, eps,
                       "%s: (p=%zu,N=%zu)[%lu,%lu]: %22.18g   %22.18g\n",
                       desc, p, N, i, j, Aij, mij);
        }
    }

  gsl_matrix_free(V);
  gsl_matrix_free(A);
  gsl_matrix_free(L);
  gsl_matrix_free(LT);
  gsl_vector_free(D);
  gsl_vector_free(work);

  return s;
}

static int
test_ldlt_band_decomp(gsl_rng * r)
{
  int s = 0;
  const size_t N_max = 50;
  size_t N, p;

  for (N = 1; N <= N_max; ++N)
    {
      gsl_matrix * m = gsl_matrix_alloc(N, N);

      for (p = 0; p < GSL_MIN(N, 10); ++p)
        {
          create_posdef_band_matrix(p, m, r);

#if 0
          /*XXX*/
          {
            gsl_vector_view v;
            gsl_matrix_set_zero(m);
            double val;
            size_t i;

            v = gsl_matrix_subdiagonal(m, 0);
            val = 162.0;
            for (i = 0; i < v.vector.size; ++i, val += 2.0)
              gsl_vector_set(&v.vector, i, val);

            v = gsl_matrix_subdiagonal(m, 1);
            val = 1.0;
            for (i = 0; i < v.vector.size; ++i, val += 1.0)
              gsl_vector_set(&v.vector, i, val);

            v = gsl_matrix_subdiagonal(m, 2);
            val = 7.0;
            for (i = 0; i < v.vector.size; ++i, val -= 1.0)
              gsl_vector_set(&v.vector, i, val);

            v = gsl_matrix_subdiagonal(m, 3);
            val = -1.0;
            for (i = 0; i < v.vector.size; ++i, val -= 1.0)
              gsl_vector_set(&v.vector, i, val);

            v = gsl_matrix_subdiagonal(m, 4);
            val = 1.0;
            for (i = 0; i < v.vector.size; ++i, val += 1.0)
              gsl_vector_set(&v.vector, i, val);

            gsl_matrix_transpose_tricpy('L', 0, m, m);
          }
#endif

          s += test_ldlt_band_decomp_eps(p, m, 1.0e2 * N * GSL_DBL_EPSILON, "ldlt_band_decomp random");
        }

      gsl_matrix_free(m);
    }

  return s;
}

int
test_ldlt_band_solve_eps(const size_t p, const gsl_matrix * m, const gsl_vector * rhs,
                         const gsl_vector * sol, const double eps, const char * desc)
{
  int s = 0;
  size_t i, N = m->size1;
  gsl_matrix * u  = gsl_matrix_alloc(N, p + 1);
  gsl_vector * x = gsl_vector_alloc(N);

  /* convert m to packed banded format */
  symm2band_matrix(p, m, u);

  s += gsl_linalg_ldlt_band_decomp(u, x);
  s += gsl_linalg_ldlt_band_solve(u, rhs, x);

  for (i = 0; i < N; i++)
    {
      double xi = gsl_vector_get(x, i);
      double yi = gsl_vector_get(sol, i);

      gsl_test_rel(xi, yi, eps,
                   "%s: p=%zu N=%zu [%lu]: %22.18g   %22.18g\n",
                   desc, p, N, i, xi, yi);
    }

  gsl_vector_free(x);
  gsl_matrix_free(u);

  return s;
}

static int
test_ldlt_band_solve(gsl_rng * r)
{
  int s = 0;
  const size_t N_max = 50;
  size_t N, p;

  for (N = 1; N <= N_max; ++N)
    {
      gsl_matrix * m = gsl_matrix_alloc(N, N);
      gsl_vector * rhs = gsl_vector_alloc(N);
      gsl_vector * sol = gsl_vector_alloc(N);

      for (p = 0; p < GSL_MIN(N, 10); ++p)
        {
          create_posdef_band_matrix(p, m, r);
          create_random_vector(sol, r);
          gsl_blas_dsymv(CblasLower, 1.0, m, sol, 0.0, rhs);

          test_ldlt_band_solve_eps(p, m, rhs, sol, 64.0 * N * GSL_DBL_EPSILON, "ldlt_band_solve random");
        }

      gsl_matrix_free(m);
      gsl_vector_free(rhs);
      gsl_vector_free(sol);
    }

  return s;
}
