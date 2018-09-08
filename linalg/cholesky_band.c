/* linalg/cholesky_band.c
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

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_cblas.h>

/*
gsl_linalg_cholesky_band_decomp()
  Cholesky decomposition of a square symmetric positive definite banded
matrix

Inputs: A - matrix in banded format, N-by-ndiag where N is the size of
            the matrix and ndiag is the number of nonzero diagonals.

Notes:
1) The main diagonal is stored in the first column of A; the first subdiagonal
in the second column and so on.
*/

int
gsl_linalg_cholesky_band_decomp(gsl_matrix * A)
{
  const size_t N = A->size1;
  const size_t ndiag = A->size2; /* number of diagonals in band, including main diagonal */
  const int kld = GSL_MAX(1, (int) ndiag - 1);
  size_t j;

  for (j = 0; j < N; ++j)
    {
      double ajj = gsl_matrix_get(A, j, 0);
      size_t kn;

      if (ajj <= 0.0)
        {
          GSL_ERROR("matrix is not positive definite", GSL_EDOM);
        }

      ajj = sqrt(ajj);
      gsl_matrix_set(A, j, 0, ajj);

      kn = GSL_MIN(ndiag - 1, N - j - 1);
      if (kn > 0)
        {
          gsl_vector_view v = gsl_matrix_subrow(A, j, 1, kn);
          gsl_matrix_view m = gsl_matrix_submatrix(A, j + 1, 0, kn, kn);

          gsl_blas_dscal(1.0 / ajj, &v.vector);

          m.matrix.tda = kld;
          gsl_blas_dsyr(CblasUpper, -1.0, &v.vector, &m.matrix);
        }
    }

  return GSL_SUCCESS;
}

int
gsl_linalg_cholesky_band_solve (const gsl_matrix * LLT,
                                const gsl_vector * b,
                                gsl_vector * x)
{
  if (LLT->size1 != b->size)
    {
      GSL_ERROR ("matrix size must match b size", GSL_EBADLEN);
    }
  else if (LLT->size1 != x->size)
    {
      GSL_ERROR ("matrix size must match solution size", GSL_EBADLEN);
    }
  else
    {
      int status;

      /* copy x <- b */
      gsl_vector_memcpy (x, b);

      status = gsl_linalg_cholesky_band_svx(LLT, x);

      return status;
    }
}

int
gsl_linalg_cholesky_band_svx (const gsl_matrix * LLT, gsl_vector * x)
{
  if (LLT->size1 != x->size)
    {
      GSL_ERROR ("matrix size must match solution size", GSL_EBADLEN);
    }
  else
    {
      /* solve for c using forward-substitution, L c = b */
      cblas_dtbsv(CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit,
                  (int) LLT->size1, (int) (LLT->size2 - 1), LLT->data, LLT->tda,
                  x->data, x->stride);

      /* perform back-substitution, L^T x = c */
      cblas_dtbsv(CblasColMajor, CblasLower, CblasTrans, CblasNonUnit,
                  (int) LLT->size1, (int) (LLT->size2 - 1), LLT->data, LLT->tda,
                  x->data, x->stride);

      return GSL_SUCCESS;
    }
}

int
gsl_linalg_cholesky_band_unpack (const gsl_matrix * LLT, gsl_matrix * L)
{
  const size_t N = LLT->size1;

  if (N != L->size1)
    {
      GSL_ERROR("L matrix does not match LLT dimensions", GSL_EBADLEN);
    }
  else if (L->size1 != L->size2)
    {
      GSL_ERROR("L matrix is not square", GSL_ENOTSQR);
    }
  else
    {
      const size_t p = LLT->size2 - 1; /* lower bandwidth */
      size_t i;

      for (i = 0; i < p + 1; ++i)
        {
          gsl_vector_const_view v = gsl_matrix_const_subcolumn(LLT, i, 0, N - i);
          gsl_vector_view w = gsl_matrix_subdiagonal(L, i);

          gsl_vector_memcpy(&w.vector, &v.vector);
        }

      return GSL_SUCCESS;
    }
}
