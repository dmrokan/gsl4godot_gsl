/* linalg/ldlt_band.c
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

/* L D L^T decomposition of a symmetric banded positive semi-definite matrix */

#include <config.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

static _gsl_vector_view symband_row_before_diag(gsl_matrix * AB, const size_t i);

#define NULL_VECTOR {0, 0, 0, 0, 0}
#define NULL_VECTOR_VIEW {{0, 0, 0, 0, 0}}

/*
symband_row_before_diag()
  Returns a vector view of the non-zero elements of row i not
including the diagonal element for a matrix in symmetric banded format.
The row structure looks like this:

AB = [  *      *      *   a_{i1} * ... * ]
     [             a_{i2}    *   * ... * ]
     [  *   a_{i3}    *      *   * ... * ]
     [ ...     *      *      *   * ... * ]
*/

static _gsl_vector_view
symband_row_before_diag(gsl_matrix * AB, const size_t i)
{
  _gsl_vector_view view = NULL_VECTOR_VIEW;

  if (i >= AB->size1)
    {
      GSL_ERROR_VAL ("row index is out of range", GSL_EINVAL, view);
    }
  else
    {
      const size_t p = AB->size2 - 1;         /* lower bandwidth */
      const size_t nz = (i <= p) ? 0 : i - p; /* number of zero entries at beginning of row i */
      gsl_vector v = NULL_VECTOR;

      v.data = AB->data + (nz * p + i);
      v.size = i - nz;
      v.stride = p;
      v.block = AB->block;
      v.owner = 0;

      view.vector = v;
      return view;
    }
}

/*
gsl_linalg_ldlt_band_decomp()
  Perform L D L^T decomposition of a symmetric banded positive
semi-definite matrix using lower triangle

Inputs: A    - (input) symmetric banded, positive semi-definite matrix
               (output) column 1 contains diagonal matrix D
                        columns 2:end contain lower triangle L factor
        work - workspace, size N

Return: success/error

Notes:
1) Based on algorithm 4.1.1 of Golub and Van Loan, Matrix Computations (4th ed).
*/

int
gsl_linalg_ldlt_band_decomp (gsl_matrix * A, gsl_vector * work)
{
  const size_t N = A->size1;
  const size_t ndiag = A->size2;

  if (ndiag > N)
    {
      GSL_ERROR ("invalid matrix dimensions", GSL_EBADLEN);
    }
  else if (work->size != N)
    {
      GSL_ERROR ("workspace does not match matrix", GSL_EBADLEN);
    }
  else
    {
      const size_t p = ndiag - 1; /* lower bandwidth */
      size_t i, j;
      double a00;
      gsl_vector_view v;

      /* check for quick return */
      if (N == 1 || ndiag == 1)
        return GSL_SUCCESS;

      /* special case first column */
      a00 = gsl_matrix_get(A, 0, 0);
      if (a00 == 0.0)
        {
          GSL_ERROR ("matrix is singular", GSL_EDOM);
        }

      v = gsl_matrix_subrow(A, 0, 1, p);
      gsl_blas_dscal(1.0 / a00, &v.vector);

      for (j = 1; j < N; ++j)
        {
          const size_t nz = (j <= p) ? 0 : j - p; /* number of zero entries at beginning of row j */
          double ajj = gsl_matrix_get(A, j, 0);   /* A(j,j) */
          double dval;
          gsl_vector_view w;

          v = symband_row_before_diag(A, j);
          w = gsl_vector_subvector(work, 0, v.vector.size);

          for (i = 0; i < v.vector.size; ++i)
            {
              double aii = gsl_matrix_get(A, i + nz, 0); /* A(i+nz,i+nz) */
              double aji = gsl_vector_get(&v.vector, i); /* A(j,i+nz) */
              gsl_vector_set(&w.vector, i, aji * aii);
            }

          gsl_blas_ddot(&v.vector, &w.vector, &dval);
          ajj -= dval;

          if (ajj == 0.0)
            {
              GSL_ERROR ("matrix is singular", GSL_EDOM);
            }

          gsl_matrix_set(A, j, 0, ajj);

          if (j < N - 1)
            {
              const size_t nleft = GSL_MIN(p, N - j - 1); /* number of non-zero entries beneath diagonal of column j */
              const double ajjinv = 1.0 / ajj;
              size_t idx;

              /* A(j+1:end, j) */
              v = gsl_matrix_subrow(A, j, 1, nleft);

              /* the following loop updates the subcolumn A(j+1:end,j),
               * represented by the block vector v below.
               *
               * [ already_updated |   0    |       0         ]
               * [ already_updated | A(j,j) |       0         ]
               * [        M        |   v    | not_yet_updated ]
               *
               * The vector v is updated as follows:
               *
               * v := 1/A(j,j) * (v - M*work)
               *
               * where:
               *
               * M    = A(j+1:N, 1:j-1)
               * work = A(j,1:j-1) .* diag(1:j-1)  ( computed above )
               *
               * We don't use a DGEMV operation due to the sparse structure
               * of M and work, so we loop through the rows of M and compute
               * dot products of the non-zero elements of M with the non-zero
               * elements of work
               */

              for (idx = 0; idx < nleft; ++idx)
                {
                  size_t k = idx + j + 1;                              /* current row index */
                  size_t nzk = (k <= p) ? 0 : k - p;                   /* number of zero entries at beginning of row k */
                  double * ptr = gsl_matrix_ptr(A, j, k - j);          /* pointer to A(k,j) entry */
                  gsl_vector_view mk0 = symband_row_before_diag(A, k); /* non-zero entries in row k not including diagonal */
                  gsl_vector_view mk = gsl_vector_subvector(&mk0.vector, 0, mk0.vector.size - idx - 1);

                  if (mk.vector.size > 0)
                    {
                      w = gsl_vector_subvector(work, nzk - nz, mk.vector.size);
                      gsl_blas_ddot(&mk.vector, &w.vector, &dval);
                      *ptr = ajjinv * (*ptr - dval);
                    }
                  else
                    {
                      *ptr *= ajjinv;
                    }
                }
            }
        }

      return GSL_SUCCESS;
    }
}

int
gsl_linalg_ldlt_band_solve (const gsl_matrix * LDLT,
                            const gsl_vector * b,
                            gsl_vector * x)
{
  if (LDLT->size1 != b->size)
    {
      GSL_ERROR ("matrix size must match b size", GSL_EBADLEN);
    }
  else if (LDLT->size1 != x->size)
    {
      GSL_ERROR ("matrix size must match solution size", GSL_EBADLEN);
    }
  else
    {
      int status;

      /* copy x <- b */
      gsl_vector_memcpy (x, b);

      status = gsl_linalg_ldlt_band_svx(LDLT, x);

      return status;
    }
}

int
gsl_linalg_ldlt_band_svx (const gsl_matrix * LDLT, gsl_vector * x)
{
  if (LDLT->size1 != x->size)
    {
      GSL_ERROR ("matrix size must match solution size", GSL_EBADLEN);
    }
  else
    {
      gsl_vector_const_view diag = gsl_matrix_const_column(LDLT, 0);

      /* solve for z using forward-substitution, L z = b */
      cblas_dtbsv(CblasColMajor, CblasLower, CblasNoTrans, CblasUnit,
                  (int) LDLT->size1, (int) (LDLT->size2 - 1), LDLT->data, LDLT->tda,
                  x->data, x->stride);

      /* solve for y, D y = z */
      gsl_vector_div(x, &diag.vector);

      /* perform back-substitution, L^T x = y */
      cblas_dtbsv(CblasColMajor, CblasLower, CblasTrans, CblasUnit,
                  (int) LDLT->size1, (int) (LDLT->size2 - 1), LDLT->data, LDLT->tda,
                  x->data, x->stride);

      return GSL_SUCCESS;
    }
}

/*
gsl_linalg_ldlt_band_unpack()
  Unpack symmetric banded format matrix LDLT into
larger matrix L and diagonal vector D
*/

int
gsl_linalg_ldlt_band_unpack (const gsl_matrix * LDLT, gsl_matrix * L, gsl_vector * D)
{
  const size_t N = LDLT->size1;

  if (N != L->size1)
    {
      GSL_ERROR("L matrix does not match LDLT dimensions", GSL_EBADLEN);
    }
  else if (L->size1 != L->size2)
    {
      GSL_ERROR("L matrix is not square", GSL_ENOTSQR);
    }
  else if (N != D->size)
    {
      GSL_ERROR("D vector does not match LDLT dimensions", GSL_EBADLEN);
    }
  else
    {
      const size_t p = LDLT->size2 - 1; /* lower bandwidth */
      gsl_vector_const_view diag = gsl_matrix_const_column(LDLT, 0);
      gsl_vector_view diagL = gsl_matrix_diagonal(L);
      size_t i;

      /* copy diagonal entries */
      gsl_vector_memcpy(D, &diag.vector);

      /* copy subdiagonals into L */
      for (i = 1; i <= p; ++i)
        {
          gsl_vector_const_view v = gsl_matrix_const_subcolumn(LDLT, i, 0, N - i);
          gsl_vector_view w = gsl_matrix_subdiagonal(L, i);
          gsl_vector_memcpy(&w.vector, &v.vector);
        }

      /* set main diagonal of L */
      gsl_vector_set_all(&diagL.vector, 1.0);

      /* zero out remaining subdiagonals */
      for (i = p + 1; i < N; ++i)
        {
          gsl_vector_view w = gsl_matrix_subdiagonal(L, i);
          gsl_vector_set_zero(&w.vector);
        }

      return GSL_SUCCESS;
    }
}
