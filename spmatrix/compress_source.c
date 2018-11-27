/* spmatrix/compress_source.c
 * 
 * Copyright (C) 2012, 2013, 2014, 2015, 2016, 2017, 2018 Patrick Alken
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

/*
gsl_spmatrix_csc()
  Create a sparse matrix in compressed column format

Inputs: T - sparse matrix in triplet format

Return: pointer to new matrix (should be freed when finished with it)
*/

TYPE (gsl_spmatrix) *
FUNCTION (gsl_spmatrix, csc) (const TYPE (gsl_spmatrix) * T)
{
  if (!GSL_SPMATRIX_ISCOO(T))
    {
      GSL_ERROR_NULL("matrix must be in triplet/COO format", GSL_EINVAL);
    }
  else
    {
      const int *Tj = T->p; /* column indices of triplet matrix */
      int *Cp;              /* column pointers of compressed column matrix */
      int *w;               /* copy of column pointers */
      TYPE (gsl_spmatrix) * m;
      size_t n, r;

      m = FUNCTION (gsl_spmatrix, alloc_nzmax) (T->size1, T->size2, T->nz, GSL_SPMATRIX_CSC);
      if (!m)
        return NULL;

      Cp = m->p;

      /* initialize column pointers to 0 */
      for (n = 0; n < m->size2 + 1; ++n)
        Cp[n] = 0;

      /*
       * compute the number of elements in each column:
       * Cp[j] = # non-zero elements in column j
       */
      for (n = 0; n < T->nz; ++n)
        Cp[Tj[n]]++;

      /* compute column pointers: p[j] = p[j-1] + nnz[j-1] */
      gsl_spmatrix_cumsum(m->size2, Cp);

      /* make a copy of the column pointers */
      w = m->work.work_int;
      for (n = 0; n < m->size2; ++n)
        w[n] = Cp[n];

      /* transfer data from triplet format to CSC */
      for (n = 0; n < T->nz; ++n)
        {
          int k = w[Tj[n]]++;
          m->i[k] = T->i[n];

          for (r = 0; r < MULTIPLICITY; ++r)
            m->data[MULTIPLICITY * k + r] = T->data[MULTIPLICITY * n + r];
        }

      m->nz = T->nz;

      return m;
    }
}

/* XXX deprecated function */
TYPE (gsl_spmatrix) *
FUNCTION (gsl_spmatrix, compcol) (const TYPE (gsl_spmatrix) * T)
{
  return FUNCTION (gsl_spmatrix, csc) (T);
}

/* XXX deprecated function */
TYPE (gsl_spmatrix) *
FUNCTION (gsl_spmatrix, ccs) (const TYPE (gsl_spmatrix) * T)
{
  return FUNCTION (gsl_spmatrix, csc) (T);
}

/*
gsl_spmatrix_crs()
  Create a sparse matrix in compressed row format

Inputs: T - sparse matrix in triplet format

Return: pointer to new matrix (should be freed when finished with it)
*/

TYPE (gsl_spmatrix) *
FUNCTION (gsl_spmatrix, csr) (const TYPE (gsl_spmatrix) * T)
{
  if (!GSL_SPMATRIX_ISCOO(T))
    {
      GSL_ERROR_NULL("matrix must be in triplet/COO format", GSL_EINVAL);
    }
  else
    {
      const int *Ti = T->i; /* row indices of triplet matrix */
      int *Cp;              /* row pointers of compressed row matrix */
      int *w;               /* copy of column pointers */
      TYPE (gsl_spmatrix) * m;
      size_t n, r;

      m = FUNCTION (gsl_spmatrix, alloc_nzmax) (T->size1, T->size2, T->nz, GSL_SPMATRIX_CSR);
      if (!m)
        return NULL;

      Cp = m->p;

      /* initialize row pointers to 0 */
      for (n = 0; n < m->size1 + 1; ++n)
        Cp[n] = 0;

      /*
       * compute the number of elements in each row:
       * Cp[i] = # non-zero elements in row i
       */
      for (n = 0; n < T->nz; ++n)
        Cp[Ti[n]]++;

      /* compute row pointers: p[i] = p[i-1] + nnz[i-1] */
      gsl_spmatrix_cumsum(m->size1, Cp);

      /* make a copy of the row pointers */
      w = m->work.work_int;
      for (n = 0; n < m->size1; ++n)
        w[n] = Cp[n];

      /* transfer data from triplet format to CSR */
      for (n = 0; n < T->nz; ++n)
        {
          int k = w[Ti[n]]++;
          m->i[k] = T->p[n];

          for (r = 0; r < MULTIPLICITY; ++r)
            m->data[MULTIPLICITY * k + r] = T->data[MULTIPLICITY * n + r];
        }

      m->nz = T->nz;

      return m;
    }
}

/* XXX deprecated function */
TYPE (gsl_spmatrix) *
FUNCTION (gsl_spmatrix, crs) (const TYPE (gsl_spmatrix) * T)
{
  return FUNCTION (gsl_spmatrix, csr) (T);
}

TYPE (gsl_spmatrix) *
FUNCTION (gsl_spmatrix, compress) (const TYPE (gsl_spmatrix) * T, const int sptype)
{
  if (sptype == GSL_SPMATRIX_CSC)
    {
      return FUNCTION (gsl_spmatrix, csc) (T);
    }
  else if (sptype == GSL_SPMATRIX_CSR)
    {
      return FUNCTION (gsl_spmatrix, csr) (T);
    }
  else if (sptype == GSL_SPMATRIX_COO)
    {
      /* make a copy of T */
      TYPE (gsl_spmatrix) * A = FUNCTION (gsl_spmatrix, alloc_nzmax) (T->size1, T->size2, T->nz, sptype);
      FUNCTION (gsl_spmatrix, memcpy) (A, T);
      return A;
    }
  else
    {
      GSL_ERROR_NULL ("unknown sparse matrix format", GSL_EINVAL);
    }
}
