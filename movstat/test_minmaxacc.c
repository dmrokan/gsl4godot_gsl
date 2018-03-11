/* movstat/test_minmaxacc.c
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

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_movstat.h>

static void
test_minmaxacc_n(const double tol, const int n, const int k)
{
  gsl_rng *rng_p = gsl_rng_alloc(gsl_rng_default);
  gsl_movstat_minmaxacc_workspace *w = gsl_movstat_minmaxacc_alloc(k);
  gsl_vector *x = gsl_vector_alloc(n);
  gsl_vector_view wv;
  int i;

  /* generate random integer vector in [-100,100] */
  random_vector_int(-100, 100, x, rng_p);

  for (i = 0; i < n; ++i)
    {
      double xi = gsl_vector_get(x, i);
      double min, max;
      double min_exact, max_exact;
      int idx = i - k + 1;
      size_t wlen;

      gsl_movstat_minmaxacc_insert(xi, w);
      gsl_movstat_minmaxacc_minmax(&min, &max, w);

      if (idx < 0)
        {
          wlen = i + 1;
          wv = gsl_vector_subvector(x, 0, wlen);
        }
      else
        {
          wlen = k;
          wv = gsl_vector_subvector(x, idx, wlen);
        }

      gsl_vector_minmax(&wv.vector, &min_exact, &max_exact);

      gsl_test_rel(min, min_exact, tol, "test_minmaxacc_n min n=%zu k=%zu i=%zu", n, k, i);
      gsl_test_rel(max, max_exact, tol, "test_minmaxacc_n max n=%zu k=%zu i=%zu", n, k, i);
    }

  gsl_movstat_minmaxacc_free(w);
  gsl_rng_free(rng_p);
  gsl_vector_free(x);
}

static void
test_minmaxacc(void)
{
  test_minmaxacc_n(GSL_DBL_EPSILON, 10, 3);
  test_minmaxacc_n(GSL_DBL_EPSILON, 100, 3);
  test_minmaxacc_n(GSL_DBL_EPSILON, 500000, 40);
  test_minmaxacc_n(GSL_DBL_EPSILON, 100000, 1001);
}
