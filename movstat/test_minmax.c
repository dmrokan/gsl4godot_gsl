/* movstat/test_minmax.c
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

/* compute filtered data by explicitely constructing window, sorting it and finding min/max */
int
slow_minmax(const gsl_movstat_end_t etype, const gsl_vector * x, gsl_vector * y_min, gsl_vector * y_max,
            const int H, const int J)
{
  const int n = (int) x->size;
  const int K = H + J + 1;
  double *window = malloc(K * sizeof(double));
  int i;

  for (i = 0; i < n; ++i)
    {
      int wsize;
      gsl_vector_view v;
      double min, max;

      wsize = test_window(etype, i, H, J, x, window);
      v = gsl_vector_view_array(window, wsize);
      gsl_vector_minmax(&v.vector, &min, &max);

      gsl_vector_set(y_min, i, min);
      gsl_vector_set(y_max, i, max);
    }

  free(window);

  return GSL_SUCCESS;
}

/* test alternating sequence [a,b,a,b,...] input */
static void
test_minmax_alt(const double tol, const size_t n, const int H, const int J,
                const gsl_movstat_end_t endtype)
{
  const double a = 5.0;
  const double b = -5.0;
  gsl_vector * x = gsl_vector_alloc(n);
  gsl_vector * y_min = gsl_vector_alloc(n);
  gsl_vector * y_max = gsl_vector_alloc(n);
  gsl_movstat_workspace * w = gsl_movstat_alloc2(H, J);
  size_t i;

  for (i = 0; i < n; ++i)
    {
      if (i % 2 == 0)
        gsl_vector_set(x, i, a);
      else
        gsl_vector_set(x, i, b);
    }

  /* compute moving min/max */
  gsl_movstat_minmax(endtype, x, y_min, y_max, w);

  for (i = 0; i < n; ++i)
    {
      double ymini = gsl_vector_get(y_min, i);
      double ymaxi = gsl_vector_get(y_max, i);

      gsl_test_rel(ymini, b, tol, "test_minmax_alt: min[%zu] endtype=%d n=%zu H=%d J=%d", i, endtype, n, H, J);
      gsl_test_rel(ymaxi, a, tol, "test_minmax_alt: max[%zu] endtype=%d n=%zu H=%d J=%d", i, endtype, n, H, J);
    }

  gsl_vector_free(x);
  gsl_vector_free(y_min);
  gsl_vector_free(y_max);
  gsl_movstat_free(w);
}

/* test noisy sine wave input */
static void
test_minmax_sine(const double tol, const size_t n, const int H, const int J,
                 const gsl_movstat_end_t endtype, gsl_rng * rng_p)
{
  gsl_vector * x = gsl_vector_alloc(n);
  gsl_vector * y_min = gsl_vector_alloc(n);
  gsl_vector * y_max = gsl_vector_alloc(n);
  gsl_vector * z_min = gsl_vector_alloc(n);
  gsl_vector * z_max = gsl_vector_alloc(n);
  gsl_movstat_workspace * w = gsl_movstat_alloc2(H, J);
  size_t i;

  /* construct noisy sine signal */
  test_noisy_sine(0.5, x, rng_p);

  /* compute moving min/max with different algorithms */
  gsl_movstat_minmax(endtype, x, y_min, y_max, w);
  slow_minmax(endtype, x, z_min, z_max, H, J);

  for (i = 0; i < n; ++i)
    {
      double ymini = gsl_vector_get(y_min, i);
      double ymaxi = gsl_vector_get(y_max, i);
      double zmini = gsl_vector_get(z_min, i);
      double zmaxi = gsl_vector_get(z_max, i);

      gsl_test_rel(ymini, zmini, tol, "test_minmax_sine: min endtype=%d n=%zu H=%d J=%d", endtype, n, H, J);
      gsl_test_rel(ymaxi, zmaxi, tol, "test_minmax_sine: max endtype=%d n=%zu H=%d J=%d", endtype, n, H, J);
    }

  gsl_vector_free(x);
  gsl_vector_free(y_min);
  gsl_vector_free(y_max);
  gsl_vector_free(z_min);
  gsl_vector_free(z_max);
  gsl_movstat_free(w);
}

static void
test_minmax(void)
{
  gsl_rng *rng_p = gsl_rng_alloc(gsl_rng_default);

  /* alternating input */

  test_minmax_alt(GSL_DBL_EPSILON, 1000, 7, 7, GSL_MOVSTAT_END_PADZERO);
  test_minmax_alt(GSL_DBL_EPSILON, 1000, 5, 2, GSL_MOVSTAT_END_PADZERO);
  test_minmax_alt(GSL_DBL_EPSILON, 500, 1, 3, GSL_MOVSTAT_END_PADZERO);

  /* noisy sine wave input */

  test_minmax_sine(GSL_DBL_EPSILON, 1000, 5, 7, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_minmax_sine(GSL_DBL_EPSILON, 2000, 0, 2, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_minmax_sine(GSL_DBL_EPSILON, 500, 3, 0, GSL_MOVSTAT_END_PADZERO, rng_p);

  test_minmax_sine(GSL_DBL_EPSILON, 500, 5, 7, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_minmax_sine(GSL_DBL_EPSILON, 1000, 10, 20, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_minmax_sine(GSL_DBL_EPSILON, 1000, 3, 3, GSL_MOVSTAT_END_PADVALUE, rng_p);

  test_minmax_sine(GSL_DBL_EPSILON, 500, 5, 7, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_minmax_sine(GSL_DBL_EPSILON, 1000, 10, 20, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_minmax_sine(GSL_DBL_EPSILON, 1000, 3, 3, GSL_MOVSTAT_END_TRUNCATE, rng_p);

  gsl_rng_free(rng_p);
}
