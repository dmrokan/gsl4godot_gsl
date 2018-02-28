/* movstat/test_median.c
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

/* compute filtered data by explicitely constructing window, sorting it and finding median */
int
slow_median(const gsl_movstat_end_t etype, const gsl_vector * x, gsl_vector * y,
            const int H, const int J)
{
  const int n = (int) x->size;
  const int K = H + J + 1;
  double *window = malloc(K * sizeof(double));
  int i, j;

  for (i = 0; i < n; ++i)
    {
      double yi;
      int idx1, idx2, wsize;

      if (etype == GSL_MOVSTAT_END_TRUNCATE)
        {
          idx1 = GSL_MAX(i - H, 0);
          idx2 = GSL_MIN(i + J, n - 1);
        }
      else
        {
          idx1 = i - H;
          idx2 = i + J;
        }

      wsize = idx2 - idx1 + 1;

      /* fill sliding window */
      for (j = idx1; j <= idx2; ++j)
        {
          int widx = j - idx1;

          if (j < 0)
            {
              /* initial condition */
              if (etype == GSL_MOVSTAT_END_PADZERO)
                window[widx] = 0.0;
              else if (etype == GSL_MOVSTAT_END_PADVALUE)
                window[widx] = gsl_vector_get(x, 0);
            }
          else if (j >= (int) n)
            {
              if (etype == GSL_MOVSTAT_END_PADZERO)
                window[widx] = 0.0;
              else if (etype == GSL_MOVSTAT_END_PADVALUE)
                window[widx] = gsl_vector_get(x, n - 1);
            }
          else
            {
              window[widx] = gsl_vector_get(x, j);
            }
        }

      yi = median_find(wsize, window);
      gsl_vector_set(y, i, yi);
    }

  free(window);

  return GSL_SUCCESS;
}

/* test root sequence */
static void
test_median_root(const double tol, const size_t n, const size_t k, const gsl_movstat_end_t etype)
{
  gsl_movstat_median_workspace *w = gsl_movstat_median_alloc(k);
  gsl_vector *x = gsl_vector_alloc(n);
  gsl_vector *y = gsl_vector_alloc(n);
  char buf[2048];
  size_t i;

  /* test a root sequence (square input): x = [zero one zero] */
  gsl_vector_set_all(x, 0.0);

  for (i = n / 3; i <= n / 2; ++i)
    gsl_vector_set(x, i, 1.0);

  /* compute y = median(x) and test y = x */
  gsl_movstat_median(etype, x, y, w);

  sprintf(buf, "n=%zu k=%zu SMF root sequence", n, k);
  compare_vectors(tol, y, x, buf);

  gsl_vector_free(x);
  gsl_vector_free(y);
  gsl_movstat_median_free(w);
}

/* symmetric window */
static void
test_median_symmetric(const double tol, const size_t n, const size_t K,
                      const gsl_movstat_end_t etype, gsl_rng *rng_p)
{
  gsl_movstat_median_workspace *w = gsl_movstat_median_alloc(K);
  gsl_vector *x = gsl_vector_alloc(n);
  gsl_vector *y = gsl_vector_alloc(n);
  gsl_vector *z = gsl_vector_alloc(n);
  char buf[2048];

  /* test median filter with random input */

  random_vector(x, rng_p);

  /* y = median(x) */
  gsl_movstat_median(etype, x, y, w);

  /* z = median(x) with slow algorithm */
  slow_median(etype, x, z, w->H, w->J);

  /* test y = z */
  sprintf(buf, "n=%zu k=%zu endtype=%u median symmetric random slow test", n, K, etype);
  compare_vectors(tol, z, y, buf);

  /* z = median(x) in-place */
  gsl_vector_memcpy(z, x);
  gsl_movstat_median(etype, z, z, w);

  sprintf(buf, "n=%zu k=%zu endtype=%u median symmetric random in-place", n, K, etype);
  compare_vectors(tol, z, y, buf);

  gsl_vector_free(x);
  gsl_vector_free(y);
  gsl_vector_free(z);
  gsl_movstat_median_free(w);
}

/* non-symmetric window */
void
test_median_nonsymmetric(const gsl_movstat_end_t etype)
{
  const double tol = GSL_DBL_EPSILON;
  const size_t n = 20;
  const double data[] = { 4.38, 3.81, 7.65, 7.95, 1.86, 4.89, 4.45, 6.46, 0.09, 7.54,
                          2.76, 6.79, 6.55, 1.62, 1.18, 4.98, 9.59, 3.40, 5.85, 2.23 };
  const double expected_y1[] = { 0.0, 1.905, 2.835, 4.095, 4.415, 4.67, 4.67, 5.675, 4.67, 4.67,
                                 5.675, 5.455, 4.61, 3.87, 5.765, 4.19, 5.415, 4.19, 2.815, 2.815 }; /* H=5, J=2 */
  const double expected_y2[] = { 4.38, 4.38, 4.89, 4.89, 4.89, 4.45, 4.89, 4.45, 6.46, 6.55,
                                 6.55, 2.76, 4.98, 4.98, 3.40, 4.98, 4.98, 3.40, 2.23, 0.0 };        /* H=1, J=3 */
  const double expected_y3[] = { 0.0, 0.0, 3.81, 4.38, 4.38, 4.89, 4.89, 4.89, 4.45, 4.89,
                                 4.45, 6.46, 6.55, 6.55, 2.76, 4.98, 4.98, 3.40, 4.98, 4.98 };       /* H=4, J=0 */
  gsl_movstat_median_workspace *w;
  gsl_vector *y = gsl_vector_alloc(n);
  gsl_vector_const_view x = gsl_vector_const_view_array(data, n);
  gsl_vector_const_view y1 = gsl_vector_const_view_array(expected_y1, n);
  gsl_vector_const_view y2 = gsl_vector_const_view_array(expected_y2, n);
  gsl_vector_const_view y3 = gsl_vector_const_view_array(expected_y3, n);
  size_t H, J;

  /* test 1: H = 5, J = 2 */
  H = 5;
  J = 2;
  w = gsl_movstat_median_alloc2(H, J);

  /* y = median(x) */
  gsl_movstat_median(etype, &x.vector, y, w);
  compare_vectors(tol, y, &y1.vector, "SMF nonsymmetric H=5 J=2");

  gsl_movstat_median_free(w);

  /* test 2: H = 1, J = 3 */
  H = 1;
  J = 3;
  w = gsl_movstat_median_alloc2(H, J);

  /* y = median(x) */
  gsl_movstat_median(etype, &x.vector, y, w);
  compare_vectors(tol, y, &y2.vector, "SMF nonsymmetric H=1 J=3");

  gsl_movstat_median_free(w);

  /* test 3: H = 4, J = 0 */
  H = 4;
  J = 0;
  w = gsl_movstat_median_alloc2(H, J);

  /* y = median(x) */
  gsl_movstat_median(etype, &x.vector, y, w);
  compare_vectors(tol, y, &y3.vector, "SMF nonsymmetric H=4 J=0");

  gsl_movstat_median_free(w);

  gsl_vector_free(y);
}

/* nonsymmetric window */
static void
test_median_nonsymmetric2(const double tol, const size_t n, const size_t H, const size_t J,
                          const gsl_movstat_end_t etype, gsl_rng *rng_p)
{
  gsl_movstat_median_workspace *w = gsl_movstat_median_alloc2(H, J);
  gsl_vector *x = gsl_vector_alloc(n);
  gsl_vector *y = gsl_vector_alloc(n);
  gsl_vector *z = gsl_vector_alloc(n);
  char buf[2048];

  /* test median filter with random input */

  random_vector(x, rng_p);

  /* y = median(x) */
  gsl_movstat_median(etype, x, y, w);

  /* z = median(x) with slow algorithm */
  slow_median(etype, x, z, w->H, w->J);

  /* test y = z */
  sprintf(buf, "n=%zu H=%zu J=%zu endtype=%u median nonsymmetric random slow test", n, H, J, etype);
  compare_vectors(tol, z, y, buf);

  /* z = median(x) in-place */
  gsl_vector_memcpy(z, x);
  gsl_movstat_median(etype, z, z, w);

  sprintf(buf, "n=%zu H=%zu J=%zu endtype=%u median nonsymmetric random in-place", n, H, J, etype);
  compare_vectors(tol, z, y, buf);

  gsl_vector_free(x);
  gsl_vector_free(y);
  gsl_vector_free(z);
  gsl_movstat_median_free(w);
}

void
test_median(void)
{
  gsl_rng *rng_p = gsl_rng_alloc(gsl_rng_default);

  test_median_root(GSL_DBL_EPSILON, 1000, 3, GSL_MOVSTAT_END_PADZERO);
  test_median_root(GSL_DBL_EPSILON, 200, 15, GSL_MOVSTAT_END_PADVALUE);
  test_median_root(GSL_DBL_EPSILON, 100, 5, GSL_MOVSTAT_END_TRUNCATE);

  test_median_symmetric(GSL_DBL_EPSILON, 1000, 3, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_median_symmetric(GSL_DBL_EPSILON, 100, 301, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_median_symmetric(GSL_DBL_EPSILON, 5000, 17, GSL_MOVSTAT_END_PADZERO, rng_p);

  test_median_symmetric(GSL_DBL_EPSILON, 2000, 31, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_median_symmetric(GSL_DBL_EPSILON, 200, 501, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_median_symmetric(GSL_DBL_EPSILON, 1000, 501, GSL_MOVSTAT_END_PADVALUE, rng_p);

#if 0
  test_median_symmetric(GSL_DBL_EPSILON, 400, 11, GSL_MOVSTAT_END_TRUNCATE, rng_p);
  test_median_symmetric(GSL_DBL_EPSILON, 200, 5, GSL_MOVSTAT_END_TRUNCATE, rng_p);
#endif

  test_median_nonsymmetric2(GSL_DBL_EPSILON, 200, 5, 4, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_median_nonsymmetric2(GSL_DBL_EPSILON, 2000, 7, 10, GSL_MOVSTAT_END_PADZERO, rng_p);
  test_median_nonsymmetric2(GSL_DBL_EPSILON, 500, 0, 1, GSL_MOVSTAT_END_PADZERO, rng_p);

  test_median_nonsymmetric2(GSL_DBL_EPSILON, 200, 1, 2, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_median_nonsymmetric2(GSL_DBL_EPSILON, 1000, 20, 40, GSL_MOVSTAT_END_PADVALUE, rng_p);
  test_median_nonsymmetric2(GSL_DBL_EPSILON, 500, 2, 30, GSL_MOVSTAT_END_PADVALUE, rng_p);

  test_median_nonsymmetric2(GSL_DBL_EPSILON, 200, 3, 1, GSL_MOVSTAT_END_TRUNCATE, rng_p);

  gsl_rng_free(rng_p);
}
