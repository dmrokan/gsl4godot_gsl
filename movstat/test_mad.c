/* movstat/test_mad.c
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

static void
test_mad_symmetric(void)
{
  const double tol = GSL_DBL_EPSILON;

  /* test movmad on a dataset with a 3-point sliding window */
  {
    const size_t n = 8;
    const size_t k = 3;
    const double data[] = { -1.0, 5.7, 3.4, 1.1, 9.5, -23.7, -5.6, 0.2 };
    const double expected_median[] = { 0.0, 3.4, 3.4, 3.4, 1.1, -5.6, -5.6, 0.0 };
    const double expected_MAD[] = { 1.0, 2.3, 2.3, 2.3, 8.4, 15.1, 5.8, 0.2 };
    gsl_vector_const_view x = gsl_vector_const_view_array(data, n);
    gsl_movstat_mad_workspace *w = gsl_movstat_mad_alloc(k);
    gsl_vector *xmedian = gsl_vector_alloc(n);
    gsl_vector *xmad = gsl_vector_alloc(n);
    size_t i;

    gsl_movstat_mad(GSL_MOVSTAT_END_PADZERO, &x.vector, xmedian, xmad, w);

    for (i = 0; i < n; ++i)
      {
        double xmedi = gsl_vector_get(xmedian, i);
        double xmadi = gsl_vector_get(xmad, i);

        gsl_test_rel(xmedi, expected_median[i], tol, "movmad symm[1] median i=%zu", i);
        gsl_test_rel(xmadi, expected_MAD[i], tol, "movmad symm[1] MAD i=%zu", i);
      }

    gsl_vector_free(xmedian);
    gsl_vector_free(xmad);
    gsl_movstat_mad_free(w);
  }

  /* test movmad on the dataset x = [1:100] with alternating signs with a 7-point window */
  {
    const size_t n = 100;
    const size_t k = 7;
    gsl_vector *x = gsl_vector_alloc(n);
    gsl_vector *xmedian = gsl_vector_alloc(n);
    gsl_vector *xmad = gsl_vector_alloc(n);
    gsl_movstat_mad_workspace *w = gsl_movstat_mad_alloc(k);
    size_t i;
    double s = 1.0;

    /* x = [1 -2 3 -4 ... 99 -100] */
    for (i = 0; i < n; ++i)
      {
        gsl_vector_set(x, i, s * (i + 1.0));
        s *= -1.0;
      }

    gsl_movstat_mad(GSL_MOVSTAT_END_PADZERO, x, xmedian, xmad, w);

    /* test median results (compared with MATLAB medfilt1) */

    gsl_test_abs(gsl_vector_get(xmedian, 0), 0.0, tol, "movmad symm[2] median i=%zu", 0);
    gsl_test_abs(gsl_vector_get(xmedian, 1), 0.0, tol, "movmad symm[2] median i=%zu", 1);
    gsl_test_abs(gsl_vector_get(xmedian, 2), 0.0, tol, "movmad symm[2] median i=%zu", 2);

    gsl_test_abs(gsl_vector_get(xmedian, n - 1), 0.0, tol, "movmad symm[2] median i=%zu", n - 1);
    gsl_test_abs(gsl_vector_get(xmedian, n - 2), 0.0, tol, "movmad symm[2] median i=%zu", n - 2);
    gsl_test_abs(gsl_vector_get(xmedian, n - 3), 0.0, tol, "movmad symm[2] median i=%zu", n - 3);

    s = 1.0;
    for (i = 3; i < n - 3; ++i)
      {
        double xmedi = gsl_vector_get(xmedian, i);
        gsl_test_rel(xmedi, s * (i - 2.0), tol, "movmad symm[2] median i=%zu", i);
        s *= -1.0;
      }

    /* test MAD results (compared with MATLAB movmad) */

    gsl_test_abs(gsl_vector_get(xmad, 0), 1.0, tol, "movmad symm[2] MAD i=%zu", 0);
    gsl_test_abs(gsl_vector_get(xmad, 1), 2.0, tol, "movmad symm[2] MAD i=%zu", 1);
    gsl_test_abs(gsl_vector_get(xmad, 2), 3.0, tol, "movmad symm[2] MAD i=%zu", 2);
    gsl_test_abs(gsl_vector_get(xmad, 3), 4.0, tol, "movmad symm[2] MAD i=%zu", 3);
    gsl_test_abs(gsl_vector_get(xmad, 4), 5.0, tol, "movmad symm[2] MAD i=%zu", 4);

    gsl_test_abs(gsl_vector_get(xmad, n - 1), 97.0, tol, "movmad symm[2] MAD i=%zu", n - 1);
    gsl_test_abs(gsl_vector_get(xmad, n - 2), 97.0, tol, "movmad symm[2] MAD i=%zu", n - 2);
    gsl_test_abs(gsl_vector_get(xmad, n - 3), 97.0, tol, "movmad symm[2] MAD i=%zu", n - 3);

    for (i = 5; i < n - 3; ++i)
      {
        double xmadi = gsl_vector_get(xmad, i);
        gsl_test_rel(xmadi, 6.0, tol, "movmad symm[2] MAD i=%zu", i);
      }

    gsl_vector_free(x);
    gsl_vector_free(xmedian);
    gsl_vector_free(xmad);
    gsl_movstat_mad_free(w);
  }
}

static void
test_mad_nonsymmetric(void)
{
  const double tol = 1.0e-12;
  const size_t n = 20;
  const size_t H = 5;
  const size_t J = 2;
  const double data[] = { 4.38, 3.81, 7.65, 7.95, 1.86, 4.89, 4.45, 6.46, 0.09, 7.54,
                          2.76, 6.79, 6.55, 1.62, 1.18, 4.98, 9.59, 3.40, 5.85, 2.23 };
  const double expected_median[] = { 0.0, 1.905, 2.835, 4.095, 4.415, 4.67, 4.67, 5.675, 4.67, 4.67,
                                     5.675, 5.455, 4.61, 3.87, 5.765, 4.19, 5.415, 4.19, 2.815, 2.815 };
  const double expected_MAD[] = { 0.0, 1.905, 2.835, 2.895, 1.58, 1.325, 2.30, 1.92, 2.36, 2.015,
                                  1.17, 1.71, 2.555, 2.685, 2.39, 2.465, 1.695, 2.16, 1.90, 2.49 };
  gsl_vector_const_view x = gsl_vector_const_view_array(data, n);
  gsl_movstat_mad_workspace *w = gsl_movstat_mad_alloc2(H, J);
  gsl_vector *xmedian = gsl_vector_alloc(n);
  gsl_vector *xmad = gsl_vector_alloc(n);
  size_t i;

  gsl_movstat_mad(GSL_MOVSTAT_END_PADZERO, &x.vector, xmedian, xmad, w);

  for (i = 0; i < n; ++i)
    {
      double xmedi = gsl_vector_get(xmedian, i);
      double xmadi = gsl_vector_get(xmad, i);

      gsl_test_rel(xmedi, expected_median[i], tol, "movmad[3] median i=%zu", i);
      gsl_test_rel(xmadi, expected_MAD[i], tol, "movmad[3] MAD i=%zu", i);
    }

  gsl_vector_free(xmedian);
  gsl_vector_free(xmad);
  gsl_movstat_mad_free(w);
}

static void
test_mad(void)
{
  test_mad_symmetric();
  test_mad_nonsymmetric();
}
