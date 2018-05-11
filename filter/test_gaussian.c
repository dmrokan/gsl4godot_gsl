/* filter/test_gaussian.c
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
#include <gsl/gsl_randist.h>

static void
test_gaussian_deriv(const double alpha, const size_t n, const size_t K)
{
#if 0
  const double f_low = 1.0;
  const double f_high = 50.0;
  const double gamma = 2.0 * M_PI / (n - 1.0);
  const double dt = 1.0 / (n - 1.0);
  gsl_vector *x = gsl_vector_alloc(n);
  gsl_vector *dx = gsl_vector_alloc(n);
  gsl_vector *y1 = gsl_vector_alloc(n);
  gsl_vector *y2 = gsl_vector_alloc(n);
  gsl_filter_gaussian_workspace *w = gsl_filter_gaussian_alloc(K);
  size_t i;

  /* make input signal composed of two sine waves at different frequencies */
  for (i = 0; i < n; ++i)
    {
      double xi = sin(gamma * f_low * i) + sin(gamma * f_high * i);
      double dxi = gamma * f_low * cos(gamma * f_low * i) +
                   gamma * f_high * cos(gamma * f_high * i);

      gsl_vector_set(x, i, xi);
      gsl_vector_set(dx, i, dxi);
    }

  /* compute y1 = G * dx(t)/dt */
  gsl_filter_gaussian(alpha, 0, dx, y1, w);

  /* compute y2 = dG/dt * x(t) */
  gsl_filter_gaussian(alpha, 1, x, y2, w);

  for (i = 0; i < n; ++i)
    {
      printf("%zu %.12e %.12e %.12e %.12e\n",
             i,
             gsl_vector_get(x, i),
             gsl_vector_get(dx, i),
             gsl_vector_get(y1, i),
             gsl_vector_get(y2, i));
    }

  gsl_vector_free(x);
  gsl_vector_free(dx);
  gsl_vector_free(y1);
  gsl_vector_free(y2);
  gsl_filter_gaussian_free(w);
#endif
}

static void
test_gaussian(void)
{
  test_gaussian_deriv(2.5, 1000, 15);
}
