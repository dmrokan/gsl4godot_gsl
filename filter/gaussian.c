/* filter/gaussian.c
 *
 * Gaussian smoothing filters
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
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_filter.h>
#include <gsl/gsl_poly.h>

/* maximum derivative order allowed for Gaussian filter */
#define GSL_FILTER_GAUSSIAN_MAX_ORDER     10

/*
gsl_filter_gaussian_alloc()
  Allocate a workspace for Gaussian filtering.

Inputs: K - number of samples in window; if even, it is rounded up to
            the next odd, to have a symmetric window

Return: pointer to workspace
*/

gsl_filter_gaussian_workspace *
gsl_filter_gaussian_alloc(const size_t K)
{
  gsl_filter_gaussian_workspace *w;

  w = calloc(1, sizeof(gsl_filter_gaussian_workspace));
  if (w == 0)
    {
      GSL_ERROR_NULL ("failed to allocate space for workspace", GSL_ENOMEM);
    }

  w->H = K / 2;
  w->K = 2 * w->H + 1;

  w->kernel = malloc(w->K * sizeof(double));
  if (w->kernel == 0)
    {
      gsl_filter_gaussian_free(w);
      GSL_ERROR_NULL ("failed to allocate space for kernel", GSL_ENOMEM);
      return NULL;
    }

  return w;
}

void
gsl_filter_gaussian_free(gsl_filter_gaussian_workspace * w)
{
  if (w->kernel)
    free(w->kernel);

  free(w);
}

/*
gsl_filter_gaussian()
  Apply a Gaussian filter to an input vector:

G_{sigma}(x) = exp [ -x^2 / (2 sigma^2) ]

Inputs: alpha - number of standard deviations to include in Gaussian kernel
        order - derivative order of Gaussian
        x     - input vector, size n
        y     - (output) filtered vector, size n
        w     - workspace

Notes:
1) If alpha = 3, then the Gaussian kernel will be a Gaussian of +/- 3 standard deviations
*/

int
gsl_filter_gaussian(const gsl_filter_end_t endtype, const double alpha, const size_t order, const gsl_vector * x,
                    gsl_vector * y, gsl_filter_gaussian_workspace * w)
{
  if (x->size != y->size)
    {
      GSL_ERROR("input and output vectors must have same length", GSL_EBADLEN);
    }
  else if (alpha <= 0.0)
    {
      GSL_ERROR("alpha must be positive", GSL_EDOM);
    }
  else
    {
      const size_t n = x->size;
      gsl_vector_view kernel = gsl_vector_view_array(w->kernel, w->K);
      double *window = malloc(w->K * sizeof(double));
      size_t i;

      /* construct Gaussian kernel of length K */
      gsl_filter_gaussian_kernel(alpha, order, 1, &kernel.vector);

      for (i = 0; i < n; ++i)
        {
          size_t wsize = gsl_movstat_fill(endtype, x, i, w->H, w->H, window);
          double sum = 0.0;
          size_t j;

          for (j = 0; j < wsize; ++j)
            sum += window[j] * w->kernel[wsize - j - 1];

          gsl_vector_set(y, i, sum);
        }

      free(window);

      return GSL_SUCCESS;
    }
}

/*
gsl_filter_gaussian_kernel()
  Construct Gaussian kernel with given sigma and order

Inputs: alpha     - number of standard deviations to include in window
        order     - kernel order (0 = gaussian, 1 = first derivative, ...)
        normalize - normalize so sum(G) = 1
        kernel    - (output) Gaussian kernel

Return: success/error

Notes:
1) If alpha = 3, then the output kernel will contain a Gaussian with +/- 3 standard deviations
*/

int
gsl_filter_gaussian_kernel(const double alpha, const size_t order, const int normalize, gsl_vector * kernel)
{
  const size_t N = kernel->size;

  if (alpha <= 0.0)
    {
      GSL_ERROR("alpha must be positive", GSL_EDOM);
    }
  else if (order > GSL_FILTER_GAUSSIAN_MAX_ORDER)
    {
      GSL_ERROR("derivative order is too large", GSL_EDOM);
    }
  else
    {
      const double half = 0.5 * (N - 1.0); /* (N - 1) / 2 */
      double sum = 0.0;
      size_t i;

      /* check for quick return */
      if (N == 1)
        {
          if (order == 0)
            gsl_vector_set(kernel, 0, 1.0);
          else
            gsl_vector_set(kernel, 0, 0.0);

          return GSL_SUCCESS;
        }

      for (i = 0; i < N; ++i)
        {
          double xi = ((double)i - half) / half;
          double yi = alpha * xi;
          double gi = exp(-0.5 * yi * yi);

          gsl_vector_set(kernel, i, gi);
          sum += gi;
        }

      /* normalize so sum(kernel) = 1 */
      if (normalize)
        gsl_vector_scale(kernel, 1.0 / sum);

      if (order > 0)
        {
          const double beta = -0.5 * alpha * alpha;
          double q[GSL_FILTER_GAUSSIAN_MAX_ORDER + 1];
          size_t k;

          /*
           * Need to calculate derivatives of the Gaussian window; define
           *
           * w(n) = C * exp [ p(n) ]
           *
           * p(n) = beta * n^2
           * beta = -1/2 * ( alpha / ((N-1)/2) )^2
           *
           * Then:
           *
           * d^k/dn^k w(n) = q_k(n) * w(n)
           *
           * where q_k(n) is a degree-k polynomial in n, which satisfies:
           *
           * q_k(n) = d/dn q_{k-1}(n) + q_{k-1}(n) * dp(n)/dn
           * q_0(n) = 1 / half^{order}
           */

          /* initialize q_0(n) = 1 / half^{order} */
          q[0] = 1.0 / gsl_pow_uint(half, order);
          for (i = 1; i <= GSL_FILTER_GAUSSIAN_MAX_ORDER; ++i)
            q[i] = 0.0;

          /* loop through derivative orders and calculate q_k(n) for k = 1,...,order */
          for (k = 1; k <= order; ++k)
            {
              double qm1 = q[0];

              q[0] = q[1];
              for (i = 1; i <= k; ++i)
                {
                  double tmp = q[i];
                  q[i] = (i + 1.0) * q[i + 1] + /* d/dn q_{k-1} */
                         2.0 * beta * qm1;      /* q_{k-1}(n) p'(n) */
                  qm1 = tmp;
                }
            }

          /* now set w(n) := q(n) * w(n) */
          for (i = 0; i < N; ++i)
            {
              double xi = ((double)i - half) / half;
              double qn = gsl_poly_eval(q, order + 1, xi);
              double *wn = gsl_vector_ptr(kernel, i);

              *wn *= qn;
            }
        }

      return GSL_SUCCESS;
    }
}
