/* filter/rmedian.c
 *
 * Contains routines related to the recursive median filter
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
#include <gsl/gsl_errno.h>

#include "mmacc.c"

gsl_filter_rmedian_workspace *
gsl_filter_rmedian_alloc(const size_t K)
{
  gsl_filter_rmedian_workspace *w;
  size_t state_size;

  w = calloc(1, sizeof(gsl_filter_rmedian_workspace));
  if (w == 0)
    {
      GSL_ERROR_NULL ("failed to allocate space for workspace", GSL_ENOMEM);
    }

  w->H = K / 2;
  w->K = 2*w->H + 1;

  state_size = mmacc_size(w->H + 1);

  w->state = malloc(state_size);
  if (w->state == NULL)
    {
      gsl_filter_rmedian_free(w);
      GSL_ERROR_NULL ("failed to allocate space for mmacc state", GSL_ENOMEM);
    }

  return w;
}

void
gsl_filter_rmedian_free(gsl_filter_rmedian_workspace * w)
{
  if (w->state)
    free(w->state);

  free(w);
}

int
gsl_filter_rmedian(const gsl_vector * x, gsl_vector * y, gsl_filter_rmedian_workspace * w)
{
  const size_t n = x->size;

  if (n != y->size)
    {
      GSL_ERROR("input and output vectors must have same length", GSL_EBADLEN);
    }
  else
    {
      double yprev = 0.0;
      double xmin, xmax, yval;
      size_t i;

      mmacc_init(w->H + 1, w->state);

      for (i = 0; i < n; ++i)
        {
          double xi = gsl_vector_get(x, i);

          mmacc_insert(xi, w->state);

          if (i >= w->H)
            {
              xmin = mmacc_min(w->state);
              xmax = mmacc_max(w->state);

              /* y_{i-H} = median [ yprev, xmin, xmax ] */
              if (yprev <= xmin)
                yval = xmin;
              else if (yprev <= xmax)
                yval = yprev;
              else
                yval = xmax;

              gsl_vector_set(y, i - w->H, yval);
              yprev = yval;
            }
        }

      /* now handle the last H elements of y */
      for (i = 0; i < w->H; ++i)
        {
          int idx = (int) n - (int) w->H + (int) i;

          /* zero pad input vector */
          mmacc_insert(0.0, w->state);

          if (idx >= 0)
            {
              xmin = mmacc_min(w->state);
              xmax = mmacc_max(w->state);

              /* y_{n-H+i} = median [ yprev, xmin, xmax ] */
              if (yprev <= xmin)
                yval = xmin;
              else if (yprev <= xmax)
                yval = yprev;
              else
                yval = xmax;

              gsl_vector_set(y, idx, yval);
              yprev = yval;
            }
        }

      return GSL_SUCCESS;
    }
}

int
gsl_filter_rmedian2(const gsl_vector * x, gsl_vector * y, gsl_filter_rmedian_workspace * w)
{
  const size_t n = x->size;

  if (n != y->size)
    {
      GSL_ERROR("input and output vectors must have same length", GSL_EBADLEN);
    }
  else
    {
      double yprev = 0.0;
      size_t i;

      for (i = 0; i < n; ++i)
        {
          double xi = gsl_vector_get(x, i);
          double xmax = xi;
          double xmin = xi;
          double yi;
          size_t j;

          /* find min/max [ x_{i+1}, ..., x_{i+H} ] */
          for (j = i + 1; j <= GSL_MIN(i + w->H, n - 1); ++j)
            {
              double xj = gsl_vector_get(x, j);
              if (xj > xmax)
                xmax = xj;
              if (xj < xmin)
                xmin = xj;
            }

          /* yi = median [ yprev, xmin, xmax ] */
          if (yprev <= xmin)
            yi = xmin;
          else if (yprev <= xmax)
            yi = yprev;
          else
            yi = xmax;

          gsl_vector_set(y, i, yi);
          yprev = yi;
        }

      return GSL_SUCCESS;
    }
}
