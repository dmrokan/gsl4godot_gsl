/* movstat/movminmax.c
 *
 * Routines related to a moving window min/max
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
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_movstat.h>

#include "movstat_common.c"
#include "mmacc.c"

/*
gsl_movstat_minmax()
  Apply minmax filter to input vector

Inputs: endtype - end point handling criteria
        x       - input vector, size n
        y_min   - output vector of minimum values, size n
        y_max   - output vector of maximum values, size n
        w       - workspace
*/

int
gsl_movstat_minmax(const gsl_movstat_end_t endtype, const gsl_vector * x, gsl_vector * y_min, gsl_vector * y_max, gsl_movstat_workspace * w)
{
  if (x->size != y_min->size)
    {
      GSL_ERROR("input and y_min vectors must have same length", GSL_EBADLEN);
    }
  else if (x->size != y_max->size)
    {
      GSL_ERROR("input and y_max vectors must have same length", GSL_EBADLEN);
    }
  else
    {
      const int n = (int) x->size;
      const int H = w->H; /* number of samples to left of current sample */
      const int J = w->J; /* number of samples to right of current sample */
      int i;
      double x1 = 0.0;    /* pad values for data edges */
      double xN = 0.0;

      /* initialize sum accumulator */
      mmacc_init(w->K, w->state);

      /* pad initial window if necessary */
      if (endtype != GSL_MOVSTAT_END_TRUNCATE)
        {
          if (endtype == GSL_MOVSTAT_END_PADZERO)
            {
              x1 = 0.0;
              xN = 0.0;
            }
          else if (endtype == GSL_MOVSTAT_END_PADVALUE)
            {
              x1 = gsl_vector_get(x, 0);
              xN = gsl_vector_get(x, n - 1);
            }

          /* pad initial windows with H values */
          for (i = 0; i < H; ++i)
            mmacc_insert(x1, w->state);
        }

      /* process input vector and fill y(0:n - J - 1) */
      for (i = 0; i < n; ++i)
        {
          double xi = gsl_vector_get(x, i);
          int idx = i - J;

          mmacc_insert(xi, w->state);

          if (idx >= 0)
            {
              gsl_vector_set(y_min, idx, mmacc_min(w->state));
              gsl_vector_set(y_max, idx, mmacc_max(w->state));
            }
        }

      if (endtype == GSL_MOVSTAT_END_TRUNCATE)
        {
          /* need to fill y(n-J:n-1) using shrinking windows */
          int idx1 = GSL_MAX(n - J, 0);
          int idx2 = n - 1;

          for (i = idx1; i <= idx2; ++i)
            {
              if (i - H > 0)
                {
                  /* delete oldest window sample as we move closer to edge */
                  mmacc_delete(w->state);
                }

              /* yi = acc_get [ work(i:K-2) ] */
              gsl_vector_set(y_min, i, mmacc_min(w->state));
              gsl_vector_set(y_max, i, mmacc_max(w->state));
            }
        }
      else
        {
          /* pad final windows and fill y(n-J:n-1) */
          for (i = 0; i < J; ++i)
            {
              int idx = n - J + i;

              mmacc_insert(xN, w->state);

              if (idx >= 0)
                {
                  gsl_vector_set(y_min, idx, mmacc_min(w->state));
                  gsl_vector_set(y_max, idx, mmacc_max(w->state));
                }
            }
        }

      return GSL_SUCCESS;
    }
}

/*
gsl_movstat_min()
  Apply minimum filter to input vector

Inputs: endtype - end point handling criteria
        x       - input vector, size n
        y       - output vector of minimum values, size n
        w       - workspace
*/

int
gsl_movstat_min(const gsl_movstat_end_t endtype, const gsl_vector * x, gsl_vector * y, gsl_movstat_workspace * w)
{
  int status = movstat_apply(gsl_movstat_accum_min, endtype, x, y, w);
  return status;
}

/*
gsl_movstat_max()
  Apply maximum filter to input vector

Inputs: endtype - end point handling criteria
        x       - input vector, size n
        y       - output vector of maximum values, size n
        w       - workspace
*/

int
gsl_movstat_max(const gsl_movstat_end_t endtype, const gsl_vector * x, gsl_vector * y, gsl_movstat_workspace * w)
{
  int status = movstat_apply(gsl_movstat_accum_max, endtype, x, y, w);
  return status;
}
