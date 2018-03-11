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
      int idx1, idx2;

      if (endtype == GSL_MOVSTAT_END_TRUNCATE)
        {
          /* reset minmaxacc since it might contain old samples from previous calls */
          gsl_movstat_minmaxacc_reset(w->minmaxacc_workspace_p);

          /* save last K - 1 samples of x for later (needed for in-place input/output) */
          idx1 = GSL_MAX(n - J - H, 0);
          idx2 = n - 1;
          for (i = idx1; i <= idx2; ++i)
            w->work[i - idx1] = gsl_vector_get(x, i);
        }
      else
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
            gsl_movstat_minmaxacc_insert(x1, w->minmaxacc_workspace_p);
        }

      /* process input vector and fill y(0:n - J - 1) */
      for (i = 0; i < n; ++i)
        {
          double xi = gsl_vector_get(x, i);
          int idx = i - J;

          gsl_movstat_minmaxacc_insert(xi, w->minmaxacc_workspace_p);

          if (idx >= 0)
            {
              double *y_min_ptr = gsl_vector_ptr(y_min, idx);
              double *y_max_ptr = gsl_vector_ptr(y_max, idx);

              gsl_movstat_minmaxacc_minmax(y_min_ptr, y_max_ptr, w->minmaxacc_workspace_p);
            }
        }

      if (endtype == GSL_MOVSTAT_END_TRUNCATE)
        {
          int wsize = n - GSL_MAX(n - J - H, 0); /* size of work array */

          /* need to fill y(n-J:n-1) using shrinking windows */
          idx1 = GSL_MAX(n - J, 0);
          idx2 = n - 1;

          for (i = idx1; i <= idx2; ++i)
            {
              double *y_min_ptr = gsl_vector_ptr(y_min, i);
              double *y_max_ptr = gsl_vector_ptr(y_max, i);
              int nsamp = n - GSL_MAX(i - H, 0); /* number of samples in this window */
              int j;

              gsl_movstat_minmaxacc_reset(w->minmaxacc_workspace_p);

              for (j = wsize - nsamp; j < wsize; ++j)
                gsl_movstat_minmaxacc_insert(w->work[j], w->minmaxacc_workspace_p);

              /* [ymini,ymaxi] = minmax [ work(i:K-2) ] */
              gsl_movstat_minmaxacc_minmax(y_min_ptr, y_max_ptr, w->minmaxacc_workspace_p);
            }
        }
      else
        {
          /* pad final windows and fill y(n-J:n-1) */
          for (i = 0; i < J; ++i)
            {
              int idx = n - J + i;

              gsl_movstat_minmaxacc_insert(xN, w->minmaxacc_workspace_p);

              if (idx >= 0)
                {
                  double *y_min_ptr = gsl_vector_ptr(y_min, idx);
                  double *y_max_ptr = gsl_vector_ptr(y_max, idx);

                  gsl_movstat_minmaxacc_minmax(y_min_ptr, y_max_ptr, w->minmaxacc_workspace_p);
                }
            }
        }

      return GSL_SUCCESS;
    }
}
