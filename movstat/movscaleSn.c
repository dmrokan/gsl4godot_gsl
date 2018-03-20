/* movstat/movscaleSn.c
 *
 * Compute moving "S_n" statistic from Croux and Rousseeuw, 1992
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
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_movstat.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>

#include "movstat_common.c"

/*
gsl_movstat_scaleSn()
  Calculate moving S_n statistic for input vector

Inputs: endtype - how to handle end points
        x       - input vector, size n
        xscale  - (output) vector of "S_n" statistics, size n
                  xscale_i = (S_n)_i for i-th window:
        w       - workspace
*/

int
gsl_movstat_scaleSn(const gsl_movstat_end_t endtype, const gsl_vector * x,
                    gsl_vector * xscale, gsl_movstat_workspace * w)
{
  if (x->size != xscale->size)
    {
      GSL_ERROR("x and xscale vectors must have same length", GSL_EBADLEN);
    }
  else
    {
      const int n = (int) x->size;
      const int H = (int) w->H; /* number of samples to left of current sample */
      const int J = (int) w->J; /* number of samples to right of current sample */
      double *window = w->work;
      int window_size;
      int i;

      for (i = 0; i < n; ++i)
        {
          double *xscalei = gsl_vector_ptr(xscale, i);
          gsl_vector_view v;

          /* fill window centered on x_i */
          window_size = movstat_fill_window(endtype, i, H, J, x, window);

          /* compute S_n for this window */
          v = gsl_vector_view_array(window, window_size);
          *xscalei = gsl_movstat_full_scaleSn(&v.vector);
        }

      return GSL_SUCCESS;
    }
}

/*
gsl_movstat_full_scaleSn()
  Calculate S_n statistic for input vector

Inputs: x - input vector, size n

Return: S_n(x)
*/

/*XXX*/
double
gsl_movstat_full_scaleSn(const gsl_vector * x)
{
  const int n = (int) x->size;
  double *work1 = malloc(n * sizeof(double));
  double *work2 = malloc(n * sizeof(double));
  double Sn;
  int j, k;

  for (j = 0; j < n; ++j)
    {
      double xj = gsl_vector_get(x, j);

      for (k = 0; k < n; ++k)
        {
          double xk = gsl_vector_get(x, k);
          work1[k] = fabs(xj - xk);
        }

      /* find med_k | x_j - x_k | */
      gsl_sort(work1, 1, n);
      work2[j] = gsl_stats_median_from_sorted_data(work1, 1, n);
    }

  /* find med_j { med_k | x_j - x_k | } */
  gsl_sort(work2, 1, n);
  Sn = gsl_stats_median_from_sorted_data(work2, 1, n);

  free(work1);
  free(work2);

  return Sn;
}
