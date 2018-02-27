/* movstat/movmad.c
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

/*
gsl_movstat_mad_alloc()
  Allocate a workspace for computing a moving median absolute deviation (MAD).
The workspace is set up to calculate a moving MAD with a given window size.
The window is defined by two integers H and J, where H is the number of
samples prior to the current sample, and J is the number of samples after
the current sample x_i:

W_i^{H,J} = {x_{i-H},...,x_i,...,x_{i+J}}

The total window size is K = H + J + 1

Inputs: K - total window size (H = J = K/2)

Return: pointer to workspace

Notes:
1) if K is even, the window size is rounded up to the next odd integer
to ensure equal number of samples before and after the current sample
*/

gsl_movstat_mad_workspace *
gsl_movstat_mad_alloc(const size_t K)
{
  const size_t H = K / 2;
  return gsl_movstat_mad_alloc2(H, H);
}
 
/*
gsl_movstat_mad_alloc2()
  Allocate a workspace for computing a moving median absolute deviation (MAD).
The workspace is set up to calculate a moving MAD with a given window size.
The window is defined by two integers H and J, where H is the number of
samples prior to the current sample, and J is the number of samples after
the current sample:

W_i^{H,J} = {x_{i-H},...,x_i,...,x_{i+J}}

The total window size is K = H + J + 1

Inputs: H - number of samples prior to current sample
        J - number of samples after current sample

Return: pointer to workspace
*/

gsl_movstat_mad_workspace *
gsl_movstat_mad_alloc2(const size_t H, const size_t J)
{
  gsl_movstat_mad_workspace *w;

  w = calloc(1, sizeof(gsl_movstat_mad_workspace));
  if (w == 0)
    {
      GSL_ERROR_NULL ("failed to allocate space for workspace", GSL_ENOMEM);
    }

  w->H = H;
  w->J = J;
  w->K = H + J + 1;

  w->median_workspace_p = gsl_movstat_median_alloc2(H, J);
  if (w->median_workspace_p == 0)
    {
      gsl_movstat_mad_free(w);
      GSL_ERROR_NULL ("failed to allocate space for median workspace", GSL_ENOMEM);
    }

  w->medacc_workspace_p = gsl_movstat_medacc_alloc(w->K);
  if (w->medacc_workspace_p == 0)
    {
      gsl_movstat_mad_free(w);
      GSL_ERROR_NULL ("failed to allocate space for mediator workspace", GSL_ENOMEM);
    }

  return w;
}

void
gsl_movstat_mad_free(gsl_movstat_mad_workspace * w)
{
  if (w->median_workspace_p)
    gsl_movstat_median_free(w->median_workspace_p);

  if (w->medacc_workspace_p)
    gsl_movstat_medacc_free(w->medacc_workspace_p);

  free(w);
}

/*
gsl_movstat_mad()
  Apply a moving MAD to an input vector

Inputs: x       - input vector, size n
        xmedian - (output) vector of median values of x, size n
                  xmedian_i = median of window centered on x_i
        xmad    - (output) vector of estimated standard deviations of x, size n
                  xmad_i = MAD of i-th window: median(|x_i - xmedian_i|)
        w       - workspace
*/

int
gsl_movstat_mad(const gsl_vector * x, gsl_vector * xmedian, gsl_vector * xmad,
                gsl_movstat_mad_workspace * w)
{
  const size_t n = x->size;
  const size_t k = w->K;

  if (n != xmedian->size)
    {
      GSL_ERROR("x and xmedian vectors must have same length", GSL_EBADLEN);
    }
  else if (n != xmad->size)
    {
      GSL_ERROR("x and xmad vectors must have same length", GSL_EBADLEN);
    }
  else if (k > n)
    {
      GSL_ERROR("window size must be less than input vector length", GSL_EBADLEN);
    }
  else
    {
      const int H = (int) w->H; /* number of samples to left of current sample */
      size_t i, j;

      /* first calculate median values of each window in x */
      gsl_movstat_median(GSL_MOVSTAT_EDGE_PADZERO, x, xmedian, w->median_workspace_p);

      /* loop over windows */
      for (i = 0; i < n; ++i)
        {
          double xmedi = gsl_vector_get(xmedian, i);
          int idx_base = (int) i - H; /* index of first sample in window i */

          /* loop over samples inside this window */
          for (j = 0; j < k; ++j)
            {
              int idx = idx_base + (int) j;

              /* check for samples near edge */
              if (idx >= 0 && idx < (int) n)
                {
                  double xj = gsl_vector_get(x, idx);
                  double absdev = fabs(xj - xmedi); /* absolute deviation for this sample inside window i */

                  /* add absolute deviation for this sample to median accumulator */
                  gsl_movstat_medacc_insert(absdev, w->medacc_workspace_p);
                }
              else
                {
                  /* zero pad window at edges xj = 0, absdev = |xj - xmedi| */
                  double absdev = fabs(xmedi); /* absolute deviation for this sample inside window i */
                  gsl_movstat_medacc_insert(absdev, w->medacc_workspace_p);
                }
            }

          gsl_vector_set(xmad, i, gsl_movstat_medacc_median(w->medacc_workspace_p));
        }

      return GSL_SUCCESS;
    }
}
