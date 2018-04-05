/* movstat/movstat_common.c
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

#ifndef __GSL_MOVSTAT_COMMON_C__
#define __GSL_MOVSTAT_COMMON_C__
 
#include <config.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_movstat.h>

/*
movstat_fill_window()
  Fill window for sample 'idx' from x using given end conditions

Inputs: endtype - how to handle end points
        idx     - index of center sample in window
        H       - number of samples left of center to include
        J       - number of samples right of center to include
        x       - input vector
        window  - (output) window of samples centered on x_{idx}

Return: number of samples in window
*/

static int
movstat_fill_window(const gsl_movstat_end_t endtype, const int idx, const int H, const int J,
                    const gsl_vector * x, double * window)
{
  const int n = x->size;
  int idx1, idx2, j;
  int wsize;

  if (endtype == GSL_MOVSTAT_END_TRUNCATE)
    {
      idx1 = GSL_MAX(idx - H, 0);
      idx2 = GSL_MIN(idx + J, n - 1);
    }
  else
    {
      idx1 = idx - H;
      idx2 = idx + J;
    }

  wsize = idx2 - idx1 + 1;

  /* fill sliding window */
  for (j = idx1; j <= idx2; ++j)
    {
      int widx = j - idx1;

      if (j < 0)
        {
          /* initial condition */
          if (endtype == GSL_MOVSTAT_END_PADZERO)
            window[widx] = 0.0;
          else if (endtype == GSL_MOVSTAT_END_PADVALUE)
            window[widx] = gsl_vector_get(x, 0);
        }
      else if (j >= n)
        {
          if (endtype == GSL_MOVSTAT_END_PADZERO)
            window[widx] = 0.0;
          else if (endtype == GSL_MOVSTAT_END_PADVALUE)
            window[widx] = gsl_vector_get(x, n - 1);
        }
      else
        {
          window[widx] = gsl_vector_get(x, j);
        }
    }

  wsize = idx2 - idx1 + 1;

  return wsize;
}

/*
movstat_apply()
  Apply moving window statistic to input vector. This is a generalized
routine to handle window endpoints and apply a given accumulator to
the input vector.

Inputs: endtype    - end point handling criteria
        x          - input vector, size n
        y          - output vector, size n
        acc_init   - initialize accumulator
        acc_insert - add a single sample to accumulator
        acc_delete - delete oldest sample from accumulator
        acc_get    - get current accumulated value
        w          - workspace

Notes:
1) It is allowed to have x = y for in-place moving statistics
*/

static int
movstat_apply(const gsl_movstat_end_t endtype,
              const gsl_vector * x,
              gsl_vector * y,
              int (*acc_init)(const size_t n, void * vstate),
              int (*acc_insert)(const double x, void * vstate),
              int (*acc_delete)(void * vstate),
              double (*acc_get)(const void * vstate),
              gsl_movstat_workspace * w)
{
  if (x->size != y->size)
    {
      GSL_ERROR("input and output vectors must have same length", GSL_EBADLEN);
    }
  else
    {
      const int n = (int) x->size;
      const int H = w->H; /* number of samples to left of current sample */
      const int J = w->J; /* number of samples to right of current sample */
      int i;
      double x1 = 0.0;    /* pad values for data edges */
      double xN = 0.0;

      /* initialize accumulator */
      (*acc_init)(w->K, w->state);

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
            (*acc_insert)(x1, w->state);
        }
      else if (acc_delete == NULL) /* FIXME XXX */
        {
          /* save last K - 1 samples of x for later (needed for in-place input/output) */
          int idx1 = GSL_MAX(n - J - H, 0);
          int idx2 = n - 1;

          for (i = idx1; i <= idx2; ++i)
            w->work[i - idx1] = gsl_vector_get(x, i);
        }

      /* process input vector and fill y(0:n - J - 1) */
      for (i = 0; i < n; ++i)
        {
          double xi = gsl_vector_get(x, i);
          int idx = i - J;

          (*acc_insert)(xi, w->state);

          if (idx >= 0)
            gsl_vector_set(y, idx, (*acc_get)(w->state));
        }

      if (endtype == GSL_MOVSTAT_END_TRUNCATE)
        {
          /* need to fill y(n-J:n-1) using shrinking windows */
          int idx1 = GSL_MAX(n - J, 0);
          int idx2 = n - 1;

          if (acc_delete == NULL)
            {
              int wsize = n - GSL_MAX(n - J - H, 0); /* size of work array */

              for (i = idx1; i <= idx2; ++i)
                {
                  int nsamp = n - GSL_MAX(i - H, 0); /* number of samples in this window */
                  int j;

                  (*acc_init)(w->K, w->state);

                  for (j = wsize - nsamp; j < wsize; ++j)
                    (*acc_insert)(w->work[j], w->state);

                  /* yi = acc_get [ work(i:K-2) ] */
                  gsl_vector_set(y, i, (*acc_get)(w->state));
                }
            }
          else
            {
              for (i = idx1; i <= idx2; ++i)
                {
                  if (i - H > 0)
                    {
                      /* delete oldest window sample as we move closer to edge */
                      (*acc_delete)(w->state);
                    }

                  /* yi = acc_get [ work(i:K-2) ] */
                  gsl_vector_set(y, i, (*acc_get)(w->state));
                }
            }
        }
      else
        {
          /* pad final windows and fill y(n-J:n-1) */
          for (i = 0; i < J; ++i)
            {
              int idx = n - J + i;

              (*acc_insert)(xN, w->state);

              if (idx >= 0)
                gsl_vector_set(y, idx, (*acc_get)(w->state));
            }
        }

      return GSL_SUCCESS;
    }
}

#endif /* __GSL_MOVSTAT_COMMON_C__ */
